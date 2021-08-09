#setup
rm(list=ls())
library(clusso)
library(MASS)
source("clustack.R")
source("helperfuncs.R")

#get arg from script
args<- commandArgs(TRUE)
startsim <- as.numeric(args[1])
setsims <- read.csv("setsims.csv")
rrsims <- read.csv("rrsims.csv")
utmJapan <- read.csv("utmJapan.csv")
padder <- function(maxcols, mincols, tabn){
    #browser()
    tabncols <- ncol(tabn)
    if (is.null(tabncols)){
        tabncols <- length(tabn)
    }
    pad <- rep(NA,maxcols-tabncols)
    if(!is.vector(tabn)){
        tabnrow <- nrow(tabn)
        tabnout <- cbind(tabn, matrix(rep(pad, tabnrow), nrow=tabnrow))
    }else{
        tabnout <- c(tabn, pad)
    }
    return(tabnout)
    
}
#################################################################################################
#################################################################################################
#SIM
#################################################################################################
#################################################################################################
#set global params
maxclust <- 15
nsimstep <- 1000
nsim <- 100

x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000

#set global
rMax <- 20 
Time <- 5
tim <- 3:5


#
#setsims <- setsims[startsim:(startsim+4),]


#Loop through batches of 10
for (m in startsim:(startsim+4)){
#for (m in 1:(1+4)){
    print(m)
#for (m in startsim:(startsim+9)){
#for (m in startsim:(startsim+1)){
    #print(m)
    #set sim
    rr <- as.matrix(rrsims[m,])
    YSIM <- as.vector(as.matrix(setsims[m, 1:1040]))
    Ex <- as.vector(as.matrix(setsims[m, 1041:2080]))
    cent <- setsims[m,2081]
    rad <- setsims[m, 2082]
    risk <- setsims[m, 2083]
    theta <- setsims[m, 2084]
    mod <- setsims[m, 2085]
    
    #create set of potential clusters based on distances
    potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
    n_uniq <- length(unique(potentialclusters$center))
    numCenters <- n_uniq
    #create giant sparse design matrix (single potential clusters)
    sparsematrix <- spacetimeMat(potentialclusters, numCenters, Time) 
    clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
    n <- length(x)
    
    outExp <- t(sparsematrix)%*%Ex
    outObs <- t(sparsematrix)%*%YSIM
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY LOCATION
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by location for each sim (apply)
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- glm(YSIM ~ as.factor(rep(c("1","2","3","4","5"), each=length(Ex)/Time)) + offset(log(Ex)),
                                                     family=quasipoisson)
        overdisp.est <- overdisp(offset_reg, overdispfloor = TRUE)
        
    }
    sim_superclust_loc.time <- system.time(sim_superclust_loc <- detectclusters(sparsematrix, Ex, YSIM, numCenters, Time, maxclust,
                                                                                byloc = TRUE, mod="poisson", 
                                                                                overdisp.est = overdisp.est))[[3]]
    print("finished stacking: by LOC")
    if(mod=="space"){
        sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                        "risk", risk, "_", "theta", as.character(theta),"_spaceonly_",m)
    } else{
        sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                        "risk", risk, "_", "theta", as.character(theta),
                        as.numeric(paste(tim, collapse = "")),"_",m)
    }
    print(filename <- paste0(sim.i,"_superclustLOC_",m,".RData"))
    #save .RData
    #save(sim_superclust_loc, file=filename)
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    id.bic_loc <- as.vector(unlist(sim_superclust_loc$selection.bic))
    id.aic_loc <- as.vector(unlist(sim_superclust_loc$selection.aic))
    #####################################################################################
    #####################################################################################
    ##################################
    #FP and POWER
    ##################################
    #Which PCs overlap true cluster?
    rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
    clusteroverlap<- rrbin_cluster%*%sparsematrix
    rrbin_inside <- ifelse(sparsematrix%*%t(clusteroverlap)!=0,1,0)
    
    ##############################
    if (sim_superclust_loc$selection.aic==0 ){
        incluster.aic <- NULL
        outcluster.aic <- NULL
        
    } else{
        #AIC
        ident.aic <- matrix(sim_superclust_loc$wtMAT[, 1:sim_superclust_loc$selection.aic], ncol=sim_superclust_loc$selection.aic)
        mat0.aic <- sapply(1:sim_superclust_loc$selection.aic, function(i) t(ident.aic[,i])%*%t(sparsematrix))
        matSum.aic <- rowSums(sapply(1:sim_superclust_loc$selection.aic, function(i) ifelse(mat0.aic[[i]]>=0.1,1,0)))
        mat.aic <- ifelse(matSum.aic!=0,1,0)
        
        #1) Did it find anything INSIDE the cluster?
        #start here
        incluster0.aic <- mat.aic%*%matrix(rrbin_inside, ncol=1)
        incluster.aic <- ifelse(incluster0.aic!=0,1,0)
        #2) Did it find anything OUTSIDE the cluster?
        rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
        #this should be everything that doesn't touch the cluster
        outcluster0.aic <- mat.aic%*%matrix(rrbin_outside, ncol=1)
        outcluster.aic <- ifelse(outcluster0.aic!=0,1,0)
    }  
    if (sim_superclust_loc$selection.bic==0 ){
        incluster.bic <- NULL
        outcluster.bic <- NULL
    } else{
        #BIC
        ident.bic <- matrix(sim_superclust_loc$wtMAT[, 1:sim_superclust_loc$selection.bic], ncol=sim_superclust_loc$selection.bic)
        mat0.bic <- sapply(1:sim_superclust_loc$selection.bic, function(i) t(ident.bic[,i])%*%t(sparsematrix))
        matSum.bic <- rowSums(sapply(1:sim_superclust_loc$selection.bic, function(i) ifelse(mat0.bic[[i]]>=0.1,1,0)))
        mat.bic <- ifelse(matSum.bic!=0,1,0)
        #1) Did it find anything INSIDE the cluster?
        #start here
        incluster0.bic <- mat.bic%*%matrix(rrbin_inside, ncol=1)
        incluster.bic <- ifelse(incluster0.bic!=0,1,0)
        #2) Did it find anything OUTSIDE the cluster?
        rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
        #this should be everything that doesn't touch the cluster
        outcluster0.bic <- mat.bic%*%matrix(rrbin_outside, ncol=1)
        outcluster.bic <- ifelse(outcluster0.bic!=0,1,0)
    }

    
    ##################################
    #Add sim results to table
    ##################################
    maxwLambda <- max(sim_superclust_loc$selection.bic, sim_superclust_loc$selection.aic)
    if(length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.bic,]))==0 & length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.aic,]))==0){
        wLambda.bic <- rep(1, 1040)
    } else {
        if(length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.bic,]))==0){
            wLambda.bic <- rep(1, 1040)   
            pad <-rep(NA, 1040*(maxwLambda-(sim_superclust_loc$selection.bic+1)), ncol=maxwLambda-(sim_superclust_loc$selection.bic+1))
        } else {
            wLambda.bic <- sim_superclust_loc$wLambda[1:sim_superclust_loc$selection.bic,]    
            pad <-rep(NA, 1040*(maxwLambda-sim_superclust_loc$selection.bic), ncol=maxwLambda-sim_superclust_loc$selection.bic)
        }
        wLambda.bic <- cbind(wLambda.bic, pad)
        
    }
    if(length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.aic,]))==0){
        wLambda.aic <- rep(1, 1040)
    } else {
        wLambda.aic <- sim_superclust_loc$wLambda[1:sim_superclust_loc$selection.aic,]
        pad <-rep(NA, 1040*(maxwLambda-sim_superclust_loc$selection.aic), ncol=maxwLambda-sim_superclust_loc$selection.aic)
        wLambda.aic <- cbind(wLambda.aic, pad)
    }
    
    
    
    tabn.loc <- rbind(c(IC="BIC",rad, risk, cent, theta,
                        timeperiod=as.numeric(paste(tim, collapse="")),
                        mod=mod, type="NA",time=sim_superclust_loc.time, method="loc", 
                        incluster.bic =  ifelse(length(incluster.bic)==0,0, incluster.bic),
                        outcluster.bic =  ifelse(length(outcluster.bic)==0,0, outcluster.bic),
                        iter=m,
                        numclusters = sim_superclust_loc$selection.bic,
                        wLambda.bic),
                      c(IC="AIC",rad, risk, cent, theta,
                        timeperiod=as.numeric(paste(tim, collapse="")),
                        mod=mod, type="NA",time=sim_superclust_loc.time, method="loc",
                        incluster.aic =  ifelse(length(incluster.aic)==0,0, incluster.aic),
                        outcluster.aic = ifelse(length(outcluster.aic)==0,0, outcluster.aic),
                        iter=m,
                        numclusters = sim_superclust_loc$selection.bic,
                        wLambda.aic))
    
    
    print("Finished clustack LOC")
    ###############################################################################################################################################
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY Potential Cluster
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by PC for each sim (apply)
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- glm(YSIM ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                  each=length(Ex)/Time)) + offset(log(Ex)),
                                             family=quasipoisson)
        overdisp.est <- overdisp(offset_reg, overdispfloor = TRUE)
        
    }
    
    sim_superclust_pc.time <- system.time(sim_superclust_pc<- detectclusters(sparsematrix, Ex, YSIM,numCenters, Time, maxclust,
                                                                             byloc = FALSE, model="poisson",
                                                                             overdisp.est = overdisp.est))[[3]]
    print(filename <- paste0(sim.i,"_superclustPC_",m,".RData"))
    #save(sim_superclust_pc, file=filename)
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    id.bic_pc <- as.vector(unlist(sim_superclust_pc$selection.bic))
    id.aic_pc <- as.vector(unlist(sim_superclust_pc$selection.aic))
    #####################################################################################
    #####################################################################################
    ##################################
    #FP and POWER
    ##################################
    #Which PCs overlap true cluster?
    rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
    clusteroverlap<- rrbin_cluster%*%sparsematrix
    rrbin_inside <- ifelse(sparsematrix%*%t(clusteroverlap)!=0,1,0)
    
    if (sim_superclust_pc$selection.aic==0 ){
        incluster.aic <- NULL
        outcluster.aic <- NULL
        #print("a")
        
    } else{
        #print("b")
        #AIC
        ident.aic <- matrix(sim_superclust_pc$wtMAT[, 1:sim_superclust_pc$selection.aic], ncol=sim_superclust_pc$selection.aic)
        mat0.aic <- sapply(1:sim_superclust_pc$selection.aic, function(i) t(ident.aic[,i])%*%t(sparsematrix))
        matSum.aic <- rowSums(sapply(1:sim_superclust_pc$selection.aic, function(i) ifelse(mat0.aic[[i]]>=0.1,1,0)))
        mat.aic <- ifelse(matSum.aic!=0,1,0)
        
        #1) Did it find anything INSIDE the cluster?
        #start here
        incluster0.aic <- mat.aic%*%matrix(rrbin_inside, ncol=1)
        incluster.aic <- ifelse(incluster0.aic!=0,1,0)
        #2) Did it find anything OUTSIDE the cluster?
        rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
        #this should be everything that doesn't touch the cluster
        outcluster0.aic <- mat.aic%*%matrix(rrbin_outside, ncol=1)
        outcluster.aic <- ifelse(outcluster0.aic!=0,1,0)
    }  
    if (sim_superclust_pc$selection.bic==0 ){
        incluster.bic <- NULL
        outcluster.bic <- NULL
    } else{
        #BIC
        ident.bic <- matrix(sim_superclust_pc$wtMAT[, 1:sim_superclust_pc$selection.bic], ncol=sim_superclust_pc$selection.bic)
        mat0.bic <- sapply(1:sim_superclust_pc$selection.bic, function(i) t(ident.bic[,i])%*%t(sparsematrix))
        matSum.bic <- rowSums(sapply(1:sim_superclust_pc$selection.bic, function(i) ifelse(mat0.bic[[i]]>=0.1,1,0)))
        mat.bic <- ifelse(matSum.bic!=0,1,0)
        #1) Did it find anything INSIDE the cluster?
        #start here
        incluster0.bic <- mat.bic%*%matrix(rrbin_inside, ncol=1)
        incluster.bic <- ifelse(incluster0.bic!=0,1,0)
        #2) Did it find anything OUTSIDE the cluster?
        rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
        #this should be everything that doesn't touch the cluster
        outcluster0.bic <- mat.bic%*%matrix(rrbin_outside, ncol=1)
        outcluster.bic <- ifelse(outcluster0.bic!=0,1,0)
    }
    ##################################
    #Add sim results to table
    ##################################
    ##################################
    maxwLambda <- max(sim_superclust_pc$selection.bic, sim_superclust_pc$selection.aic)
    if(length(as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.bic,]))==0 & length(as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.aic,]))==0){
        wLambda.bic <- rep(1, 1040)
    } else {
        if(length(as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.bic,]))==0){
            wLambda.bic <- rep(1, 1040)   
            pad <-rep(NA, 1040*(maxwLambda-(sim_superclust_pc$selection.bic+1)), ncol=maxwLambda-(sim_superclust_pc$selection.bic+1))
        } else {
            wLambda.bic <- sim_superclust_pc$wLambda[1:sim_superclust_pc$selection.bic,]    
            pad <-rep(NA, 1040*(maxwLambda-sim_superclust_pc$selection.bic), ncol=maxwLambda-sim_superclust_pc$selection.bic)
        }
        wLambda.bic <- cbind(wLambda.bic, pad)
        
    }
    if(length(as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.aic,]))==0){
        wLambda.aic <- rep(1, 1040)
    } else {
        wLambda.aic <- sim_superclust_pc$wLambda[1:sim_superclust_pc$selection.aic,]
        pad <-rep(NA, 1040*(maxwLambda-sim_superclust_pc$selection.aic), ncol=maxwLambda-sim_superclust_pc$selection.aic)
        wLambda.aic <- cbind(wLambda.aic, pad)
    }
    
    
    
    tabn.pc <- rbind(c(IC="BIC",rad, risk, cent, theta,
                       timeperiod=as.numeric(paste(tim, collapse="")),
                       mod=mod, type="NA",time=sim_superclust_pc.time, method="pc",
                       incluster.bic = ifelse(length(incluster.bic)==0,0, incluster.bic),
                       outcluster.bic = ifelse(length(outcluster.bic)==0,0, outcluster.bic),
                       iter = m,
                       numclusters = sim_superclust_loc$selection.bic,
                       wLambda.bic),
                     
                     c(IC="AIC",rad, risk, cent, theta,
                       timeperiod=as.numeric(paste(tim, collapse="")),
                       mod=mod,type="NA",time=sim_superclust_pc.time, method="pc",
                       incluster.aic = ifelse(length(incluster.aic)==0,0, incluster.aic),
                       outcluster.aic = ifelse(length(outcluster.aic)==0,0, outcluster.aic),
                       iter = m,
                       numclusters = sim_superclust_loc$selection.aic,
                       wLambda.aic))
    
    
    
    print("Finished clustack PC")
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY CLUSSO
    #####################################################################################
    #####################################################################################
    #####################################################################################
    YSIM1 <- as.vector(matrix(YSIM, ncol=Time, byrow = FALSE))
    
    ##############
    E01 <- as.vector(matrix(Ex, ncol=Time, byrow = FALSE))
    truth1 <- as.vector(matrix(YSIM, ncol=Time, byrow=FALSE))
    period1 <- as.vector(matrix(as.factor(rep(1:5, each=208)), byrow=FALSE))
    # init <- setVectors(japanbreastcancer$period, japanbreastcancer$expdeath, japanbreastcancer$death,
    #                    covars=NULL, Time=Time)
    # 
    # 
    # E01 <- as.vector(matrix(init$E0, ncol=Time, byrow=FALSE))
    # truth1 <- as.vector(matrix(init$Y.vec, ncol=Time, byrow = FALSE))
    # period1 <- as.vector(matrix(init$Year, ncol=Time, byrow = FALSE))
    id <- rep(1:208, times = 5)
    #create list of dataframes
    jbcSIM <- cbind.data.frame(expected = E01,observed = YSIM1,period = period1,id = id)
    jbcSIM <- jbcSIM[order(jbcSIM$id, jbcSIM$period),]
    jbcSIM$period <- as.factor(jbcSIM$period)
    #create list of dataframes
    
    #run clusso
    if(mod=="space"){
        sim_clusso.time <- system.time(sim_clusso <- clusso(df=jbcSIM,
                                                            expected = expected,
                                                            observed = observed,
                                                            timeperiod = period,
                                                            covars=FALSE,
                                                            id = id,
                                                            x = x,
                                                            y = y,
                                                            rMax = rMax,
                                                            utm=TRUE,
                                                            analysis = "space",
                                                            model = "poisson",
                                                            maxclust = maxclust))[[3]]
    } else {
        sim_clusso.time <- system.time(sim_clusso <- clusso(df=jbcSIM,
                                                            expected = expected,
                                                            observed = observed,
                                                            timeperiod = period,
                                                            covars=FALSE,
                                                            id = id,
                                                            x = x,
                                                            y = y,
                                                            rMax = rMax,
                                                            utm=TRUE,
                                                            analysis = "spacetime",
                                                            model = "poisson",
                                                            maxclust = maxclust))[[3]]
    }
    
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    vec <- rep(0, 208*Time)
    #position <- list(vec)[rep(1, nsim)]
    if(mod=="space"){
        #SPACE
        #background rates
        ##Quasi-P
        bgRate_i.bic.qp <- sapply(1:Time,
                                  function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.qp.s$E.qbic,ncol=Time)[,j])))))
        bgRate_i.aic.qp <- sapply(1:Time,
                                  function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.qp.s$E.qaic,ncol=Time)[,j])))))
        
        
        bgRate.bic.qp <- rep(bgRate_i.bic.qp, each = 208)
        bgRate.aic.qp <- rep(bgRate_i.aic.qp, each = 208)
        
        #detect
        ix.bic.qp <- which(abs(as.vector(sim_clusso$lassoresult.qp.s$E.qbic) - bgRate.bic.qp)>=10^-3)
        ix.aic.qp <- which(abs(as.vector(sim_clusso$lassoresult.qp.s$E.qaic) - bgRate.aic.qp)>=10^-3)
        
        ##Poisson
        bgRate_i.bic.p <- sapply(1:Time,
                                 function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.p.s$E.qbic,ncol=Time)[,j])))))
        bgRate_i.aic.p <- sapply(1:Time,
                                 function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.p.s$E.qaic,ncol=Time)[,j])))))
        
        bgRate.bic.p <- rep(bgRate_i.bic.p, each = 208)
        bgRate.aic.p <- rep(bgRate_i.aic.p, each = 208)
        
        #detect
        ix.bic.p <- which(abs(as.vector(sim_clusso$lassoresult.p.s$E.qbic) - bgRate.bic.p)>=10^-3)
        ix.aic.p <- which(abs(as.vector(sim_clusso$lassoresult.p.s$E.qaic) - bgRate.aic.p)>=10^-3)
        
        ############################
        ##FP and Power
        ############################
        #QP
        listpow.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                    sim_clusso$lassoresult.qp.s$selections$select.qbic,rr,
                                                    risk,nsim=1,Time, numCenters, pow=TRUE)
        listpow.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                     sim_clusso$lassoresult.qp.s$selections$select.qaic,rr,
                                                     risk,nsim=1,Time, numCenters, pow=TRUE)
        
        listfp.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                   sim_clusso$lassoresult.qp.s$selections$select.qbic,rr,
                                                   risk,nsim=1,Time, numCenters, pow=FALSE)
        listfp.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                    sim_clusso$lassoresult.qp.s$selections$select.qaic,rr,
                                                    risk,nsim=1,Time, numCenters, pow=FALSE)
        
        
        ###Poisson
        #QP
        listpow.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                   sim_clusso$lassoresult.qp.s$selections$select.qbic,rr,
                                                   risk,nsim=1,Time, numCenters, pow=TRUE)
        listpow.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                    sim_clusso$lassoresult.qp.s$selections$select.qaic,rr,
                                                    risk,nsim=1,Time, numCenters, pow=TRUE)
        listfp.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                  sim_clusso$lassoresult.p.s$selections$select.qbic,rr,
                                                  risk,nsim=1,Time, numCenters, pow=FALSE)
        listfp.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                   sim_clusso$lassoresult.p.s$selections$select.qaic,rr,
                                                   risk,nsim=1,Time, numCenters, pow=FALSE)
        
    } else{
        #SPACETIME
        #background rates
        ##Quasi-P
        bgRate_i.bic.qp <- sapply(1:Time,
                                  function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.qp.st$E.qbic,ncol=Time)[,j])))))
        bgRate_i.aic.qp <- sapply(1:Time,
                                  function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.qp.st$E.qaic,ncol=Time)[,j])))))
        
        
        bgRate.bic.qp <- rep(bgRate_i.bic.qp, each = 208)
        bgRate.aic.qp <- rep(bgRate_i.aic.qp, each = 208)
        
        #detect
        ix.bic.qp <- which(abs(as.vector(sim_clusso$lassoresult.qp.st$E.qbic) - bgRate.bic.qp)>=10^-3)
        ix.aic.qp <- which(abs(as.vector(sim_clusso$lassoresult.qp.st$E.qaic) - bgRate.aic.qp)>=10^-3)
        
        ##Poisson
        bgRate_i.bic.p <-sapply(1:Time,
                                function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.p.st$E.qbic,ncol=Time)[,j])))))
        bgRate_i.aic.p <- sapply(1:Time,
                                 function(j) as.numeric(names(which.max(table(matrix(sim_clusso$lassoresult.p.st$E.qaic,ncol=Time)[,j])))))
        
        bgRate.bic.p <- rep(bgRate_i.bic.p, each = 208)
        bgRate.aic.p <- rep(bgRate_i.aic.p, each = 208)
        
        #detect
        ix.bic.p <- which(abs(as.vector(sim_clusso$lassoresult.p.st$E.qbic) - bgRate.bic.p)>=10^-3)
        ix.aic.p <- which(abs(as.vector(sim_clusso$lassoresult.p.st$E.qaic) - bgRate.aic.p)>=10^-3)
        
        ############################
        ##FP and Power
        ############################
        #QP
        listpow.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                    sim_clusso$lassoresult.qp.st$selections$select.qbic,rr,
                                                    risk,nsim=1,Time, numCenters, pow=TRUE)
        listpow.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                     sim_clusso$lassoresult.qp.st$selections$select.qaic,rr,
                                                     risk,nsim=1,Time, numCenters, pow=TRUE)
        
        listfp.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                   sim_clusso$lassoresult.qp.st$selections$select.qbic,rr,
                                                   risk,nsim=1,Time, numCenters, pow=FALSE)
        listfp.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                    sim_clusso$lassoresult.qp.st$selections$select.qaic,rr,
                                                    risk,nsim=1,Time, numCenters, pow=FALSE)
        
        ###Poisson
        #QP
        listpow.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                   sim_clusso$lassoresult.qp.st$selections$select.qbic,rr,
                                                   risk,nsim=1,Time, numCenters, pow=TRUE)
        listpow.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                    sim_clusso$lassoresult.qp.st$selections$select.qaic,rr,
                                                    risk,nsim=1,Time, numCenters, pow=TRUE)
        listfp.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                  sim_clusso$lassoresult.p.st$selections$select.qbic,rr,
                                                  risk,nsim=1,Time, numCenters, pow=FALSE)
        listfp.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                   sim_clusso$lassoresult.p.st$selections$select.qaic,rr,
                                                   risk,nsim=1,Time, numCenters, pow=FALSE)
        
    }
    
    ##################################
    #Add sim results to table
    ##################################
    if(mod == "space"){
        rrest.bic.qp <- sim_clusso$lassoresult.qp.s$E.qbic
        rrest.aic.qp <- sim_clusso$lassoresult.qp.s$E.qaic
        rrest.bic.p <- sim_clusso$lassoresult.p.s$E.qbic
        rrest.aic.p <- sim_clusso$lassoresult.p.s$E.qaic
    } else{
        rrest.bic.qp <- sim_clusso$lassoresult.qp.st$E.qbic
        rrest.aic.qp <- sim_clusso$lassoresult.qp.st$E.qaic
        rrest.bic.p <- sim_clusso$lassoresult.p.st$E.qbic
        rrest.aic.p <- sim_clusso$lassoresult.p.st$E.qaic
    }
    
    if(mod == "spacetime"){
        tabn.clusso <- rbind(c(IC="BIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod,  type = "QP",time=sim_clusso.time, method="clusso",
                               incluster.bic = as.vector(listpow.bic.qp),
                               outcluster.bic = as.vector(listfp.bic.qp),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.qp.st$numclust.qbic,
                               rrest.bic.qp ),
                             c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod,type = "QP",time=sim_clusso.time, method="clusso",
                               incluster.aic = as.vector(listpow.aic.qp),
                               outcluster.aic = as.vector(listfp.aic.qp),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.qp.st$numclust.qaic,
                               rrest.aic.qp ),
                             c(IC="BIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                               incluster.bic = as.vector(listpow.bic.p),
                               outcluster.bic = as.vector(listfp.bic.p),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.p.st$numclust.qbic,
                               rrest.bic.p ),
                             c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                               incluster.aic = as.vector(listpow.aic.p),
                               outcluster.aic = as.vector(listfp.aic.p),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.p.st$numclust.qaic,
                               rrest.aic.p))
        
    } else{
        tabn.clusso <- rbind(c(IC="BIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod,  type = "QP",time=sim_clusso.time, method="clusso",
                               incluster.bic = as.vector(listpow.bic.qp),
                               outcluster.bic = as.vector(listfp.bic.qp),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.qp.s$numclust.qbic,
                               rrest.bic.qp ),
                             c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod,type = "QP",time=sim_clusso.time, method="clusso",
                               incluster.aic = as.vector(listpow.aic.qp),
                               outcluster.aic = as.vector(listfp.aic.qp),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.qp.s$numclust.qaic,
                               rrest.aic.qp ),
                             c(IC="BIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                               incluster.bic = as.vector(listpow.bic.p),
                               outcluster.bic = as.vector(listfp.bic.p),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.p.s$numclust.qbic,
                               rrest.bic.p ),
                             c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                               incluster.aic = as.vector(listpow.aic.p),
                               outcluster.aic = as.vector(listfp.aic.p),
                               iter = m,
                               numclusters = sim_clusso$lassoresult.p.s$numclust.qaic,
                               rrest.aic.p))
    }
   
    
    print("Finished clusso")
    
    ###################################################################################
    ###################################################################################
    ###################################################################################
    #STEPWISE SPATIAL SCAN VARIANT
    ###################################################################################
    ###################################################################################
    ###################################################################################
    sim_stepscan.time <- system.time(sim_stepscan <- stepscan(YSIM, Ex, Time=5, sparsematrix, nsim=nsimstep, maxclust=maxclust))[[3]]
    
    simid_stepscan <- paste0(sim.i,"stepscan")
    
    ########################################################################
    #STEPWISE SCAN ANALYSIS
    ########################################################################
    #power/FP
    ##########################################
    numclustersid <- which(sim_stepscan$pvals>0.05)-1
    ixids <- step_clusterix(sparsematrix, sim_stepscan, numclustersid=numclustersid)
    #pow/fp
    outpow.stepscan<- spatscanfs_prob_clusteroverlap(sim_stepscan,ixids, numclustersid ,sparsematrix, rr, risk,pow=TRUE)
    outfp.stepscan <- spatscanfs_prob_clusteroverlap(sim_stepscan,ixids, numclustersid ,sparsematrix, rr, risk,pow=FALSE)
    #cell detection
    if(numclustersid==0){
        rrest.mc <- rep(1,1040)    
    } else{
        rrest.mc <- sim_stepscan$clusters[, max(which(sim_stepscan$pvals<0.05))]
    }
    
    
    tabn.stepscan <- c(IC="MC",rad, risk, cent, theta,
                       timeperiod=as.numeric(paste(tim, collapse="")),
                       mod=mod,type="NA", time =  sim_stepscan.time ,method = "stepscan",
                       incluster.mc = as.vector(outpow.stepscan),
                       outcluster.mc = as.vector(outfp.stepscan),
                       iter = m,
                       numclusters = numclustersid,
                       rrest.mc)
    print("Finished stepwise scan")
    
    ###################################################################################
    ###################################################################################
    ###################################################################################
    #FORWARD STAGEWISE
    ###################################################################################
    ###################################################################################
    ###################################################################################
    clusnum=NULL
    start=end=NULL
    for(i in 1:Time)
        for(j in 0:(Time-i)){
            start=c(start,i)
            end=c(end,i+j)
        }
    
    for(i in 1:Time) #for each of the time intervals
        for(j in 1:(Time+1-i))  #for each j in each time interval plus the next time interval minues1-5  of the time interval
            clusnum=c(clusnum,j*clusters$n)
    #clusters contains data of the center, x, y (from utm), radius, n, and last
    #number of potential clusters
    sd=sqrt((clusnum-(clusnum)^2/n/Time)/(n*Time-1))     # s.d of predictors
    #this has 134,440 items; for each element of 'clusnum'; clusters$n has 8,960 elements
    ###### Parameters #######
    delta = 0.001           ## stagewise step size
    max.ndelta = 15000      ## max number of deltas
    n = length(x)    ## number of cells=N
    
    ind = (1:n)[!duplicated(cbind(x,y))] ## extract unique elements
    
    tmp = as.matrix(dist(cbind(x,y)))[ind,]
    
    last=apply(tmp,1,function(x,r) order(x)[1:sum(x<=r)],r=rMax)
    nc = unlist2(lapply(last,length))
    last = unlist2(last)
    
    r=unlist2(apply(tmp,1, function(x,r) { sort(x[x<=r]) },r=rMax))
    
    clusters=data.frame(center=rep(ind,nc),
                        x=x[rep(ind,nc)],y=y[rep(ind,nc)],
                        r=r, 
                        n=unlist(lapply(nc,seq)),
                        last=last)
    
    
    sim_stage.time <- system.time(sim_stage <- cluster_model(delta,YSIM,Ex,sd, Time))[[3]]
    simid_stage <- paste0(sim.i,"fstagewise")
    
    ########################################################################
    #FORWARD STAGEWISE ANALYSIS
    ########################################################################
    ##########################################
    #power/FP
    ##########################################
    #pow/fp
    outfp.stage <- forwardstage_prob_clusteroverlap(sim_stage,sparsematrix, rr, risk,pow=FALSE)
    outfp.stage.bic <- outfp.stage$fp.bic
    outfp.stage.aic <- outfp.stage$fp.aic
    
    outpow.stage <- forwardstage_prob_clusteroverlap(sim_stage,sparsematrix, rr, risk,pow=TRUE)
    outpow.stage.bic <- outpow.stage$power.bic
    outpow.stage.aic <- outpow.stage$power.aic
    
    
    tabn.fstagewise <- rbind(c(IC="BIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod,  type="NA", time = sim_stage.time  ,method = "fstagewise",
                               incluster.bic = as.vector(outpow.stage.bic),
                               outcluster.bic = as.vector(outfp.stage.bic),
                               iter = m,
                               numclustersid = sim_stage$K[sim_stage$n.bic],
                               as.vector(sim_stage$RRbic)),
                             c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type="NA", time =  sim_stage.time ,method = "fstagewise" ,
                               incluster.bic = as.vector(outpow.stage.bic),
                               outcluster.bic = as.vector(outfp.stage.bic),
                               iter = m,
                               numclustersid = sim_stage$K[sim_stage$n.aic],
                               as.vector(sim_stage$RRaic)))
    
    print("Finished forward stagewise")
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    ##WRITE TO CSV
    #####################################################################################
    #####################################################################################
    #####################################################################################
    maxcols <- max(ncol(tabn.loc), ncol(tabn.pc), ncol(tabn.clusso), ncol(tabn.fstagewise), ncol(tabn.stepscan))
    mincols <- min(ncol(tabn.loc), ncol(tabn.pc), ncol(tabn.clusso), ncol(tabn.fstagewise), ncol(tabn.stepscan))
    
    tabn.locpad <- padder(maxcols, mincols,tabn.loc)
    tabn.pcpad <- padder(maxcols, mincols,tabn.pc)
    tabn.clussopad <- padder(maxcols, mincols,tabn.clusso)
    tabn.stepscanpad <- padder(maxcols, mincols,tabn.stepscan)
    tabn.fstagewisepad <- padder(maxcols, mincols,tabn.fstagewise)
    
    out <- rbind(tabn.locpad, tabn.pcpad, tabn.clussopad, tabn.stepscanpad, tabn.fstagewisepad)
    
    
    #out <- rbind(tabn.loc, cbind(tabn.pc, rbind(pad2, pad2)), tabn.clusso, c(tabn.stepscan, pad2), tabn.fstagewise) 
    if(mod=="space"){
        out <- cbind(stsmodel=rep("space", times=nrow(out)), out)
    }else{
        out <- cbind(stsmodel=rep("ST", times=nrow(out)), out)
    }
    write.csv(out, file=paste0("theta",theta,"_singlecluster","_iter",m,"_startsim",startsim,".csv"), row.names = TRUE)    
    
}
