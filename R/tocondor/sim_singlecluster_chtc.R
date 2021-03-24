#setup
library(clusso)
library(MASS)
source("clustack.R")
source("helperfuncs.R")

#get arg from script
args<- commandArgs(TRUE)
startsim <- args[1]
setsims <- read.csv("setsims.csv")
rrsims <- read.csv("rrsims.csv")
utmJapan <- read.csv("utmJapan.csv")
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

#Loop through batches of 5
for (m in startsim:(startsim+4)){
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
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
        
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
    if(id.bic_loc!=0){
        outbic.loc <- calcbounds(id.aic_loc, IC="bic", sim_superclust_loc, 
                                 byloc = TRUE,Ex, YSIM, target="cluster", 
                                 sparsemat = sparsematrix)
        #outbic.loc <- calcbounds(id.bic_loc, IC="bic", sim_superclust_loc, bylocation = TRUE,outExp)
    } else {
        print("No clusters identified: BIC")
    }
    if(id.aic_loc!=0){
        
        if(id.bic_loc==id.aic_loc){
            outaic.loc <- outbic.loc
        } else {
            outaic.loc <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc, 
                                     byloc = TRUE,Ex, YSIM, target="cluster", 
                                     sparsemat = sparsematrix)
            #outaic.loc <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc, bylocation = TRUE,outExp)
        }
    } else {
        print("No clusters identified: AIC")
    }
    ##################################
    #FP and POWER
    ##################################
    #Which PCs overlap true cluster?
    rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
    clusteroverlap<- rrbin_cluster%*%sparsematrix
    rrbin_inside <- ifelse(sparsematrix%*%t(clusteroverlap)!=0,1,0)
    
    ##############################
    #BIC
    ident.bic <- matrix(sim_superclust_loc$wtMAT[, sim_superclust_loc$selection.bic], ncol=sim_superclust_loc$selection.bic)
    mat0.bic <- sapply(1:sim_superclust_loc$selection.bic, function(i) t(ident.bic[,i])%*%t(sparsematrix))
    matSum.bic <- rowSums(sapply(1:sim_superclust_loc$selection.bic, function(i) ifelse(mat0.bic[[i]]>=0.1,1,0)))
    mat.bic <- ifelse(matSum.bic!=0,1,0)
    
    #AIC
    ident.aic <- matrix(sim_superclust_loc$wtMAT[, 1:sim_superclust_loc$selection.aic], ncol=sim_superclust_loc$selection.aic)
    mat0.aic <- sapply(1:sim_superclust_loc$selection.aic, function(i) t(ident.aic[,i])%*%t(sparsematrix))
    matSum.aic <- rowSums(sapply(1:sim_superclust_loc$selection.aic, function(i) ifelse(mat0.aic[[i]]>=0.1,1,0)))
    mat.aic <- ifelse(matSum.aic!=0,1,0)
    
    #1) Did it find anything INSIDE the cluster?
    #start here
    incluster0.bic <- mat.bic%*%matrix(rrbin_inside, ncol=1)
    incluster0.aic <- mat.aic%*%matrix(rrbin_inside, ncol=1)
    
    incluster.bic <- ifelse(incluster.bic!=0,1,0)
    incluster.aic <- ifelse(incluster.aic!=0,1,0)
    #2) Did it find anything OUTSIDE the cluster?
    rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
    #this should be everything that doesn't touch the cluster
    outcluster0.bic <- mat.bic%*%matrix(rrbin_outside, ncol=1)
    outcluster0.aic <- mat.aic%*%matrix(rrbin_outside, ncol=1)
    outcluster.bic <- ifelse(outcluster0.bic!=0,1,0)
    outcluster.aic <- ifelse(outcluster0.aic!=0,1,0)
    
    ##################################
    #Add sim results to table
    ##################################
    maxwLambda <- max(sim_superclust_loc$selection.bic, sim_superclust_loc$selection.aic)
    if(length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.bic,]))==0){
        wLambda.bic <- rep(1, 1040)
    } else {
        wLambda.bic <- sim_superclust_loc$wLambda[1:sim_superclust_loc$selection.bic,]
        pad <-rep(NA, 1040*(maxwLambda-sim_superclust_loc$selection.bic), ncol=maxwLambda-sim_superclust_loc$selection.bic)
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
                        wLambda.bic),
                      c(IC="AIC",rad, risk, cent, theta,
                        timeperiod=as.numeric(paste(tim, collapse="")),
                        mod=mod, type="NA",time=sim_superclust_loc.time, method="loc",
                        incluster.aic =  ifelse(length(incluster.aic)==0,0, incluster.aic),
                        outcluster.aic = ifelse(length(outcluster.aic)==0,0, outcluster.aic),
                        iter=m,
                        wLambda.aic))
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(id.bic_loc==0) & isTRUE(id.aic_loc==0)){
        print("No clusters identified")
        bounds.loc <- as.data.frame(matrix(rep(NA*41), ncol=41))
        
    } else{
        if(id.bic_loc==0){
            fill.bic <- matrix(rep(NA,41*(length(which(id.bic_loc!=0)))),
                               nrow=length(which(id.bic_loc!=0)))
        } else {
            fill.bic <- NULL
        }
        if(id.aic_loc==0){
            fill.aic <- matrix(rep(NA,41*(length(which(id.aic_loc!=0)))),
                               nrow=length(which(id.aic_loc!=0)))
        } else {
            fill.aic <- NULL
        }
        #######################
        bounds.aic <-  matrix(rep(NA, 41*sim_superclust_loc$selection.aic),nrow=sim_superclust_loc$selection.aic)
        bounds.bic <- matrix(rep(NA, 41*sim_superclust_loc$selection.bic),nrow=sim_superclust_loc$selection.bic)
        for(i in 1:sim_superclust_loc$selection.bic){
                bounds.bic[i,] <- cbind(matrix(unlist(outbic.loc[[i]]$outnonma$nonma.theta),
                                              byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outbuck.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outmaw2.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outmata.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc[[i]]$outmataTlog.theta), byrow=TRUE, ncol=3),
                                       fill.bic,
                                       
                                       outbic.loc[[i]]$outnonma.time,
                                       outbic.loc[[i]]$outnonmaTlog.time,
                                       outbic.loc[[i]]$outnonma_asymp.time,
                                       outbic.loc[[i]]$outnonma_asympTlog.time,
                                       outbic.loc[[i]]$outbuck.theta.time,
                                       outbic.loc[[i]]$outbuckTlog.theta.time,
                                       outbic.loc[[i]]$outmaw2.theta.time,
                                       outbic.loc[[i]]$outmaw2Tlog.theta.time,
                                       outbic.loc[[i]]$outmata.theta.time,
                                       outbic.loc[[i]]$outmataTlog.theta.time, 
                                       "BIC")

        }
        for(i in 1:sim_superclust_loc$selection.aic){
            bounds.aic[i,] <- cbind(matrix(unlist(outaic.loc[[i]]$outnonma$nonma.theta),
                                                        byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outbuck.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outmaw2.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outmata.theta), byrow=TRUE, ncol=3),
                                                 matrix(unlist(outaic.loc[[i]]$outmataTlog.theta), byrow=TRUE, ncol=3),
                                    fill.aic,
                                    
                                    outaic.loc[[i]]$outnonma.time,
                                    outaic.loc[[i]]$outnonmaTlog.time,
                                    outaic.loc[[i]]$outnonma_asymp.time,
                                    outaic.loc[[i]]$outnonma_asympTlog.time,
                                    outaic.loc[[i]]$outbuck.theta.time,
                                    outaic.loc[[i]]$outbuckTlog.theta.time,
                                    outaic.loc[[i]]$outmaw2.theta.time,
                                    outaic.loc[[i]]$outmaw2Tlog.theta.time,
                                    outaic.loc[[i]]$outmata.theta.time,
                                    outaic.loc[[i]]$outmataTlog.theta.time,
                                    "AIC")
            
        }
        
        bounds.loc <- rbind.data.frame(bounds.bic, bounds.aic)
    }
    
    bounds.loc$risk <- risk
    bounds.loc$radius <- rad
    
    bounds.loc$iter <- m
    bounds.loc$theta <- theta
    bounds.loc$cent <- cent
    
    bounds.loc$select_orig.bic <- sim_superclust_loc$selection.bic
    bounds.loc$select_orig.aic <- sim_superclust_loc$selection.aic
    bounds.loc$simID <- 1:nrow(bounds.loc)
    bounds.loc$method <- "LOC"
    
    names(bounds.loc) <- c("nonma.LB", "clusterMA.1", "nonma.UB",
                           "nonmaTlog.LB", "clusterMA.2", "nonmaTlog.UB",
                           "nonma_asymp.LB", "clusterMA.3", "nonma_asymp.UB",
                           "nonma_asympTlog.LB", "clusterMA.4", "nonma_asympTlog.UB",
                           "buck.LB", "clusterMA.5", "buck.UB",
                           "bucklog.LB", "clusterMA.6", "bucklog.UB",
                           "maw2.LB", "clusterMA.7", "maw2.UB",
                           "maw2log.LB", "clusterMA.8", "maw2log.UB",
                           "mata.LB", "clusterMA.9", "mata.UB",
                           "matalog.LB", "clusterMA.10", "matalog.UB",
                           
                           "nonma.time", "nonmaTlog.time", 
                           "nonma_asymp.time", "nonma_asympTlog.time", 
                           "buck.time","bucklog.time",
                           "maw2.time", "maw2log.time",
                           "mata.time", "matalog.time",
                           "IC",
                           "risk","rad","iter","theta", "cent",
                           "select_orig.bic", "select_orig.aic", "simID", "method")
    
    
    #table.bounds.loc <- rbind(table.bounds.loc, bounds.loc)
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
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
        
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
    if(id.bic_pc!=0){
        outbic.pc <- calcbounds(id.bic_pc, IC="bic", sim_superclust_pc, bylocation = FALSE,outExp)
    } else {
        print("No clusters identified: BIC")
    }
    if(id.aic_pc!=0){
        
        if(id.bic_pc==id.aic_pc){
            outaic.pc <- outbic.pc
        } else {
            outaic.pc <- calcbounds(id.aic_pc, IC="aic", sim_superclust_pc, bylocation = FALSE,outExp)
        }
    } else {
        print("No clusters identified: AIC")
    }
    ##################################
    #FP and POWER
    ##################################
    #Which PCs overlap true cluster?
    rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
    clusteroverlap<- rrbin_cluster%*%sparsematrix
    rrbin_inside <- ifelse(sparsematrix%*%t(clusteroverlap)!=0,1,0)
    
    ##############################
    ##############################
    # #Which PCs overlap true cluster?
    #what was identified in each sim by IC
    ident.bic <-sim_superclust_pc$wtMAT[, sim_superclust_pc$selection.bic]
    mat.bic <-t(ident.bic)%*%t(sparsematrix)
    
    ident.aic <- sim_superclust_pc$wtMAT[, sim_superclust_pc$selection.aic]
    mat.aic <- t(ident.aic)%*%t(sparsematrix)
    
    
    
    #1) Did it find anything INSIDE the cluster?
    incluster.bic <- mat.bic%*%matrix(rrbin_inside, ncol=1)
    incluster.aic <- mat.aic%*%matrix(rrbin_inside, ncol=1)
    
    
    #2) Did it find anything OUTSIDE the cluster?
    #this should be everything that doesn't touch the cluster
    outcluster.bic <- mat.bic%*%matrix(rrbin_outside, ncol=1)
    outcluster.aic <- mat.aic%*%matrix(rrbin_outside, ncol=1)
    
    
    ##################################
    #Add sim results to table
    ##################################
    ##################################
    if(length(as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.bic,]))==0){
        wLambda.bic <- rep(1, 1040)
    } else {
        wLambda.bic <- as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.bic,])
    }
    if(length(as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.aic,]))==0){
        wLambda.aic <- rep(1, 1040)
    } else {
        wLambda.aic <- as.vector(sim_superclust_pc$wLambda[sim_superclust_pc$selection.aic,])
    }
    
    
    tabn.pc <- rbind(c(IC="BIC",rad, risk, cent, theta,
                       timeperiod=as.numeric(paste(tim, collapse="")),
                       mod=mod, type="NA",time=sim_superclust_pc.time, method="pc",
                       incluster.bic = ifelse(length(incluster.bic@x)==0,0, incluster.bic@x),
                       outcluster.bic = ifelse(length(outcluster.bic@x)==0,0, outcluster.bic@x),
                       iter = m,
                       wLambda.bic),
                     
                     c(IC="AIC",rad, risk, cent, theta,
                       timeperiod=as.numeric(paste(tim, collapse="")),
                       mod=mod,type="NA",time=sim_superclust_pc.time, method="pc",
                       incluster.aic = ifelse(length(incluster.aic@x)==0,0, incluster.aic@x),
                       outcluster.aic = ifelse(length(outcluster.aic@x)==0,0, outcluster.aic@x),
                       iter = m,
                       wLambda.aic))
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(id.bic_pc==0) & isTRUE(id.aic_pc==0)){
        print("No clusters identified")
        bounds.pc <- as.data.frame(matrix(rep(NA*41), ncol=41))
    } else{
        #how many empty rows?
        if(id.bic_pc==0){
            fill.bic <- matrix(rep(NA,41*(length(which(id.bic_pc!=0)))),
                               nrow=length(which(id.bic_pc!=0)))
        } else {
            fill.bic <- NULL
        }
        if(id.aic_pc==0){
            fill.aic <- matrix(rep(NA,41*(length(which(id.aic_pc!=0)))),
                               nrow=length(which(id.aic_pc!=0)))
        } else {
            fill.aic <- NULL
        }
        #######################
        bounds.aic <-  matrix(rep(NA, 41*sim_superclust_pc$selection.aic),nrow=sim_superclust_pc$selection.aic)
        bounds.bic <- matrix(rep(NA, 41*sim_superclust_pc$selection.bic),nrow=sim_superclust_pc$selection.bic)
        for(i in 1:sim_superclust_pc$selection.bic){
            bounds.bic[i,] <- cbind(matrix(unlist(outbic.pc[[i]]$outnonma$nonma.theta),
                                           byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outbuck.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outmaw2.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outmata.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outbic.pc[[i]]$outmataTlog.theta), byrow=TRUE, ncol=3),
                                    fill.bic,
                                    
                                    outbic.pc[[i]]$outnonma.time,
                                    outbic.pc[[i]]$outnonmaTlog.time,
                                    outbic.pc[[i]]$outnonma_asymp.time,
                                    outbic.pc[[i]]$outnonma_asympTlog.time,
                                    outbic.pc[[i]]$outbuck.theta.time,
                                    outbic.pc[[i]]$outbuckTlog.theta.time,
                                    outbic.pc[[i]]$outmaw2.theta.time,
                                    outbic.pc[[i]]$outmaw2Tlog.theta.time,
                                    outbic.pc[[i]]$outmata.theta.time,
                                    outbic.pc[[i]]$outmataTlog.theta.time, 
                                    "BIC")
            
        }
        for(i in 1:sim_superclust_pc$selection.aic){
            bounds.aic[i,] <- cbind(matrix(unlist(outaic.pc[[i]]$outnonma$nonma.theta),
                                           byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outbuck.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outmaw2.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outmata.theta), byrow=TRUE, ncol=3),
                                    matrix(unlist(outaic.pc[[i]]$outmataTlog.theta), byrow=TRUE, ncol=3),
                                    fill.aic,
                                    
                                    outaic.pc[[i]]$outnonma.time,
                                    outaic.pc[[i]]$outnonmaTlog.time,
                                    outaic.pc[[i]]$outnonma_asymp.time,
                                    outaic.pc[[i]]$outnonma_asympTlog.time,
                                    outaic.pc[[i]]$outbuck.theta.time,
                                    outaic.pc[[i]]$outbuckTlog.theta.time,
                                    outaic.pc[[i]]$outmaw2.theta.time,
                                    outaic.pc[[i]]$outmaw2Tlog.theta.time,
                                    outaic.pc[[i]]$outmata.theta.time,
                                    outaic.pc[[i]]$outmataTlog.theta.time,
                                    "AIC")
            
        }
        
        bounds.pc <- rbind.data.frame(bounds.bic, bounds.aic)
        
     
    }
    
    bounds.pc$risk <- risk
    bounds.pc$radius <- rad
    
    bounds.pc$iter <- m
    bounds.pc$theta <- theta

    #ADD
    bounds.pc$pop <- cent
    
    bounds.pc$select_orig.bic <- sim_superclust_pc$selection.bic
    bounds.pc$select_orig.aic <- sim_superclust_pc$selection.aic
    bounds.pc$simID <- 1:nrow(bounds.pc)
    bounds.pc$method <- "PC"
    
    names(bounds.loc) <- c("nonma.LB", "clusterMA.1", "nonma.UB",
                           "nonmaTlog.LB", "clusterMA.2", "nonmaTlog.UB",
                           "nonma_asymp.LB", "clusterMA.3", "nonma_asymp.UB",
                           "nonma_asympTlog.LB", "clusterMA.4", "nonma_asympTlog.UB",
                           "buck.LB", "clusterMA.5", "buck.UB",
                           "bucklog.LB", "clusterMA.6", "bucklog.UB",
                           "maw2.LB", "clusterMA.7", "maw2.UB",
                           "maw2log.LB", "clusterMA.8", "maw2log.UB",
                           "mata.LB", "clusterMA.9", "mata.UB",
                           "matalog.LB", "clusterMA.10", "matalog.UB",
                           
                           "nonma.time", "nonmaTlog.time", 
                           "nonma_asymp.time", "nonma_asympTlog.time", 
                           "buck.time","bucklog.time",
                           "maw2.time", "maw2log.time",
                           "mata.time", "matalog.time",
                           "IC",
                           "risk","rad","iter","theta", "cent",
                           "select_orig.bic", "select_orig.aic", "simID", "method")
    
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
    
    
    tabn.clusso <- rbind(c(IC="BIC",rad, risk, cent, theta,
                           timeperiod=as.numeric(paste(tim, collapse="")),
                           mod=mod,  type = "QP",time=sim_clusso.time, method="clusso",
                           incluster.bic = as.vector(listpow.bic.qp),
                           outcluster.bic = as.vector(listfp.bic.qp),
                           iter = m,
                           rrest.bic.qp ),
                         c(IC="AIC",rad, risk, cent, theta,
                           timeperiod=as.numeric(paste(tim, collapse="")),
                           mod=mod,type = "QP",time=sim_clusso.time, method="clusso",
                           incluster.aic = as.vector(listpow.aic.qp),
                           outcluster.aic = as.vector(listfp.aic.qp),
                           iter = m,
                           rrest.aic.qp ),
                         c(IC="BIC",rad, risk, cent, theta,
                           timeperiod=as.numeric(paste(tim, collapse="")),
                           mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                           incluster.bic = as.vector(listpow.bic.p),
                           outcluster.bic = as.vector(listfp.bic.p),
                           iter = m,
                           rrest.bic.p ),
                         c(IC="AIC",rad, risk, cent, theta,
                           timeperiod=as.numeric(paste(tim, collapse="")),
                           mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                           incluster.aic = as.vector(listpow.aic.p),
                           outcluster.aic = as.vector(listfp.aic.p),
                           iter = m,
                           rrest.aic.p))
    
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
                               as.vector(sim_stage$RRbic)),
                             c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type="NA", time =  sim_stage.time ,method = "fstagewise" ,
                               incluster.bic = as.vector(outpow.stage.bic),
                               outcluster.bic = as.vector(outfp.stage.bic),
                               iter = m,
                               as.vector(sim_stage$RRaic)))
    
    print("Finished forward stagewise")
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    ##WRITE TO CSV
    #####################################################################################
    #####################################################################################
    #####################################################################################
    out <- rbind(tabn.loc, tabn.pc, tabn.clusso, tabn.stepscan, tabn.fstagewise) 
    out_bounds <- rbind(bounds.loc, bounds.pc)
    if(mod=="space"){
        out <- cbind(stsmodel=rep("space", times=nrow(out)), out)
        out_bounds <- cbind(stsmodel=rep("space", times=nrow(out)), out)
    }else{
        out <- cbind(stsmodel=rep("ST", times=nrow(out)), out)
        out_bounds <- cbind(stsmodel=rep("ST", times=nrow(out)), out)
    }
    write.csv(out, file=paste0("theta",theta,"_singlecluster","_iter",m,"_startsim",startsim,".csv"), row.names = TRUE)    
    write.csv(out_bounds, file=paste0("theta",theta,"_singlecluster_bounds", "_iter",m,"_startsim",startsim,".csv"), row.names = TRUE)
    
}
