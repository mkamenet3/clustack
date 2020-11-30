rm(list=ls())
set.seed(20200921)
library(clusso)
#source("clustack_20191222.R")
source("clustack.R")
##############################
#NULL
##############################

#0) Setup
#dframe1 <- read.csv("clusso-newpenalty/clusso/clusso/data/jap.breast.F.9.10.11.csv")
#dframe2 <- read.csv("clusso-newpenalty/clusso/clusso/data/utmJapan.csv")
#dframe1 <- read.csv("../../../../clusso-newpenalty/clusso/data/jap.breast.F.9.10.11.csv")
#dframe2 <- read.csv("../../../../clusso-newpenalty/clusso/data/utmJapan.csv")
dframe1 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

#dframe.poly2 <- read.csv("clusso-newpenalty/clusso/clusso/data/japan_poly2.csv")
#dframe.poly2 <- read.csv("../../../../clusso-newpenalty/clusso/data/japan_poly2.csv")
dframe.poly2 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
#dframe.prefect2 <- read.csv("clusso-newpenalty/clusso/clusso/data/japan_prefect2.csv")
#dframe.prefect2 <- read.csv("../../../../clusso-newpenalty/clusso/data/japan_prefect2.csv")
dframe.prefect2 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]
japanbreastcancer <- dframe3
japanbreastcancer$period <- japanbreastcancer$year
japanbreastcancer$year <- NULL

cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath

x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]



#set global
rMax <- 20 
Time <- 5
#create set of potential clusters based on distances
potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparsematrix <- spacetimeMat(potentialclusters, numCenters, Time) 

##################################################
######args definitions###########
########Single Cluster:
#args1 = dispersion (Inf (pure poisson), 1000 (moderate overdispersion), 2 (extreme overdispersion))
#args2 = nsims
##################################################

#Arguments passed will only control theta parameter
##R CMD BATCH "--args arg1 arg2" myscript.R -> 
args<- commandArgs(TRUE)

#################################
#TESTER SIM
#################################
#quick function to recode
reval <- function(probs, ix){
    probs[ix] <-1
    return(probs)
}
#nsims
#nsim <- 2#100
cent=1
rad=18
tim = c(1:5)
risk =1
thetas = c(Inf, 1000)
overdispfloor = TRUE
nullmod <- TRUE
# center=1
# radius=18
# tim = c(1:5)
# risk.ratio =1
# thetas = c(Inf, 1000, 2)
# overdispfloor = TRUE
# nullmod <- TRUE

#arguments passed
theta <- as.numeric(args[1]) 
nsim <- as.numeric(args[2])
maxclust <- 15


table.detection.loc.st <- NULL
table.detection.pc.st <- NULL
table.detection.clusso.st <- NULL
table.detection.loc.space <- NULL
table.detection.pc.space <- NULL
table.detection.clusso.space <- NULL
table.bounds.loc.st <- NULL
table.bounds.pc.st <- NULL
table.bounds.loc.space <- NULL
table.bounds.pc.space <- NULL


eps <- 3
path.figures <- "../../../figures/OUTDEC2020/"
path.tables <- "../../../results/OUTDEC2020/"
# path.figures <- "."
# path.tables <- "."


#################################################################################################
#################################################################################################
#SPACETIME
#################################################################################################
#################################################################################################
# Start the clock!
ptm <- proc.time()
for(theta in thetas){
    print(paste0("Params:", "center: ",cent,"; radius: ",
                 rad, "; Timeperiods: ", as.numeric(paste(tim, collapse = "")),
                 "; RR: ", risk, "; Theta :", theta))
    simid <- 1:nsim
    #put the cluster in
    clusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
    n <- length(x)
    init <- clusso::setVectors(japanbreastcancer$period,
                               japanbreastcancer$expdeath, japanbreastcancer$death,
                               covars=NULL, Time=Time)
    E1 <- init$E0
    #Ex <- clusso::scale(init, Time)
    Yx <- init$Y.vec
    vectors <- list(Period = init$Year, Ex=init$E0, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)
    n_uniq <- length(unique(clusters$center))
    numCenters <- n_uniq
    tmp <- clusters[clusters$center==cent,]
    cluster <- tmp[(tmp$r <= rad),]
    rr = matrix(1, nrow=n, ncol=Time)
    rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
    E1 <- as.vector(rr)*init$E0
    message(paste("Running model for periods",tim[1],"through", tail(tim, n=1)))
    #simulate data here
    if (is.infinite(theta)){
        print("Inf theta: pure poisson")
        YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
    } else {
        print("OverD: NB")
        YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
    }
    
    Ex <- clusso::scale_sim(YSIM, init, nsim, Time)
    
    outExp <- lapply(1:nsim, function(i) t(sparsematrix)%*%Ex[[i]])
    outObs <- lapply(1:nsim, function(i) t(sparsematrix)%*%YSIM[[i]])
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY LOC
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by location for each sim (apply)
    #run superclust by location for each sim (apply)
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- lapply(1:nsim, function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
                                                     family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
    }
    sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                    numCenters, Time, maxclust,
                                                                    bylocation = TRUE, model="poisson",
                                                                    overdisp.est = overdisp.est))
    print("finished stacking: by LOC")
    sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                    "risk", risk, "_", "theta", as.character(theta),
                    as.numeric(paste(tim, collapse = "")))
    print(filename <- paste0(sim.i,"_superclustLOC",".RData"))
    #save .RData
    save(sim_superclust_loc, file=filename)
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #NO forceidentify, save
    #####################################################################################
    #####################################################################################
    #####################################################################################
    id.bic <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)))
    id.aic <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)))
    #####################################################################################
    #####################################################################################
    if(any(id.aic!=0)){
        outaic.loc_st <- calcbounds(id.aic, IC="aic", sim_superclust_loc)
    }
    if(any(id.bic!=0)){
        outbic.loc_st <- calcbounds(id.bic, IC="bic", sim_superclust_loc)
    }
    

    #####################################################################################
    #BY BIC
    #####################################################################################
    
    
    
    
    
    if(any(id.bic!=0)){
        #do all diagnostics
        idix.bic <- which(id.bic!=0)
        #prep
        wslarge <- lapply(1:length(idix.bic), function(j) sim_superclust_loc[[idix.bic[[j]]]]$wtMAT[,sim_superclust_loc[[idix.bic[[j]]]]$selection.bic])
        clusterRR_uniqlarge <- lapply(1:length(idix.bic), function(j) sapply(1:nrow(sim_superclust_loc[[idix.bic[[j]]]]$Lambda_dense), 
                                                                 function(k) unique(sim_superclust_loc[[idix.bic[[j]]]]$Lambda_dense[k,]))) 
        clusterRR_ilarge <- lapply(1:length(idix.bic), function(i) rep(NA, 66870))
        clusterRR_uniq_ilarge <- lapply(1:length(idix.bic), function(j) as.matrix(do.call(rbind, clusterRR_uniqlarge[[idix.bic[[j]]]]), ncol=2))
        clusterRR_ilarge <- lapply(1:length(idix.bic), function(j) selectuniqRR(clusterRR_uniq_ilarge[[idix.bic[[j]]]]))
        
        print("BIC")
        #Perform
        outnonma.bic <- nonma(sim_superclust_loc, clusterRR_ilarge, wslarge, idix.bic, IC="bic")
        outnonma_asymp.bic <- nonma_asymp(sim_superclust_loc, clusterRR_ilarge, wslarge, idix.bic, IC="bic")
        print("nonma finished")
        outbuck.theta.bic.time <- system.time(outbuck.theta.bic <- lapply(1:nsim, function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                             thetaa = cluster_thetaa[[i]], 
                                                                                                             w_q=wslarge[[i]], 
                                                                                                             sparsematrix=t(sparsematrix), 
                                                                                                             outExp[[i]],overdisp.est = NULL)))
        outbuck.theta.bic.time <- system.time(outbuck.theta.bic <- lapply(1:nsim, function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]],
                                                                                                             thetaa = cluster_thetaa[[i]], 
                                                                                                             w_q=wslarge[[i]], 
                                                                                                             sparsematrix=t(sparsematrix), 
                                                                                                             outExp[[i]],overdisp.est = NULL)))
        outbuckTlog.theta.bic.time <- system.time(outbuckTlog.theta.bic <- lapply(1:nsim, function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]],
                                                                                                                     thetaa =cluster_thetaa[[i]], 
                                                                                                                     w_q=wslarge[[i]], 
                                                                                                                     sparsematrix=t(sparsematrix),
                                                                                                                     outExp[[i]],
                                                                                                                     overdisp.est = NULL, 
                                                                                                                     transform=TRUE)))
        print("buckland finished")
        outmaw2.theta.bic.time <- system.time(outmaw2.theta.bic <- lapply(1:nsim, function(i) maw2(thetai=clusterRR_ilarge[[i]], 
                                                                                                   thetaa = cluster_thetaa[[i]], 
                                                                                                   w_q=wslarge[[i]],
                                                                                                   sparsematrix=t(sparsematrix ), 
                                                                                                   outExp[[i]],overdisp.est = NULL)))
        
        outmaw2Tlog.theta.bic.time <- system.time(outmaw2Tlog.theta.bic  <- lapply(1:nsim, function(i) maw2(thetai=clusterRR_ilarge[[i]], 
                                                                                                            thetaa = cluster_thetaa[[i]], 
                                                                                                            w_q=wslarge[[i]], 
                                                                                                            sparsematrix=t(sparsematrix), 
                                                                                                            outExp[[i]], 
                                                                                                            overdisp.est = NULL,
                                                                                                            transform=TRUE)))
        print("maw2 finished")
        outmata.theta.bic.time <- system.time(outmata.theta.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                         thetaa = cluster_thetaa[[i]], 
                                                                                                         w_q=wslarge[[i]], 
                                                                                                         sparsematrix=t(sparsematrix ), 
                                                                                                         outExp = outExp[[i]],
                                                                                                         overdisp.est = NULL,
                                                                                                         transform="none")))
        outmataT.theta.bic.time <- system.time(outmataT.theta.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                           thetaa = cluster_thetaa[[i]], 
                                                                                                           w_q=wslarge[[i]], 
                                                                                                           sparsematrix=t(sparsematrix ), 
                                                                                                           outExp = outExp[[i]],
                                                                                                           overdisp.est = NULL, 
                                                                                                           transform="sqrt")))
        outmataTlog.theta.bic.time <- system.time(outmataTlog.theta.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                                 thetaa = cluster_thetaa[[i]], 
                                                                                                                 w_q=wslarge[[i]], 
                                                                                                                 sparsematrix=t(sparsematrix ),
                                                                                                                 outExp = outExp[[i]], 
                                                                                                                 overdisp.est = NULL,
                                                                                                                 transform="log")))
        print("mata finished")
        
    } else {
        print("No clusters identified: BIC")
    }
    
    
    #####################################################################################
    #BY AIC
    #####################################################################################
    if(any(id.aic!=0)){
    
        
    } else {
        print("No clusters identified: AIC")
        print("AIC identifies same as BIC")
    #     nonma.aic <- lapply(1:nsim, function(i) rep(NA,3))
    #     nonma_asymp.aic <- lapply(1:nsim, function(i) rep(NA,3))
    #     outmaw1.aic <- lapply(1:nsim, function(i) rep(NA,3))
    #     outmaw2.aic <- lapply(1:nsim, function(i) rep(NA,3))
    #     outmata.aic <- lapply(1:nsim, function(i) rep(NA,3))
    #     outmataT.aic <- lapply(1:nsim, function(i) rep(NA,3))
    #     outmataTlog.aic <- lapply(1:nsim, function(i) rep(NA,3))
    # }
    }
  
    #####################################################################################        
    outfp.bic <- sum(ifelse(unlist(id.bic)!=0,1,0))/nsim
    outfp.aic <- sum(ifelse(unlist(id.aic)!=0,1,0))/nsim
    #####################################################################################
    
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    ##BIC
    selects <- sapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="bic", flav="loc")
    ##AIC
    selects <- sapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="aic", flav="loc")
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", pow=outpow.bic, fp = outfp.bic),
                      cbind(IC="AIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", pow=outpow.aic, fp = outfp.aic),
                      cbind(IC="AICc",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", pow=outpow.aicc, fp = outfp.aicc))
    table.detection.loc.st <- rbind(table.detection.loc.st, tabn.loc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    bounds.loc <- cbind.data.frame(matrix(unlist(nonma.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.aic), byrow=TRUE, ncol=3))
    bounds.loc$risk <- risk
    bounds.loc$radius <- rad
    bounds.loc$select_orig.bic <- unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic))
    bounds.loc$select_orig.aic <- unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic))
    bounds.loc$simID <- 1:nrow(bounds.loc)
    table.bounds.loc.st <- rbind(table.bounds.loc.st, bounds.loc)

    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY Potential Cluster
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by PC for each sim (apply)
    #run superclust by location for each sim (apply)
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- lapply(1:nsim, function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
                                                     family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
    }
    sim_superclust_pc<- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                  numCenters, Time, maxclust,
                                                                  bylocation = FALSE, model="poisson",
                                                                  overdisp.est = overdisp.est))
    print(filename <- paste0(sim.i,"_superclustPC",".RData"))
    #save .RData
    save(sim_superclust_pc, file=filename)
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #Use forceidentify
    #####################################################################################
    #####################################################################################
    #####################################################################################
    
    #####################################################################################
    #####################################################################################
    #BY BIC
    #####################################################################################
    #####################################################################################
    wslarge <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wtMAT[,sim_superclust_pc[[i]]$selection.bic])
    clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc[[i]]$Lambda_dense), 
                                                             function(k) unique(sim_superclust_pc[[i]]$Lambda_dense[k,]))) 
    
    clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
    clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
    clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
    ##################################################
    #NON-MA VARIANCE
    ##################################################
    clusterRRlarge <- lapply(1:nsim, 
                             function(i) unique(sim_superclust_pc[[i]]$Lambda_dense[sim_superclust_pc[[i]]$maxpcs[sim_superclust_pc[[i]]$selection.bic],])[2])
    se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(clusterRRlarge[[i]]/(Time*n)))
    nonma.bic <- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                  clusterMA = clusterRRlarge[[i]],
                                                  ub=clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]]))
    cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
    
    #asymptotic
    se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/(sum(YSIM[[i]][ix]))))
    nonma_asymp.bic <- lapply(1:nsim, function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                        clusterMA = clusterRRlarge[[i]],
                                                        ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
    
    ##################################################
    #Buckland 1997
    ##################################################
    outbuck.bic <- lapply(1:nsim, 
                          function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                     w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
    ##################################################
    outmaw1.bic <- lapply(1:nsim, 
                          function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                           w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #MAW2 (B&A pg. 345)
    ##################################################
    outmaw2.bic <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                   w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #Turek-Fletcher MATA Bounds (for non-normal data)
    ##################################################
    outmata.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="none"))
    ##################################################
    #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
    ##################################################
    outmataT.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                          w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="sqrt"))
    ##################################################
    #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
    ##################################################
    outmataTlog.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="log"))
    
    #####################################################################################
    #####################################################################################
    #BY AIC
    #####################################################################################
    #####################################################################################
    if (any(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)) != unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic)))){
        print("AIC identifies different from BIC")
        wslarge <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wtMAT[,sim_superclust_pc[[i]]$selection.aic])
        clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc[[i]]$Lambda_dense), 
                                                                 function(k) unique(sim_superclust_pc[[i]]$Lambda_dense[k,]))) 
        
        clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
        clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
        clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
        ##################################################
        #NON-MA VARIANCE
        ##################################################
        clusterRRlarge <- lapply(1:nsim, 
                                 function(i) unique(sim_superclust_pc[[i]]$Lambda_dense[sim_superclust_pc[[i]]$maxpcs[sim_superclust_pc[[i]]$selection.aic],])[2])
        se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(clusterRRlarge[[i]]/(Time*n)))
        nonma.aic <- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                      clusterMA = clusterRRlarge[[i]],
                                                      ub=clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]]))
        cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
        
        #asymptotic
        se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/(sum(YSIM[[i]][ix]))))
        nonma_asymp.aic <- lapply(1:nsim, function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                            clusterMA = clusterRRlarge[[i]],
                                                            ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
        
        ##################################################
        #Buckland 1997
        ##################################################
        outbuck.aic <- lapply(1:nsim, 
                              function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
        ##################################################
        outmaw1.aic <- lapply(1:nsim, 
                              function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                               w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #MAW2 (B&A pg. 345)
        ##################################################
        outmaw2.aic <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                       w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #Turek-Fletcher MATA Bounds (for non-normal data)
        ##################################################
        outmata.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="none"))
        ##################################################
        #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
        ##################################################
        outmataT.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                              w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="sqrt"))
        ##################################################
        #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
        ##################################################
        outmataTlog.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                 w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="log"))
    } else {
        print("AIC identifies same as BIC")
        nonma.aic <- lapply(1:nsim, function(i) rep(NA,3))
        nonma_asymp.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmaw1.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmaw2.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmata.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmataT.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmataTlog.aic <- lapply(1:nsim, function(i) rep(NA,3))
    }
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #POWER AND FALSE POSITIVE RATE
    #####################################################################################
    #####################################################################################
    #####################################################################################
    ##################################
    #################################
    #DIAGNOSTICS: #calc power and FB rate
    #################################
    # #Which PCs overlap true cluster?
    #what was identified in each sim by IC
    ident.bic <- lapply(1:nsim, function(i)
        round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],eps))
    ident.aic <- lapply(1:nsim, function(i)
        round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],eps))
    ident.aicc <- lapply(1:nsim, function(i)
        round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aicc,],eps))
    #1) Did it find anything INSIDE the cluster?
    incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
    incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
    incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
    #calc power
    pow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
    pow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
    pow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
    outpow.bic <- paste0(pow.bic*100, "%")
    outpow.aic <- paste0(pow.aic*100, "%")
    outpow.aicc <- paste0(pow.aicc*100, "%")
    #2) Did it find anything OUTSIDE the cluster?
    #rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
    #this should be everything that doesn't touch the cluster
    outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
    outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
    outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
    #calc FP rate
    fp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
    fp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
    fp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim
    outfp.bic <- paste0(fp.bic*100, "%")
    outfp.aic <- paste0(fp.aic*100, "%")
    outfp.aicc <- paste0(fp.aicc*100, "%")
    # ##################################
    # #plot probability maps
    # ##################################
    # #create empties
    # vec <- rep(0, 208 * Time)
    # position.bic <- list(vec)[rep(1, nsim)]
    # position.aic <- list(vec)[rep(1, nsim)]
    # position.aicc <- list(vec)[rep(1, nsim)]
    # #recode identified cells as 1's, all other zeros
    # ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
    # ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1) ==1))
    # ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
    # #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
    # simindicator.bic <- mapply(reval, position.bic, ix.bic)
    # simindicator.aic <- mapply(reval, position.aic, ix.aic)
    # simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
    # #find probability of detection for each location in time
    # probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
    # probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
    # probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
    # #map probability detections by IC to grey scale
    # colprob <- colormapping(list(probs.bic,
    #                              probs.aic,
    #                              probs.aicc,
    #                              as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
    
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    ##BIC
    selects <- sapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="bic", flav="pc")
    ##AIC
    selects <- sapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="aic", flav="pc")

    ##################################
    #Add sim results to table
    ##################################
    tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", pow=outpow.bic, fp = outfp.bic),
                     cbind(IC="AIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", pow=outpow.aic, fp = outfp.aic),
                     cbind(IC="AICc",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", pow=outpow.aicc, fp = outfp.aicc))
    table.detection.pc.st <- rbind(table.detection.pc.st, tabn.pc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    bounds.pc <- cbind.data.frame(matrix(unlist(nonma.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.aic), byrow=TRUE, ncol=3))
    bounds.pc$risk <- risk
    bounds.pc$radius <- rad
    bounds.pc$select_orig.bic <- unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic))
    bounds.pc$select_orig.aic <- unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic))
    bounds.pc$simID <- 1:nrow(bounds.pc)
    
    table.bounds.pc.st <- rbind(table.bounds.pc.st, bounds.pc)

    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY CLUSSO
    #####################################################################################
    #####################################################################################
    #####################################################################################
    YSIM1 <- lapply(1:nsim, function(i) as.vector(matrix(YSIM[[i]], ncol=Time, byrow = FALSE)))
    E01 <- as.vector(matrix(init$E0, ncol=Time, byrow=FALSE))
    truth1 <- as.vector(matrix(init$Y.vec, ncol=Time, byrow = FALSE))
    period1 <- as.vector(matrix(init$Year, ncol=Time, byrow = FALSE))
    id <- rep(1:208, times = 5)
    #create list of dataframes
    jbcSIM <- lapply(1:nsim, function(i) cbind.data.frame(expected = E01,
                                                          observed = YSIM1[[i]],
                                                          period = period1,
                                                          id = id))
    jbcSIM <- lapply(1:nsim, function(i) jbcSIM[[i]][order(jbcSIM[[i]]$id, jbcSIM[[i]]$period),])
    #create list of dataframes

    #run clusso
    sim_clusso <- lapply(1:nsim, function(i) clusso::clusso(df=jbcSIM[[i]],
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
                                                            maxclust = maxclust))
    print(filename <- paste0(sim.i,"_clusso",".RData"))
    #save .RData
    save(sim_clusso, file=filename)
    #system(paste0("gzip ", filename))
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################

    vec <- rep(0, 208*Time)
    position <- list(vec)[rep(1, nsim)]
    #background rates
    ##Quasi-P
    bgRate_i.bic.qp <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.st$E.qbic,ncol=Time)[,j]))))))
    bgRate_i.aic.qp <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.st$E.qaic,ncol=Time)[,j]))))))
    bgRate_i.aicc.qp <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.st$E.qaicc,ncol=Time)[,j]))))))

    bgRate.bic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.bic.qp[[i]], each = 208))
    bgRate.aic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aic.qp[[i]], each = 208))
    bgRate.aicc.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aicc.qp[[i]], each = 208))

    #detect
    ix.bic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qbic) - bgRate.bic.qp[[i]])>=10^-3))
    ix.aic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qaic) - bgRate.aic.qp[[i]])>=10^-3))
    ix.aicc.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qaicc) - bgRate.aicc.qp[[i]])>=10^-3))

    simindicator.bic.qp <- mapply(reval, position, ix.bic.qp)
    simindicator.aic.qp <- mapply(reval, position, ix.aic.qp)
    simindicator.aicc.qp <- mapply(reval, position, ix.aicc.qp)

    probs.bic.qp <- Matrix::rowSums(simindicator.bic.qp)/nsim
    probs.aic.qp <- Matrix::rowSums(simindicator.aic.qp)/nsim
    probs.aicc.qp <- Matrix::rowSums(simindicator.aicc.qp)/nsim

    ##Poisson
    bgRate_i.bic.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qbic,ncol=Time)[,j]))))))
    bgRate_i.aic.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qaic,ncol=Time)[,j]))))))
    bgRate_i.aicc.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qaicc,ncol=Time)[,j]))))))

    bgRate.bic.p <- lapply(1:nsim, function(i) rep(bgRate_i.bic.p[[i]], each = 208))
    bgRate.aic.p <- lapply(1:nsim, function(i) rep(bgRate_i.aic.p[[i]], each = 208))
    bgRate.aicc.p <- lapply(1:nsim, function(i) rep(bgRate_i.aicc.p[[i]], each = 208))

    #detect
    ix.bic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qbic) - bgRate.bic.p[[i]])>=10^-3))
    ix.aic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qaic) - bgRate.aic.p[[i]])>=10^-3))
    ix.aicc.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qaicc) - bgRate.aicc.p[[i]])>=10^-3))

    simindicator.bic.p <- mapply(reval, position, ix.bic.p)
    simindicator.aic.p <- mapply(reval, position, ix.aic.p)
    simindicator.aicc.p <- mapply(reval, position, ix.aicc.p)

    probs.bic.p <- Matrix::rowSums(simindicator.bic.p)/nsim
    probs.aic.p <- Matrix::rowSums(simindicator.aic.p)/nsim
    probs.aicc.p <- Matrix::rowSums(simindicator.aicc.p)/nsim

    # #map probability detections by IC to grey scale
    # colprob.qp <- colormapping(list(probs.bic.qp,
    #                                 probs.aic.qp,
    #                                 probs.aicc.qp,
    #                                 as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob.qp,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_qp_clusso.pdf"))
    # 
    # colprob.p <- colormapping(list(probs.bic.p,
    #                                probs.aic.p,
    #                                probs.aicc.p,
    #                                as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob.p,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_p_clusso.pdf"))

    ##################################
    #POWER/FP Rate
    ##################################
    ##Power
    ###QP
    listpow.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                           sim_clusso[[i]]$lassoresult.qp.st$selections$select.qbic,rr,
                                                                           risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                            sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaic,rr,
                                                                            risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aicc.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                             sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaicc,rr,
                                                                             risk,nsim,Time, numCenters, pow=TRUE))
    outpow.bic.qp <- paste0(sum(unlist(listpow.bic.qp))/nsim*100, "%")
    outpow.aic.qp <- paste0(sum(unlist(listpow.aic.qp))/nsim*100, "%")
    outpow.aicc.qp <- paste0(sum(unlist(listpow.aicc.qp))/nsim*100, "%")
    ###Poisson
    listpow.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                          sim_clusso[[i]]$lassoresult.p.st$selections$select.qbic,rr,
                                                                          risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                           sim_clusso[[i]]$lassoresult.p.st$selections$select.qaic,rr,
                                                                           risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aicc.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                            sim_clusso[[i]]$lassoresult.p.st$selections$select.qaicc,rr,
                                                                            risk,nsim,Time, numCenters, pow=TRUE))
    outpow.bic.p <- paste0(sum(unlist(listpow.bic.p))/nsim*100, "%")
    outpow.aic.p <- paste0(sum(unlist(listpow.aic.p))/nsim*100, "%")
    outpow.aicc.p <- paste0(sum(unlist(listpow.aicc.p))/nsim*100, "%")

    ##FP
    listfp.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                          sim_clusso[[i]]$lassoresult.qp.st$selections$select.qbic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                           sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaic,rr,
                                                                           risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aicc.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                            sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaicc,rr,
                                                                            risk,nsim,Time, numCenters, pow=FALSE))
    outfp.bic.qp <- paste0(sum(unlist(listfp.bic.qp))/nsim*100, "%")
    outfp.aic.qp <- paste0(sum(unlist(listfp.aic.qp))/nsim*100, "%")
    outfp.aicc.qp <- paste0(sum(unlist(listfp.aicc.qp))/nsim*100, "%")
    ###Poisson
    listfp.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                         sim_clusso[[i]]$lassoresult.p.st$selections$select.qbic,rr,
                                                                         risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                          sim_clusso[[i]]$lassoresult.p.st$selections$select.qaic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aicc.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                           sim_clusso[[i]]$lassoresult.p.st$selections$select.qaicc,rr,
                                                                           risk,nsim,Time, numCenters, pow=FALSE))
    outfp.bic.p <- paste0(sum(unlist(listfp.bic.p))/nsim*100, "%")
    outfp.aic.p <- paste0(sum(unlist(listfp.aic.p))/nsim*100, "%")
    outfp.aicc.p <- paste0(sum(unlist(listfp.aicc.p))/nsim*100, "%")

    ##################################
    #RR maps
    ##################################
    ##################################
    #Add sim results to table
    ##################################
    tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.bic.qp, fp = outfp.bic.qp, type = "QP"),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aic.qp, fp = outfp.aic.qp, type = "QP"),
                         cbind(IC="AICc",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aicc.qp, fp = outfp.aicc.qp, type = "QP"),
                         cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.bic.p, fp = outfp.bic.p, type = "Pois"),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aic.p, fp = outfp.aic.p, type = "Pois"),
                         cbind(IC="AICc",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aicc.p, fp = outfp.aicc.p, type = "Pois"))
    table.detection.clusso.st <- rbind(table.detection.clusso.st, tabn.clusso)
}

# Stop the clock
proc.time() - ptm
#####################################################################################
#####################################################################################
#####################################################################################
##WRITE TO CSV
#####################################################################################
#####################################################################################
#####################################################################################
#superclust by loc
print(table.detection.loc.st)
write.csv(table.detection.loc.st, file=paste0(path.tables,"null_singlecluster_loc_ST.csv"), row.names=TRUE)
#bound by loc
print(table.bounds.loc.st)
write.csv(table.bounds.loc.st, file=paste0(path.tables,"null_singlecluster_loc_ST_bounds.csv"), row.names=TRUE)

#superclust by pc
print(table.detection.pc.st)
write.csv(table.detection.pc.st, file=paste0(path.tables,"null_singlecluster_pc_ST.csv"), row.names=TRUE)
#bounds by pc
print(table.bounds.pc.st)
write.csv(table.bounds.pc.st, file=paste0(path.tables, "null_singlecluster_pc_ST_bounds.csv"), row.names = TRUE)


#clusso
print(table.detection.clusso.st)
write.csv(table.detection.clusso.st, file=paste0(path.tables,"null_singlecluster_clusso_ST.csv"), row.names=TRUE)




#################################################################################################
#################################################################################################
#SPACE-ONLY
#################################################################################################
#################################################################################################
# Start the clock!
ptm <- proc.time()
for(theta in thetas){
    print(paste0("Params:", "center: ",cent,"; radius: ",
                 rad, "; Timeperiods: ", "spaceonly", 
                 "; RR: ", risk, "; Theta :", theta))
    simid <- 1:nsim
    #put the cluster in
    clusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
    n <- length(x)
    init <- clusso::setVectors(japanbreastcancer$period, 
                               japanbreastcancer$expdeath, japanbreastcancer$death,
                               covars=NULL, Time=Time)
    E1 <- init$E0
    #Ex <- clusso::scale(init, Time)
    Yx <- init$Y.vec
    vectors <- list(Period = init$Year, Ex=init$E0, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)    
    n_uniq <- length(unique(clusters$center))
    numCenters <- n_uniq
    tmp <- clusters[clusters$center==cent,]
    cluster <- tmp[(tmp$r <= rad),]
    
    rr = matrix(1, nrow=n, ncol=Time)
    #rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
    allTime <- 1:Time
    rr[cluster$last, allTime[1]:tail(allTime, n=1)] <- risk
    E1 <- as.vector(rr)*init$E0
    message(paste("Running model for spaceonly"))
    #simulate data here
    if (is.infinite(theta)){
        print("Inf theta: pure poisson")    
        YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
    } else {
        print("OverD: NB")
        YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))    
    }
    Ex <- clusso::scale_sim(YSIM, init, nsim, Time)
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY LOCATION
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by location for each sim (apply)
    #run superclust by location for each sim (apply)
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- lapply(1:nsim, function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
                                                     family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
    }
    sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                    numCenters, Time, maxclust,
                                                                    bylocation = TRUE, model="poisson",
                                                                    overdisp.est = overdisp.est))
    sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                    "risk", risk, "_", "theta", as.character(theta),"_spaceonly")
    print(filename <- paste0(sim.i,"_superclustLOC",".RData"))
    #save .RData
    save(sim_superclust_loc, file=filename)
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #Use forceidentify
    #####################################################################################
    #####################################################################################
    #####################################################################################
    
    #####################################################################################
    #####################################################################################
    #BY BIC
    #####################################################################################
    #####################################################################################
    wslarge <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wtMAT[,sim_superclust_loc[[i]]$selection.bic])
    clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_loc[[i]]$Lambda_dense), 
                                                             function(k) unique(sim_superclust_loc[[i]]$Lambda_dense[k,]))) 
    
    clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
    clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
    clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
    ##################################################
    #NON-MA VARIANCE
    #TODO 2020-09-22
    ##################################################
    ######################
    # #MA by maxloc
    # ix_locs<- rep(0,1040);ix_locs[ix] <-1
    # ix_locs_mat <- matrix(aa, nrow = 1)%*%sparsematrix 
    # ix_locs_mods <- which(ix_locs_mat!=0)
    cluster_thetaa_locs <- lapply(1:nsim, function(i) 1)#lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]][ix_locs_mods]*wslarge[[i]][ix_locs_mods]))
    cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
    ######################
    clusterRRlarge <- lapply(1:nsim, 
                             function(i) unique(sim_superclust_loc[[i]]$Lambda_dense[sim_superclust_loc[[i]]$maxpcs[sim_superclust_loc[[i]]$selection.bic],])[2])
    se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(cluster_thetaa_locs[[i]]/(Time*n)))
    
    nonma.bic <- lapply(1:nsim, function(i) cbind(lb=cluster_thetaa_locs[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                  clusterMA = cluster_thetaa_locs[[i]],
                                                  ub=cluster_thetaa_locs[[i]]+1.96*se_clusterRRlarge[[i]]))
    #asymptotic
    se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(cluster_thetaa_locs[[i]]/(sum(YSIM[[i]]))))
    nonma_asymp.bic <- lapply(1:nsim, function(i) cbind(lbasymp=cluster_thetaa_locs[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                        clusterMA = cluster_thetaa_locs[[i]],
                                                        ubasymp=cluster_thetaa_locs[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
    
    ##################################################
    #Buckland 1997
    ##################################################
    outbuck.bic <- lapply(1:nsim, 
                          function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                     w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
    ##################################################
    outmaw1.bic <- lapply(1:nsim, 
                          function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                           w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #MAW2 (B&A pg. 345)
    ##################################################
    outmaw2.bic <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                   w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #Turek-Fletcher MATA Bounds (for non-normal data)
    ##################################################
    outmata.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="none"))
    ##################################################
    #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
    ##################################################
    outmataT.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                          w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="sqrt"))
    ##################################################
    #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
    ##################################################
    outmataTlog.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="log"))
    
    #####################################################################################
    #####################################################################################
    #BY AIC
    #####################################################################################
    #####################################################################################
    if (any(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)) != unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)))){
        print("AIC identifies different from BIC")
        wslarge <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wtMAT[,sim_superclust_loc[[i]]$selection.aic])
        clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_loc[[i]]$Lambda_dense), 
                                                                 function(k) unique(sim_superclust_loc[[i]]$Lambda_dense[k,]))) 
        
        clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
        clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
        clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
        ##################################################
        #NON-MA VARIANCE
        #TODO 2020-09-22
        ##################################################
        ######################
        #MA by maxloc
        #ix_locs<- rep(0,1040);ix_locs[ix] <-1
        #ix_locs_mat <- matrix(aa, nrow = 1)%*%sparsematrix 
        #ix_locs_mods <- which(ix_locs_mat!=0)
        cluster_thetaa_locs <- lapply(1:nsim, function(i) 1)#lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]][ix_locs_mods]*wslarge[[i]][ix_locs_mods]))
        cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
        ######################
        clusterRRlarge <- lapply(1:nsim, 
                                 function(i) unique(sim_superclust_loc[[i]]$Lambda_dense[sim_superclust_loc[[i]]$maxpcs[sim_superclust_loc[[i]]$selection.aic],])[2])
        se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(cluster_thetaa_locs[[i]]/(Time*n)))
        
        nonma.aic <- lapply(1:nsim, function(i) cbind(lb=cluster_thetaa_locs[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                      clusterMA = cluster_thetaa_locs[[i]],
                                                      ub=cluster_thetaa_locs[[i]]+1.96*se_clusterRRlarge[[i]]))
        #asymptotic
        se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(cluster_thetaa_locs[[i]]/(sum(YSIM[[i]]))))
        nonma_asymp.aic <- lapply(1:nsim, function(i) cbind(lbasymp=cluster_thetaa_locs[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                            clusterMA = cluster_thetaa_locs[[i]],
                                                            ubasymp=cluster_thetaa_locs[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
        
        ##################################################
        #Buckland 1997
        ##################################################
        outbuck.aic <- lapply(1:nsim, 
                              function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
        ##################################################
        outmaw1.aic <- lapply(1:nsim, 
                              function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                               w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #MAW2 (B&A pg. 345)
        ##################################################
        outmaw2.aic <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                       w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #Turek-Fletcher MATA Bounds (for non-normal data)
        ##################################################
        outmata.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="none"))
        ##################################################
        #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
        ##################################################
        outmataT.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                              w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="sqrt"))
        ##################################################
        #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
        ##################################################
        outmataTlog.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                 w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="log"))
    } else {
        print("AIC identifies same as BIC")
        nonma.aic <- lapply(1:nsim, function(i) rep(NA,3))
        nonma_asymp.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmaw1.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmaw2.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmata.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmataT.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmataTlog.aic <- lapply(1:nsim, function(i) rep(NA,3))
    }
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #POWER AND FALSE POSITIVE RATE
    #####################################################################################
    #####################################################################################
    #####################################################################################
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    #Which PCs overlap true cluster?
    rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
    clusteroverlap<- rrbin_cluster%*%sparsematrix
    rrbin_inside <- ifelse(sparsematrix%*%t(clusteroverlap)!=0,1,0)
    #what was identified in each sim by IC
    ident.bic <- lapply(1:nsim, function(i) 
        round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,],eps))
    ident.aic <- lapply(1:nsim, function(i)
        round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,],eps))
    ident.aicc <- lapply(1:nsim, function(i) 
        round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aicc,],eps))
    #1) Did it find anything INSIDE the cluster?
    incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
    incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
    incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
    #calc power
    pow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
    pow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
    pow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
    outpow.bic <- paste0(pow.bic*100, "%")
    outpow.aic <- paste0(pow.aic*100, "%")
    outpow.aicc <- paste0(pow.aicc*100, "%")
    #2) Did it find anything OUTSIDE the cluster?
    rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0) 
    #this should be everything that doesn't touch the cluster
    outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
    outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
    outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
    #calc FP rate
    fp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
    fp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
    fp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim
    outfp.bic <- paste0(fp.bic*100, "%")
    outfp.aic <- paste0(fp.aic*100, "%")
    outfp.aicc <- paste0(fp.aicc*100, "%")
    # ##################################
    # #plot probability maps
    # ##################################
    # #create empties
    # vec <- rep(0, 208 * Time)
    # position.bic <- list(vec)[rep(1, nsim)]
    # position.aic <- list(vec)[rep(1, nsim)]
    # position.aicc <- list(vec)[rep(1, nsim)]
    # #recode identified cells as 1's, all other zeros
    # ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
    # ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1)==1))
    # ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
    # #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
    # simindicator.bic <- mapply(reval, position.bic, ix.bic)
    # simindicator.aic <- mapply(reval, position.aic, ix.aic)
    # simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
    # #find probability of detection for each location in time
    # probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
    # probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
    # probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
    # #map probability detections by IC to grey scale
    # colprob <- colormapping(list(probs.bic,
    #                              probs.aic,
    #                              probs.aicc,
    #                              as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_loc.pdf"))
    
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    ##BIC
    selects <- sapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="bic", flav="loc")
    ##AIC
    selects <- sapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="aic", flav="loc")

    ##################################
    #Add sim results to table
    ##################################
    tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", pow=outpow.bic, fp = outfp.bic),
                      cbind(IC="AIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", pow=outpow.aic, fp = outfp.aic),
                      cbind(IC="AICc",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", pow=outpow.aicc, fp = outfp.aicc))
    table.detection.loc.space <- rbind(table.detection.loc.space, tabn.loc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    bounds.loc <- cbind.data.frame(matrix(unlist(nonma.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.aic), byrow=TRUE, ncol=3))
    bounds.loc$risk <- risk
    bounds.loc$radius <- rad
    bounds.loc$select_orig.bic <- unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic))
    bounds.loc$select_orig.aic <- unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic_orig))
    bounds.loc$simID <- 1:nrow(bounds.loc)
    
    table.bounds.loc.space <- rbind(table.bounds.loc.space, bounds.loc)
    
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
        offset_reg <- lapply(1:nsim, function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
                                                     family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
    }
    sim_superclust_pc<- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                  numCenters, Time, maxclust,
                                                                  bylocation = FALSE, model="poisson",
                                                                  overdisp.est = overdisp.est))
    print(filename <- paste0(sim.i,"_superclustPC",".RData"))
    #save .RData
    save(sim_superclust_pc, file=filename)
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #Use forceidentify
    #####################################################################################
    #####################################################################################
    #####################################################################################
    
    #####################################################################################
    #####################################################################################
    #BY BIC
    #####################################################################################
    #####################################################################################
    wslarge <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wtMAT[,sim_superclust_pc[[i]]$selection.bic])
    clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc[[i]]$Lambda_dense), 
                                                             function(k) unique(sim_superclust_pc[[i]]$Lambda_dense[k,]))) 
    
    clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
    clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
    clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
    ##################################################
    #NON-MA VARIANCE
    ##################################################
    clusterRRlarge <- lapply(1:nsim, 
                             function(i) unique(sim_superclust_pc[[i]]$Lambda_dense[sim_superclust_pc[[i]]$maxpcs[sim_superclust_pc[[i]]$selection.bic],])[2])
    se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(clusterRRlarge[[i]]/(Time*n)))
    cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
    
    nonma.bic <- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                  clusterMA = clusterRRlarge[[i]],
                                                  ub=clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]]))
    #asymptotic
    se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/(sum(YSIM[[i]][ix]))))
    nonma_asymp.bic <- lapply(1:nsim, function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                        clusterMA = clusterRRlarge[[i]],
                                                        ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
    
    ##################################################
    #Buckland 1997
    ##################################################
    outbuck.bic <- lapply(1:nsim, 
                          function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                     w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
    ##################################################
    outmaw1.bic <- lapply(1:nsim, 
                          function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                           w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #MAW2 (B&A pg. 345)
    ##################################################
    outmaw2.bic <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                   w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
    ##################################################
    #Turek-Fletcher MATA Bounds (for non-normal data)
    ##################################################
    outmata.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="none"))
    ##################################################
    #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
    ##################################################
    outmataT.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                          w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="sqrt"))
    ##################################################
    #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
    ##################################################
    outmataTlog.bic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="log"))
    
    #####################################################################################
    #####################################################################################
    #BY AIC
    #####################################################################################
    #####################################################################################
    if (any(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)) != unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic)))){
        print("AIC identifies different from BIC")
        wslarge <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wtMAT[,sim_superclust_pc[[i]]$selection.aic])
        clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc[[i]]$Lambda_dense), 
                                                                 function(k) unique(sim_superclust_pc[[i]]$Lambda_dense[k,]))) 
        
        clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
        clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
        clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
        ##################################################
        #NON-MA VARIANCE
        ##################################################
        clusterRRlarge <- lapply(1:nsim, 
                                 function(i) unique(sim_superclust_pc[[i]]$Lambda_dense[sim_superclust_pc[[i]]$maxpcs[sim_superclust_pc[[i]]$selection.aic],])[2])
        se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(clusterRRlarge[[i]]/(Time*n)))
        cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
        
        nonma.aic <- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                      clusterMA = clusterRRlarge[[i]],
                                                      ub=clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]]))
        #asymptotic
        se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/(sum(YSIM[[i]][ix]))))
        nonma_asymp.aic <- lapply(1:nsim, function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                            clusterMA = clusterRRlarge[[i]],
                                                            ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
        
        ##################################################
        #Buckland 1997
        ##################################################
        outbuck.aic <- lapply(1:nsim, 
                              function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
        ##################################################
        outmaw1.aic <- lapply(1:nsim, 
                              function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                               w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #MAW2 (B&A pg. 345)
        ##################################################
        outmaw2.aic <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                       w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL))
        ##################################################
        #Turek-Fletcher MATA Bounds (for non-normal data)
        ##################################################
        outmata.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="none"))
        ##################################################
        #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
        ##################################################
        outmataT.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                              w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="sqrt"))
        ##################################################
        #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
        ##################################################
        outmataTlog.aic <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                 w_q=wslarge[[i]], sparsematrix=t(sparsematrix ), overdisp.est = NULL, transform="log"))
    } else {
        print("AIC identifies same as BIC")
        nonma.aic <- lapply(1:nsim, function(i) rep(NA,3))
        nonma_asymp.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmaw1.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmaw2.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmata.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmataT.aic <- lapply(1:nsim, function(i) rep(NA,3))
        outmataTlog.aic <- lapply(1:nsim, function(i) rep(NA,3))
    }
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #POWER AND FALSE POSITIVE RATE
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #################################
    #DIAGNOSTICS: #calc power and FB rate
    #################################
    # #Which PCs overlap true cluster?
    #what was identified in each sim by IC
    ident.bic <- lapply(1:nsim, function(i) 
        round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],eps))
    ident.aic <- lapply(1:nsim, function(i) 
        round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],eps))
    ident.aicc <- lapply(1:nsim, function(i)
        round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aicc,],eps))
    #1) Did it find anything INSIDE the cluster?
    incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
    incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
    incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
    #calc power
    pow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
    pow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
    pow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
    outpow.bic <- paste0(pow.bic*100, "%")
    outpow.aic <- paste0(pow.aic*100, "%")
    outpow.aicc <- paste0(pow.aicc*100, "%")
    #2) Did it find anything OUTSIDE the cluster?
    #rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
    #this should be everything that doesn't touch the cluster
    outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
    outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
    outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
    #calc FP rate
    fp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
    fp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
    fp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim
    outfp.bic <- paste0(fp.bic*100, "%")
    outfp.aic <- paste0(fp.aic*100, "%")
    outfp.aicc <- paste0(fp.aicc*100, "%")
    # ##################################
    # #plot probability maps
    # ##################################
    # #create empties
    # vec <- rep(0, 208 * Time)
    # position.bic <- list(vec)[rep(1, nsim)]
    # position.aic <- list(vec)[rep(1, nsim)]
    # position.aicc <- list(vec)[rep(1, nsim)]
    # #recode identified cells as 1's, all other zeros
    # ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
    # ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1) ==1))
    # ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
    # #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
    # simindicator.bic <- mapply(reval, position.bic, ix.bic)
    # simindicator.aic <- mapply(reval, position.aic, ix.aic)
    # simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
    # #find probability of detection for each location in time
    # probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
    # probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
    # probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
    # #map probability detections by IC to grey scale
    # colprob <- colormapping(list(probs.bic,
    #                              probs.aic,
    #                              probs.aicc,
    #                              as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    ##BIC
    selects <- sapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="bic", flav="pc")
    ##AIC
    selects <- sapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic)
    meanrr <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wLambda[selects[i],])
    names(meanrr) <- paste0("l",1:nsim)
    meanrr[which(selects==0)] <- NULL
    meanrr.df <- do.call(rbind, meanrr)
    meanrrs <- apply(meanrr.df,2, mean)
    ric <- matrix(meanrrs, ncol = Time)
    plotmeanrr_stack(ric, Time, sim.i,ic="aic", flav="pc")

    # ##################################
    #Add sim results to table
    ##################################
    tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", pow=outpow.bic, fp = outfp.bic),
                     cbind(IC="AIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", pow=outpow.aic, fp = outfp.aic),
                     cbind(IC="AICc",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", pow=outpow.aicc, fp = outfp.aicc))
    table.detection.pc.space <- rbind(table.detection.pc.space, tabn.pc)
    ##################################
    #Add clustackbounds to table
    ##################################
    bounds.pc <- cbind.data.frame(matrix(unlist(nonma.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.bic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(nonma_asymp.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw1.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmaw2.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmata.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataT.aic), byrow=TRUE, ncol=3),
                                   matrix(unlist(outmataTlog.aic), byrow=TRUE, ncol=3))
    bounds.pc$risk <- risk
    bounds.pc$radius <- rad
    bounds.pc$select_orig.bic <- unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)
    bounds.pc$select_orig.aic <- unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic))
    bounds.pc$simID <- 1:nrow(bounds.loc)
    table.bounds.pc.space <- rbind(table.bounds.pc.space, bounds.pc)
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY CLUSSO
    #####################################################################################
    #####################################################################################
    #####################################################################################
    YSIM1 <- lapply(1:nsim, function(i) as.vector(matrix(YSIM[[i]], ncol=Time, byrow = FALSE)))
    E01 <- as.vector(matrix(init$E0, ncol=Time, byrow=FALSE))
    truth1 <- as.vector(matrix(init$Y.vec, ncol=Time, byrow = FALSE))
    period1 <- as.vector(matrix(init$Year, ncol=Time, byrow = FALSE))
    id <- rep(1:208, times = 5)
    #create list of dataframes
    jbcSIM <- lapply(1:nsim, function(i) cbind.data.frame(expected = E01,
                                                          observed = YSIM1[[i]],
                                                          period = period1,
                                                          id = id))
    jbcSIM <- lapply(1:nsim, function(i) jbcSIM[[i]][order(jbcSIM[[i]]$id, jbcSIM[[i]]$period),])
    #create list of dataframes
    
    #run clusso
    sim_clusso <- lapply(1:nsim, function(i) clusso::clusso(df=jbcSIM[[i]],
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
                                                            maxclust = maxclust))
    print(filename <- paste0(sim.i,"_clusso",".RData"))
    #save .RData
    save(sim_clusso, file=filename)
    #system(paste0("gzip ", filename))
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    
    vec <- rep(0, 208*Time)
    position <- list(vec)[rep(1, nsim)]
    #background rates
    ##Quasi-P
    bgRate_i.bic.qp <- lapply(1:nsim, function(i) 
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.s$E.qbic,ncol=Time)[,j]))))))
    bgRate_i.aic.qp <- lapply(1:nsim, function(i) 
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.s$E.qaic,ncol=Time)[,j]))))))
    bgRate_i.aicc.qp <- lapply(1:nsim, function(i) 
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.s$E.qaicc,ncol=Time)[,j]))))))
    
    bgRate.bic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.bic.qp[[i]], each = 208))
    bgRate.aic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aic.qp[[i]], each = 208))
    bgRate.aicc.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aicc.qp[[i]], each = 208))
    
    #detect 
    ix.bic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.s$E.qbic) - bgRate.bic.qp[[i]])>=10^-3))
    ix.aic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.s$E.qaic) - bgRate.aic.qp[[i]])>=10^-3))
    ix.aicc.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.s$E.qaicc) - bgRate.aicc.qp[[i]])>=10^-3))
    
    simindicator.bic.qp <- mapply(reval, position, ix.bic.qp)
    simindicator.aic.qp <- mapply(reval, position, ix.aic.qp)
    simindicator.aicc.qp <- mapply(reval, position, ix.aicc.qp)
    
    probs.bic.qp <- Matrix::rowSums(simindicator.bic.qp)/nsim
    probs.aic.qp <- Matrix::rowSums(simindicator.aic.qp)/nsim
    probs.aicc.qp <- Matrix::rowSums(simindicator.aicc.qp)/nsim
    
    ##Poisson
    bgRate_i.bic.p <- lapply(1:nsim, function(i) 
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.s$E.qbic,ncol=Time)[,j]))))))
    bgRate_i.aic.p <- lapply(1:nsim, function(i) 
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.s$E.qaic,ncol=Time)[,j]))))))
    bgRate_i.aicc.p <- lapply(1:nsim, function(i) 
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.s$E.qaicc,ncol=Time)[,j]))))))
    
    bgRate.bic.p <- lapply(1:nsim, function(i) rep(bgRate_i.bic.p[[i]], each = 208))
    bgRate.aic.p <- lapply(1:nsim, function(i) rep(bgRate_i.aic.p[[i]], each = 208))
    bgRate.aicc.p <- lapply(1:nsim, function(i) rep(bgRate_i.aicc.p[[i]], each = 208))
    
    #detect 
    ix.bic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.s$E.qbic) - bgRate.bic.p[[i]])>=10^-3))
    ix.aic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.s$E.qaic) - bgRate.aic.p[[i]])>=10^-3))
    ix.aicc.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.s$E.qaicc) - bgRate.aicc.p[[i]])>=10^-3))
    
    simindicator.bic.p <- mapply(reval, position, ix.bic.p)
    simindicator.aic.p <- mapply(reval, position, ix.aic.p)
    simindicator.aicc.p <- mapply(reval, position, ix.aicc.p)
    
    probs.bic.p <- Matrix::rowSums(simindicator.bic.p)/nsim
    probs.aic.p <- Matrix::rowSums(simindicator.aic.p)/nsim
    probs.aicc.p <- Matrix::rowSums(simindicator.aicc.p)/nsim
    
    # #map probability detections by IC to grey scale
    # colprob.qp <- colormapping(list(probs.bic.qp,
    #                                 probs.aic.qp,
    #                                 probs.aicc.qp,
    #                                 as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob.qp,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_qp_clusso.pdf"))
    # 
    # colprob.p <- colormapping(list(probs.bic.p,
    #                                probs.aic.p,
    #                                probs.aicc.p,
    #                                as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob.p,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_p_clusso.pdf"))
    
    ##################################
    #POWER/FP Rate
    ##################################
    ##Power
    ###QP
    listpow.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                           sim_clusso[[i]]$lassoresult.qp.s$selections$select.qbic,rr,
                                                                           risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                            sim_clusso[[i]]$lassoresult.qp.s$selections$select.qaic,rr, 
                                                                            risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aicc.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                             sim_clusso[[i]]$lassoresult.qp.s$selections$select.qaicc,rr, 
                                                                             risk,nsim,Time, numCenters, pow=TRUE))
    outpow.bic.qp <- paste0(sum(unlist(listpow.bic.qp))/nsim*100, "%")
    outpow.aic.qp <- paste0(sum(unlist(listpow.aic.qp))/nsim*100, "%")
    outpow.aicc.qp <- paste0(sum(unlist(listpow.aicc.qp))/nsim*100, "%")
    ###Poisson
    listpow.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                          sim_clusso[[i]]$lassoresult.p.s$selections$select.qbic,rr,
                                                                          risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                           sim_clusso[[i]]$lassoresult.p.s$selections$select.qaic,rr, 
                                                                           risk,nsim,Time, numCenters, pow=TRUE))
    listpow.aicc.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                            sim_clusso[[i]]$lassoresult.p.s$selections$select.qaicc,rr, 
                                                                            risk,nsim,Time, numCenters, pow=TRUE))
    outpow.bic.p <- paste0(sum(unlist(listpow.bic.p))/nsim*100, "%")
    outpow.aic.p <- paste0(sum(unlist(listpow.aic.p))/nsim*100, "%")
    outpow.aicc.p <- paste0(sum(unlist(listpow.aicc.p))/nsim*100, "%")
    
    ##FP
    listfp.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                          sim_clusso[[i]]$lassoresult.qp.s$selections$select.qbic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                           sim_clusso[[i]]$lassoresult.qp.s$selections$select.qaic,rr, 
                                                                           risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aicc.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                            sim_clusso[[i]]$lassoresult.qp.s$selections$select.qaicc,rr, 
                                                                            risk,nsim,Time, numCenters, pow=FALSE))
    outfp.bic.qp <- paste0(sum(unlist(listfp.bic.qp))/nsim*100, "%")
    outfp.aic.qp <- paste0(sum(unlist(listfp.aic.qp))/nsim*100, "%")
    outfp.aicc.qp <- paste0(sum(unlist(listfp.aicc.qp))/nsim*100, "%")
    ###Poisson
    listfp.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                         sim_clusso[[i]]$lassoresult.p.s$selections$select.qbic,rr,
                                                                         risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                          sim_clusso[[i]]$lassoresult.p.s$selections$select.qaic,rr, 
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aicc.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                           sim_clusso[[i]]$lassoresult.p.s$selections$select.qaicc,rr, 
                                                                           risk,nsim,Time, numCenters, pow=FALSE))
    outfp.bic.p <- paste0(sum(unlist(listfp.bic.p))/nsim*100, "%")
    outfp.aic.p <- paste0(sum(unlist(listfp.aic.p))/nsim*100, "%")
    outfp.aicc.p <- paste0(sum(unlist(listfp.aicc.p))/nsim*100, "%")

    ##################################
    #Add sim results to table
    ##################################
    tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.bic.qp, fp = outfp.bic.qp, type = "QP"),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aic.qp, fp = outfp.aic.qp, type = "QP"),
                         cbind(IC="AICc",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aicc.qp, fp = outfp.aicc.qp, type = "QP"),
                         cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.bic.p, fp = outfp.bic.p, type = "Pois"),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aic.p, fp = outfp.aic.p, type = "Pois"),
                         cbind(IC="AICc",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", pow=outpow.aicc.p, fp = outfp.aicc.p, type = "Pois"))
    table.detection.clusso.space <- rbind(table.detection.clusso.space, tabn.clusso)
}   


# Stop the clock
proc.time() - ptm
#####################################################################################
#####################################################################################
#####################################################################################
##WRITE TO CSV
#####################################################################################
#####################################################################################
#####################################################################################
#superclust by loc
print(table.detection.loc.space)
write.csv(table.detection.loc.space, file=paste0(path.tables,"null_singlecluster_loc_space.csv"), row.names=TRUE)
#bound by loc
print(table.bounds.loc.space)
write.csv(table.bounds.loc.space, file=paste0(path.tables,"null_singlecluster_loc_space_bounds.csv"), row.names=TRUE)


#superclust by loc
print(table.detection.pc.space)
write.csv(table.detection.pc.space, file=paste0(path.tables,"null_singlecluster_pc_space.csv"), row.names=TRUE)
#bounds by pc
print(table.bounds.pc.space)
write.csv(table.bounds.pc.space, file=paste0(path.tables, "null_singlecluster_pc_space_bounds.csv"), row.names = TRUE)


#clusso
print(table.detection.clusso.space)
write.csv(table.detection.clusso.space, file=paste0(path.tables,"null_singlecluster_clusso_space.csv"), row.names=TRUE)



