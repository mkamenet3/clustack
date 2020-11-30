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
    sim_superclust_loc.time <- system.time(sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                    numCenters, Time, maxclust,
                                                                    bylocation = TRUE, model="poisson",
                                                                    overdisp.est = overdisp.est)))[[3]]
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
    id.bic_loc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)))
    id.aic_loc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)))
    #####################################################################################
    #####################################################################################
    if(any(id.bic_loc!=0)){
        outbic.loc_st <- calcbounds(id.bic_loc, IC="bic", sim_superclust_loc)
    } else {
        print("No clusters identified: BIC")
    }
    if(any(id.aic_loc!=0)){
        outaic.loc_st <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc)
    } else {
        print("No clusters identified: AIC")
    }
    
    #####################################################################################        
    outfp.bic_loc <- sum(ifelse(unlist(id.bic_loc)!=0,1,0))/nsim
    outfp.aic_loc <- sum(ifelse(unlist(id.aic_loc)!=0,1,0))/nsim
    #####################################################################################
    
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    #no identification - return 1 for all cells
    create_plotmeanrr_stack(sim_superclust_loc, IC="bic", flav="loc", Time=Time, nsim = nsim, sim.i)
    create_plotmeanrr_stack(sim_superclust_loc, IC="aic", flav="loc", Time=Time, nsim = nsim, sim.i)
   
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", fp = outfp.bic_loc,
                            time=sim_superclust_loc.time),
                      cbind(IC="AIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", fp = outfp.aic_loc,
                            time=sim_superclust_loc.time))
    table.detection.loc.st <- rbind(table.detection.loc.st, tabn.loc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(all(id.bic_loc==0)) & isTRUE(all(id.aic_loc==0))){
        bounds.loc <- cbind.data.frame(matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3))
    } else{
        bounds.loc <- cbind.data.frame(matrix(unlist(outbic.loc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outmata.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.loc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                       
                                       matrix(unlist(outaic.loc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outmata.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.loc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                       
                                       rep(outbic.loc_st$outnonma.time, nsim),
                                       rep(outbic.loc_st$outnonma_asymp.time, nsim),
                                       rep(outbic.loc_st$outbuck.theta.time, nsim),
                                       rep(outbic.loc_st$outbuckTlog.theta.time, nsim),
                                       rep(outbic.loc_st$outmaw2.theta.time, nsim),
                                       rep(outbic.loc_st$outmaw2Tlog.theta.time, nsim),
                                       rep(outbic.loc_st$outmata.theta.time, nsim),
                                       rep(outbic.loc_st$outmataT.theta.time, nsim),
                                       rep(outbic.loc_st$outmataTlog.theta.time, nsim),
                                       
                                       rep(outaic.loc_st$outnonma.time, nsim),
                                       rep(outaic.loc_st$outnonma_asymp.time, nsim),
                                       rep(outaic.loc_st$outbuck.theta.time, nsim),
                                       rep(outaic.loc_st$outbuckTlog.theta.time, nsim),
                                       rep(outaic.loc_st$outmaw2.theta.time, nsim),
                                       rep(outaic.loc_st$outmaw2Tlog.theta.time, nsim),
                                       rep(outaic.loc_st$outmata.theta.time, nsim),
                                       rep(outaic.loc_st$outmataT.theta.time, nsim),
                                       rep(outaic.loc_st$outmataTlog.theta.time, nsim))
    }

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
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- lapply(1:nsim, function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
                                                     family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
    }
    sim_superclust_pc.time <- system.time(sim_superclust_pc<- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                  numCenters, Time, maxclust,
                                                                  bylocation = FALSE, model="poisson",
                                                                  overdisp.est = overdisp.est)))[[3]]
    print("finished stacking: by PC")
    print(filename <- paste0(sim.i,"_superclustPC",".RData"))
    #save .RData
    save(sim_superclust_pc, file=filename)
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #NO forceidentify, save
    #####################################################################################
    #####################################################################################
    #####################################################################################
    id.bic_pc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)))
    id.aic_pc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic)))
    #####################################################################################
    #####################################################################################
    if(any(id.bic_pc!=0)){
        outbic.loc_st <- calcbounds(id.bic_pc, IC="bic", sim_superclust_pc)
    } else {
        print("No clusters identified: BIC")
    }
    if(any(id.aic_pc!=0)){
        outaic.loc_st <- calcbounds(id.aic_pc, IC="aic", sim_superclust_pc)
    } else {
        print("No clusters identified: AIC")
    }
    
    #####################################################################################        
    outfp.bic_pc <- sum(ifelse(unlist(id.bic_pc)!=0,1,0))/nsim
    outfp.aic_pc <- sum(ifelse(unlist(id.aic_pc)!=0,1,0))/nsim
    #####################################################################################
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    #no identification - return 1 for all cells
    create_plotmeanrr_stack(sim_superclust_pc, IC="bic", flav="pc", Time=Time, nsim = nsim, sim.i)
    create_plotmeanrr_stack(sim_superclust_pc, IC="aic", flav="pc", Time=Time, nsim = nsim, sim.i)
    
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", fp = outfp.bic_pc,
                           time=sim_superclust_pc.time),
                     cbind(IC="AIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", fp = outfp.aic_pc,
                           time=sim_superclust_pc.time ))
    table.detection.pc.st <- rbind(table.detection.pc.st, tabn.pc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(all(id.bic_pc==0)) & isTRUE(all(id.aic_pc==0))){
        bounds.pc <- cbind.data.frame(matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3))
    } else{
        bounds.pc <- cbind.data.frame(matrix(unlist(outbic.pc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outmata.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbic.pc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                       
                                       matrix(unlist(outaic.pc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outmata.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                       matrix(unlist(outaic.pc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                      
                                      rep(outbic.pc_st$outnonma.time, nsim),
                                      rep(outbic.pc_st$outnonma_asymp.time, nsim),
                                      rep(outbic.pc_st$outbuck.theta.time, nsim),
                                      rep(outbic.pc_st$outbuckTlog.theta.time, nsim),
                                      rep(outbic.pc_st$outmaw2.theta.time, nsim),
                                      rep(outbic.pc_st$outmaw2Tlog.theta.time, nsim),
                                      rep(outbic.pc_st$outmata.theta.time, nsim),
                                      rep(outbic.pc_st$outmataT.theta.time, nsim),
                                      rep(outbic.pc_st$outmataTlog.theta.time, nsim),
                                      
                                      rep(outaic.pc_st$outnonma.time, nsim),
                                      rep(outaic.pc_st$outnonma_asymp.time, nsim),
                                      rep(outaic.pc_st$outbuck.theta.time, nsim),
                                      rep(outaic.pc_st$outbuckTlog.theta.time, nsim),
                                      rep(outaic.pc_st$outmaw2.theta.time, nsim),
                                      rep(outaic.pc_st$outmaw2Tlog.theta.time, nsim),
                                      rep(outaic.pc_st$outmata.theta.time, nsim),
                                      rep(outaic.pc_st$outmataT.theta.time, nsim),
                                      rep(outaic.pc_st$outmataTlog.theta.time, nsim))
    }
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
    sim_clusso.time <- system.time(sim_clusso <- lapply(1:nsim, function(i) clusso::clusso(df=jbcSIM[[i]],
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
                                                            maxclust = maxclust)))[[3]]
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


    bgRate.bic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.bic.qp[[i]], each = 208))
    bgRate.aic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aic.qp[[i]], each = 208))

    #detect
    ix.bic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qbic) - bgRate.bic.qp[[i]])>=10^-3))
    ix.aic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qaic) - bgRate.aic.qp[[i]])>=10^-3))

    # simindicator.bic.qp <- mapply(reval, position, ix.bic.qp)
    # simindicator.aic.qp <- mapply(reval, position, ix.aic.qp)
# 
# 
#     probs.bic.qp <- Matrix::rowSums(simindicator.bic.qp)/nsim
#     probs.aic.qp <- Matrix::rowSums(simindicator.aic.qp)/nsim


    ##Poisson
    bgRate_i.bic.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qbic,ncol=Time)[,j]))))))
    bgRate_i.aic.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qaic,ncol=Time)[,j]))))))

    bgRate.bic.p <- lapply(1:nsim, function(i) rep(bgRate_i.bic.p[[i]], each = 208))
    bgRate.aic.p <- lapply(1:nsim, function(i) rep(bgRate_i.aic.p[[i]], each = 208))

    #detect
    ix.bic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qbic) - bgRate.bic.p[[i]])>=10^-3))
    ix.aic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qaic) - bgRate.aic.p[[i]])>=10^-3))

    # simindicator.bic.p <- mapply(reval, position, ix.bic.p)
    # simindicator.aic.p <- mapply(reval, position, ix.aic.p)
    # 
    # probs.bic.p <- Matrix::rowSums(simindicator.bic.p)/nsim
    # probs.aic.p <- Matrix::rowSums(simindicator.aic.p)/nsim

    
    ##FP
    listfp.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                          sim_clusso[[i]]$lassoresult.qp.st$selections$select.qbic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                           sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaic,rr,
                                                                           risk,nsim,Time, numCenters, pow=FALSE))

    outfp.bic.qp <- sum(unlist(listfp.bic.qp))/nsim
    outfp.aic.qp <- sum(unlist(listfp.aic.qp))/nsim

    ###Poisson
    listfp.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                         sim_clusso[[i]]$lassoresult.p.st$selections$select.qbic,rr,
                                                                         risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                          sim_clusso[[i]]$lassoresult.p.st$selections$select.qaic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    outfp.bic.p <- sum(unlist(listfp.bic.p))/nsim
    outfp.aic.p <- sum(unlist(listfp.aic.p))/nsim


    ##################################
    #RR maps
    ##################################
    ##################################
    #Add sim results to table
    ##################################
    tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST",  fp = outfp.bic.qp, type = "QP",
                               time=sim_clusso.time),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", fp = outfp.aic.qp, type = "QP",
                               time=sim_clusso.time),
                         cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", fp = outfp.bic.p, type = "Pois",
                               time=sim_clusso.time),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", fp = outfp.aic.p, type = "Pois",
                               time=sim_clusso.time))
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
    sim_superclust_loc.time <- system.time(sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                    numCenters, Time, maxclust,
                                                                    bylocation = TRUE, model="poisson",
                                                                    overdisp.est = overdisp.est)))[[3]]
    sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                    "risk", risk, "_", "theta", as.character(theta),"_spaceonly")
    print(filename <- paste0(sim.i,"_superclustLOC",".RData"))
    #save .RData
    save(sim_superclust_loc, file=filename)
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #NO forceidentify
    #####################################################################################
    #####################################################################################
    #####################################################################################
    id.bic_loc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)))
    id.aic_loc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)))
    #####################################################################################
    #####################################################################################
    if(any(id.bic_loc!=0)){
        outbic.loc_st <- calcbounds(id.bic_loc, IC="bic", sim_superclust_loc)
    } else {
        print("No clusters identified: BIC")
    }
    if(any(id.aic_loc!=0)){
        outaic.loc_st <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc)
    } else {
        print("No clusters identified: AIC")
    }
    
    #####################################################################################        
    outfp.bic_loc <- sum(ifelse(unlist(id.bic_loc)!=0,1,0))/nsim
    outfp.aic_loc <- sum(ifelse(unlist(id.aic_loc)!=0,1,0))/nsim
    #####################################################################################
    create_plotmeanrr_stack(sim_superclust_loc, IC="bic", flav="loc", Time=Time, nsim = nsim, sim.i)
    create_plotmeanrr_stack(sim_superclust_loc, IC="aic", flav="loc", Time=Time, nsim = nsim, sim.i)
    
    
    
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", fp = outfp.bic_loc,
                            time=sim_superclust_loc.time),
                      cbind(IC="AIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod="ST", fp = outfp.aic_loc,
                            time=sim_superclust_loc.time))
      table.detection.loc.space <- rbind(table.detection.loc.space, tabn.loc)
    
      #################################
      #Add clustackbounds to table
      ##################################
      if(isTRUE(all(id.bic_loc==0)) & isTRUE(all(id.aic_loc==0))){
          bounds.loc <- cbind.data.frame(matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                         matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3))
      } else{
          bounds.loc <- cbind.data.frame(matrix(unlist(outbic.loc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outmata.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outbic.loc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                         
                                         matrix(unlist(outaic.loc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outmata.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                         matrix(unlist(outaic.loc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                         
                                         rep(outbic.loc_st$outnonma.time, nsim),
                                         rep(outbic.loc_st$outnonma_asymp.time, nsim),
                                         rep(outbic.loc_st$outbuck.theta.time, nsim),
                                         rep(outbic.loc_st$outbuckTlog.theta.time, nsim),
                                         rep(outbic.loc_st$outmaw2.theta.time, nsim),
                                         rep(outbic.loc_st$outmaw2Tlog.theta.time, nsim),
                                         rep(outbic.loc_st$outmata.theta.time, nsim),
                                         rep(outbic.loc_st$outmataT.theta.time, nsim),
                                         rep(outbic.loc_st$outmataTlog.theta.time, nsim),
                                         
                                         rep(outaic.loc_st$outnonma.time, nsim),
                                         rep(outaic.loc_st$outnonma_asymp.time, nsim),
                                         rep(outaic.loc_st$outbuck.theta.time, nsim),
                                         rep(outaic.loc_st$outbuckTlog.theta.time, nsim),
                                         rep(outaic.loc_st$outmaw2.theta.time, nsim),
                                         rep(outaic.loc_st$outmaw2Tlog.theta.time, nsim),
                                         rep(outaic.loc_st$outmata.theta.time, nsim),
                                         rep(outaic.loc_st$outmataT.theta.time, nsim),
                                         rep(outaic.loc_st$outmataTlog.theta.time, nsim))
      }
      
      
      
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
    sim_superclust_pc.time <- system.time(sim_superclust_pc<- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                  numCenters, Time, maxclust,
                                                                  bylocation = FALSE, model="poisson",
                                                                  overdisp.est = overdisp.est)))[[3]]
    print(filename <- paste0(sim.i,"_superclustPC",".RData"))
    #save .RData
    save(sim_superclust_pc, file=filename)
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTACK BOUNDS
    #no forceidentify
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #####################################################################################
    id.bic_pc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic)))
    id.aic_pc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.aic)))
    
    #####################################################################################
    #####################################################################################
    if(any(id.bic_pc!=0)){
        outbic.loc_st <- calcbounds(id.bic_pc, IC="bic", sim_superclust_pc)
    } else {
        print("No clusters identified: BIC")
    }
    if(any(id.aic_pc!=0)){
        outaic.loc_st <- calcbounds(id.aic_pc, IC="aic", sim_superclust_pc)
    } else {
        print("No clusters identified: AIC")
    }
    
    #####################################################################################        
    outfp.bic_pc <- sum(ifelse(unlist(id.bic_pc)!=0,1,0))/nsim
    outfp.aic_pc <- sum(ifelse(unlist(id.aic_pc)!=0,1,0))/nsim
    #####################################################################################
    ##################################
    #RR maps
    ##################################
    #plot mean RR across sims
    #no identification - return 1 for all cells
    create_plotmeanrr_stack(sim_superclust_pc, IC="bic", flav="pc", Time=Time, nsim = nsim, sim.i)
    create_plotmeanrr_stack(sim_superclust_pc, IC="aic", flav="pc", Time=Time, nsim = nsim, sim.i)
    
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", fp = outfp.bic_pc,
                           time=sim_superclust_pc.time),
                     cbind(IC="AIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod="ST", fp = outfp.aic_pc,
                           time=sim_superclust_pc.time ))

    table.detection.pc.space <- rbind(table.detection.pc.space, tabn.pc)
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(all(id.bic_pc==0)) & isTRUE(all(id.aic_pc==0))){
        bounds.pc <- cbind.data.frame(matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3))
    } else{
        bounds.pc <- cbind.data.frame(matrix(unlist(outbic.pc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outmata.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outbic.pc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                      
                                      matrix(unlist(outaic.pc_st$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outbuck.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outmaw2.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outmata.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outmataT.theta), byrow=TRUE, ncol=3),
                                      matrix(unlist(outaic.pc_st$outmataTlog.theta), byrow=TRUE, ncol=3),
                                      
                                      rep(outbic.pc_st$outnonma.time, nsim),
                                      rep(outbic.pc_st$outnonma_asymp.time, nsim),
                                      rep(outbic.pc_st$outbuck.theta.time, nsim),
                                      rep(outbic.pc_st$outbuckTlog.theta.time, nsim),
                                      rep(outbic.pc_st$outmaw2.theta.time, nsim),
                                      rep(outbic.pc_st$outmaw2Tlog.theta.time, nsim),
                                      rep(outbic.pc_st$outmata.theta.time, nsim),
                                      rep(outbic.pc_st$outmataT.theta.time, nsim),
                                      rep(outbic.pc_st$outmataTlog.theta.time, nsim),
                                      
                                      rep(outaic.pc_st$outnonma.time, nsim),
                                      rep(outaic.pc_st$outnonma_asymp.time, nsim),
                                      rep(outaic.pc_st$outbuck.theta.time, nsim),
                                      rep(outaic.pc_st$outbuckTlog.theta.time, nsim),
                                      rep(outaic.pc_st$outmaw2.theta.time, nsim),
                                      rep(outaic.pc_st$outmaw2Tlog.theta.time, nsim),
                                      rep(outaic.pc_st$outmata.theta.time, nsim),
                                      rep(outaic.pc_st$outmataT.theta.time, nsim),
                                      rep(outaic.pc_st$outmataTlog.theta.time, nsim))
    }
    
    bounds.pc$risk <- risk
    bounds.pc$radius <- rad
    bounds.pc$select_orig.bic <- unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic))
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
    sim_clusso.time <- system.time(sim_clusso <- lapply(1:nsim, function(i) clusso::clusso(df=jbcSIM[[i]],
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
                                                            maxclust = maxclust)))[[3]]
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
    
    
    bgRate.bic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.bic.qp[[i]], each = 208))
    bgRate.aic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aic.qp[[i]], each = 208))
    
    #detect
    ix.bic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qbic) - bgRate.bic.qp[[i]])>=10^-3))
    ix.aic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.st$E.qaic) - bgRate.aic.qp[[i]])>=10^-3))
    
    # simindicator.bic.qp <- mapply(reval, position, ix.bic.qp)
    # simindicator.aic.qp <- mapply(reval, position, ix.aic.qp)
    # 
    # 
    #     probs.bic.qp <- Matrix::rowSums(simindicator.bic.qp)/nsim
    #     probs.aic.qp <- Matrix::rowSums(simindicator.aic.qp)/nsim
    
    
    ##Poisson
    bgRate_i.bic.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qbic,ncol=Time)[,j]))))))
    bgRate_i.aic.p <- lapply(1:nsim, function(i)
        sapply(1:Time,
               function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.st$E.qaic,ncol=Time)[,j]))))))
    
    bgRate.bic.p <- lapply(1:nsim, function(i) rep(bgRate_i.bic.p[[i]], each = 208))
    bgRate.aic.p <- lapply(1:nsim, function(i) rep(bgRate_i.aic.p[[i]], each = 208))
    
    #detect
    ix.bic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qbic) - bgRate.bic.p[[i]])>=10^-3))
    ix.aic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qaic) - bgRate.aic.p[[i]])>=10^-3))
    
    # simindicator.bic.p <- mapply(reval, position, ix.bic.p)
    # simindicator.aic.p <- mapply(reval, position, ix.aic.p)
    # 
    # probs.bic.p <- Matrix::rowSums(simindicator.bic.p)/nsim
    # probs.aic.p <- Matrix::rowSums(simindicator.aic.p)/nsim
    
    
    ##FP
    listfp.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                          sim_clusso[[i]]$lassoresult.qp.st$selections$select.qbic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.st,
                                                                           sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaic,rr,
                                                                           risk,nsim,Time, numCenters, pow=FALSE))
    
    outfp.bic.qp <- sum(unlist(listfp.bic.qp))/nsim
    outfp.aic.qp <- sum(unlist(listfp.aic.qp))/nsim
    
    ###Poisson
    listfp.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                         sim_clusso[[i]]$lassoresult.p.st$selections$select.qbic,rr,
                                                                         risk,nsim,Time, numCenters, pow=FALSE))
    listfp.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.st,
                                                                          sim_clusso[[i]]$lassoresult.p.st$selections$select.qaic,rr,
                                                                          risk,nsim,Time, numCenters, pow=FALSE))
    outfp.bic.p <- sum(unlist(listfp.bic.p))/nsim
    outfp.aic.p <- sum(unlist(listfp.aic.p))/nsim
    
    
    

    ##################################
    #Add sim results to table
    ##################################
    tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST",  fp = outfp.bic.qp, type = "QP",
                               time=sim_clusso.time),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", fp = outfp.aic.qp, type = "QP",
                               time=sim_clusso.time),
                         cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", fp = outfp.bic.p, type = "Pois",
                               time=sim_clusso.time),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod="ST", fp = outfp.aic.p, type = "Pois",
                               time=sim_clusso.time))
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



