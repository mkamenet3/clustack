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
dframe1 <- read.csv("../../../../clusso-newpenalty/clusso/data/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("../../../../clusso-newpenalty/clusso/data/utmJapan.csv")
#dframe1 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/jap.breast.F.9.10.11.csv")
#dframe2 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

#dframe.poly2 <- read.csv("clusso-newpenalty/clusso/clusso/data/japan_poly2.csv")
dframe.poly2 <- read.csv("../../../../clusso-newpenalty/clusso/data/japan_poly2.csv")
#dframe.poly2 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japan_poly2.csv")
#japan.poly2 <- dframe.poly2[,2:3]
#dframe.prefect2 <- read.csv("clusso-newpenalty/clusso/clusso/data/japan_prefect2.csv")
dframe.prefect2 <- read.csv("../../../../clusso-newpenalty/clusso/data/japan_prefect2.csv")
#dframe.prefect2 <- read.csv("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japan_prefect2.csv")
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
#thetas = c(Inf, 1000)
overdispfloor = TRUE
nullmod <- TRUE
maxclust <- 15
nsimstep <- 1000 #for stepwise monte carlo

# center=1
# radius=18
# tim = c(1:5)
# risk.ratio =1
# nsim =2
# nsimstep <- 10 #for stepwise monte carlo
# theta = Inf
# model <- "ST" #or "space"
# thetas = c(Inf, 1000, 2)
# overdispfloor = TRUE
# nullmod <- TRUE

#arguments passed
#arguments passed
theta1 <- as.numeric(args[1])
theta2 <- as.numeric(args[2])
thetas <- c(theta1,theta2)
print(thetas)
nsim <- as.numeric(args[3])
print(nsim)
model <- as.character(args[4]) #space or spacetime
print(model)



table.detection.loc <- NULL
table.detection.pc <- NULL
table.detection.clusso <- NULL
table.detection.stepscan <- NULL
table.detection.fstage <- NULL

table.bounds.loc <- NULL
table.bounds.pc <- NULL

eps <- 3
path.figures <- "../../../figures/OUTDEC2020/"
path.tables <- "../../../results/OUTDEC2020/"
# path.figures <- "./"
# path.tables <- "./"


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
    if(model=="space"){
        allTime <- 1:Time
        rr[cluster$last, allTime[1]:tail(allTime, n=1)] <- risk
        message(paste("Running model for spaceonly"))
    } else{
        rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
        message(paste("Running model for periods",tim[1],"through", tail(tim, n=1)))
    }
    E1 <- as.vector(rr)*init$E0
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
    if(model=="space"){
        sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                        "risk", risk, "_", "theta", as.character(theta),"_spaceonly")
    } else{
        sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                        "risk", risk, "_", "theta", as.character(theta),
                        as.numeric(paste(tim, collapse = "")))
    }
    print(filename <- paste0(sim.i,"_superclustLOC",".RData"))
    #save .RData
    #save(sim_superclust_loc, file=filename)
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
 
         if(all(id.bic_loc==id.aic_loc)){
             outaic.loc_st <- outbic.loc_st
         } else {
             outaic.loc_st <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc)
         }
    } else {
        print("No clusters identified: AIC")
    }
    
    #####################################################################################        
    outfp.bic_loc <- sum(ifelse(unlist(id.bic_loc)!=0,1,0))/nsim
    outfp.aic_loc <- sum(ifelse(unlist(id.aic_loc)!=0,1,0))/nsim
    #####################################################################################
    
    ##################################
    #FPR Maps
    ##################################
    create_plotFPR_stack(sim_superclust_loc,IC="bic", flav="loc",Time=Time, nsim=nsim, sim.i)
    create_plotFPR_stack(sim_superclust_loc,IC="aic", flav="loc",Time=Time, nsim=nsim, sim.i)

    ##################################
    #Add sim results to table
    ##################################
    tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod=model, fp = outfp.bic_loc,
                            type="NA",
                            time=sim_superclust_loc.time,
                            method="LOC"),
                      cbind(IC="AIC",rad, risk, cent, theta,
                            time=as.numeric(paste(tim, collapse="")),
                            mod=model, fp = outfp.aic_loc,
                            type="NA",
                            time=sim_superclust_loc.time,
                            method="LOC"))
    table.detection.loc <- rbind(table.detection.loc, tabn.loc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(all(id.bic_loc==0)) & isTRUE(all(id.aic_loc==0))){
        print("No clusters identified")
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
    bounds.loc$method <- "LOC"
    
    names(bounds.loc) <- c("nonma.bic.LB", "clusterMA.bic", "nonma.bic.UB",
                       "nonma_asymp.bic.LB", "clusterMA.bic.1", "nonma_asymp.bic.UB",
                       "buck.bic.LB", "clusterMA.bic.2", "buck.bic.UB",
                       "bucklog.bic.LB", "clusterMA.bic.3", "bucklog.bic.UB",
                       "maw2.bic.LB", "clusterMA.bic.6", "maw2.bic.UB",
                       "maw2log.bic.LB", "clusterMA.bic.7", "maw2log.bic.UB",
                       "mata.bic.LB", "clusterMA.bic.8", "mata.bic.UB",
                       "matasqrt.bic.LB", "clusterMA.bic.9", "matasqrt.bic.UB",
                       "matalog.bic.LB", "clusterMA.bic.10", "matalog.bic.UB",
                       
                       "nonma.aic.LB", "clusterMA.aic", "nonma.aic.UB",
                       "nonma_asymp.aic.LB", "clusterMA.aic.1", "nonma_asymp.aic.UB",
                       "buck.aic.LB", "clusterMA.aic.2", "buck.aic.UB",
                       "bucklog.aic.LB", "clusterMA.aic.3", "bucklog.aic.UB",
                       "maw2.aic.LB", "clusterMA.aic.6", "maw2.aic.UB",
                       "maw2log.aic.LB", "clusterMA.aic.7", "maw2log.aic.UB",
                       "mata.aic.LB", "clusterMA.aic.8", "mata.aic.UB",
                       "matasqrt.aic.LB", "clusterMA.aic.9", "matasqrt.aic.UB",
                       "matalog.aic.LB", "clusterMA.aic.10", "matalog.aic.UB",
                       
                       "nonma.bic.time", "nonma_asymp.bic.time", "buck.bic.time",
                       "bucklog.bic.time","maw2.bic.time", "maw2log.bic.time",
                       "mata.bic.time", "matasqrt.bic.time", "matalog.bic.time",
                       
                       "nonma.aic.time", "nonma_asymp.aic.time", "buck.aic.time",
                       "bucklog.aic.time","maw2.aic.time", "maw2log.aic.time",
                       "mata.aic.time", "matasqrt.aic.time", "matalog.aic.time",
                       "risk","rad","select_orig.bic", "select_orig.aic", "simID", "method")
    
    table.bounds.loc <- rbind(table.bounds.loc, bounds.loc)

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
    #save(sim_superclust_pc, file=filename)
    
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
        outbic.pc_st <- calcbounds(id.bic_pc, IC="bic", sim_superclust_pc)
    } else {
        print("No clusters identified: BIC")
    }
    if(any(id.aic_pc!=0)){
        
        if(all(id.bic_pc==id.aic_pc)){
            outaic.pc_st <- outbic.pc_st
        } else {
            outaic.pc_st <- calcbounds(id.aic_pc, IC="aic", sim_superclust_pc)
        }
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
    # create_plotmeanrr_stack(sim_superclust_pc, IC="bic", flav="pc", Time=Time, nsim = nsim, sim.i)
    # create_plotmeanrr_stack(sim_superclust_pc, IC="aic", flav="pc", Time=Time, nsim = nsim, sim.i)
    ##################################
    #FPR Maps
    ##################################
    create_plotFPR_stack(sim_superclust_pc,IC="bic", flav="pc",Time=Time, nsim=nsim, sim.i)
    create_plotFPR_stack(sim_superclust_pc,IC="aic", flav="pc",Time=Time, nsim=nsim, sim.i)
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod=model, fp = outfp.bic_pc,
                           type="NA",
                           time=sim_superclust_pc.time,
                           method="PC"),
                     cbind(IC="AIC",rad, risk, cent, theta,
                           time=as.numeric(paste(tim, collapse="")),
                           mod=model, fp = outfp.aic_pc,
                           type="NA",
                           time=sim_superclust_pc.time,
                           method="PC"))
    table.detection.pc <- rbind(table.detection.pc, tabn.pc)
    
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
    bounds.pc$method <- "PC"
    
    names(bounds.pc) <- c("nonma.bic.LB", "clusterMA.bic", "nonma.bic.UB",
                           "nonma_asymp.bic.LB", "clusterMA.bic.1", "nonma_asymp.bic.UB",
                           "buck.bic.LB", "clusterMA.bic.2", "buck.bic.UB",
                           "bucklog.bic.LB", "clusterMA.bic.3", "bucklog.bic.UB",
                           "maw2.bic.LB", "clusterMA.bic.6", "maw2.bic.UB",
                           "maw2log.bic.LB", "clusterMA.bic.7", "maw2log.bic.UB",
                           "mata.bic.LB", "clusterMA.bic.8", "mata.bic.UB",
                           "matasqrt.bic.LB", "clusterMA.bic.9", "matasqrt.bic.UB",
                           "matalog.bic.LB", "clusterMA.bic.10", "matalog.bic.UB",
                           
                           "nonma.aic.LB", "clusterMA.aic", "nonma.aic.UB",
                           "nonma_asymp.aic.LB", "clusterMA.aic.1", "nonma_asymp.aic.UB",
                           "buck.aic.LB", "clusterMA.aic.2", "buck.aic.UB",
                           "bucklog.aic.LB", "clusterMA.aic.3", "bucklog.aic.UB",
                           "maw2.aic.LB", "clusterMA.aic.6", "maw2.aic.UB",
                           "maw2log.aic.LB", "clusterMA.aic.7", "maw2log.aic.UB",
                           "mata.aic.LB", "clusterMA.aic.8", "mata.aic.UB",
                           "matasqrt.aic.LB", "clusterMA.aic.9", "matasqrt.aic.UB",
                           "matalog.aic.LB", "clusterMA.aic.10", "matalog.aic.UB",
                           
                           "nonma.bic.time", "nonma_asymp.bic.time", "buck.bic.time",
                           "bucklog.bic.time","maw2.bic.time", "maw2log.bic.time",
                           "mata.bic.time", "matasqrt.bic.time", "matalog.bic.time",
                           
                           "nonma.aic.time", "nonma_asymp.aic.time", "buck.aic.time",
                           "bucklog.aic.time","maw2.aic.time", "maw2log.aic.time",
                           "mata.aic.time", "matasqrt.aic.time", "matalog.aic.time",
                           "risk","rad","select_orig.bic", "select_orig.aic", "simID", "method")
    
    table.bounds.pc <- rbind(table.bounds.pc, bounds.pc)

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
    if(model=="space"){
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
    } else {
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
    }

    #print(filename <- paste0(sim.i,"_clusso",".RData"))
    #save .RData
    #save(sim_clusso, file=filename)
    #system(paste0("gzip ", filename))
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    vec <- rep(0, 208*Time)
    position <- list(vec)[rep(1, nsim)]
    if(model=="space"){
        #background rates
        ##Quasi-P
        bgRate_i.bic.qp <- lapply(1:nsim, function(i)
            sapply(1:Time,
                   function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.s$E.qbic,ncol=Time)[,j]))))))
        bgRate_i.aic.qp <- lapply(1:nsim, function(i)
            sapply(1:Time,
                   function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.qp.s$E.qaic,ncol=Time)[,j]))))))
        
        
        bgRate.bic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.bic.qp[[i]], each = 208))
        bgRate.aic.qp <- lapply(1:nsim, function(i) rep(bgRate_i.aic.qp[[i]], each = 208))
        
        #detect
        ix.bic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.s$E.qbic) - bgRate.bic.qp[[i]])>=10^-3))
        ix.aic.qp <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.qp.s$E.qaic) - bgRate.aic.qp[[i]])>=10^-3))
        
        ##Poisson
        bgRate_i.bic.p <- lapply(1:nsim, function(i)
            sapply(1:Time,
                   function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.s$E.qbic,ncol=Time)[,j]))))))
        bgRate_i.aic.p <- lapply(1:nsim, function(i)
            sapply(1:Time,
                   function(j) as.numeric(names(which.max(table(matrix(sim_clusso[[i]]$lassoresult.p.s$E.qaic,ncol=Time)[,j]))))))
        
        bgRate.bic.p <- lapply(1:nsim, function(i) rep(bgRate_i.bic.p[[i]], each = 208))
        bgRate.aic.p <- lapply(1:nsim, function(i) rep(bgRate_i.aic.p[[i]], each = 208))
        
        #detect
        ix.bic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.s$E.qbic) - bgRate.bic.p[[i]])>=10^-3))
        ix.aic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.s$E.qaic) - bgRate.aic.p[[i]])>=10^-3))
        
        ##FP
        listfp.bic.qp<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                              sim_clusso[[i]]$lassoresult.qp.s$selections$select.qbic,rr,
                                                                              risk,nsim,Time, numCenters, pow=FALSE))
        listfp.aic.qp <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.qp.s,
                                                                               sim_clusso[[i]]$lassoresult.qp.s$selections$select.qaic,rr,
                                                                               risk,nsim,Time, numCenters, pow=FALSE))
        
        outfp.bic.qp <- sum(unlist(listfp.bic.qp))/nsim
        outfp.aic.qp <- sum(unlist(listfp.aic.qp))/nsim
        
        ###Poisson
        listfp.bic.p<- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                             sim_clusso[[i]]$lassoresult.p.s$selections$select.qbic,rr,
                                                                             risk,nsim,Time, numCenters, pow=FALSE))
        listfp.aic.p <- lapply(1:nsim, function(i) clusso_prob_clusteroverlap(sparsematrix,sim_clusso[[i]]$lassoresult.p.s,
                                                                              sim_clusso[[i]]$lassoresult.p.s$selections$select.qaic,rr,
                                                                              risk,nsim,Time, numCenters, pow=FALSE))
    } else{
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
    }
    
    outfp.bic.p <- sum(unlist(listfp.bic.p))/nsim
    outfp.aic.p <- sum(unlist(listfp.aic.p))/nsim


    ##################################
    #Prob maps
    ##################################
    ##AIC - QP
    simindicator.aic.qp <- mapply(reval, position, ix.aic.qp)
    probs.aic.qp <- Matrix::rowSums(simindicator.aic.qp)/nsim
    plotmeanrr_stack(matrix(probs.aic.qp,ncol=Time), Time=Time, sim.i=sim.i, ic="AIC_QP", flav="clusso", greys=TRUE)
    ##BIC - QP
    simindicator.bic.qp <- mapply(reval, position, ix.bic.qp)
    probs.bic.qp <- Matrix::rowSums(simindicator.bic.qp)/nsim
    plotmeanrr_stack(matrix(probs.bic.qp,ncol=Time), Time=Time, sim.i=sim.i, ic="BIC_QP", flav="clusso", greys=TRUE)
    
    ##AIC - P
    simindicator.aic.p <- mapply(reval, position, ix.aic.p)
    probs.aic.p <- Matrix::rowSums(simindicator.aic.p)/nsim
    plotmeanrr_stack(matrix(probs.aic.p,ncol=Time), Time=Time, sim.i=sim.i, ic="AIC_P", flav="clusso", greys=TRUE)
    ##BIC - P
    simindicator.bic.p <- mapply(reval, position, ix.bic.p)
    probs.bic.p <- Matrix::rowSums(simindicator.bic.p)/nsim
    plotmeanrr_stack(matrix(probs.bic.p,ncol=Time), Time=Time, sim.i=sim.i, ic="BIC_P", flav="clusso", greys=TRUE)
    
    ##################################
    #Add sim results to table
    ##################################
    tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod=model,  fp = outfp.bic.qp, type = "QP",
                               time=sim_clusso.time,
                               method="clusso"),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod=model, fp = outfp.aic.qp, type = "QP",
                               time=sim_clusso.time,
                               method="clusso"),
                         cbind(IC="BIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod=model, fp = outfp.bic.p, type = "Pois",
                               time=sim_clusso.time,
                               method="clusso"),
                         cbind(IC="AIC",rad, risk, cent, theta,
                               time=as.numeric(paste(tim, collapse="")),
                               mod=model, fp = outfp.aic.p, type = "Pois",
                               time=sim_clusso.time,
                               method="clusso"))
    table.detection.clusso <- rbind(table.detection.clusso, tabn.clusso)
    
    ###################################################################################
    ###################################################################################
    ###################################################################################
    #STEPWISE SPATIAL SCAN VARIANT
    ###################################################################################
    ###################################################################################
    ###################################################################################
    sim_stepscan.time <- system.time(sim_stepscan <- lapply(1:nsim, function(i) stepscan(YSIM[[i]], Ex[[i]], Time=5, sparsematrix, nsim=nsimstep, maxclust=maxclust)))[[3]]
    
    simid_stepscan <- paste0(sim.i,"stepscan")
    #save(simid_stepscan, file = paste0(simid_stepscan,".RData"))
    ########################################################################
    #STEPWISE SCAN ANALYSIS
    ########################################################################
    #power/FP
    ##########################################
    numclustersid <- lapply(1:nsim, function(i) which(sim_stepscan[[i]]$pvals>0.05)-1)
    ixids <- lapply(1:nsim, function(i) step_clusterix(sparsematrix, sim_stepscan[[i]], numclustersid=numclustersid[[i]]))
    #pow/fp
    #outpow.stepscan<- spatscanfs_prob_clusteroverlap(sim_stepscan,ixids, numclustersid ,sparsematrix, rr, risk,pow=TRUE, nsim)
    outfp.stepscan <- spatscanfs_prob_clusteroverlap(sim_stepscan,ixids, numclustersid ,sparsematrix, rr, risk,pow=FALSE,nsim)
    #cell detection
    vec <- rep(0, 208*Time)
    position <- list(vec)[rep(1, nsim)]
    simindicator.stepscan <- mapply(reval, position, ixids)
    probs.stepscan <- Matrix::rowSums(simindicator.stepscan)/nsim
    #plot
    plotmeanrr_stack(matrix(probs.stepscan,ncol=Time), Time=Time, sim.i=sim.i, ic="MC", flav="stepscan", greys=TRUE)
    print("Finished stepwise scan")
    
    tabn.stepscan <- cbind(IC="MC",rad, risk, cent, theta,
          time="345",
          mod=model, fp = outfp.stepscan, type="NA", time =  sim_stepscan.time ,method = "stepscan" )
    table.detection.stepscan <- rbind(table.detection.stepscan, tabn.stepscan)
    
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
    
    
    sim_stage.time <- system.time(sim_stage <- lapply(1:nsim, function(i) cluster_model(delta,YSIM[[i]],Ex[[i]],sd, Time)))[[3]]
    simid_stage <- paste0(sim.i,"fstagewise")
    #save(simid_stage, file = paste0(simid_stage,".RData"))
    
    ########################################################################
    #FORWARD STAGEWISE ANALYSIS
    ########################################################################
    ##########################################
    #power/FP
    ##########################################
    #pow/fp
    outfp.stage <- forwardstage_prob_clusteroverlap(sim_stage,sparsematrix, rr, risk,pow=FALSE,nsim)
    outfp.stage.bic <- outfp.stage$fp.bic
    outfp.stage.aic <- outfp.stage$fp.aic
    
    #cell detection
    vec <- rep(0, 208*Time)
    position <- list(vec)[rep(1, nsim)]
    
    probs.stage.bic  <- Matrix::rowSums(matrix(unlist(outfp.stage$betaselect_bin.bic), ncol=nsim, byrow = FALSE))/nsim
    plotmeanrr_stack(matrix(probs.stage.bic,ncol=Time), Time=Time, sim.i=sim.i, ic="BIC", flav="fstagewise", greys=TRUE)
    
    probs.stage.aic  <- Matrix::rowSums(matrix(unlist(outfp.stage$betaselect_bin.aic), ncol=nsim, byrow = FALSE))/nsim
    plotmeanrr_stack(matrix(probs.stage.aic,ncol=Time), Time=Time, sim.i=sim.i, ic="AIC", flav="fstagewise", greys=TRUE)
    
    
    tabn.fstagewise <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                           time="345",
                           mod=model, fp = outfp.stage.bic, type="NA", time = sim_stage.time  ,method = "fstagewise" ),
                           cbind(IC="AIC",rad, risk, cent, theta,
                                 time="345",
                                 mod=model, fp = outfp.stage.aic, type="NA", time =  sim_stage.time ,method = "fstagewise" ))
    table.detection.fstage <- rbind(table.detection.fstage, tabn.fstagewise )
    print("Finished forward stagewise")
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
if(model=="space"){
    out_space <- rbind(table.detection.loc, table.detection.pc, table.detection.clusso, table.detection.stepscan, table.detection.fstage)
    print(out_space)
    write.csv(out_space, file=paste0(path.tables, "null_singlecluster_space.csv"), row.names = TRUE)
    
    out_bounds_space <- rbind(table.bounds.loc ,table.bounds.pc)
    print(out_bounds_space)
    write.csv(out_bounds_space, file=paste0(path.tables, "null_singlecluster_bounds_space.csv"), row.names = TRUE)
    
}else{
    out_st <- rbind(table.detection.loc, table.detection.pc, table.detection.clusso, table.detection.stepscan, table.detection.fstage)
    print(out_st)
    write.csv(out_st, file=paste0(path.tables, "null_singlecluster_ST.csv"), row.names = TRUE)
    
    out_bounds_st <- rbind(table.bounds.loc ,table.bounds.pc)
    print(out_bounds_st)
    write.csv(out_bounds_st, file=paste0(path.tables, "null_singlecluster_bounds_ST.csv"), row.names = TRUE)
}

