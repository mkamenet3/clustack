rm(list=ls())
set.seed(20190326)
library(clusso)
library(MASS)
#source("clustack_20191222.R")
source("clustack.R") 



##################################################
#SINGLE CLUSTER SIMS
##################################################

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
japan.poly2 <- dframe.poly2[,2:3]
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




##################################################
######args definitions###########
########Single Cluster:
#args1 = theta1 (Inf (pure poisson))
#args2 = theta2 (60 (moderate overdispersion))
#args3 = nsims
#args4 = "space" or "spacetime"
##################################################
#Arguments passed will only control theta parameter
##R CMD BATCH "--args arg1 arg2" myscript.R -> 
args<- commandArgs(TRUE)

#quick function to recode
reval <- function(probs, ix){
    probs[ix] <-1
    return(probs)
}
#arguments passed
maxclust <- 15
nsimstep <- 1000
theta1 <- as.numeric(args[1])
theta2 <- as.numeric(args[2])
thetas <- c(theta1, theta2)
print(thetas)
nsim <- as.numeric(args[3])
model <- as.character(args[4])
print(model)
#scenarios
centers <- c(150,35)
radii <- c(9,11,18)
timeperiod <- c(3:5)
tim <- c(3:5)
risk.ratios <- c(1.1, 1.5, 2, 0.5)

#test
cent <- 150
rad <- 9
risk<- c(1.5)
tim <- c(3:5)
theta <- Inf
thetas <- c(Inf, 60)
nsim <- 5
model <- "spacetime"
models <- c("spacetime", "space")
nsimstep <- 10
path.figures <- "./"
path.tables <- "./"


eps <- 3
path.figures <- "../../../figures/OUTDEC2020/"
path.tables <- "../../../results/OUTDEC2020/"


#################################################################################################
#################################################################################################
#Set-Up for Condor
#################################################################################################
#################################################################################################
nsettings <-96
#2085  = 1040 (Yx) + 1040 (Ex) + cent + rad +risk +theta +mod (each row)
setsims <- matrix(rep(NA, 2085*nsettings*nsim), nrow=nsettings*nsim)
rrsims <- matrix(rep(NA, 1040*nsettings*nsim), nrow=nsettings*nsim)
#set global
rMax <- 20 
Time <- 5
tim <- 3:5

start <- 0
for (theta in thetas){
    for (cent in centers){
        for (rad in radii){
            for (risk in risk.ratios){
                for (mod in models){
                print(paste0("Params:", "center: ",cent,"; radius: ",
                             rad, "; Timeperiods: ", as.numeric(paste(tim, collapse = "")),
                             "; RR: ", risk, "; Theta :", theta, "; Model :", mod))
                simid <- 1:nsim
                #put the cluster in
                clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
                n <- length(x)
                init <- setVectors(japanbreastcancer$period,
                                   japanbreastcancer$expdeath, japanbreastcancer$death,
                                   covars=NULL, Time=Time)
                E1 <- init$E0
                #Ex <- scale(init, Time)
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
                #if (theta=="Inf"){
                if (is.infinite(theta)){
                    print("Inf theta: pure poisson")
                    YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
                } else {
                    print("OverD: NB")
                    YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
                }
               Ex <- scale_sim(YSIM, init, nsim, Time)
               #setsims <- lapply(1:nsim, function(i) list(YSIM[[i]], Ex[[i]], cent, rad, tim, risk, theta, mod))
               for (i in 1:nsim){
                   print(start)
                   print(i)
                   setsims[i+start,] <- c(YSIM[[i]], Ex[[i]],cent, rad, risk, theta, mod)
                   rrsims[i+start,] <- as.vector(rr)
               }
                    #print(start)
               start <- start+nsim

                }
                
            }
        }
    }
    
}

write.csv(setsims, file="setsims.csv", row.names = FALSE)
write.csv(rrsims, file="rrsims.csv", row.names = FALSE)

setsims <- read.csv("setsims.csv")
rrsims <- read.csv("rrsims.csv")
#setsims <- 
#################################################################################################
#################################################################################################
#SIM
#################################################################################################
#################################################################################################
# Start the clock!
#for(i in 1:nrow(setsims)){
#for(i in 1:nrow(setsims)){
    #set sim
    rr <- as.matrix(rrsims[i,])
    YSIM <- as.vector(as.matrix(setsims[i, 1:1040]))
    Ex <- as.vector(as.matrix(setsims[i, 1041:2080]))
    cent <- setsims[i,2081]
    rad <- setsims[i, 2082]
    risk <- setsims[i, 2083]
    theta <- setsims[i, 2084]
    mod <- setsims[i, 2085]
    #set global
    rMax <- 20 
    Time <- 5
    tim <- 3:5
    #create set of potential clusters based on distances
    potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
    n_uniq <- length(unique(potentialclusters$center))
    numCenters <- n_uniq
    #create giant sparse design matrix (single potential clusters)
    sparsematrix <- spacetimeMat(potentialclusters, numCenters, Time) 
    
    
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
        offset_reg <- lapply(1:nsim, function(i) glm(YSIM ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex)/Time)) + offset(log(Ex)),
                                                     family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
        
    }
    sim_superclust_loc.time <- system.time(sim_superclust_loc <- detectclusters(sparsematrix, Ex, YSIM, numCenters, Time, maxclust,
                                                                                                           bylocation = TRUE, mod="poisson", 
                                                                                                           overdisp.est = overdisp.est))[[3]]
    print("finished stacking: by LOC")
    if(mod=="space"){
        sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                        "risk", risk, "_", "theta", as.character(theta),"_spaceonly_",i)
    } else{
        sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                        "risk", risk, "_", "theta", as.character(theta),
                        as.numeric(paste(tim, collapse = "")),"_",i)
    }
    print(filename <- paste0(sim.i,"_superclustLOC_",i,".RData"))
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
        outbic.loc <- calcbounds(id.bic_loc, IC="bic", sim_superclust_loc, bylocation = TRUE)
    } else {
        print("No clusters identified: BIC")
    }
    if(id.aic_loc!=0){
        
        if(id.bic_loc==id.aic_loc){
            outaic.loc <- outbic.loc
        } else {
            outaic.loc <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc, bylocation = TRUE)
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
    ident.bic <- sim_superclust_loc$wtMAT[, sim_superclust_loc$selection.bic] 
    mat.bic <- t(ident.bic)%*%t(sparsematrix)

    ident.aic <- sim_superclust_loc$wtMAT[, sim_superclust_loc$selection.aic] 
    mat.aic <- t(ident.bic)%*%t(sparsematrix)
    
    #1) Did it find anything INSIDE the cluster?
    #start here
    incluster.bic <- mat.bic%*%matrix(rrbin_inside, ncol=1)
    incluster.aic <- mat.aic%*%matrix(rrbin_inside, ncol=1)
    
    #calc power
    #outpow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
    #outpow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
    #2) Did it find anything OUTSIDE the cluster?
    rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
    #this should be everything that doesn't touch the cluster
    outcluster.bic <- mat.bic%*%matrix(rrbin_outside, ncol=1)
    outcluster.aic <- mat.aic%*%matrix(rrbin_outside, ncol=1)
    
    #calc FP rate
    #outfp.bic <- sum(ifelse(unlist(outcluster.bic)==0,0,1))/nsim
    #outfp.aic <- sum(ifelse(unlist(outcluster.aic)==0,0,0))/nsim
    ##################################
    #plot probability maps
    ##################################
    #sim.bic <- matrix(unlist(mat.bic), ncol = nsim, byrow = FALSE)
    #sim.aic <- matrix(unlist(mat.aic), ncol = nsim, byrow = FALSE)
    #probs.bic <- Matrix::rowSums(sim.bic)/nsim
    #probs.aic <- Matrix::rowSums(sim.aic)/nsim
    #probs.aicc <- probs.aic #need to fix colormapping to allow for optional
    
    #map probability detections by IC to grey scale
    #colprob <- colormapping(list(probs.bic,
    #                             probs.aic,
    #                             probs.aicc,
    #                             as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    #probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_loc_overlap.pdf"))
    ##################################
    #Mean RR maps
    ##################################
    #create_plotmeanrr_stack(sim_superclust_loc, IC="bic", flav="loc", Time=Time, nsim = nsim, sim.i, greys=FALSE)
    #create_plotmeanrr_stack(sim_superclust_loc, IC="aic", flav="loc", Time=Time, nsim = nsim, sim.i, greys=FALSE)
    ##################################
    #Add sim results to table
    ##################################
    if(length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.bic,]))==0){
        wLambda.bic <- rep(1, 1040)
    } else {
        wLambda.bic <- as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.bic,])
    }
    if(length(as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.aic,]))==0){
        wLambda.aic <- rep(1, 1040)
    } else {
        wLambda.aic <- as.vector(sim_superclust_loc$wLambda[sim_superclust_loc$selection.aic,])
    }
           
           
    
    tabn.loc <- rbind(c(IC="BIC",rad, risk, cent, theta,
                            timeperiod=as.numeric(paste(tim, collapse="")),
                            mod=mod, type="NA",time=sim_superclust_loc.time, method="loc", 
                            incluster.bic =  ifelse(length(incluster.bic@x)==0,0, incluster.bic@x),
                            #incluster.aic = as.vector(incluster.aic),
                            outcluster.bic =  ifelse(length(outcluster.bic@x)==0,0, outcluster.bic@x),
                            iter=i,
                            #outcluster.aic = as.vector(outcluster.aic)
                            wLambda.bic),
                      c(IC="AIC",rad, risk, cent, theta,
                            timeperiod=as.numeric(paste(tim, collapse="")),
                            mod=mod, type="NA",time=sim_superclust_loc.time, method="loc",
                            #incluster.bic = as.vector(incluster.bic),
                            incluster.aic =  ifelse(length(incluster.aic@x)==0,0, incluster.aic@x),
                            #outcluster.bic = as.vector(outcluster.bic),
                            outcluster.aic = ifelse(length(outcluster.aic@x)==0,0, outcluster.aic@x),
                            iter=i,
                            wLambda.aic))

    #table.detection.loc <- rbind(table.detection.loc, tabn.loc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(id.bic_loc==0) & isTRUE(id.aic_loc==0)){
        print("No clusters identified")
        bounds.loc <- cbind.data.frame(matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                       
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1),
                                       matrix(rep(NA), byrow=TRUE, ncol=1))
        
    } else{
        if(id.bic_loc==0){
            fill.bic <- matrix(rep(NA,30*(length(which(id.bic_loc!=0)))),
                               nrow=length(which(id.bic_loc!=0)))
        } else {
            fill.bic <- NULL
        }
        if(id.aic_loc==0){
            fill.aic <- matrix(rep(NA,30*(length(which(id.aic_loc!=0)))),
                               nrow=length(which(id.aic_loc!=0)))
        } else {
            fill.aic <- NULL
        }
        print(str(fill.bic))
        print(str(fill.aic))
        bounds.loc <- cbind.data.frame(rbind.data.frame(cbind(matrix(unlist(outbic.loc$outnonma$nonma.theta), 
                                                                     byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outbuck.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outmaw2.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outmata.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outbic.loc$outmataTlog.theta), byrow=TRUE, ncol=3)),
                                                        fill.bic),
                                       
                                       rbind.data.frame(cbind(matrix(unlist(outaic.loc$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outbuck.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outmaw2.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outmata.theta), byrow=TRUE, ncol=3),
                                                              matrix(unlist(outaic.loc$outmataTlog.theta), byrow=TRUE, ncol=3)),
                                                        fill.aic),
                                       
                                       rep(outbic.loc$outnonma.time),
                                       rep(outbic.loc$outnonmaTlog.time),
                                       rep(outbic.loc$outnonma_asymp.time),
                                       rep(outbic.loc$outnonma_asympTlog.time),
                                       rep(outbic.loc$outbuck.theta.time),
                                       rep(outbic.loc$outbuckTlog.theta.time),
                                       rep(outbic.loc$outmaw2.theta.time),
                                       rep(outbic.loc$outmaw2Tlog.theta.time),
                                       rep(outbic.loc$outmata.theta.time),
                                       rep(outbic.loc$outmataTlog.theta.time),
                                       
                                       rep(outaic.loc$outnonma.time),
                                       rep(outaic.loc$outnonmaTlog.time),
                                       rep(outaic.loc$outnonma_asymp.time),
                                       rep(outaic.loc$outnonma_asympTlog.time),
                                       rep(outaic.loc$outbuck.theta.time),
                                       rep(outaic.loc$outbuckTlog.theta.time),
                                       rep(outaic.loc$outmaw2.theta.time),
                                       rep(outaic.loc$outmaw2Tlog.theta.time),
                                       rep(outaic.loc$outmata.theta.time),
                                       rep(outaic.loc$outmataTlog.theta.time))
    }
    
    bounds.loc$risk <- risk
    bounds.loc$radius <- rad
    
    bounds.loc$iter <- i
    bounds.loc$theta <- theta
    bounds.loc$cent <- cent
    
    bounds.loc$select_orig.bic <- sim_superclust_loc$selection.bic
    bounds.loc$select_orig.aic <- sim_superclust_loc$selection.aic
    bounds.loc$simID <- 1:nrow(bounds.loc)
    bounds.loc$method <- "LOC"
    
    names(bounds.loc) <- c("nonma.bic.LB", "clusterMA.bic.1", "nonma.bic.UB",
                           "nonmaTlog.bic.LB", "clusterMA.bic.2", "nonmaTlog.bic.UB",
                           "nonma_asymp.bic.LB", "clusterMA.bic.3", "nonma_asymp.bic.UB",
                           "nonma_asympTlog.bic.LB", "clusterMA.bic.4", "nonma_asympTlog.bic.UB",
                           "buck.bic.LB", "clusterMA.bic.5", "buck.bic.UB",
                           "bucklog.bic.LB", "clusterMA.bic.6", "bucklog.bic.UB",
                           "maw2.bic.LB", "clusterMA.bic.7", "maw2.bic.UB",
                           "maw2log.bic.LB", "clusterMA.bic.8", "maw2log.bic.UB",
                           "mata.bic.LB", "clusterMA.bic.9", "mata.bic.UB",
                           "matalog.bic.LB", "clusterMA.bic.10", "matalog.bic.UB",
                           
                           "nonma.aic.LB", "clusterMA.aic.1", "nonma.aic.UB",
                           "nonmaTlog.aic.LB", "clusterMA.aic.2", "nonmaTlog.aic.UB",
                           "nonma_asymp.aic.LB", "clusterMA.aic.3", "nonma_asymp.aic.UB",
                           "nonma_asympTlog.aic.LB", "clusterMA.aic.4", "nonma_asympTlog.aic.UB",
                           "buck.aic.LB", "clusterMA.aic.5", "buck.aic.UB",
                           "bucklog.aic.LB", "clusterMA.aic.6", "bucklog.aic.UB",
                           "maw2.aic.LB", "clusterMA.aic.7", "maw2.aic.UB",
                           "maw2log.aic.LB", "clusterMA.aic.8", "maw2log.aic.UB",
                           "mata.aic.LB", "clusterMA.aic.9", "mata.aic.UB",
                           "matalog.aic.LB", "clusterMA.aic.10", "matalog.aic.UB",
                           
                           "nonma.bic.time", "nonmaTlog.bic.time",
                           "nonma_asymp.bic.time", "nonma_asympTlog.bic.time", 
                           "buck.bic.time","bucklog.bic.time",
                           "maw2.bic.time", "maw2log.bic.time",
                           "mata.bic.time", "matalog.bic.time",
                           
                           "nonma.aic.time", "nonmaTlog.aic.time", 
                           "nonma_asymp.aic.time", "nonma_asympTlog.aic.time", 
                           "buck.aic.time","bucklog.aic.time",
                           "maw2.aic.time", "maw2log.aic.time",
                           "mata.aic.time", "matalog.aic.time",
                           "risk","rad","iter","theta", "cent",
                           "select_orig.bic", "select_orig.aic", "simID", "method")

    
    #table.bounds.loc <- rbind(table.bounds.loc, bounds.loc)
    
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
        offset_reg <- lapply(1:nsim, 
                             function(i) glm(YSIM ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                       each=length(Ex)/Time)) + offset(log(Ex)),
                                             family=quasipoisson))
        overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
        
    }
    
    sim_superclust_pc.time <- system.time(sim_superclust_pc<- detectclusters(sparsematrix, Ex, YSIM,numCenters, Time, maxclust,
                                                                                                        bylocation = FALSE, model="poisson",
                                                                                                        overdisp.est = overdisp.est))[[3]]
    print(filename <- paste0(sim.i,"_superclustPC_",i,".RData"))
    #save(sim_superclust_pc, file=filename)
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    id.bic_pc <- as.vector(unlist(sim_superclust_pc$selection.bic))
    id.aic_pc <- as.vector(unlist(sim_superclust_pc$selection.aic))
    #####################################################################################
    #####################################################################################
    if(id.bic_pc!=0){
        outbic.pc <- calcbounds(id.bic_pc, IC="bic", sim_superclust_pc, bylocation = FALSE)
    } else {
        print("No clusters identified: BIC")
    }
    if(id.aic_pc!=0){
        
        if(id.bic_pc==id.aic_pc){
            outaic.pc <- outbic.pc
        } else {
            outaic.pc <- calcbounds(id.aic_pc, IC="aic", sim_superclust_pc, bylocation = FALSE)
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
    
    
    
    #calc power
    #outpow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
    #outpow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
    
    #2) Did it find anything OUTSIDE the cluster?
    #this should be everything that doesn't touch the cluster
    outcluster.bic <- mat.bic%*%matrix(rrbin_outside, ncol=1)
    outcluster.aic <- mat.aic%*%matrix(rrbin_outside, ncol=1)
    
    
    #calc FP rate
    #outfp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
    #outfp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
    
    # ##################################
    # #plot probability maps
    # ##################################
    # sim.bic <- matrix(unlist(mat.bic), ncol = nsim, byrow = FALSE)
    # sim.aic <- matrix(unlist(mat.aic), ncol = nsim, byrow = FALSE)
    # probs.bic <- Matrix::rowSums(sim.bic)/nsim
    # probs.aic <- Matrix::rowSums(sim.aic)/nsim
    # probs.aicc <- probs.aic #need to fix colormapping to allow for optional
    # 
    # #map probability detections by IC to grey scale
    # colprob <- colormapping(list(probs.bic,
    #                              probs.aic,
    #                              probs.aicc,
    #                              as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
    # #plot map with probability detection by each IC
    # probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
    # 
    # ##################################
    # #Mean RR maps
    # ##################################
    # create_plotmeanrr_stack(sim_superclust_pc, IC="bic", flav="pc", Time=Time, nsim = nsim, sim.i, greys=FALSE)
    # create_plotmeanrr_stack(sim_superclust_pc, IC="aic", flav="pc", Time=Time, nsim = nsim, sim.i, greys=FALSE)
    
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
                       iter = i,
                       wLambda.bic),
                     
                     c(IC="AIC",rad, risk, cent, theta,
                       timeperiod=as.numeric(paste(tim, collapse="")),
                       mod=mod,type="NA",time=sim_superclust_pc.time, method="pc",
                       incluster.aic = ifelse(length(incluster.aic@x)==0,0, incluster.aic@x),
                       outcluster.aic = ifelse(length(outcluster.aic@x)==0,0, outcluster.aic@x),
                       iter = i,
                       wLambda.aic))
    #table.detection.pc <- rbind(table.detection.pc, tabn.pc)
    
    ##################################
    #Add clustackbounds to table
    ##################################
    if(isTRUE(id.bic_pc==0) & isTRUE(id.aic_pc==0)){
        print("No clusters identified")
        bounds.pc <- cbind.data.frame(matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      matrix(rep(NA,3), byrow=TRUE, ncol=3),
                                      
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1),
                                      matrix(rep(NA), byrow=TRUE, ncol=1))
    } else{
        #how many empty rows?
        if(id.bic_pc==0){
            fill.bic <- matrix(rep(NA,30*(length(which(id.bic_pc!=0)))),
                               nrow=length(which(id.bic_pc!=0)))
        } else {
            fill.bic <- NULL
        }
        if(id.aic_pc==0){
            fill.aic <- matrix(rep(NA,30*(length(which(id.aic_pc!=0)))),
                               nrow=length(which(id.aic_pc!=0)))
        } else {
            fill.aic <- NULL
        }
        
        
        bounds.pc <- cbind.data.frame(rbind.data.frame(cbind(matrix(unlist(outbic.pc$outnonma$nonma.theta), 
                                                                    byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outbuck.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outmaw2.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outmata.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outbic.pc$outmataTlog.theta), byrow=TRUE, ncol=3)),
                                                       fill.bic),
                                      
                                      rbind.data.frame(cbind(matrix(unlist(outaic.pc$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outbuck.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outmaw2.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outmata.theta), byrow=TRUE, ncol=3),
                                                             matrix(unlist(outaic.pc$outmataTlog.theta), byrow=TRUE, ncol=3)),
                                                       fill.aic),
                                      
                                      rep(outbic.pc$outnonma.time),
                                      rep(outbic.pc$outnonmaTlog.time),
                                      rep(outbic.pc$outnonma_asymp.time),
                                      rep(outbic.pc$outnonma_asympTlog.time),
                                      rep(outbic.pc$outbuck.theta.time),
                                      rep(outbic.pc$outbuckTlog.theta.time),
                                      rep(outbic.pc$outmaw2.theta.time),
                                      rep(outbic.pc$outmaw2Tlog.theta.time),
                                      rep(outbic.pc$outmata.theta.time),
                                      rep(outbic.pc$outmataTlog.theta.time),
                                      
                                      rep(outaic.pc$outnonma.time),
                                      rep(outaic.pc$outnonmaTlog.time),
                                      rep(outaic.pc$outnonma_asymp.time),
                                      rep(outaic.pc$outnonma_asympTlog.time),
                                      rep(outaic.pc$outbuck.theta.time),
                                      rep(outaic.pc$outbuckTlog.theta.time),
                                      rep(outaic.pc$outmaw2.theta.time),
                                      rep(outaic.pc$outmaw2Tlog.theta.time),
                                      rep(outaic.pc$outmata.theta.time),
                                      rep(outaic.pc$outmataTlog.theta.time))
    }
    
    bounds.pc$risk <- risk
    bounds.pc$radius <- rad
    
    bounds.pc$iter <- i
    bounds.pc$theta <- theta
    bounds.pc$cent <- cent
    
    bounds.pc$select_orig.bic <- sim_superclust_pc$selection.bic
    bounds.pc$select_orig.aic <- sim_superclust_pc$selection.aic
    bounds.pc$simID <- 1:nrow(bounds.pc)
    bounds.pc$method <- "PC"
    
    names(bounds.pc) <- c("nonma.bic.LB", "clusterMA.bic.1", "nonma.bic.UB",
                          "nonmaTlog.bic.LB", "clusterMA.bic.2", "nonmaTlog.bic.UB",
                          "nonma_asymp.bic.LB", "clusterMA.bic.3", "nonma_asymp.bic.UB",
                          "nonma_asympTlog.bic.LB", "clusterMA.bic.4", "nonma_asympTlog.bic.UB",
                          "buck.bic.LB", "clusterMA.bic.5", "buck.bic.UB",
                          "bucklog.bic.LB", "clusterMA.bic.6", "bucklog.bic.UB",
                          "maw2.bic.LB", "clusterMA.bic.7", "maw2.bic.UB",
                          "maw2log.bic.LB", "clusterMA.bic.8", "maw2log.bic.UB",
                          "mata.bic.LB", "clusterMA.bic.9", "mata.bic.UB",
                          "matalog.bic.LB", "clusterMA.bic.10", "matalog.bic.UB",
                          
                          "nonma.aic.LB", "clusterMA.aic.1", "nonma.aic.UB",
                          "nonmaTlog.aic.LB", "clusterMA.aic.2", "nonmaTlog.aic.UB",
                          "nonma_asymp.aic.LB", "clusterMA.aic.3", "nonma_asymp.aic.UB",
                          "nonma_asympTlog.aic.LB", "clusterMA.aic.4", "nonma_asympTlog.aic.UB",
                          "buck.aic.LB", "clusterMA.aic.5", "buck.aic.UB",
                          "bucklog.aic.LB", "clusterMA.aic.6", "bucklog.aic.UB",
                          "maw2.aic.LB", "clusterMA.aic.7", "maw2.aic.UB",
                          "maw2log.aic.LB", "clusterMA.aic.8", "maw2log.aic.UB",
                          "mata.aic.LB", "clusterMA.aic.9", "mata.aic.UB",
                          "matalog.aic.LB", "clusterMA.aic.10", "matalog.aic.UB",
                          
                          "nonma.bic.time", "nonmaTlog.bic.time",
                          "nonma_asymp.bic.time", "nonma_asympTlog.bic.time", 
                          "buck.bic.time","bucklog.bic.time",
                          "maw2.bic.time", "maw2log.bic.time",
                          "mata.bic.time", "matalog.bic.time",
                          
                          "nonma.aic.time", "nonmaTlog.aic.time", 
                          "nonma_asymp.aic.time", "nonma_asympTlog.aic.time", 
                          "buck.aic.time","bucklog.aic.time",
                          "maw2.aic.time", "maw2log.aic.time",
                          "mata.aic.time", "matalog.aic.time",
                          "risk","rad","iter","theta", "cent",
                          "select_orig.bic", "select_orig.aic", "simID", "method")
    
    #table.bounds.pc <- rbind(table.bounds.pc, bounds.pc)
    
    
    
    
    
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY CLUSSO
    #####################################################################################
    #####################################################################################
    #####################################################################################
    YSIM1 <- as.vector(matrix(YSIM, ncol=Time, byrow = FALSE))
    E01 <- as.vector(matrix(init$E0, ncol=Time, byrow=FALSE))
    truth1 <- as.vector(matrix(init$Y.vec, ncol=Time, byrow = FALSE))
    period1 <- as.vector(matrix(init$Year, ncol=Time, byrow = FALSE))
    id <- rep(1:208, times = 5)
    #create list of dataframes
    jbcSIM <- cbind.data.frame(expected = E01,observed = YSIM1,period = period1,id = id)
    jbcSIM <- jbcSIM[order(jbcSIM$id, jbcSIM$period),]
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
    #print(filename <- paste0(sim.i,"_clusso",".RData"))
    #save .RData
    #save(sim_clusso, file=filename)
    ##################################
    #DIAGNOSTICS: #calc power and FB rate
    ##################################
    vec <- rep(0, 208*Time)
    position <- list(vec)[rep(1, nsim)]
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
                                                                               risk,nsim,Time, numCenters, pow=TRUE)
        listpow.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                                                sim_clusso$lassoresult.qp.s$selections$select.qaic,rr,
                                                                                risk,nsim,Time, numCenters, pow=TRUE)
        
        listfp.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                                              sim_clusso$lassoresult.qp.s$selections$select.qbic,rr,
                                                                              risk,nsim,Time, numCenters, pow=FALSE)
        listfp.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.s,
                                                                               sim_clusso$lassoresult.qp.s$selections$select.qaic,rr,
                                                                               risk,nsim,Time, numCenters, pow=FALSE)
        
        #outfp.bic.qp <- sum(unlist(listfp.bic.qp))/nsim
        #outfp.aic.qp <- sum(unlist(listfp.aic.qp))/nsim
        #outpow.bic.qp <- sum(unlist(listpow.bic.qp))/nsim
        #outpow.aic.qp <- sum(unlist(listpow.aic.qp))/nsim
        
        ###Poisson
        #QP
        listpow.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                                              sim_clusso$lassoresult.qp.s$selections$select.qbic,rr,
                                                                              risk,nsim,Time, numCenters, pow=TRUE)
        listpow.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                                               sim_clusso$lassoresult.qp.s$selections$select.qaic,rr,
                                                                               risk,nsim,Time, numCenters, pow=TRUE)
        listfp.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                                             sim_clusso$lassoresult.p.s$selections$select.qbic,rr,
                                                                             risk,nsim,Time, numCenters, pow=FALSE)
        listfp.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.s,
                                                                              sim_clusso$lassoresult.p.s$selections$select.qaic,rr,
                                                                              risk,nsim,Time, numCenters, pow=FALSE)
        #outfp.bic.p <- sum(unlist(listfp.bic.p))/nsim
        #outfp.aic.p <- sum(unlist(listfp.aic.p))/nsim
        #outpow.bic.p <- sum(unlist(listpow.bic.p))/nsim
        #outpow.aic.p <- sum(unlist(listpow.aic.p))/nsim
        
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
        ix.bic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qbic) - bgRate.bic.p[[i]])>=10^-3))
        ix.aic.p <- lapply(1:nsim, function(i) which(abs(as.vector(sim_clusso[[i]]$lassoresult.p.st$E.qaic) - bgRate.aic.p[[i]])>=10^-3))
        
        ############################
        ##FP and Power
        ############################
        #QP
        listpow.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                                               sim_clusso$lassoresult.qp.st$selections$select.qbic,rr,
                                                                               risk,nsim,Time, numCenters, pow=TRUE)
        listpow.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                                                sim_clusso$lassoresult.qp.st$selections$select.qaic,rr,
                                                                                risk,nsim,Time, numCenters, pow=TRUE)
        
        listfp.bic.qp<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                                              sim_clusso$lassoresult.qp.st$selections$select.qbic,rr,
                                                                              risk,nsim,Time, numCenters, pow=FALSE)
        listfp.aic.qp <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.qp.st,
                                                                               sim_clusso$lassoresult.qp.st$selections$select.qaic,rr,
                                                                               risk,nsim,Time, numCenters, pow=FALSE)
        
        #outfp.bic.qp <- sum(unlist(listfp.bic.qp))/nsim
        #outfp.aic.qp <- sum(unlist(listfp.aic.qp))/nsim
        #outpow.bic.qp <- sum(unlist(listpow.bic.qp))/nsim
        #outpow.aic.qp <- sum(unlist(listpow.aic.qp))/nsim
        
        ###Poisson
        #QP
        listpow.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                                              sim_clusso$lassoresult.qp.st$selections$select.qbic,rr,
                                                                              risk,nsim,Time, numCenters, pow=TRUE)
        listpow.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                                               sim_clusso$lassoresult.qp.st$selections$select.qaic,rr,
                                                                               risk,nsim,Time, numCenters, pow=TRUE)
        listfp.bic.p<- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                                             sim_clusso$lassoresult.p.st$selections$select.qbic,rr,
                                                                             risk,nsim,Time, numCenters, pow=FALSE)
        listfp.aic.p <- clusso_prob_clusteroverlap(sparsematrix,sim_clusso$lassoresult.p.st,
                                                                              sim_clusso$lassoresult.p.st$selections$select.qaic,rr,
                                                                              risk,nsim,Time, numCenters, pow=FALSE)
        #outfp.bic.p <- sum(unlist(listfp.bic.p))/nsim
        #outfp.aic.p <- sum(unlist(listfp.aic.p))/nsim
        #outpow.bic.p <- sum(unlist(listpow.bic.p))/nsim
        #outpow.aic.p <- sum(unlist(listpow.aic.p))/nsim
        
    }
    ##################################
    #Prob maps
    ##################################
    # ##AIC - QP
    # simindicator.aic.qp <- mapply(reval, position, ix.aic.qp)
    # probs.aic.qp <- Matrix::rowSums(simindicator.aic.qp)/nsim
    # plotmeanrr_stack(matrix(probs.aic.qp,ncol=Time), Time=Time, sim.i=sim.i, ic="AIC_QP", flav="clusso", greys=TRUE)
    # ##BIC - QP
    # simindicator.bic.qp <- mapply(reval, position, ix.bic.qp)
    # probs.bic.qp <- Matrix::rowSums(simindicator.bic.qp)/nsim
    # plotmeanrr_stack(matrix(probs.bic.qp,ncol=Time), Time=Time, sim.i=sim.i, ic="BIC_QP", flav="clusso", greys=TRUE)
    # 
    # ##AIC - P
    # simindicator.aic.p <- mapply(reval, position, ix.aic.p)
    # probs.aic.p <- Matrix::rowSums(simindicator.aic.p)/nsim
    # plotmeanrr_stack(matrix(probs.aic.p,ncol=Time), Time=Time, sim.i=sim.i, ic="AIC_P", flav="clusso", greys=TRUE)
    # ##BIC - P
    # simindicator.bic.p <- mapply(reval, position, ix.bic.p)
    # probs.bic.p <- Matrix::rowSums(simindicator.bic.p)/nsim
    # plotmeanrr_stack(matrix(probs.bic.p,ncol=Time), Time=Time, sim.i=sim.i, ic="BIC_P", flav="clusso", greys=TRUE)
    
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
                               iter = i,
                           rrest.bic.qp ),
                         c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod,type = "QP",time=sim_clusso.time, method="clusso",
                               incluster.aic = as.vector(listpow.aic.qp),
                               outcluster.aic = as.vector(listfp.aic.qp),
                               iter = i,
                               rrest.aic.qp ),
                         c(IC="BIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                               incluster.bic = as.vector(listpow.bic.p),
                               outcluster.bic = as.vector(listfp.bic.p),
                               iter = i,
                               rrest.bic.p ),
                         c(IC="AIC",rad, risk, cent, theta,
                               timeperiod=as.numeric(paste(tim, collapse="")),
                               mod=mod, type = "Pois",time=sim_clusso.time, method="clusso",
                               incluster.aic = as.vector(listpow.aic.p),
                               outcluster.aic = as.vector(listfp.aic.p),
                               iter = i,
                               rrest.aic.p))
    #table.detection.clusso <- rbind(table.detection.clusso, tabn.clusso)
    
    
    
    ###################################################################################
    ###################################################################################
    ###################################################################################
    #STEPWISE SPATIAL SCAN VARIANT
    ###################################################################################
    ###################################################################################
    ###################################################################################
    sim_stepscan.time <- system.time(sim_stepscan <- stepscan(YSIM, Ex, Time=5, sparsematrix, nsim=nsimstep, maxclust=maxclust))[[3]]
    
    simid_stepscan <- paste0(sim.i,"stepscan")
    #save(simid_stepscan, file = paste0(simid_stepscan,".RData"))
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
    # vec <- rep(0, 208*Time)
    # position <- list(vec)[rep(1, nsim)]
    # simindicator.stepscan <- mapply(reval, position, ixids)
    # probs.stepscan <- Matrix::rowSums(simindicator.stepscan)/nsim
    # #plot
    # plotmeanrr_stack(matrix(probs.stepscan,ncol=Time), Time=Time, sim.i=sim.i, ic="MC", flav="stepscan", greys=TRUE)
    # print("Finished stepwise scan")
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
                       iter = i,
                       rrest.mc)
                       
    #table.detection.stepscan <- rbind(table.detection.stepscan, tabn.stepscan)
    
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
    #save(simid_stage, file = paste0(simid_stage,".RData"))
    
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
    
    #cell detection
    # vec <- rep(0, 208*Time)
    # position <- list(vec)[rep(1, nsim)]
    # 
    # probs.stage.bic  <- Matrix::rowSums(matrix(unlist(outfp.stage$betaselect_bin.bic), ncol=nsim, byrow = FALSE))/nsim
    # plotmeanrr_stack(matrix(probs.stage.bic,ncol=Time), Time=Time, sim.i=sim.i, ic="BIC", flav="fstagewise", greys=TRUE)
    # 
    # probs.stage.aic  <- Matrix::rowSums(matrix(unlist(outfp.stage$betaselect_bin.aic), ncol=nsim, byrow = FALSE))/nsim
    # plotmeanrr_stack(matrix(probs.stage.aic,ncol=Time), Time=Time, sim.i=sim.i, ic="AIC", flav="fstagewise", greys=TRUE)
    
    
    tabn.fstagewise <- rbind(c(IC="BIC",rad, risk, cent, theta,
                                   timeperiod=as.numeric(paste(tim, collapse="")),
                                   mod=mod,  type="NA", time = sim_stage.time  ,method = "fstagewise",
                                   incluster.bic = as.vector(outpow.stage.bic),
                               outcluster.bic = as.vector(outfp.stage.bic),
                               iter =i ,
                               as.vector(sim_stage$RRbic)),
                             c(IC="AIC",rad, risk, cent, theta,
                                   timeperiod=as.numeric(paste(tim, collapse="")),
                                   mod=mod, type="NA", time =  sim_stage.time ,method = "fstagewise" ,
                               incluster.bic = as.vector(outpow.stage.bic),
                               outcluster.bic = as.vector(outfp.stage.bic),
                               iter = i,
                               as.vector(sim_stage$RRaic)))
    #table.detection.fstage <- rbind(table.detection.fstage, tabn.fstagewise )
    print("Finished forward stagewise")
    
    
#}


#####################################################################################
#####################################################################################
#####################################################################################
##WRITE TO CSV
#####################################################################################
#####################################################################################
#####################################################################################
if(model=="space"){
    #out_space <- rbind(table.detection.loc, table.detection.pc, table.detection.clusso, table.detection.stepscan, table.detection.fstage)
    out_space <- rbind(tabn.loc, tabn.pc, tabn.clusso, tabn.stepscan, tabn.fstagewise)
    print(out_space)
    write.csv(out_space, file=paste0(path.tables,"theta",theta,"_singlecluster_space.csv"), row.names = TRUE)
    
    out_bounds_space <- rbind(tabn.loc ,tabn.pc)
    print(out_bounds_space)
    write.csv(out_bounds_space, file=paste0(path.tables,"theta",theta,"_singlecluster_bounds_space.csv"), row.names = TRUE)
    
}else{
    #out_st <- rbind(table.detection.loc, table.detection.pc, table.detection.clusso, table.detection.stepscan, table.detection.fstage)
    out_st <- rbind(tabn.loc, tabn.pc, tabn.clusso, tabn.stepscan, tabn.fstagewise)
    print(out_st)
    write.csv(out_st, file=paste0(path.tables,"theta",theta,"_singlecluster_ST.csv"), row.names = TRUE)
    
    out_bounds_st <- rbind(table.bounds.loc ,table.bounds.pc)
    print(out_bounds_st)
    write.csv(out_bounds_st, file=paste0(path.tables,"theta",theta,"_singlecluster_bounds_ST.csv"), row.names = TRUE)
}



#rm(list=ls())