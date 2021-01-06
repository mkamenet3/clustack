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
# theta1 <- as.numeric(args[1]) 
# theta2 <- as.numeric(args[2]) 
# thetas <- c(theta1, theta2)
# print(thetas)
# nsim <- as.numeric(args[3])
# model <- as.character(args[4])
# print(model)
# #scenarios
# centers <- c(150,35)
# radii <- c(9,11,18)
# timeperiod <- c(3:5)
# tim <- c(3:5)
# risk.ratios <- c(1, 1.1, 1.5, 2)

#test
cent <- 150
rad <- 11
risk <- 1.5
tim <- c(3:5)
theta <- Inf
nsim <- 2
model <- "spacetime"



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



#################################################################################################
#################################################################################################
#SIM
#################################################################################################
#################################################################################################
# Start the clock!
ptm <- proc.time()
for (cent in centers){
    for (rad in radii){
        for (risk in risk.ratios){
            print(paste0("Params:", "center: ",cent,"; radius: ",
                         rad, "; Timeperiods: ", as.numeric(paste(tim, collapse = "")),
                         "; RR: ", risk, "; Theta :", theta))
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
            
            outExp <- lapply(1:nsim, function(i) t(sparsematrix)%*%Ex[[i]])
            outObs <- lapply(1:nsim, function(i) t(sparsematrix)%*%YSIM[[i]])
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
            #print(filename <- paste0(sim.i,"_superclustLOC",".RData"))
            #save .RData
            #save(sim_superclust_loc, file=filename)
            ##################################
            #DIAGNOSTICS: #calc power and FB rate
            ##################################
            id.bic_loc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic)))
            id.aic_loc <- as.vector(unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic)))
            #####################################################################################
            #####################################################################################
            if(any(id.bic_loc!=0)){
                outbic.loc <- calcbounds(id.bic_loc, IC="bic", sim_superclust_loc)
            } else {
                print("No clusters identified: BIC")
            }
            if(any(id.aic_loc!=0)){
                
                if(all(id.bic_loc==id.aic_loc)){
                    outaic.loc <- outbic.loc
                } else {
                    outaic.loc <- calcbounds(id.aic_loc, IC="aic", sim_superclust_loc)
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
            outpow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
            outpow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
            outpow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
            #2) Did it find anything OUTSIDE the cluster?
            rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
            #this should be everything that doesn't touch the cluster
            outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            #calc FP rate
            outfp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
            outfp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
            outfp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim
            ##################################
            #plot probability maps
            ##################################
            #create empties
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            #recode identified cells as 1's, all other zeros
            ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1)==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
            #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            #find probability of detection for each location in time
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            #map probability detections by IC to grey scale
            colprob <- colormapping(list(probs.bic,
                                         probs.aic,
                                         probs.aicc,
                                         as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            # #plot map with probability detection by each IC
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_loc.pdf"))
            ##################################
            #Mean RR maps
            ##################################
            create_plotmeanrr_stack(sim_superclust_loc, IC="bic", flav="loc", Time=Time, nsim = nsim, sim.i, greys=FALSE)
            create_plotmeanrr_stack(sim_superclust_loc, IC="aic", flav="loc", Time=Time, nsim = nsim, sim.i, greys=FALSE)
            ##################################
            #Add sim results to table
            ##################################
            tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                    time=as.numeric(paste(tim, collapse="")),
                                    mod=model, pow=outpow.bic, fp = outfp.bic),
                              cbind(IC="AIC",rad, risk, cent, theta,
                                    time=as.numeric(paste(tim, collapse="")),
                                    mod=model, pow=outpow.aic, fp = outfp.aic),
                              cbind(IC="AICc",rad, risk, cent, theta,
                                    time=as.numeric(paste(tim, collapse="")),
                                    mod=model, pow=outpow.aicc, fp = outfp.aicc))
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
                                               matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3),
                                               matrix(rep(NA,3*nsim), byrow=TRUE, ncol=3))
            } else{
                bounds.loc <- cbind.data.frame(matrix(unlist(outbic.loc$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outbuck.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outmaw2.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outmata.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outbic.loc$outmataTlog.theta), byrow=TRUE, ncol=3),
                                               
                                               matrix(unlist(outaic.loc$outnonma$nonma.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outnonmaTlog$nonma.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outnonma_asymp$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outnonma_asympTlog$nonma_asymp.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outbuck.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outbuckTlog.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outmaw2.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outmaw2Tlog.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outmata.theta), byrow=TRUE, ncol=3),
                                               matrix(unlist(outaic.loc$outmataTlog.theta), byrow=TRUE, ncol=3),
                                               
                                               rep(outbic.loc$outnonma.time, nsim),
                                               rep(outbic.loc$outnonmaTlog.time, nsim),
                                               rep(outbic.loc$outnonma_asymp.time, nsim),
                                               rep(outbic.loc$outnonma_asympTlog.time, nsim),
                                               rep(outbic.loc$outbuck.theta.time, nsim),
                                               rep(outbic.loc$outbuckTlog.theta.time, nsim),
                                               rep(outbic.loc$outmaw2.theta.time, nsim),
                                               rep(outbic.loc$outmaw2Tlog.theta.time, nsim),
                                               rep(outbic.loc$outmata.theta.time, nsim),
                                               rep(outbic.loc$outmataTlog.theta.time, nsim),
                                               
                                               rep(outaic.loc$outnonma.time, nsim),
                                               rep(outaic.loc$outnonmaTlog.time, nsim),
                                               rep(outaic.loc$outnonma_asymp.time, nsim),
                                               rep(outaic.loc$outnonma_asympTlog.time, nsim),
                                               rep(outaic.loc$outbuck.theta.time, nsim),
                                               rep(outaic.loc$outbuckTlog.theta.time, nsim),
                                               rep(outaic.loc$outmaw2.theta.time, nsim),
                                               rep(outaic.loc$outmaw2Tlog.theta.time, nsim),
                                               rep(outaic.loc$outmata.theta.time, nsim),
                                               rep(outaic.loc$outmataTlog.theta.time, nsim))
            }
            
            bounds.loc$risk <- risk
            bounds.loc$radius <- rad
            bounds.loc$select_orig.bic <- unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.bic))
            bounds.loc$select_orig.aic <- unlist(lapply(1:nsim, function(i) sim_superclust_loc[[i]]$selection.aic))
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
                                   "risk","rad","select_orig.bic", "select_orig.aic", "simID", "method")
            
            table.bounds.loc <- rbind(table.bounds.loc, bounds.loc)

##########################################################################################################################################################################
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
                                     function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
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
            ##################################
            #plot probability maps
            ##################################
            #create empties
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            #recode identified cells as 1's, all other zeros
            ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1) ==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
            #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            #find probability of detection for each location in time
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            #map probability detections by IC to grey scale
            colprob <- colormapping(list(probs.bic,
                                         probs.aic,
                                         probs.aicc,
                                         as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
            ##################################
            #RR maps
            ##################################
            # #plot each RR for each sim
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],
            #                                         res.aic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_pc_simID",simid[i],".pdf")))
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
            table.detection.pc.st <- rbind(table.detection.pc.st, tabn.pc)

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
            sim_clusso <- lapply(1:nsim, function(i) clusso(df=jbcSIM[[i]],
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

            #map probability detections by IC to grey scale
            colprob.qp <- colormapping(list(probs.bic.qp,
                                            probs.aic.qp,
                                            probs.aicc.qp,
                                            as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob.qp,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_qp_clusso.pdf"))

            colprob.p <- colormapping(list(probs.bic.p,
                                           probs.aic.p,
                                           probs.aicc.p,
                                           as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob.p,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_p_clusso.pdf"))

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
            # #plot each RR for each sim
            # ##QP
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_clusso[[i]]$lassoresult.qp.st$E.qbic,
            #                                         res.aic = sim_clusso[[i]]$lassoresult.qp.st$E.qaic,
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_qp_clusso_simID",simid[i],".pdf")))
            # ##P
            # ##QP
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_clusso[[i]]$lassoresult.p.st$E.qbic,
            #                                         res.aic = sim_clusso[[i]]$lassoresult.p.st$E.qaic,
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_p_clusso_simID",simid[i],".pdf")))
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
    }
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
write.csv(table.detection.loc.st, 
          file=paste0(path.tables,"theta",theta,"_singlecluster_loc_ST.csv"), row.names=TRUE)

#superclust by loc
print(table.detection.pc.st)
write.csv(table.detection.pc.st, 
          file=paste0(path.tables,"theta", theta,"_singlecluster_pc_ST.csv"), row.names=TRUE)

#clusso
print(table.detection.clusso.st)
write.csv(table.detection.clusso.st, 
          file=paste0(path.tables,"theta",theta,"_singlecluster_clusso_ST.csv"), row.names=TRUE)




#################################################################################################
#################################################################################################
#SPACE-ONLY
#################################################################################################
#################################################################################################
# Start the clock!
ptm <- proc.time()
for (cent in centers){
    for (rad in radii){
        for (risk in risk.ratios){
            print(paste0("Params:", "center: ",cent,"; radius: ",
                         rad, "; Timeperiods: ", "spaceonly", 
                         "; RR: ", risk, "; Theta :", theta))
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
            #rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
            allTime <- 1:Time
            rr[cluster$last, allTime[1]:tail(allTime, n=1)] <- risk
            E1 <- as.vector(rr)*init$E0
            message(paste("Running model for spaceonly"))
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
                offset_reg <- lapply(1:nsim, 
                                     function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
                                                                               each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
                                                     family=quasipoisson))
                overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
                
            }
            
            sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                            numCenters, Time, maxclust,
                                                                            bylocation = TRUE, model="poisson",
                                                                            overdisp.est = overdisp.est))
            sim.i <- paste0(path.figures,"sim","_", "center",cent,"_" ,"radius", rad, "_",
                            "risk", risk, "_", "theta", as.character(theta),"_spaceonly_")
            print(filename <- paste0(sim.i,"_superclustLOC",".RData"))
            #save .RData
            save(sim_superclust_loc, file=filename)
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
            ##################################
            #plot probability maps
            ##################################
            #create empties
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            #recode identified cells as 1's, all other zeros
            ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1)==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
            #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            #find probability of detection for each location in time
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            #map probability detections by IC to grey scale
            colprob <- colormapping(list(probs.bic,
                                         probs.aic,
                                         probs.aicc,
                                         as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_loc.pdf"))
            ##################################
            #RR maps
            ##################################
            # #plot each RR for each sim            
            # lapply(1:nsim, function(i) 
            #     plotmapAllIC(res.bic = sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,],
            #                  res.aic = sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,],
            #                  oracle = rr,
            #                  genpdf = TRUE,
            #                  pdfname = paste0(sim.i,"_loc_simID",simid[i],".pdf")))
            ##################################
            #Add sim results to table
            ##################################
            tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                    time="12345",
                                    mod="space", pow=outpow.bic, fp = outfp.bic),
                              cbind(IC="AIC",rad, risk, cent, theta,
                                    time="12345",
                                    mod="space", pow=outpow.aic, fp = outfp.aic),
                              cbind(IC="AICc",rad, risk, cent, theta,
                                    time="12345",
                                    mod="space", pow=outpow.aicc, fp = outfp.aicc))
            table.detection.loc.space <- rbind(table.detection.loc.space, tabn.loc)
            
            
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
                                     function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
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
            ##################################
            #plot probability maps
            ##################################
            #create empties
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            #recode identified cells as 1's, all other zeros
            ix.bic <- lapply(1:nsim, function(i) which(ifelse((length(ident.bic[[i]]!=0) & ident.bic[[i]]==1),0,1)==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse((length(ident.aic[[i]]!=0) & ident.aic[[i]]==1),0,1) ==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse((length(ident.aicc[[i]]!=0) & ident.aicc[[i]]==1),0,1) ==1))
            #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            #find probability of detection for each location in time
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            #map probability detections by IC to grey scale
            colprob <- colormapping(list(probs.bic,
                                         probs.aic,
                                         probs.aicc,
                                         as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
            ##################################
            #RR maps
            ##################################
            # #plot each RR for each sim
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],
            #                                         res.aic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_pc_simID",simid[i],".pdf")))
            ##################################
            #Add sim results to table
            ##################################
            tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                   time="12345",
                                   mod="space", pow=outpow.bic, fp = outfp.bic),
                             cbind(IC="AIC",rad, risk, cent, theta,
                                   time="12345",
                                   mod="space", pow=outpow.aic, fp = outfp.aic),
                             cbind(IC="AICc",rad, risk, cent, theta,
                                   time="12345",
                                   mod="space", pow=outpow.aicc, fp = outfp.aicc))
            table.detection.pc.space <- rbind(table.detection.pc.space, tabn.pc)
            
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
            sim_clusso <- lapply(1:nsim, function(i) clusso(df=jbcSIM[[i]],
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
            
            #map probability detections by IC to grey scale
            colprob.qp <- colormapping(list(probs.bic.qp,
                                            probs.aic.qp,
                                            probs.aicc.qp,
                                            as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob.qp,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_qp_clusso.pdf"))
            
            colprob.p <- colormapping(list(probs.bic.p,
                                           probs.aic.p,
                                           probs.aicc.p,
                                           as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob.p,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_p_clusso.pdf"))
            
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
            #RR maps
            ##################################
            # #plot each RR for each sim
            # ##QP
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_clusso[[i]]$lassoresult.qp.s$E.qbic,
            #                                         res.aic = sim_clusso[[i]]$lassoresult.qp.s$E.qaic,
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_qp_clusso_simID",simid[i],".pdf")))
            # ##P
            # ##QP
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_clusso[[i]]$lassoresult.p.s$E.qbic,
            #                                         res.aic = sim_clusso[[i]]$lassoresult.p.s$E.qaic,
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_p_clusso_simID",simid[i],".pdf")))
            ##################################
            #Add sim results to table
            ##################################
            tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                       time="12345",
                                       mod="space", pow=outpow.bic.qp, fp = outfp.bic.qp, type = "QP"),
                                 cbind(IC="AIC",rad, risk, cent, theta,
                                       time="12345",
                                       mod="space", pow=outpow.aic.qp, fp = outfp.aic.qp, type = "QP"),
                                 cbind(IC="AICc",rad, risk, cent, theta,
                                       time="12345",
                                       mod="space", pow=outpow.aicc.qp, fp = outfp.aicc.qp, type = "QP"),
                                 cbind(IC="BIC",rad, risk, cent, theta,
                                       time="12345",
                                       mod="space", pow=outpow.bic.p, fp = outfp.bic.p, type = "Pois"),
                                 cbind(IC="AIC",rad, risk, cent, theta,
                                       time="12345",
                                       mod="space", pow=outpow.aic.p, fp = outfp.aic.p, type = "Pois"),
                                 cbind(IC="AICc",rad, risk, cent, theta,
                                       time="12345",
                                       mod="space", pow=outpow.aicc.p, fp = outfp.aicc.p, type = "Pois"))
            table.detection.clusso.space <- rbind(table.detection.clusso.space, tabn.clusso)
        }   
    }
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
write.csv(table.detection.loc.space, 
          file=paste0(path.tables,"theta", theta,"_singlecluster_loc_space.csv"), row.names=TRUE)

#superclust by loc
print(table.detection.pc.space)
write.csv(table.detection.pc.space, 
          file=paste0(path.tables,"theta",theta,"_singlecluster_pc_space.csv"), row.names=TRUE)

#clusso
print(table.detection.clusso.space)
write.csv(table.detection.clusso.space, 
          file=paste0(path.tables,"theta",theta,"_singlecluster_clusso_space.csv"), row.names=TRUE)



