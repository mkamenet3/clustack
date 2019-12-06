set.seed(20190326)
library(clusso)
source("R/clustack.R")


#0) Setup
load("data/japanbreastcancer.RData")
cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath
# centroids <- japanbreastcancer
# periods <- japanbreastcancer
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]


#set global
rMax <- 20 
Time <- 5
#maxclust <- 10
locLambdas <- vector("list", maxclust)
#create set of potential clusters based on distances

#create giant sparse design matrix (single potential clusters)
sparseMAT <- spacetimeMat(potentialclusters, numCenters, Time) 

#################################
#TESTER SIM
#################################
#nsims
nsim <- 10
#theta - overdispersion
theta <- 100
maxclust <- 15


#scenarios
centers <- c(150,35)
radii <- c(9,18)
timeperiod <- c(3:5)
#timeperiods <- list(c(3:5), c(1:2))
risk.ratios <- c(1, 1.1, 1.5, 2)
table.detection <- NULL
eps <- 8

#test
cent <- 150
rad <- 18
risk <- 2
tim <- c(3:5)

for (cent in centers){
    for (rad in radii){
        for (risk in risk.ratios){
            print(paste0("Params:", "center: ",cent,"; radius: ",
                         rad, "; Timeperiods: ", as.numeric(paste(tim, collapse = "")), 
                         "; RR: ", risk, "; Theta :", theta))
            #put the cluster in
            clusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
            n <- length(x)
            init <- clusso::setVectors(japanbreastcancer$period, japanbreastcancer$expdeath, japanbreastcancer$death,covars=NULL, Time=Time)
            E1 <- init$E0
            #Ex <- clusso::scale(init, Time)
            Yx <- init$Y.vec
            vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)    
            n_uniq <- length(unique(potentialclusters$center))
            numCenters <- n_uniq
            tmp <- clusters[clusters$center==cent,]
            cluster <- tmp[(tmp$r <= rad),]
            rr = matrix(1, nrow=n, ncol=Time)
            rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
            E1 <- as.vector(rr)*init$E0
            message(paste("Running model for periods",tim[1],"through", tail(tim, n=1)))
            #simulate data here
            YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
            Ex <- clusso::scale_sim(YSIM, init, nsim, Time)
            
            #run superclust by location for each sim (apply)
            sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                     numCenters, Time, maxclust,
                                                                     bylocation = TRUE, model="poisson"))
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_superclustLOC"))
            filename <- paste0(sim_superclust_loc,".RData")
            #plot
            # plotmap(sim_superclust_loc$wLambda[sim_superclust_loc$selection.bic,],pdfname = paste0(sim.i, "_bic.pdf"),genpdf = TRUE)
            # plotmap(sim_superclust_loc$wLambda[sim_superclust_loc$selection.aic,],pdfname = paste0(sim.i, "_aic"),genpdf = TRUE)
            
            #run superclust by PC for each sim (apply)
            sim_superclust_pc <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                            numCenters, Time, maxclust,
                                                                            bylocation = FALSE, model="poisson"))
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_superclustPC"))
            filename <- paste0(sim_superclust_pc,".RData")

            #if(null==TRUE){
                #TODO
        #}
        #else{}
            #calc power and FB rate
            rrbin_inside <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
            ident.bic <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,])
            ident.aic <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,])
            ident.aicc <- lapply(1:nsim, function(i) sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aicc,])
            #1) Did it find anything INSIDE the cluster?
            incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ident.bic[[i]],ncol=1))
            incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ident.aic[[i]],ncol=1))
            incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ident.aicc[[i]],ncol=1))
            
            pow.bic <- paste0(sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim*100, "%")
            pow.aic <- paste0(sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim*100, "%")
            pow.aicc <- paste0(sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim*100, "%")
            
            
            #2) Did it find anything OUTSIDE the cluster?
            rrbin_outside <- matrix(as.vector(ifelse(rr!=1,0,1)),nrow=1)
            outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ident.bic[[i]],ncol=1))
            outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ident.aic[[i]],ncol=1))
            outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ident.aicc[[i]],ncol=1))
            
            fp.bic <- paste0(sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim*100, "%")
            fp.aic <- paste0(sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim*100, "%")
            fp.aicc <- paste0(sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim*100, "%")
            
            
            
            
             #OUTPUT
            #detections
            # print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
            #                       "risk", risk, "_", "theta", as.character(theta),
            #                       as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson"))
            # filename <- paste0("OUTJune19/",sim.i,".RData")
            # save(res, file = filename)
            # ##TABLE
            # (tabn <- rbind("******",cbind(rad,risk,cent,theta,
            #                               time=as.numeric(paste(tim, collapse = "")),
            #                               mod = "ST", rbind("QuasiPois",res$detect_out.qp.st),
            #                               rbind("Pois",res$detect_out.p.st)),
            #                cbind(rad,risk,cent,theta,time=as.numeric(paste(tim, collapse = "")),
            #                      mod = "Space",
            #                      rbind("QuasiPois",res$detect_out.qp.s),
            #                      rbind("Pois",res$detect_out.p.s))))
            # table.detection <- rbind(table.detection, tabn)
            # ##Plots
            # pdfnamerr <- paste0("OUTJune19/",sim.i,"_RR.pdf")
            # pdfnameprob <- paste0("OUTJune19/",sim.i,"_prob.pdf")
            # easyplot(japan.prefect2, japan.poly2 ,pdfnamerr, res$rrcolors, mods,
            #          space="both",probmap=FALSE, obs=NULL, rr=TRUE)
            # easyplot(japan.prefect2, japan.poly2 , pdfnameprob, res$probcolors, mods,
            #          space="both",probmap=TRUE, obs=NULL, rr=FALSE)
            
            
        }
    }
}
