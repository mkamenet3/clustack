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
#quick function to recode
reval <- function(probs, ix){
    probs[ix] <-1
    return(probs)
}
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
eps <- 3

#test
cent <- 150
rad <- 11
risk <- 1.5
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
            
            #tester
            # test_loc <-detectclusters(sparseMAT, Ex[[2]], YSIM[[2]], numCenters, Time, maxclust, bylocation = TRUE, model="poisson")
            # #BIC
            # plotmap(test_loc$wLambda[test_loc$selection.bic,],genpdf = FALSE)
            # summary(test_loc$wLambda[test_loc$selection.bic,])
            # 
            # #AIC/AICc
            # plotmap(test_loc$wLambda[test_loc$selection.aic,],genpdf = FALSE)
            # summary(test_loc$wLambda[test_loc$selection.aic,])
            
            
            ###
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_superclustLOC"))
            simid <- 1:nsim
            filename <- paste0(sim_superclust_loc,".RData")


            #if(null==TRUE){
                #TODO
        #}
        #else{}
            
            #calc power and FB rate
            rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
            clusteroverlap<- rrbin_cluster%*%sparseMAT
            rrbin_inside <- ifelse(sparseMAT%*%t(clusteroverlap)!=0,1,0)
            
            #identify what is different from background
            # bgRate_i <- lapply(1:nsim, function(i) sapply(1:Time,
            #                                               function(j) as.numeric(names(which.max(table(matrix(as.vector(select_mu[[i]]),
            #                                                                                                   ncol=Time)[,j]))))))
            # 
            # 
            
            
            
            ident.bic <- lapply(1:nsim, function(i) round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,],eps))
            ident.aic <- lapply(1:nsim, function(i) round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,],eps))
            ident.aicc <- lapply(1:nsim, function(i) round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aicc,],eps))
            
            
            
            #1) Did it find anything INSIDE the cluster?
            incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            
            pow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
            pow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
            pow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
            
            outpow.bic <- paste0(pow.bic*100, "%")
            outpow.aic <- paste0(pow.aic*100, "%")
            outpow.aicc <- paste0(pow.aicc*100, "%")
            
            
            #2) Did it find anything OUTSIDE the cluster?
            rrbin_outside <- ifelse(sparseMAT%*%t(clusteroverlap)==0,1,0) #matrix(as.vector(ifelse(rr!=1,0,1)),nrow=1)
            #this should be everything that doesn't touch the cluster
            outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            
            fp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
            fp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
            fp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim

            outfp.bic <- paste0(fp.bic*100, "%")
            outfp.aic <- paste0(fp.aic*100, "%")
            outfp.aicc <- paste0(fp.aicc*100, "%")
            
            
            ##################################
            #plot probability maps
            
            #TODO - need to calculate still what's different from the background
            ##################################
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            ix.bic <- lapply(1:nsim, function(i) which(ifelse(ident.bic[[i]]==1,0,1) ==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse(ident.aic[[i]]==1,0,1) ==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse(ident.aicc[[i]]==1,0,1) ==1))
            
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            
            
            #colormapping
            colprob <- colormapping(list(probs.bic,
                                             probs.aic,
                                             probs.aicc,
                                             as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob.pdf"))
            # pb.p.s <- get_prob(lassoresult.p.s, initial.s,
            #                    E1.s ,n, Time, nsim,thresh)
            # probcolors.p.s <- colormapping(pb.p.s$probs, Time, cv = NULL, prob = TRUE)
            
            
            
            
            #plot each RR for each sim            
            lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,],
                                                    res.aic = sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,],
                                                    oracle = rr,
                                                    genpdf = TRUE,
                                                    pdfname = paste0(sim.i,"_simID",simid[i],".pdf")))
            
            
            #run superclust by PC for each sim (apply)
            sim_superclust_pc <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                           numCenters, Time, maxclust,
                                                                           bylocation = FALSE, model="poisson"))
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_superclustPC"))
            filename <- paste0(sim_superclust_pc,".RData")
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
