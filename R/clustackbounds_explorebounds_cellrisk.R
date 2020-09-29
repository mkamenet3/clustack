rm(list=ls())

library(clusso)
library(Matrix)

set.seed(20200914)
#source("../scripts/clustack/R/clustack.R")
source("clustack.R")


###########################################################
#LOAD DATA
###########################################################
#0) Setup
load("../data/japanbreastcancer.RData")
#load("../scripts/clustack/data/japanbreastcancer.RData")
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
maxclust <-12
#maxclust <- 10
locLambdas <- vector("list", maxclust)
#create set of potential clusters based on distances
potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
n <- length(x)
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparseMAT <- spacetimeMat(potentialclusters, numCenters, Time) 


#######################################################
radii <- c(9,11,18)
cent <- 150
risks <- c(1.1,1.5,2)
ecounts <- c(5,10,50,100,250,500,1000)#5#10#50#100#250#500#1000
overdisp.est <- NULL
nsim <-50#100
tim <- 1:5
masterout <- NULL

#put the cluster in
for(risk in risks){
    for(ecount in ecounts){
        for(rad in radii){
            
            print(c(risk,ecount,rad))
            
            clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
            n <- length(x)
            
            n_uniq <- length(unique(clusters$center))
            numCenters <- n_uniq
            tmp <- clusters[clusters$center==cent,]
            cluster <- tmp[(tmp$r <= rad),]
            rr = matrix(1, nrow=n, ncol=Time)
            rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
            
            E0 <- rep(ecount, 1040)
            E1 <- rr*E0
            YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
            Ex <- lapply(1:nsim, function(i) E0)
            outExp <- lapply(1:nsim, function(i) t(sparseMAT)%*%Ex[[i]])
            outObs <- lapply(1:nsim, function(i) t(sparseMAT)%*%YSIM[[i]])
            ix <- which(rr==risk) #ix ST cells wi
            
            
            sim_superclust_pc_large <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                                 numCenters, Time, maxclust,
                                                                                 bylocation = FALSE, model="poisson",
                                                                                 overdisp.est = overdisp.est))
            
            ##################################################
            ##################################################
            #NON-MA VARIANCE
            ##################################################
            ##################################################
            
            clusterRRlarge <- lapply(1:nsim, function(i) YSIM[[i]]/Ex[[i]])
                #lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$wLambda[sim_superclust_pc_large[[i]]$selection.bic_forceid,])
            # clusterRRlarge <- lapply(1:nsim, 
            #                          function(i) unique(sim_superclust_pc_large[[i]]$Lambda_dense[sim_superclust_pc_large[[i]]$maxpcs[sim_superclust_pc_large[[i]]$selection.bic_forceid],])[2])
            
            se_clusterRRlarge <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/Ex[[i]]))
            # se_clusterRRlarge <- lapply(1:nsim, 
            #                             function(i)sqrt(clusterRRlarge[[i]]/outExp[[i]][sim_superclust_pc_large[[i]]$maxpcs[sim_superclust_pc_large[[i]]$selection.bic_forceid]]))
            
            nonma.time <- system.time(nonma<- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]],
                                                     clusterMA = clusterRRlarge[[i]],
                                                     ub = clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]])))
            # nonma.time <- system.time(nonma <- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]], 
            #                                                                     clusterMA = clusterRRlarge[[i]],
            #                                                                     ub=clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]])))
            # 
            
            #n
            se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/YSIM[[i]]))
            # se_clusterRRlarge_asymp <- lapply(1:nsim, 
            #                                   function(i) sqrt(clusterRRlarge[[i]]/outObs[[i]][sim_superclust_pc_large[[i]]$maxpcs[sim_superclust_pc_large[[i]]$selection.bic_forceid]]))
            nonma_asymp.time <- system.time(nonma_asymp <- lapply(1:nsim, 
                                                                  function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                                                    clusterMA = clusterRRlarge[[i]],
                                                                                    ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]])))
            # nonma_asymp.time <- system.time(nonma_asymp <- lapply(1:nsim, 
            #                                                       function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
            #                                                                         clusterMA = clusterRRlarge[[i]],
            #                                                                         ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]])))
            
            
            ##################################################
            ##################################################
            #Buckland 1997
            ##################################################
            ##################################################
            wslarge <- lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$wtMAT[,sim_superclust_pc_large[[i]]$selection.bic_forceid])
            clusterRR_ilarge <- lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$Lambda_dense) #lapply(1:nsim, function(i)aa)
            
            #test <- wslarge[[1]]%*%clusterRR_uniqlarge[[1]]
            # clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc_large[[i]]$Lambda_dense), 
            #                                                          function(k) unique(sim_superclust_pc_large[[i]]$Lambda_dense[k,]))) 
            ##########################
            
            
            #cluster_thetaa <- lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$wLambda[sim_superclust_pc_large[[i]]$selection.bic_forceid,])
            #model-average clusterMA
            #clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
            #clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
            #clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
            #cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
            
            
            ##########################
            # bucklandbounds2 <- function(thetai,thetaa, w_q,sparsematrix, outExp,overdisp.est, transform=NULL) {
            #     #NT <- rowSums(sparsematrix)
            #     if(!is.null(transform)){
            #         #print("log-scale")
            #         #log transform
            #         if(!is.null(overdisp.est)){
            #             varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k]*outExp[k])))
            #         } else {
            #             varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))
            #         }
            #         withintheta <- (log(thetai) - log(thetaa))^2
            #         varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))^2
            #         var_thetaa <- as.vector(varthetas_w)
            #         UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            #         LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
            #         
            #     } else {
            #         #print("No transform")
            #         #no transform
            #         if(!is.null(overdisp.est)){
            #             varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
            #         } else {
            #             varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
            #         }
            #         withintheta <- (thetai - thetaa)^2
            #         test <- sapply(1:nrow(sparsematrix), function(k) sqrt(varthetai[,k] + withintheta[k,]))
            #         test2 <- matrix(w_q, nrow = 1)%*%t(test)
            #         #varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))
            #         var_thetaa <- as.vector(varthetas_w)
            #         UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            #         LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
            #     }
            #     return(list(buckland.LB = LBa,
            #                 clusterMA = thetaa,
            #                 buckland.UB = UBa))
            # }
            
            outbuck <- lapply(1:nsim, 
                              function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa =cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparseMAT), Ex[[i]],
                                                         overdisp.est = NULL, cellrates = TRUE))
            
            # outbuck.time <- system.time(outbuck <- lapply(1:nsim, 
            #                                               function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa =cluster_thetaa[[i]], 
            #                                                                          w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp[[i]],
            #                                                                          overdisp.est = NULL)))
            # outbuck
            ##########################
            #Log-scale
            outbuckTlog.time <- system.time(outbuckTlog <- lapply(1:nsim, 
                                                                  function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa =cluster_thetaa[[i]], 
                                                                                             w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp[[i]],
                                                                                             overdisp.est = NULL, transform=TRUE, cellrates = TRUE)))
            outbuckTlog
            
            ##################################################
            ##################################################
            #(adjusted) MAW1 (B&A pg. 164)
            ##################################################
            ##################################################
            
            
            outmaw1.time <- system.time(outmaw1 <- lapply(1:nsim, 
                                                          function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                           w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp[[i]], overdisp.est = NULL)))
            outmaw1
            
            ##########################
            #Log-scale
            outmaw1Tlog.time <- system.time(outmaw1Tlog <- lapply(1:nsim, 
                                                                  function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                                   w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp[[i]], overdisp.est = NULL,
                                                                                   transform=TRUE)))
            outmaw1Tlog
            
            
            
            ##################################################
            ##################################################
            #MAW2 (B&A pg. 345)
            ##################################################
            ##################################################
            
            outmaw2.time <- system.time(outmaw2 <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                                   w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp[[i]], overdisp.est = NULL)))
            outmaw2
            
            ##########################
            #Log-scale
            outmaw2Tlog.time <- system.time(outmaw2Tlog  <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                                            w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp[[i]], 
                                                                                            overdisp.est = NULL,
                                                                                            transform=TRUE)))
            outmaw2Tlog
            
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds (for non-normal data)
            ##################################################
            ##################################################
            
            outmata.time <- system.time(outmata <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                         thetaa = cluster_thetaa[[i]], 
                                                                                         w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp = outExp[[i]],
                                                                                         overdisp.est = NULL, transform="none")))
            outmata
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
            ##################################################
            ##################################################
            
            outmataTsqrt.time <- system.time(outmataTsqrt <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                   thetaa = cluster_thetaa[[i]], 
                                                                                                   w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp = outExp[[i]],
                                                                                                   overdisp.est = NULL, transform="sqrt")))
            outmataTsqrt
            
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
            ##################################################
            ##################################################
            
            outmataTlog.time <- system.time(outmataTlog <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                 thetaa = cluster_thetaa[[i]], 
                                                                                                 w_q=wslarge[[i]], sparsematrix=t(sparseMAT), outExp = outExp[[i]],
                                                                                                 overdisp.est = NULL, transform="log")))
            outmataTlog
            
            ##################################################
            ##################################################
            #Create Master for Output
            ##################################################
            ##################################################
            
            master <- cbind.data.frame(matrix(unlist(nonma), byrow=TRUE, ncol=3),
                                       matrix(unlist(nonma_asymp), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbuck), byrow=TRUE, ncol=3),
                                       matrix(unlist(outbuckTlog), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmaw1), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmaw1Tlog), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmaw2), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmaw2Tlog), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmata), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmataTsqrt), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmataTlog), byrow=TRUE, ncol=3))
            
            
            master$risk <- risk
            master$exp <- ecount
            master$radius <- rad
            master$anyforced <- ifelse(unlist(lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$selection.bic))==0,"yes","no")
            master$simID <- 1:nrow(master)
            master$nonma.time <- rep(nonma.time[[3]], nsim)
            master$nonma_asymp.time <- rep(nonma_asymp.time[[3]], nsim)
            master$outbuck.time <- rep(outbuck.time[[3]], nsim)
            master$outbuckTlog.time <- rep(outbuckTlog.time[[3]], nsim)
            master$outmaw1.time <- rep(outmaw1.time[[3]], nsim)
            master$outmaw1Tlog.time <- rep(outmaw1Tlog.time[[3]], nsim)
            master$outmaw2.time <- rep(outmaw2.time[[3]], nsim)
            master$outmaw2Tlog.time <- rep(outmaw2Tlog.time[[3]], nsim)
            master$outmata.time <- rep(outmata.time[[3]], nsim)
            master$outmataTsqrt.time <- rep(outmataTsqrt.time[[3]], nsim)
            master$outmataTlog.time <- rep(outmataTlog.time[[3]], nsim)
            
            
            names(master) <- c("nonma.LB", "clusterMA", "nonma.UB",
                               "nonma_asymp.LB", "clusterMA.1", "nonma_asymp.UB",
                               "buck.LB", "clusterMA.2", "buck.UB",
                               "bucklog.LB", "clusterMA.3", "bucklog.UB",
                               "maw1.LB", "clusterMA.4", "maw1.UB",
                               "maw1log.LB", "clusterMA.5", "maw1log.UB",
                               "maw2.LB", "clusterMA.6", "maw2.UB",
                               "maw2log.LB", "clusterMA.7", "maw2log.UB",
                               "mata.LB", "clusterMA.8", "mata.UB",
                               "matasqrt.LB", "clusterMA.9", "matasqrt.UB",
                               "matalog.LB", "clusterMA.10", "matalog.UB", 
                               "risk", "ecount", "rad", "anyforced","simID", 
                               "nonma.time", "nonma_asymp.time", "outbuck.time",
                               "outbucklog.time","outmaw1.time","outmaw1log.time",
                               "outmaw2.time","outmaw2log.time", "outmata.time",
                               "outmataTsqrt.time","outmataTlog.time")
            masterout <- rbind(masterout, master)
            
        }
    }
    
}
write.csv(masterout, file="masterout_nsim50_cellrisk.csv")
rm(list=ls())
