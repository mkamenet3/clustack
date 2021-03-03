###########################################################
###########################################################
###########################################################
#Cell Risk Sim
###########################################################
###########################################################
###########################################################
#Notes
#################
#2020-10-14: sampling cell scheme (because doing all 1040 cells takes too long
#esp for mata-based bounds)
#1) cell in cluster center
#2) cell close to cluster center
#3) cell on border, but inside cluster
#4) cell on border, but outside cluster
#5) cell far from cluster
##Do this only for a single period
##2021-03-01: add in by loc; increase nsim to 100
###########################################################


rm(list=ls())
library(clusso)
library(Matrix)

set.seed(20200914)
source("clustack.R")


###########################################################
#LOAD DATA
###########################################################
#0) Setup
#load("../data/japanbreastcancer.RData")
load("../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japanbreastcancer.RData")

#select period 3 only for space only (period 7938)
japanbreastcancerp3 <- droplevels(subset(japanbreastcancer, period=="7938"))



#load("../scripts/clustack/data/japanbreastcancer.RData")
cases <- japanbreastcancerp3$death
expected <- japanbreastcancerp3$expdeath
# centroids <- japanbreastcancer
# periods <- japanbreastcancer
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]


#set global
rMax <- 20 
Time <- 1
maxclust <-15
#maxclust <- 10
locLambdas <- vector("list", maxclust)
#create set of potential clusters based on distances
potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
n <- length(x)
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparseMAT <- spacetimeMat(potentialclusters, numCenters, Time) 



cleanlist <- function(outlist, nsim,bounds=FALSE){
    if(bounds==TRUE){
        list1  <- lapply(1:nsim, function(x) matrix(unlist(outlist[[x]]), ncol=3))
        listcombine <-  do.call(rbind,lapply(list1,matrix,ncol=3,byrow=FALSE))
    } else {
        list1 <- lapply(1:nsim, function(x) matrix(as.vector(outlist[[x]]), ncol=3))
        listcombine <-  do.call(rbind,lapply(list1,matrix,ncol=3,byrow=FALSE))   
    }
    return(listcombine)
}

#######################################################
#######################################################


#######################################################
radii <- c(9,11,18)
cent <- 150
risks <- c(1.1,1.5,2)
ecounts <- c(5,10,50,100,250,1000)#5#10#50#100#250#500#1000
overdisp.est <- NULL
nsim <- 100#2#5#50#100
tim <- 1#:5
out.loc <- NULL
out.pc <- NULL

#put the cluster in
for(risk in risks){
    for(ecount in ecounts){
        for(rad in radii){
            
            print(c(risk,ecount,rad))
            ##################################################################
            potentialclusters$tosample <- NA
            potentialclusters$tosample[potentialclusters$center==cent & potentialclusters$n==1] <- "cell1"
            potentialclusters$tosample[potentialclusters$center==cent & potentialclusters$n==2] <- "cell2"
            #on the border
            maxid <- max(potentialclusters$n[potentialclusters$center==cent & potentialclusters$r< rad])
            potentialclusters$tosample[potentialclusters$center==cent & potentialclusters$n==maxid] <- "cell3"
            potentialclusters$tosample[potentialclusters$center==cent & potentialclusters$n==(maxid+1)] <- "cell4"
            
            #cell completely outside: select cell 80
            #which(1:208 %in% potentialclusters$last[potentialclusters$center==150] == FALSE)
            potentialclusters$tosample[potentialclusters$center==80 & potentialclusters$n==1] <- "cell5"
            
            cellsix <- potentialclusters$last[!is.na(potentialclusters$tosample)]
            #################################################################
            clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
            n <- length(x)
            
            n_uniq <- length(unique(clusters$center))
            numCenters <- n_uniq
            tmp <- clusters[clusters$center==cent,]
            cluster <- tmp[(tmp$r <= rad),]
            rr = as.vector(matrix(1, nrow=n, ncol=Time))
            rr[cluster$last] <- risk
            
            E0 <- rep(ecount, 208)
            E1 <- rr*E0
            YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
            Ex <- lapply(1:nsim, function(i) E0)
            outExp <- lapply(1:nsim, function(i) t(sparseMAT)%*%Ex[[i]])
            outObs <- lapply(1:nsim, function(i) t(sparseMAT)%*%YSIM[[i]])
            ix <- which(rr==risk) #ix ST cells wi
            
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            #BY PC
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            sim_superclust_pc <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                                 numCenters, Time, maxclust,
                                                                                 byloc = FALSE, model="poisson",
                                                                                 overdisp.est = overdisp.est))
            
            
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            #BY LOC
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                           numCenters, Time, maxclust,
                                                                           byloc= TRUE, model="poisson",
                                                                           overdisp.est = overdisp.est))
            
            
            
            
            id.bic_pc <- as.vector(unlist(sim_superclust_pc$selection.bic))
            id.aic_pc <- as.vector(unlist(sim_superclust_pc$selection.aic))
            
            outbic.pc <- calcbounds(id.bic_loc, IC="bic", sim_superclust_pc, bylocation = FALSE, cellrates=TRUE, cellsix = cellsix)
            
            print("finished stacking by PC")
            ##################################################
            ##################################################
            #NON-MA VARIANCE
            ##################################################
            ##################################################
            clusterRRlarge <- lapply(1:nsim, function(i) YSIM[[i]][cellsix]/Ex[[i]][cellsix])
            se_clusterRRlarge <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/Ex[[i]][cellsix]))
           
            nonma.time <- system.time(nonma<- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]],
                                                     clusterMA = clusterRRlarge[[i]],
                                                     ub = clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]])))
           
            se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/YSIM[[i]][cellsix]))
            nonma_asymp.time <- system.time(nonma_asymp <- lapply(1:nsim, 
                                                                  function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                                                    clusterMA = clusterRRlarge[[i]],
                                                                                    ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]])))
            
            print("nonma finished")
            
            ##################################################
            ##################################################
            #Buckland 1997
            ##################################################
            ##################################################
            wslarge <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wtMAT[,sim_superclust_pc[[i]]$selection.bic_forceid])
            clusterRR_ilarge <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$Lambda_dense) 
            cluster_thetaa <- lapply(1:nsim, function(i) sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic_forceid,])
            
            
            outbuck.time <- system.time(outbuck <- lapply(1:nsim, 
                              function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]][,cellsix],
                                                         thetaa =cluster_thetaa[[i]][cellsix], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparseMAT), Ex[[i]][cellsix],
                                                         overdisp.est = NULL, cellrates = TRUE)))
            
            #str(outbuck)
            ##########################
            #Log-scale
            outbuckTlog.time <- system.time(outbuckTlog <- lapply(1:nsim, 
                                                                  function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]][,cellsix],
                                                                                             thetaa =cluster_thetaa[[i]][cellsix], 
                                                                                             w_q=wslarge[[i]], 
                                                                                             sparsematrix=t(sparseMAT), Ex[[i]][cellsix],
                                                                                             overdisp.est = NULL, transform=TRUE,
                                                                                             cellrates = TRUE)))
            #str(outbuckTlog)
            print("buckland finished")
           
            ##################################################
            ##################################################
            #MAW2 (B&A pg. 345)
            ##################################################
            ##################################################
            
            outmaw2.time <- system.time(outmaw2 <- lapply(1:nsim, function(i) maw2(thetai=clusterRR_ilarge[[i]][,cellsix], 
                                                                                   thetaa = cluster_thetaa[[i]][cellsix], 
                                                                                   w_q=wslarge[[i]], 
                                                                                   sparsematrix=t(sparseMAT), 
                                                                                   outExp = Ex[[i]][cellsix], overdisp.est = NULL,
                                                                                   cellrates = TRUE)))
            #str(outmaw2)
            
            ##########################
            #Log-scale
            outmaw2Tlog.time <- system.time(outmaw2Tlog  <- lapply(1:nsim, function(i) maw2(thetai=clusterRR_ilarge[[i]][,cellsix],
                                                                                            thetaa = cluster_thetaa[[i]][cellsix], 
                                                                                            w_q=wslarge[[i]], 
                                                                                            sparsematrix=t(sparseMAT), 
                                                                                            outExp = Ex[[i]][cellsix], 
                                                                                            overdisp.est = NULL,
                                                                                            transform=TRUE,
                                                                                            cellrates = TRUE)))
            #str(outmaw2Tlog)
            print("maw2 finished")
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds (for non-normal data)
            ##################################################
            ##################################################
            
            outmata.time <- system.time(outmata <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]][,cellsix], 
                                                                                         thetaa = cluster_thetaa[[i]][cellsix], 
                                                                                         w_q=wslarge[[i]], sparsematrix=t(sparseMAT),
                                                                                         outExp = Ex[[i]][cellsix],
                                                                                         overdisp.est = NULL, 
                                                                                         transform="none",
                                                                                         cellrates=TRUE)))
            #str(outmata)
            print("outmata finished")
            
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
            ##################################################
            ##################################################
            
            outmataTlog.time <- system.time(outmataTlog <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]][,cellsix], 
                                                                                                 thetaa = cluster_thetaa[[i]][cellsix], 
                                                                                                 w_q=wslarge[[i]], 
                                                                                                 sparsematrix=t(sparseMAT), 
                                                                                                 outExp = Ex[[i]][cellsix],
                                                                                                 overdisp.est = NULL, 
                                                                                                 transform="log",
                                                                                                 cellrates=TRUE)))
            #str(outmataTlog)
            print("outmatalog finished")
            ##################################################
            ##################################################
            #Create Master for Output
            ##################################################
            ##################################################
            master <- cbind.data.frame(cleanlist(nonma,nsim, bounds = FALSE),
                                       cleanlist(nonma_asymp,nsim, bounds = FALSE),
                                       cleanlist(outbuck,nsim, bounds = TRUE),
                                       cleanlist(outbuckTlog,nsim, bounds = TRUE),
                                       cleanlist(outmaw2,nsim, bounds = TRUE),
                                       cleanlist(outmaw2Tlog,nsim, bounds = TRUE),
                                       cleanlist(outmata,nsim, bounds = TRUE),
                                       cleanlist(outmataTlog,nsim, bounds = TRUE))

            master$risk <- risk
            master$exp <- ecount
            master$radius <- rad
            master$anyforced <- ifelse(unlist(lapply(1:nsim, function(i) sim_superclust_pc[[i]]$selection.bic))==0,"yes","no")
            master$simID <- rep(1:nsim, each=length(cellsix))#1:nrow(master)
            master$nonma.time <- rep(nonma.time[[3]], nsim*length(cellsix))
            master$nonma_asymp.time <- rep(nonma_asymp.time[[3]], nsim*length(cellsix))
            master$outbuck.time <- rep(outbuck.time[[3]], nsim*length(cellsix))
            master$outbuckTlog.time <- rep(outbuckTlog.time[[3]], nsim*length(cellsix))
            master$outmaw2.time <- rep(outmaw2.time[[3]], nsim*length(cellsix))
            master$outmaw2Tlog.time <- rep(outmaw2Tlog.time[[3]], nsim*length(cellsix))
            master$outmata.time <- rep(outmata.time[[3]], nsim*length(cellsix))
            master$outmataTlog.time <- rep(outmataTlog.time[[3]], nsim*length(cellsix))
            master$cellid <- rep(cellsix,nsim)
            
            
            names(master) <- c("nonma.LB", "clusterMA", "nonma.UB",
                               "nonma_asymp.LB", "clusterMA.1", "nonma_asymp.UB",
                               "buck.LB", "clusterMA.2", "buck.UB",
                               "bucklog.LB", "clusterMA.3", "bucklog.UB",
                               "maw2.LB", "clusterMA.6", "maw2.UB",
                               "maw2log.LB", "clusterMA.7", "maw2log.UB",
                               "mata.LB", "clusterMA.8", "mata.UB",
                               "matalog.LB", "clusterMA.9", "matalog.UB", 
                               "risk", "ecount", "rad", "anyforced","simID", 
                               "nonma.time", "nonma_asymp.time", "outbuck.time",
                               "outbucklog.time",
                               "outmaw2.time","outmaw2log.time", "outmata.time",
                              "outmataTlog.time", "cellid")
            out.pc <- rbind(out.pc, master)
            print(paste0(c(risk,ecount,rad), " finished"))
            
        }
    }
    
}
#write.csv(masterout, file="masterout_nsim50_cellrisk.csv")

write.csv(masterout, file="masterout_nsim100_cellrisk.csv")
rm(list=ls())



#######################################################
#######################################################
#Figure out how to identify cells to sample
#let's do this for center 1 for now


# #let's draw this
# ixsample <- which(!is.na(potentialclusters$tosample))
# 
# #test <- 
# 
# id <- rep(0,nrow(potentialclusters))
# id[ixsample] <-1
# ixcellsamples <- as.vector(sparseMAT %*% matrix(id, ncol=1))
# ixcellsamples <- ifelse(ixcellsamples!=0,1,0)
# cellsix <- which(ixcellsamples!=0)
# 
#get plotting info:
# dframe2 <- read.csv("../data/utmJapan.csv")
# dframe.poly2 <- read.csv("../data/japan_poly2.csv")
# japan.poly2 <- dframe.poly2[,2:3]
# dframe.prefect2 <- read.csv("../data/japan_prefect2.csv")
# japan.prefect2 <- dframe.prefect2[,2:5]
# x <- utmJapan$utmx/1000
# y <- utmJapan$utmy/1000
# japan.poly2 <- dframe.poly2[,2:3]
# japan.prefect2 <- dframe.prefect2[,2:5]
# 
# colors1 <-rep(NA,208)
# colors1[potentialclusters$last[!is.na(potentialclusters$tosample)]] <- "black"
# #adjustcolor("red",alpha.f=0.5)
# #colors1[which(ixcellsamples!=0)] <- "black"
# 
# 
# par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=colors1,border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,paste0("Cells to Sample"),cex=1.00)
# 
# 
# tmp <- potentialclusters[potentialclusters$center==cent,]
# cluster <- tmp[(tmp$r <= rad),]
# rr = matrix(1, nrow=n, ncol=Time)
# rr[cluster$last, 1] <- 2
# 
# colors2 <- rep("white", 208)
# colors2[which(as.vector(rr)!=1)] <- "red"
# 
# par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=colors2,border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# polygon(japan.poly2,col=adjustcolor(colors1, alpha=0.5),border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,paste0("Cluster+Sampled Cells"),cex=1.00)
# 
# 
# colors2 <- rep("white", 208)
# colors2[which(as.vector(rr)!=1)] <- "red"
# 
# par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=colors2,border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# polygon(japan.poly2,col=colors2,border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,paste0("Cluster"),cex=1.00)

# 
