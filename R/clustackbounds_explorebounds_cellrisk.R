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
source("helperfuncs.R")


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



# cleanlist <- function(outlist, nsim,bounds=FALSE){
#     if(bounds==TRUE){
#         list1  <- lapply(1:nsim, function(x) matrix(unlist(outlist[[x]]), ncol=3))
#         listcombine <-  do.call(rbind,lapply(list1,matrix,ncol=3,byrow=FALSE))
#     } else {
#         list1 <- lapply(1:nsim, function(x) matrix(as.vector(outlist[[x]]), ncol=3))
#         listcombine <-  do.call(rbind,lapply(list1,matrix,ncol=3,byrow=FALSE))   
#     }
#     return(listcombine)
# }

cleanlist.clusters <- function(outlist, nsim){
    #nonma
    outnonma <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma$nonma.theta))))
    outnonma.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma$nonma.theta.time))))
    outnonmaTlog <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonmaTlog$nonma.theta))))
    outnonmaTlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonmaTlog$nonma.theta.time))))
    #nonma asymp
    outnonma_asymp <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asymp$nonma_asymp.theta))))
    outnonma_asymp.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asymp$nonma_asymp.theta.time))))
    outnonmaTlog_asymp <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asympTlog$nonma_asymp.theta))))
    outnonmaTlog_asymp.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asympTlog$nonma_asymp.theta.time))))
    #Buckland
    outbuck <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outbuck.theta))))
    outbuck.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outbuck.theta.time))))
    outbuckTlog <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outbuckTlog.theta))))
    outbuckTlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outbuckTlog.theta.time))))
    #Burnham & Anderson
    outmaw2 <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmaw2.theta))))
    outmaw2.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmaw2.theta.time))))
    outmaw2Tlog <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmaw2Tlog.theta))))
    outmaw2Tlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmaw2Tlog.theta.time))))
    #MATA
    outmata <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmata.theta))))
    outmata.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmata.theta.time))))
    outmataTlog <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmataTlog.theta))))
    outmataTlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmataTlog.theta.time))))
    
    res <- cbind.data.frame(outnonma, outnonmaTlog,
                 outnonma_asymp, outnonmaTlog_asymp,
                 outbuck, outbuckTlog,
                 outmaw2, outmaw2Tlog,
                 outmata, outmataTlog,
                 
                 outnonma.time, outnonmaTlog.time,
                 outnonma_asymp.time, outnonmaTlog_asymp.time,
                 outbuck.time, outbuckTlog.time,
                 outmaw2.time, outmaw2Tlog.time,
                 outmata.time, outmataTlog.time)
    names(res) <- c("nonma.LB", "clusterMA.1", "nonma.UB",
                    "nonmaTlog.LB", "clusterMA.2", "nonmaTlog.UB",
                    "nonma_asymp.LB", "clusterMA.3", "nonma_asymp.UB",
                    "nonmaTlog_asymp.LB", "clusterMA.4", "nonmaTlog_asymp.UB",
                    "buckland.LB", "clusterMA.5", "buckland.UB",
                    "bucklandTlog.LB", "clusterMA.6", "bucklandTlog.UB",
                    "maw2.LB", "clusterMA.7", "maw2.UB",
                    "maw2Tlog.LB", "clusterMA.8", "maw2Tlog.UB",
                    "mata.LB", "clusterMA.9", "mata.UB",
                    "mataTlog.LB", "clusterMA.10", "mataTlog.UB",
                    "nonma.time", "nonmaTlog.time", 
                    "nonma_asymp.time", "nonmaTlog.time",
                    "buckland.time", "bucklandTlog.time",
                    "maw2.time", "maw2Tlog.time",
                    "mata.time", "mataTlog.time")
    return(res)
    
}

cleanlist.cells <- function(outlist, nsim, cellsix){
    #nonma
    outnonma <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma$nonma.theta))))
    outnonma.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma$nonma.theta.time))))
    outnonmaTlog <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonmaTlog$nonma.theta))))
    outnonmaTlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonmaTlog$nonma.theta.time))))
    #nonma asymp
    outnonma_asymp <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asymp$nonma_asymp.theta))))
    outnonma_asymp.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asymp$nonma_asymp.theta.time))))
    outnonmaTlog_asymp <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asympTlog$nonma_asymp.theta))))
    outnonmaTlog_asymp.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outnonma_asympTlog$nonma_asymp.theta.time))))
    #Buckland
    outbuck <- do.call(rbind, lapply(1:nsim, function(i) matrix(unlist((outlist[[i]]$outbuck.theta)), ncol=3, byrow=FALSE)))
    outbuck.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outbuck.theta.time))))
    outbuckTlog <- do.call(rbind, lapply(1:nsim, function(i) matrix(unlist((outlist[[i]]$outbuckTlog.theta)), ncol=3, byrow=FALSE)))
    outbuckTlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outbuckTlog.theta.time))))
    #Burnham & Anderson
    outmaw2 <- do.call(rbind, lapply(1:nsim, function(i) matrix(unlist((outlist[[i]]$outmaw2.theta)), ncol=3, byrow=FALSE)))
    outmaw2.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmaw2.theta.time))))
    outmaw2Tlog <- do.call(rbind, lapply(1:nsim, function(i) matrix(unlist((outlist[[i]]$outmaw2Tlog.theta)), ncol=3, byrow = FALSE)))
    outmaw2Tlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmaw2Tlog.theta.time))))
    #MATA
    outmata <- do.call(rbind, lapply(1:nsim, function(i) matrix(unlist((outlist[[i]]$outmata.theta)), ncol=3, byrow=FALSE)))
    outmata.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmata.theta.time))))
    outmataTlog <- do.call(rbind, lapply(1:nsim, function(i) matrix(unlist((outlist[[i]]$outmataTlog.theta)), ncol=3, byrow=FALSE)))
    outmataTlog.time <- do.call(rbind, lapply(1:nsim, function(i) unlist((outlist[[i]]$outmataTlog.theta.time))))
    
    res <- cbind.data.frame(outnonma, outnonmaTlog,
                            outnonma_asymp, outnonmaTlog_asymp,
                            outbuck, outbuckTlog,
                            outmaw2, outmaw2Tlog,
                            outmata, outmataTlog,
                            
                            outnonma.time, outnonmaTlog.time,
                            outnonma_asymp.time, outnonmaTlog_asymp.time,
                            outbuck.time, outbuckTlog.time,
                            outmaw2.time, outmaw2Tlog.time,
                            outmata.time, outmataTlog.time)
    res$cellid <- rep(cellsix, times=nsim)
    names(res) <- c("nonma.LB", "clusterMA.1", "nonma.UB",
                    "nonmaTlog.LB", "clusterMA.2", "nonmaTlog.UB",
                    "nonma_asymp.LB", "clusterMA.3", "nonma_asymp.UB",
                    "nonmaTlog_asymp.LB", "clusterMA.4", "nonmaTlog_asymp.UB",
                    "buckland.LB", "clusterMA.5", "buckland.UB",
                    "bucklandTlog.LB", "clusterMA.6", "bucklandTlog.UB",
                    "maw2.LB", "clusterMA.7", "maw2.UB",
                    "maw2Tlog.LB", "clusterMA.8", "maw2Tlog.UB",
                    "mata.LB", "clusterMA.9", "mata.UB",
                    "mataTlog.LB", "clusterMA.10", "mataTlog.UB",
                    "nonma.time", "nonmaTlog.time", 
                    "nonma_asymp.time", "nonmaTlog.time",
                    "buckland.time", "bucklandTlog.time",
                    "maw2.time", "maw2Tlog.time",
                    "mata.time", "mataTlog.time", "cellid")
    return(res)
    
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
            
            id.bic_pc <- lapply(1:nsim, function(i) as.vector(unlist(sim_superclust_pc[[i]]$selection.bic)))
            id.aic_pc <- lapply(1:nsim, function(i) as.vector(unlist(sim_superclust_pc[[i]]$selection.aic)))
            
            
            #by PC
            ##Clusters
            outcluster.pc.bic <- lapply(1:nsim, function(i) calcbounds(id.bic_pc[[i]], IC="bic", 
                                                          sim_superclust_pc[[i]], byloc = FALSE,
                                                          Ex[[i]], YSIM[[i]], target="cluster", sparsemat = sparseMAT))
            outcluster.pc.aic <- lapply(1:nsim, function(i) calcbounds(id.bic_pc[[i]], IC="aic", 
                                                          sim_superclust_pc[[i]], byloc = FALSE,
                                                          Ex[[i]], YSIM[[i]], target="cluster", sparsemat = sparseMAT))
            ##Cells
            outcells.pc.bic <- lapply(1:nsim, function(i) calcbounds(id.bic_pc[[i]], IC="bic", 
                                                          sim_superclust_pc[[i]], byloc = FALSE,
                                                          Ex[[i]], YSIM[[i]], target="cells", cellsix = cellsix,
                                                          sparsemat = sparseMAT))
            outcells.pc.aic <- lapply(1:nsim, function(i) calcbounds(id.bic_pc[[i]], IC="aic", 
                                                          sim_superclust_pc[[i]], byloc = FALSE,
                                                          Ex[[i]], YSIM[[i]], target="cells", cellsix = cellsix,
                                                          sparsemat = sparseMAT))
    
            print("finished stacking & bounds by PC")
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
            
            id.bic_loc <- lapply(1:nsim, function(i) as.vector(unlist(sim_superclust_loc[[i]]$selection.bic)))
            id.aic_loc <- lapply(1:nsim, function(i) as.vector(unlist(sim_superclust_loc[[i]]$selection.aic)))
            
            outcluster.loc.bic <- lapply(1:nsim, function(i) calcbounds(id.bic_loc[[i]], IC="bic", 
                                                          sim_superclust_loc[[i]], byloc = TRUE,
                                                          Ex[[i]], YSIM[[i]], target="cluster", sparsemat = sparseMAT))
            
            outcluster.loc.aic <- lapply(1:nsim, function(i) calcbounds(id.bic_loc[[i]], IC="aic", 
                                                          sim_superclust_loc[[i]], byloc = TRUE,
                                                          Ex[[i]], YSIM[[i]], target="cluster", sparsemat = sparseMAT))
            
            outcells.loc.bic <- lapply(1:nsim, function(i) calcbounds(id.bic_loc[[i]], IC="bic", 
                                                          sim_superclust_loc[[i]], byloc = TRUE,
                                                          Ex[[i]], YSIM[[i]], target="cells", cellsix = cellsix,
                                                          sparsemat = sparseMAT))
            outcells.loc.aic <- lapply(1:nsim, function(i) calcbounds(id.bic_loc[[i]], IC="aic", 
                                                          sim_superclust_loc[[i]], byloc = TRUE,
                                                          Ex[[i]], YSIM[[i]], target="cells", cellsix = cellsix,
                                                          sparsemat = sparseMAT))
         
            
            print("finished stacking & bounds by LOC")
           
            ##################################################
            ##################################################
            #OUTPUT: CLUSTER
            ##################################################
            ##################################################
            #PC
            ##################################################
            master.bic<- cleanlist.clusters(outcluster.pc.bic, nsim=nsim)
            master.bic$IC <- "BIC"
            master.bic$select<-  unlist(id.bic_pc)
            
            master.aic<- cleanlist.clusters(outcluster.pc.aic, nsim=nsim)
            master.aic$IC <- "AIC"
            master.aic$select <-  unlist(id.aic_pc)

            master.cluster.pc <- rbind.data.frame(master.bic, master.aic)
            master.cluster.pc$method <- "PC"
            master.cluster.pc$simID <- rep(1:nsim, times=nsim)
            ##################################################
            #LOC
            ##################################################
            master.bic<- cleanlist.clusters(outcluster.loc.bic, nsim=nsim)
            master.bic$IC <- "BIC"
            master.bic$select<-  unlist(id.bic_loc)
            
            master.aic<- cleanlist.clusters(outcluster.loc.aic, nsim=nsim)
            master.aic$IC <- "AIC"
            master.aic$select <-  unlist(id.aic_loc)
            
            master.cluster.loc <- rbind.data.frame(master.bic, master.aic)
            master.cluster.loc$method <- "LOC"
            master.cluster.loc$simID <- rep(1:nsim, times=nsim)
            #COMBINE
            master.cluster <- rbind.data.frame(master.cluster.pc, master.cluster.loc)
            master.cluster$risk <- risk
            master.cluster$exp <- ecount
            master.cluster$radius <- rad
            ###########################################################################################
            ##################################################
            ##################################################
            #OUTPUT: CELLS
            ##################################################
            ##################################################
            #PC
            ##################################################
            master.bic<- cleanlist.cells(outcells.pc.bic, nsim=nsim, cellsix = cellsix)
            master.bic$IC <- "BIC"
            master.bic$select<-  unlist(id.bic_pc)
            
            master.aic<- cleanlist.cells(outcells.pc.aic, nsim=nsim, cellsix = cellsix)
            master.aic$IC <- "AIC"
            master.aic$select <-  unlist(id.aic_pc)
            
            master.cells.pc <- rbind.data.frame(master.bic, master.aic)
            master.cells.pc$method <- "PC"
            master.cells.pc$simID <- rep(1:nsim, each=(nsim*length(cellsix)))
            ##################################################
            #LOC
            ##################################################
            master.bic<- cleanlist.cells(outcells.loc.bic, nsim=nsim, cellsix = cellsix)
            master.bic$IC <- "BIC"
            master.bic$select<-  unlist(id.bic_loc)
            
            master.aic<- cleanlist.cells(outcells.loc.aic, nsim=nsim, cellsix = cellsix)
            master.aic$IC <- "AIC"
            master.aic$select <-  unlist(id.aic_loc)
            
            master.cells.loc <- rbind.data.frame(master.bic, master.aic)
            master.cells.loc$method <- "LOC"
            master.cells.loc$simID <- rep(1:nsim, each=(nsim*length(cellsix)))
            #COMBINE
            master.cells <- rbind.data.frame(master.cells.pc, master.cells.loc)
            master.cells$risk <- risk
            master.cells$exp <- ecount
            master.cells$radius <- rad
                
        }
    }
    
}
#write.csv(masterout, file="masterout_nsim50_cellrisk.csv")

write.csv(master.cluster, file="mastercluster_nsim100.csv")
write.csv(master.cells, file="mastercells_nsim100.csv")
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
