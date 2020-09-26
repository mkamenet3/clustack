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
            #Ex <- rep(ecounts,1040)
            #YSIM <- as.vector(rr)*Ex
            YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
            #init <- clusso::setVectors(japanbreastcancer$period, japanbreastcancer$expdeath, japanbreastcancer$death,covars=NULL, Time=Time)
            # std <- lapply(1:nsim,
            #               function(i) sapply(1:Time,
            #                                  function(j) (matrix(E0, ncol = Time)[, j]) * (sum(matrix(YSIM[[i]], ncol = Time)[,j])/sum(matrix(E0, ncol = Time)[, j]))))
            # Ex <- lapply(1:nsim, function(i) as.vector(std[[i]]))
            
            # E0 <- rep(ecounts, 1040)
            # E1 <- rr*E0
            # #Ex <- rep(ecounts,1040)
            # #YSIM <- as.vector(rr)*Ex
            # YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
            # #init <- clusso::setVectors(japanbreastcancer$period, japanbreastcancer$expdeath, japanbreastcancer$death,covars=NULL, Time=Time)
            # # std <- lapply(1:nsim, 
            # #               function(i) sapply(1:Time, 
            # #                                  function(j) (matrix(E0, ncol = Time)[, j]) * (sum(matrix(YSIM[[i]], ncol = Time)[,j])/sum(matrix(E0, ncol = Time)[, j]))))
            Ex <- lapply(1:nsim, function(i) E0)#lapply(1:nsim, function(i) as.vector(std[[i]]))
            outExp <- lapply(1:nsim, function(i) t(sparseMAT)%*%Ex[[i]])
            outObs <- lapply(1:nsim, function(i) t(sparseMAT)%*%YSIM[[i]])
            ix <- which(rr==risk) #ix ST cells wi
            
            
            sim_superclust_pc_large <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
                                                                                 numCenters, Time, maxclust,
                                                                                 bylocation = FALSE, model="poisson",
                                                                                 overdisp.est = overdisp.est))
            
            
            # ##############
            # #maxlocs start
            # sim_superclust_pc_large <- lapply(1:nsim, function(i) detectclusters(sparseMAT, Ex[[i]], YSIM[[i]],
            #                                                                      numCenters, Time, maxclust,
            #                                                                      bylocation = TRUE, model="poisson", 
            #                                                                      overdisp.est = overdisp.est))
            #find all PCs that overlap max location
            wslarge <- lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$wtMAT[,sim_superclust_pc_large[[i]]$selection.bic])
            clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc_large[[i]]$Lambda_dense), 
                                                                     function(k) unique(sim_superclust_pc_large[[i]]$Lambda_dense[k,]))) 
            ###
            # #id pcs that overlap maxloc
            # ix_locs<- rep(0,1040);ix_locs[ix] <-1
            # ix_locs_mat <- matrix(aa, nrow = 1)%*%sparseMAT
            # ix_locs_mods <- which(ix_locs_mat!=0)
            # #ix_locs_mods <- ifelse(test!=0,1,0)
            # 
            # cluster_thetaa_locs <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]][ix_locs_mods]*wslarge[[i]][ix_locs_mods]))
            
            ###
            clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
            clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
            clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
            cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
            
            
            
            clusterRRlarge <- lapply(1:nsim, 
                                     function(i) unique(sim_superclust_pc_large[[i]]$Lambda_dense[sim_superclust_pc_large[[i]]$maxpcs[sim_superclust_pc_large[[i]]$selection.bic],])[2])
            
            
            #maxlocs end
            
            # 
            ##################################################
            ##################################################
            #NON-MA VARIANCE
            ##################################################
            ##################################################
            
            # clusterRRlarge <- lapply(1:nsim, 
            #                          function(i) unique(sim_superclust_pc_large[[i]]$Lambda_dense[sim_superclust_pc_large[[i]]$maxpcs[sim_superclust_pc_large[[i]]$selection.bic_orig],])[2])
            clusterRRlarge <- lapply(1:nsim, 
                                     function(i) unique(sim_superclust_pc_large[[i]]$Lambda_dense[sim_superclust_pc_large[[i]]$maxpcs[sim_superclust_pc_large[[i]]$selection.bic],])[2])
            se_clusterRRlarge <- lapply(1:nsim, function(i)sqrt(clusterRRlarge[[i]]/(Time*n)))
            nonma <- lapply(1:nsim, function(i) cbind(lb=clusterRRlarge[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                      clusterMA = clusterRRlarge[[i]],
                                                      ub=clusterRRlarge[[i]]+1.96*se_clusterRRlarge[[i]]))
            
            
            #n
            se_clusterRRlarge_asymp <- lapply(1:nsim, function(i) sqrt(clusterRRlarge[[i]]/(sum(YSIM[[i]][ix]))))
            nonma_asymp <- lapply(1:nsim, function(i) cbind(lbasymp=clusterRRlarge[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                            clusterMA = clusterRRlarge[[i]],
                                                            ubasymp=clusterRRlarge[[i]]+1.96*se_clusterRRlarge_asymp[[i]]))
            
            
            ##################################################
            ##################################################
            #Buckland 1997
            ##################################################
            ##################################################
            # wslarge <- lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$wtMAT[,sim_superclust_pc_large[[i]]$selection.bic])
            wslarge <- lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$wtMAT[,sim_superclust_pc_large[[i]]$selection.bic])
            clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc_large[[i]]$Lambda_dense), 
                                                                     function(k) unique(sim_superclust_pc_large[[i]]$Lambda_dense[k,]))) 
            # clusterRR_uniqlarge <- lapply(1:nsim, function(i) sapply(1:nrow(sim_superclust_pc_large[[i]]$Lambda_dense), 
            #                               function(k) unique(sim_superclust_pc_large[[i]]$Lambda_dense[k,]))) 
            
            # uniqRRslarge<- lapply(1:nsim, function(i) t(clusterRR_uniqlarge[[i]]))#as.matrix(do.call(rbind, clusterRR_uniqlarge), nrow=2, byrow=TRUE)
            # clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(t(clusterRR_uniqlarge[[i]])))
            
            ##########################
            #model-average clusterMA
            clusterRR_ilarge <- lapply(1:nsim, function(i) rep(NA, 66870))
            clusterRR_uniq_ilarge <- lapply(1:nsim, function(i) as.matrix(do.call(rbind, clusterRR_uniqlarge[[i]]), ncol=2))
            clusterRR_ilarge <- lapply(1:nsim, function(i) selectuniqRR(clusterRR_uniq_ilarge[[i]]))
            cluster_thetaa <- lapply(1:nsim, function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
          
            
            ##########################
            outbuck <- lapply(1:nsim, 
                              function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], thetaa =cluster_thetaa[[i]], 
                                                         w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL))
            outbuck
            
            
            
            ##################################################
            ##################################################
            #(un-adjusted) MAW1 (B&A pg. 164) = Buckland 1997
            ##################################################
            ##################################################
            
            
            outmaw1 <- lapply(1:nsim, 
                              function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                               w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL))
            outmaw1
            
            
            
            
            
            ##################################################
            ##################################################
            #MAW2 (B&A pg. 345)
            ##################################################
            ##################################################
            
            outmaw2 <- lapply(1:nsim, function(i) maw1(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                       w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL))
            outmaw2
            
            
            
            
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds (for non-normal data)
            ##################################################
            ##################################################
            
            
            # outmata <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = clusterRRlarge[[i]], 
            #w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL))
            # outmata
            outmata <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                             w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL, transform="none"))
            outmata
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds: SQRT TRANSFORMED
            ##################################################
            ##################################################
            
            # outmataT <- lapply(1:nsim, function(i) mataboundsT(thetai=clusterRR_ilarge[[i]], thetaa = clusterRRlarge[[i]], 
            #                                                    w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL))
            # outmataT
            
            outmataT <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                               w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL, transform="sqrt"))
            outmataT
            
            ##################################################
            ##################################################
            #Turek-Fletcher MATA Bounds: LOG TRANSFORMED
            ##################################################
            ##################################################
            
            # outmataTlog <- lapply(1:nsim, function(i) logmataboundsT(thetai=clusterRR_ilarge[[i]], thetaa = clusterRRlarge[[i]], 
            #                                                          w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL))
            # outmataTlog
            outmataTlog <- lapply(1:nsim, function(i) matabounds(thetai=clusterRR_ilarge[[i]], thetaa = cluster_thetaa[[i]], 
                                                                     w_q=wslarge[[i]], sparsematrix=t(sparseMAT), overdisp.est = NULL, transform="log"))
            outmataTlog
            
            
            # master_std <- cbind.data.frame(matrix(unlist(nonma), byrow=TRUE, ncol=3),
            #                     matrix(unlist(nonma_asymp), byrow=TRUE, ncol=3),
            #                     matrix(unlist(outmaw1), byrow=TRUE, ncol=3),
            #                     matrix(unlist(outmaw2), byrow=TRUE, ncol=3),
            #                     matrix(unlist(outmata), byrow=TRUE, ncol=3),
            #                     matrix(unlist(outmataT), byrow=TRUE, ncol=3),
            #                     matrix(unlist(outmataTlog), byrow=TRUE, ncol=3))
            # names(master_std) <- c("nonma.LB", "clusterMA", "nonma.UB",
            #                        "nonma_asymp.LB", "clusterMA", "nonma_asymp.UB",
            #                        "maw1.LB", "clusterMA", "maw1.UB",
            #                        "maw2.LB", "clusterMA", "maw2.UB",
            #                        "mata.LB", "clusterMA", "mata.UB",
            #                        "matasqrt.LB", "clusterMA", "matasqrt.UB",
            #                        "matalog.LB", "clusterMA", "matalog.UB")
            
            
            master <- cbind.data.frame(matrix(unlist(nonma), byrow=TRUE, ncol=3),
                                       matrix(unlist(nonma_asymp), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmaw1), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmaw2), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmata), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmataT), byrow=TRUE, ncol=3),
                                       matrix(unlist(outmataTlog), byrow=TRUE, ncol=3))
            master$risk <- risk
            master$exp <- ecount
            master$radius <- rad
            master$anyforced <- ifelse(any(unlist(lapply(1:nsim, function(i) sim_superclust_pc_large[[i]]$selection.bic_orig))==0),"yes","no")
            master$simID <- 1:nrow(master)
            
            
            names(master) <- c("nonma.LB", "clusterMA", "nonma.UB",
                               "nonma_asymp.LB", "clusterMA.1", "nonma_asymp.UB",
                               "maw1.LB", "clusterMA.2", "maw1.UB",
                               "maw2.LB", "clusterMA.3", "maw2.UB",
                               "mata.LB", "clusterMA.4", "mata.UB",
                               "matasqrt.LB", "clusterMA.5", "matasqrt.UB",
                               "matalog.LB", "clusterMA.6", "matalog.UB", 
                               "risk", "ecount", "rad", "anyforced","simID")
            #write.csv(master, file=paste0("force_nsim100_boundscompare","_risk_",risk, "_ecount_",ecount,"_rad_",rad,".csv"))
            masterout <- rbind(masterout, master)
            #write.csv(master, file="nsim10_std_boundscompare.csv")
            
            #write.csv(master, file="nsim10_nostd_boundscompare.csv")
            #write.csv(master, file="E500_nsim10_nostd_boundscompare.csv")
            #write.csv(master, file="E250_nsim10_nostd_boundscompare.csv")
            #write.csv(master, file="E100_nsim10_nostd_boundscompare.csv")
            #write.csv(master, file="E50_nsim10_nostd_boundscompare.csv")
            #write.csv(master, file="E10_nsim10_nostd_boundscompare.csv")
            #write.csv(master, file="E5_nsim10_nostd_boundscompare.csv")
        }
    }
    
}
write.csv(masterout, file="masterout_nsim50.csv")
# 
# # ####################################################################################################
# # #ANALYSIS
# ####################################################################################################
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(forcats)
# library(patchwork)
# 
# cleandat <- function(dat){
#     clean <- dat %>%
#         select(-c(clusterMA.1, clusterMA.2, clusterMA.3, clusterMA.4, clusterMA.5, clusterMA.6)) %>%
#         mutate(simiD = X) %>%
#         pivot_longer(cols = ends_with("B"),
#                      names_to="mtd",
#                      values_to="Bound") %>%
#         separate(mtd, sep="([.])", into=c("method", "LBUB"), remove=FALSE) %>%
#         select(-mtd) %>%
#         pivot_wider(names_from="LBUB",
#                     values_from="Bound") %>%
#         filter(method!="nonma_asymp") %>%
#         mutate(method = as.factor(method)) %>%
#         mutate(method= fct_relevel(method, "nonma", "maw1", "maw2", "mata", "matasqrt", "matalog"))
# 
#     return(clean)
# }
# 
# #load all csvs
# temp <- list.files(path="data",pattern="*.csv", full.names = TRUE)
# myfiles <- lapply(temp, read.csv)
# 
# myfiles.clean <- lapply(myfiles, cleandat)
# 
# compare <- do.call(rbind.data.frame, myfiles.clean)
# compare$RR <- rep(c(1.1,1.5,2), each=1260)
# compare$radius <- rep(rep(c(18,11,9), each=60), times=21)
# compare$exp <- rep(c(10,100,1000,250,5,50,500), times=180*3)
# compare$expF <- factor(compare$exp, levels=c("5", "10", "50", "100", "250","500", "1000"))
# compare$width <- compare$UB-compare$LB
# # 
# p9 <-ggplot(subset(compare, radius==9)) +
#     geom_line(aes(x=exp, y=clusterMA, group=simiD)) +
#     geom_line(aes(x=exp, y=LB, group=simiD, color=as.factor(RR))) +
#     geom_line(aes(x=exp, y=UB, group=simiD, color=as.factor(RR))) +
#     geom_line(aes(x=exp,y= RR), size=1, color="red", linetype="dashed") +
#     geom_hline(yintercept = 1, color="grey30", linetype="dashed", size=1) +
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts") +
#     ggtitle("9km Radius ")
# 
# p11 <-ggplot(subset(compare, radius==11)) +
#     geom_line(aes(x=exp, y=clusterMA, group=simiD)) +
#     geom_line(aes(x=exp, y=LB, group=simiD, color=as.factor(RR))) +
#     geom_line(aes(x=exp, y=UB, group=simiD, color=as.factor(RR))) +
#     geom_line(aes(x=exp,y= RR), size=1, color="red", linetype="dashed") +
#     geom_hline(yintercept = 1, color="grey30", linetype="dashed", size=1) +
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts")+
#     ggtitle("11km Radius ")
# 
# 
# p18 <-ggplot(subset(compare, radius==18)) +
#     geom_line(aes(x=exp, y=clusterMA, group=simiD)) +
#     geom_line(aes(x=exp, y=LB, group=simiD, color=as.factor(RR))) +
#     geom_line(aes(x=exp, y=UB, group=simiD, color=as.factor(RR))) +
#     geom_line(aes(x=exp,y= RR), size=1, color="red", linetype="dashed") +
#     geom_hline(yintercept = 1, color="grey30", linetype="dashed", size=1) +
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts")+
#     ggtitle("18km Radius ")
# 
# 
# p9 + p11 + p18
# 
# 
# #CI widths
# ggplot(compare, aes(as.factor(exp),width)) +
#     geom_boxplot(aes( fill=as.factor(RR)), outlier.shape = NA)+
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     xlab("Expected Counts")
# 
# p9w <-ggplot(subset(compare, radius==9), aes(as.factor(exp),width)) +
#     geom_boxplot(aes( fill=as.factor(RR)), outlier.shape = NA)+
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     xlab("Expected Counts")
# 
# p11w <-ggplot(subset(compare, radius==11), aes(as.factor(exp),width)) +
#     geom_boxplot(aes( fill=as.factor(RR)), outlier.shape = NA)+
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     xlab("Expected Counts")
# p18w <-ggplot(subset(compare,radius==18), aes(as.factor(exp),width)) +
#     geom_boxplot(aes( fill=as.factor(RR)), outlier.shape = NA)+
#     facet_grid(method~RR) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     xlab("Expected Counts")
# p9w + p11w + p18w
# 
# 
# ggplot(compare, aes(width)) +
#     geom_histogram(aes( fill=as.factor(RR)))+
#     #geom_line(aes(x=exp, y=width, group=simiD)) +
#     #geom_line(aes(x=exp, y=LB, group=simiD, color=as.factor(RR))) +
#     #geom_line(aes(x=exp, y=UB, group=simiD, color=as.factor(RR))) +
#     #geom_line(aes(x=exp,y= RR), size=1, color="red", linetype="dashed") +
#     #geom_hline(yintercept = 1, color="grey30", linetype="dashed", size=1) +
#     facet_grid(exp~method, scales="free_y") +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     # geom_hline(aes(yintercept = clusterMA), color="black") +
#     # geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
#     # geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
#     #geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
#     # facet_grid(exp~method) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))#+
#     #scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     #xlab("Expected Counts")
# 
# 
# 
# combo_mean <- compare %>%
#     group_by(method, exp, RR, radius) %>%
#     summarise(clusterMA_mean = mean(clusterMA),
#            LB_mean = mean(LB, na.rm=TRUE),
#            UB_mean = mean(UB, na.rm = TRUE))
# 
# 
# ggplot(subset(combo_mean, radius==11)) +
#     geom_line(aes(x=exp, y=clusterMA_mean), size=1.5) +
#     geom_line(aes(x=exp, y=LB_mean,color=as.factor(RR)), size=1.5) +
#     geom_line(aes(x=exp, y=UB_mean,  color=as.factor(RR)), size=1.5) +
#     geom_line(aes(x=exp,y= RR), size=1, color="red", linetype="dashed") +
#     geom_hline(yintercept = 1, color="grey30", linetype="dashed", size=1) +
#     facet_grid(RR~method) +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     # geom_hline(aes(yintercept = clusterMA), color="black") +
#     # geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
#     # geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
#     #geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
#     # facet_grid(exp~method) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts")
# 
# 
# 
# 
# 
# #calculate coverage probability
# 
# compare_cov <- compare %>%
#     group_by(method, exp, RR, radius) %>%
#     drop_na(LB, UB) %>%
#     # mutate(UB = ifelse(is.na(UB),0,UB),
#     #        LB = ifelse(is.na(LB),0, LB)) %>%
#     # mutate(cov_ind = ifelse(RR>LB & RR<UB,1,0)) %>%
#     summarize(coverage  = mean(cov_ind))
# 
# p9c <- ggplot(subset(compare_cov, radius==9)) +
#     geom_line(aes(x=exp, y=coverage, color=as.factor(RR)), size=1.5) +
#     facet_grid(method~RR) +
#     theme_bw() +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts") +
#     #ggtitle("NAs = 0")
#     ggtitle("9km Radius")
# 
# 
# p11c <- ggplot(subset(compare_cov, radius==11)) +
#     geom_line(aes(x=exp, y=coverage, color=as.factor(RR)), size=1.5) +
#     facet_grid(method~RR) +
#     theme_bw() +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts") +
#     #ggtitle("NAs = 0")
#     ggtitle("11km Radius")
# 
# 
# p18c <- ggplot(subset(compare_cov, radius==18)) +
#     geom_line(aes(x=exp, y=coverage, color=as.factor(RR)), size=1.5) +
#     facet_grid(method~RR) +
#     theme_bw() +
#     scale_color_manual(values=c("goldenrod", "green4", "dodgerblue")) +
#     theme(axis.text.x = element_text(angle = 45))+
#     scale_x_continuous(breaks=c(5,10,50,100,250,500,1000)) +
#     xlab("Expected Counts") +
#     #ggtitle("NAs = 0")
#     ggtitle("18km Radius")
# 
# p9c + p11c + p18c
# 
# 
# # 
# # #
# # #dat <- read.csv("E5_nsim10_nostd_boundscompare.csv")
# # # dat <- read.csv("E50_nsim10_nostd_boundscompare.csv")
# # #dat <- read.csv("E10_nsim10_nostd_boundscompare.csv")
# # #dat <- read.csv("E100_nsim10_nostd_boundscompare.csv")
# # #dat <- read.csv("E250_nsim10_nostd_boundscompare.csv")
# # #dat <- read.csv("E500_nsim10_nostd_boundscompare.csv")
# # # #dat <- read.csv("E1000_nsim10_nostd_boundscompare.csv")
# # # dat <- read.csv("nsim10_std_boundscompare.csv")
# # 
# # dat5 <- read.csv("E5_nsim10_nostd_boundscompare.csv")
# # dat50 <- read.csv("E50_nsim10_nostd_boundscompare.csv")
# # dat10 <- read.csv("E10_nsim10_nostd_boundscompare.csv")
# # dat100 <- read.csv("E100_nsim10_nostd_boundscompare.csv")
# # dat250 <- read.csv("E250_nsim10_nostd_boundscompare.csv")
# # dat500 <- read.csv("E500_nsim10_nostd_boundscompare.csv")
# # 
# # dat5c <- cleandat(dat5)
# # dat10c <- cleandat(dat10)
# # dat50c <- cleandat(dat50)
# # dat100c <- cleandat(dat100)
# # dat250c <- cleandat(dat250)
# # dat500c <- cleandat(dat500)
# # #dat <- read.csv("E1000_nsim10_nostd_boundscompare.csv")
# # # dat <- read.csv("nsim10_std_boundscompare.csv")
# # 
# # # 
# # # dat1 <- dat %>%
# # #     select(-c(clusterMA.1, clusterMA.2, clusterMA.3, clusterMA.4, clusterMA.5, clusterMA.6)) %>%
# # #     mutate(simiD = X) %>%
# # #     pivot_longer(cols = ends_with("B"),
# # #                  names_to="mtd",
# # #                  values_to="Bound") %>%
# # #     separate(mtd, sep="([.])", into=c("method", "LBUB"), remove=FALSE) %>%
# # #     select(-mtd) %>%
# # #     pivot_wider(names_from="LBUB",
# # #                 values_from="Bound") %>%
# # #     filter(method!="nonma_asymp") %>%
# # #     mutate(method = as.factor(method)) %>%
# # #     mutate(method= fct_relevel(method, "nonma", "maw1", "maw2", "mata", "matasqrt", "matalog"))
# # 
# # 
# # ggplot(dat1) +
# #     geom_hline(aes(yintercept = clusterMA), color="black") +
# #     geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
# #     geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
# #     geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
# #     facet_grid(.~method) +
# #     theme_bw() +
# #     ylim(1.5,2.5)
# # 
# # combo <- rbind.data.frame(dat5c, dat10c, dat50c,dat100c, dat250c, dat500c)
# # combo$exp <- rep(c(5,10,50,100,250,500), each=60)
# # combo$expF <- as.factor(combo$exp)
# # combo$expF <- factor(combo$expF, levels=c("5", "10", "50", "100", "250","500"))
# # 
# # 
# # ggplot(combo) +
# #     geom_hline(aes(yintercept = clusterMA), color="black") +
# #     geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
# #     geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
# #     geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
# #     facet_grid(exp~method) +
# #     theme_bw() +
# #     ylim(1.5,2.5)
# # 
# # 
# # 
# # ggplot(combo) +
# #     geom_line(aes(x=exp, y=clusterMA, group=simiD)) +
# #     facet_grid(.~method) +
# #     # geom_hline(aes(yintercept = clusterMA), color="black") +
# #     # geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
# #     # geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
# #     geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
# #     # facet_grid(exp~method) +
# #     theme_bw() +
# #     ylim(1.75,2.25)
# # 
# # 
# # ggplot(combo) +
# #     geom_line(aes(x=exp, y=clusterMA, group=simiD)) +
# #     geom_line(aes(x=exp, y=LB, group=simiD), color="blue") +
# #     geom_line(aes(x=exp, y=UB, group=simiD), color="blue") +
# #     facet_grid(.~method) +
# #     # geom_hline(aes(yintercept = clusterMA), color="black") +
# #     # geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
# #     # geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
# #     geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
# #     # facet_grid(exp~method) +
# #     theme_bw() +
# #     ylim(1.75,2.25)
# # 
# # 
# # combo_mean <- combo %>%
# #     group_by(method, exp) %>%
# #     summarise(clusterMA_mean = mean(clusterMA),
# #            LB_mean = mean(LB),
# #            UB_mean = mean(UB))
# # 
# # 
# # 
# # ggplot(combo_mean) +
# #     geom_line(aes(x=exp, y=clusterMA_mean)) +
# #     geom_line(aes(x=exp, y=LB_mean), color="blue") +
# #     geom_line(aes(x=exp, y=UB_mean), color="blue") +
# #     facet_grid(.~method) +
# #     # geom_hline(aes(yintercept = clusterMA), color="black") +
# #     # geom_hline(aes(yintercept = LB, color=as.factor(simiD))) +
# #     # geom_hline(aes(yintercept = UB, color=as.factor(simiD))) +
# #     geom_hline(yintercept = 2, size=1, color="red", linetype="dashed") +
# #     # facet_grid(exp~method) +
# #     theme_bw() +
# #     ylim(1.7,2.25)
# # 
# # 
