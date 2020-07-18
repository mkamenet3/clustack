#Calculating Clustack Bounds
#Simulation Study
#M.Kamenetsky
#v1: 2020-06-30

rm(list=ls())
library(dplyr)
library(tidyr)
library(clusso)
library(ggplot2)
library(forcats)



####################################################################
set.seed(20200630)
####################################################################
#library(clusso)
library(MASS)
library(clusso)
source("R/clustack.R")
load("../clustack/data/japanbreastcancer.RData")
cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath


#0) Setup
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]


#set global
rMax <- 20 
Time <- 5
#create set of potential clusters based on distances
potentialclusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparsematrix <- spacetimeMat(potentialclusters, numCenters, Time) 
#set maxclust
maxclust <- 15

##############################
#ARGS
##############################
nsim <- 2#10
risk = 2
cent = 150
rad = 11
tim <- c(1:5)
theta <- 60#args[1]

##############################
master <- NULL

for(s in 1:nsim){
    print(paste0("SimID: ",s))
    simid <- s
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
        #YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
        YSIM <- rpois(length(E1), lambda = E1)
    } else {
        print("OverD: NB")
        YSIM <- MASS::rnegbin(E1, theta = theta)
    }
    Ex <- unlist(scale_sim(list(YSIM), init, 1, Time))
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #SET-UP
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by location for each sim (apply)
    if (is.infinite(theta)){
        overdisp.est <- NULL
    } else {
        offset_reg <- glm(YSIM ~ as.factor(rep(c("1","2","3","4","5"),each=length(Ex)/Time)) + offset(log(Ex)),
                                                     family=quasipoisson)
        overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = TRUE)

    }
    outEx <- matrix(Ex, nrow=1)%*%sparsematrix
    outYx <- matrix(YSIM, nrow=1)%*%sparsematrix
    lambdahat <- YSIM/Ex
    
    Lambda <- as.vector(lambdahat)*sparsematrix #big Lambda matrix
    Lambda_dense <- as.matrix(Lambda)
    Lambda_dense[Lambda_dense == 0] <- 1
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY Potential Cluster
    #####################################################################################
    #####################################################################################
    #####################################################################################
    sim_superclust_pc<- detectclusters(sparsematrix, Ex, YSIM,
                                       numCenters, Time, maxclust,
                                       bylocation = FALSE, model="poisson",
                                       overdisp.est = overdisp.est)
    #BIC
    clusterlocs.bic <- clusterlocs_ident(sim_superclust_pc, selection="selection.bic", sparsematrix)
    modelthetas.bic <- extract_thetai(clusterlocs.bic$param_ix, clusterlocs.bic$w_q, Lambda_dense)
    #calculate null and model curves
    bounds.bic.log <- calcbounds(Yx, Ex,
                             modelthetas.bic$thetai, 
                             modelthetas.bic$thetaa, 
                             clusterlocs.bic$param_ix, 
                             clusterlocs.bic$w_q,
                             sparsematrix, logscale = TRUE,
                             overdisp.est=overdisp.est)
    if(sim_superclust_pc$selection.bic!=sim_superclust_pc$selection.aic){
        #AIC
        clusterlocs.aic <- clusterlocs_ident(sim_superclust_pc, selection="selection.aic", sparsematrix)
        modelthetas.aic <- extract_thetai(clusterlocs.aic$param_ix, clusterlocs.aic$w_q, Lambda_dense)
        #calculate null and model curves
        bounds.aic.log <- calcbounds(Yx, Ex,
                                 modelthetas.aic$thetai, 
                                 modelthetas.aic$thetaa, 
                                 clusterlocs.aic$param_ix, 
                                 clusterlocs.aic$w_q,
                                 sparsematrix, logscale = TRUE,
                                 overdisp.est=overdisp.est) 
    } else {
        print("BIC and AIC select the same path.")
        bounds.aic.log <- bounds.bic.log
        modelthetas.aic <- modelthetas.bic
    }
    ########################
    #Calculate MATA-Bounds
    ########################
    #BIC
    if(!is.null(overdisp.est)) {
        se.thetai.bic <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))*sqrt(overdisp.est)
    } else{
        se.thetai.bic <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))   
    }
    modelthetas_thetai.bic.log <- log(modelthetas.bic$thetai)
    mata_st.bic <- sapply(1:1040, function(x) mata_solvebounds(modelthetas_thetai.bic.log[,x],
                                                                 se.thetaii = se.thetai.bic[,x],
                                                                 w_q =  clusterlocs.bic$w_q,
                                                                 alpha = 0.025,
                                                                 tol=1e-8))
    mata_st.bic <- t(mata_st.bic)
    if(sim_superclust_pc$selection.bic!=sim_superclust_pc$selection.aic){
        #AIC
        if(!is.null(overdisp.est)) {
            se.thetai.aic <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.aic$param_ix]))*sqrt(overdisp.est)   
        } else {
            se.thetai.aic <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.aic$param_ix]))    
        }
        modelthetas_thetai.aic.log <- log(modelthetas.aic$thetai)
        mata_st.aic <- sapply(1:1040, function(x) mata_solvebounds(modelthetas_thetai.aic.log[,x],
                                                                     se.thetaii = se.thetai.aic[,x],
                                                                     w_q =  clusterlocs.aic$w_q,
                                                                     alpha = 0.025,
                                                                     tol=1e-8))
        mata_st.aic <- t(mata_st.aic)
        
    } else {
        print("BIC and AIC select the same path.")
        mata_st.aic <- mata_st.bic
    }
    #Create simrow
     row <- cbind(simID = rep(s,1040),
                 thetaa.bic= as.vector(modelthetas.bic$thetaa),
                 thetaa.aic = as.vector(modelthetas.aic$thetaa),
                 adjusted.bic.log.LB = as.vector(bounds.bic.log[[1]]),
                 adjusted.bic.log.UB = as.vector(bounds.bic.log[[2]]),
                 adjusted.aic.log.LB = as.vector(bounds.aic.log[[1]]),
                 adjusted.aic.log.UB = as.vector(bounds.aic.log[[2]]),
                 mata.bic.LB = mata_st.bic[,1],
                 mata.bic.UB = mata_st.bic[,2],
                 mata.aic.LB = mata_st.aic[,1],
                 mata.aic.UB = mata_st.aic[,2])
    master <- rbind(master, row)
}

write.csv(master, file="clustackbounds_sim_bypc.csv")


##########################################################################################
#POST ANALYSIS
##########################################################################################
master2 <- as.data.frame(master)
master2$locid <- rep(1:1040, times=nsim)
master2$locidF <- as.factor(rep(1:1040, times=2))
clusterix <- which(rr==risk)
master2 <- master2 %>%
    mutate(clusterix = ifelse(locid %in% clusterix,1,0)) 
master2 <- master2 %>%
    mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB & clusterix==1), 1, 0),
           coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB & clusterix==1), 1, 0),
           coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB & clusterix==1), 1, 0),
           coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB & clusterix==1), 1, 0))

################################
#Coverage probability
################################



ggplot(data=master2,aes(x=locid, y=thetaa.bic, group=factor(simID), color=factor(simID))) +
    theme_bw() +
    geom_ribbon(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB, fill=factor(simID)), alpha=0.1) +
    geom_line()

master2 %>%
    arrange(thetaa.bic) %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic, group=factor(simID), color=factor(simID))) +
    theme_bw() +
    geom_ribbon(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB, fill=factor(simID)), alpha=0.1) +
    geom_line()

master2 %>%
    dplyr::filter(coverage.bic==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=locidF, y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB,color=factor(simID)), alpha=0.9, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("nsim=10, estimates & bounds where RR=2 within bounds (adequate coverage)")

master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.bic.log.LB, ymax = adjusted.bic.log.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (Buckland, BIC)")
master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.aic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.aic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.aic.log.LB, ymax = adjusted.aic.log.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (Buckland, AIC)")


master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.bic.LB, ymax = mata.bic.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (MATA, BIC)")
master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.aic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.aic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.aic.LB, ymax = mata.aic.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (MATA, AIC)")

##########################################
#coverage probability
##########################################
#make a quick and dirty plot
dat <- cbind.data.frame(thetaa=as.vector(modelthetas.bic$thetaa), UBa = as.vector(bounds.bic$ma_adjusted[[2]]), 
                         LBa=as.vector(bounds.bic$ma_adjusted[[1]]))
dat <- dat %>%arrange(thetaa) %>%mutate(id = 1:nrow(.))
ggplot(data=dat) +
    geom_line(aes(x=id, y=thetaa), size=2) +
    theme_bw() +
    geom_line(aes(x=id, y=LBa), col="blue", size=1.5, alpha=0.5) +
    geom_line(aes(x=id, y=UBa), col="blue", size=1.5, alpha=0.5) +
    geom_hline(yintercept = 1, col="red", linetype="dashed", size=1) +
    geom_vline(xintercept = 823, col="red", linetype="dashed", size=1)


# #https://github.com/jlaake/RMark/blob/master/RMark/R/mata.wald.r
# thetai <- log(modelthetas.bic$thetai)
# thetaa <- log(modelthetas.bic$thetaa)
# se.thetai <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))
# w_q <- clusterlocs.bic$w_q
# thetaaa <- log(as.vector(thetaa)[827])
# thetaii <- log(thetai[,827])
# se.thetaii <- se.thetai[,827]
# 
# ########################################################
# #MATA-Intervals
# ########################################################
# thetai <- as.matrix(thetai)
# diff1 = apply(thetai,1, function(x) (as.vector(thetaa)-x))
# zquant <- sapply(1:ncol(diff1), function(x) diff1[,x]/se.thetai[x,])
# test1 <- apply(zquant, 2, function(x) pnorm(zquant))

######################################################
mataL <- uniroot(f=mata_tailareazscore, interval=c(-3, 3),
                 thetaii= thetaii,
                 se.thetaii=se.thetaii,
                 w_q=w_q, alpha=0.025, tol=1e-8)$root

mataU <- uniroot(f=mata_tailareazscore, interval=c(-3, 3),
                 thetaii= thetaii,
                 se.thetaii=se.thetaii,
                 w_q=w_q, alpha=1-0.025, tol=1e-8)$root
exp(cbind(mataL, mataU))



mata_solvebounds(thetaii = thetaii, se.thetaii = se.thetaii, w_q, alpha = 0.025, tol=1e-8)
mata_solvebounds(thetaii = thetaii, se.thetaii = se.thetaii, w_q, alpha = 0.025, tol=1e-8)
#do this for all ST locations

mata_st <- sapply(1:1040, function(x) mata_solvebounds(thetai[,x],
                                                    se.thetaii = se.thetai[,x],
                                                    w_q = w_q,
                                                    alpha = 0.025,
                                                    tol=1e-8))
mata_test <- sapply(1:10, function(x) mata_solvebounds(thetai[,x],
                                                       se.thetaii = se.thetai[,x],
                                                       w_q = w_q,
                                                       alpha = 0.025,
                                                       tol=1e-8))
mata_stt <- t(mata_st)
mata_df <- as.data.frame(mata_stt)
names(mata_df) <- c("mataLB", "mataUB")
mata_df$locid <- 1:1040

mata_df <- mata_df %>%
    mutate(clusterix = ifelse(locid %in% clusterix,1,0)) 
mata_df <- mata_df %>%
    mutate(coverage.bic = ifelse((risk < mataUB & risk > mataLB & clusterix==1), 1, 0))






