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
#Arguments passed will only control theta parameter
##R CMD BATCH "--args arg1 arg2" myscript.R -> 
args<- commandArgs(TRUE)
theta <- args[1]
nsim <- args[2]
#Fixed Params
cent <- 150
tim <- c(1:5)

#Sim-Over Params
risks <- c(1.1, 1.5, 2.0)
radii <- c(9, 11, 18)

#tester
# risk = 2
# cent = 150
# rad = 11
# tim <- c(1:5)


##############################
master <- NULL

for(rad in radii){
    for(risk in risks){
        for(s in 1:nsim){
            print(paste0("SimID: ",s,
                         "Risk: ", risk,
                         "Radius: ", rad))
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
                         risk = rep(risk, 1040),
                         radius = rep(rad, 1040),
                         theta = rep(theta, 1040),
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
    }
}



write.csv(master, file=paste0("clustackbounds_sim_bypc", "_theta", theta,".csv"))


