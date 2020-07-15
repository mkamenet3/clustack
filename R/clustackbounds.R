#Calculating Clustack Bounds
#Simulation Study
#M.Kamenetsky
#v1: 2020-06-30

rm(list=ls())
library(dplyr)
library(tidyr)
library(clusso)
library(ggplot2)
####################################################################
#HELPER FUNCTIONS
####################################################################
dpoisson_theta <- function (y, lambda, E0, ix) {
    #print(str(ix))
    ell <- rep(NA, length(lambda))
    for(i in 1:length(lambda)){
        x <- rep(1,1040)
        x[ix] <- lambda[i]
        ell[i] <- sum(y * log(x * E0 ) - (x * E0))
    }
    return(ell)
}

#identify clusters and their locations
#' @param res result of model-averaging
#' @param selection Character string indicating which information criterion to use. Options include \code{"selection.aic"}, \code{"selection.aicc"}, or \code{"selection.bic"}
#' @param  sparsematrix sparse design matrix where rows are locations and columns are potential clusters
clusterlocs_ident <- function(res, selection, sparsematrix){
    maxpc <- res$maxpcs[unlist(res[selection])]
    pcmax <- rep(0,dim(sparsematrix)[2]); pcmax[maxpc] <-1; pcmax <- matrix(pcmax,ncol=1)
    pcpcmax <- t(sparsematrix%*%pcmax)%*%sparsematrix; pcpcmax <- ifelse(pcpcmax!=0,1,0)
    #Q = which models actually have parameter of interest
    Qmods_pc_ix<- which(pcpcmax!=0)
    param_ix <- which(pcpcmax!=0)
    w2_q <- res[["wtMAT"]][,unlist(res[selection])][param_ix]#w2[param_ix]
    #pclocs should identify all  models that contrain cluster65193
    pclocs <- as.vector(matrix(pcpcmax, nrow=1)%*%t(sparsematrix))
    pclocs <- ifelse(pclocs!=0,1,0)
    return(list(param_ix = param_ix,
                pclocs = pclocs,
                w_q = w2_q))
}

#Extract all estimated theta_i for each model that overlaps max identified cluster
#' @param param_ix Indices of potential clusters (models) that overlap max identified cluster
#' @param w_q Weights for each of the identified potential clusters (models) that overlap the max identified cluster
#' @param Lambda_dense Large matrix of relative risks for each potential cluster
extract_thetai <- function(param_ix, w_q, Lambda_dense){
    #create a sparematrix
    M <- Matrix(0, nrow=dim(Lambda_dense)[2], ncol=length(param_ix), sparse=TRUE)
    for (i in 1:length(param_ix)){
        M[,i][param_ix[i]] <- 1    
    }
    thetas <- t(M)%*%t(Lambda_dense)
    thetaa <- matrix(w_q, nrow=1)%*%thetas
    return(list(thetai = thetas,
                thetaa = thetaa))
}



#calculate profile-likelihood confidence bounds and adjust for model-averaaging
#' @param null Null model likelihood
#' @param proflik Vector of profiled likelihoods for set of thetas
calcbounds <- function(Yx, Ex, thetai,thetaa, param_ix, w_q,sparsematrix, logscale) {
    #################################
    #Added 7-8-20
    #thetai <- modelthetas.bic$thetai
    #thetaa <- modelthetas.bic$thetaa@x
    #################################
    if(logscale==FALSE){
        variance <- vector(mode = "list", length = dim(thetai)[1])
        #model averaged calc
        varthetai <- apply(thetai, 2, function(x) x/outEx@x[clusterlocs.bic$param_ix])
        withintheta <- apply(thetai,1, function(x) (x-as.vector(thetaa))^2)
        var <- sapply(1:ncol(varthetas), function(k) sqrt(varthetai[k,1]+withintheta[,k]))
        varthetas_w <- var%*%matrix(w_q, ncol=1)#crossprod(matrix(w_q, ncol=1,yy2 )) #apply(yy2,2, function(x) w_q*x)
        var_thetaa <- as.vector(varthetas_w)
        UBa = exp(log(thetaa) + 1.96*sqrt(var_thetaa))
        LBa = exp(log(thetaa) - 1.96*sqrt(var_thetaa))
    } else {
        #TODO
    }

    return(ma_adjusted = c(LBa, UBa))
}


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

#test
nsim <-2
theta = 100
risk = 2
cent = 150
rad = 11
tim <- c(1:5)


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
    #if (is.infinite(theta)){
     #   print("Inf theta: pure poisson")
        #YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
        YSIM <- rpois(length(E1), lambda = E1)
    #} else {
     #   print("OverD: NB")
    #    YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
    #}
    Ex <- unlist(scale_sim(list(YSIM), init, 1, Time))
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #CLUSTER DETECTION BY LOCATION
    #####################################################################################
    #####################################################################################
    #####################################################################################
    #run superclust by location for each sim (apply)
    # if (is.infinite(theta)){
    #     overdisp.est <- NULL
    # } else {
    #     offset_reg <- lapply(1:nsim, function(i) glm(YSIM[[i]] ~ as.factor(rep(c("1","2","3","4","5"), 
    #                                                                            each=length(Ex[[i]])/Time)) + offset(log(Ex[[i]])),
    #                                                  family=quasipoisson))
    #     overdisp.est <- overdisp(offset_reg, sim = TRUE, overdispfloor = TRUE)
    #     
    # }
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
                                       overdisp.est = NULL)
    # thetas <- seq(0.1,3, length=1000)
    # null <- dpoisson_theta(Yx, 1, Ex, ix=1:dim(sparsematrix)[1])
    
    #BIC
    clusterlocs.bic <- clusterlocs_ident(sim_superclust_pc, selection="selection.bic", sparsematrix)
    modelthetas.bic <- extract_thetai(clusterlocs.bic$param_ix, clusterlocs.bic$w_q, Lambda_dense)
    #calculate null and model curves
    #out.bic <- dpoisson_theta(Yx, thetas, Ex, ix = which(clusterlocs.bic$pclocs!=0))
    bounds.bic <- calcbounds(Yx, Ex,
                             modelthetas.bic$thetai, 
                             modelthetas.bic$thetaa, 
                             clusterlocs.bic$param_ix, 
                             clusterlocs.bic$w_q,
                             sparsematrix) 
    if(sim_superclust_pc$selection.bic!=sim_superclust_pc$selection.aic){
        #AIC
        clusterlocs.aic <- clusterlocs_ident(sim_superclust_pc, selection="selection.aic", sparsematrix)
        modelthetas.aic <- extract_thetai(clusterlocs.aic$param_ix, clusterlocs.aic$w_q, Lambda_dense)
        #calculate null and model curves
        #out.aic <- dpoisson_theta(Yx, thetas, Ex, ix = which(clusterlocs.aic$pclocs!=0))
        bounds.aic <- calcbounds(Yx, Ex,
                                 modelthetas.aic$thetai, 
                                 modelthetas.aic$thetaa, 
                                 clusterlocs.aic$param_ix, 
                                 clusterlocs.aic$w_q,
                                 sparsematrix) 
    } else {
        print("BIC and AIC select the same path.")
        modelthetas.aic <- vector(mode = "list", length = 1)
        modelthetas.aic[[1]] <- rep(-999999,1040)
        names(modelthetas.aic) <- c("thetaa")  
        bounds.aic <- vector(mode = "list", length = 2)
        bounds.aic[[1]] <- rep(-999999, 1040)
        bounds.aic[[2]] <- rep(-999999, 1040)
        names(bounds.aic) <- c("profiled", "ma_adjusted")
    }
    
    row <- cbind(simID = rep(s,1040),
                 thetaa.bic= as.vector(modelthetas.bic$thetaa),
                 thetaa.aic = as.vector(modelthetas.aic$thetaa),
                 adjusted.bic.LB = as.vector(bounds.bic$ma_adjusted[[1]]),
                 adjusted.bic.UB = as.vector(bounds.bic$ma_adjusted[[2]]),
                 adjusted.aic.LB = as.vector(bounds.aic$ma_adjusted[[1]]),
                 adjusted.aic.UB = as.vector(bounds.aic$ma_adjusted[[2]]))
    
    master <- rbind(master, row)

}

write.csv(master, file="clustackbounds_sim_bypc.csv")


#make a quick and dirty plot
#dat <- cbind.data.frame(thetaa=as.vector(thetaa), UBa = as.vector(UBa), LBa=as.vector(LBa))
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

