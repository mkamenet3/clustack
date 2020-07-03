#Calculating Clustack Bounds
#Simulation Study
#M.Kamenetsky
#v1: 2020-06-30

rm(list=ls())

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
    M <- Matrix(0, nrow=dim(Lambda_dense)[1], ncol=length(param_ix), sparse=TRUE)
    for (i in 1:length(param_ix)){
        M[,i][param_ix[i]] <- 1    
    }
    aa <- t(M)%*%Lambda_dense
    thetas <- apply(aa,1, function(x) {
        if(length(unique(x))==1){
            valsout <- unique(x)
        } else {
            valsout <- unique(x)[which(unique(x)!=1)]
        }
        return(valsout)
    })
    thetaa <- sum(w_q*unlist(thetas))
    return(list(thetai = thetas,
                thetaa = thetaa))
}

# vartheta <- function(input, Yx, Ex, thetai, thetaa,sparsematrix){
#     param_ix_k <- input[[1]]
#     w_q_k <- input[[2]]
#     amod <- glm(Yx ~ sparsematrix[,param_ix_k] -1, offset=Ex, family="poisson")
#     deviances <- amod$deviance
#     conddeviances <- coef(summary(amod))[,2]
#     variance <- sum(w_q_k*sqrt((coef(summary(amod))[,2])^2 + (thetai - thetaa)^2))
#     return(variance)
# }




#calculating unconditional variance for model-averaged bounds
#' @param thetai Vector of estimated theta values for each Q model
#' @param thetaa model-averaged estimate of theta
#' bounds_adjustma <- function(thetai, thetaa, param_ix, w_q,sparsematrix){
library(data.table)
bounds_adjustma <- function(thetai, thetaa, param_ix, w_q,sparsematrix){
    variance <- vector(mode = "list", length = length(thetai))
    system.time(for (i in 1:length(thetai)){
        amod <- glm(Yx ~ sparsematrix[,param_ix[k]] -1, offset=Ex, family="poisson")
        variance[[k]] <- w_q[k]*sqrt((coef(summary(amod))[,2])^2 + (thetai[k] - thetaa)^2)
    })
    (var_theta1 <- sum(unlist(variance)))
    var_thetaa <- sum(unlist(variance))
    UBa = exp(log(thetaa) + 1.96*sqrt(var_thetaa))
    LBa = exp(log(thetaa) - 1.96*sqrt(var_thetaa))
    return(adjusted = c(LBa, UBa))
}



#calculate profile-likelihood confidence bounds and adjust for model-averaaging
#' @param null Null model likelihood
#' @param proflik Vector of profiled likelihoods for set of thetas
calcbounds <- function(null, out, Yx, Ex, thetai,thetaa, param_ix, w_q,sparsematrix, adjustedbounds) {
    variance <- vector(mode = "list", length = length(thetai))
    
    crit <- (null)-(1.92)
    changepoint <- which.max(out)
    ix_ub <- which.min(abs(out[changepoint:length(out)]-crit))
    UB <- thetas[ix_ub+changepoint-1]
    ix_lb <- which.min(abs(out[1:changepoint]-crit))
    LB <- thetas[ix_lb]
    if (adjustedbounds==TRUE){
        # system.time(mods <- lapply(1:length(thetai), function(x)  glm(Yx ~ sparsematrix[,param_ix[x]] -1, offset=Ex, family="poisson")))
        # #lapply(1:length(thetai), function(x)  glm(Yx ~ sparsematrix[,param_ix[x]] -1, offset=Ex, family="poisson"))
        # system.time(vars <- lapply(1:length(thetai), function(x) w_q[x]*sqrt((coef(summary(mods[[x]]))[,2])^2 + (thetai[x] - thetaa)^2)))
        system.time(mods <- lapply(1:20, function(x)  glm(Yx ~ sparsematrix[,param_ix[x]] -1, offset=Ex, family="poisson")))
        #lapply(1:length(thetai), function(x)  glm(Yx ~ sparsematrix[,param_ix[x]] -1, offset=Ex, family="poisson"))
        system.time(vars <- lapply(1:20, function(x) w_q[x]*sqrt((coef(summary(mods[[x]]))[,2])^2 + (thetai[x] - thetaa)^2)))
        #adjusted <- bounds_adjustma(thetai, thetaa, param_ix, w_q,sparsematrix)
        var_thetaa <- sum(unlist(vars))
        UBa = exp(log(thetaa) + 1.96*sqrt(var_thetaa))
        LBa = exp(log(thetaa) - 1.96*sqrt(var_thetaa))
        
    } else {
        adjusted <- NULL
    }
    
    return(list(profiled = c(LB, UB),
           ma_adjusted = c(LBa, UBa)))
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
nsim <-5
theta = 100
risk = 2
cent = 150
rad = 11
tim <- c(1:5)


##############################
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
#E1 <- as.vector(rr)*init$E0
message(paste("Running model for periods",tim[1],"through", tail(tim, n=1)))
#simulate data here
# YSIM <- rnegbin(E1, theta = theta)
# summary(YSIM/E1)
#Ex <- unlist(scale_sim(list(YSIM), init, nsim=1, Time))
sparseMAT <- spacetimeMat(clusters , numCenters, Time) #maxclust <- 10
sparsemat <- t(sparseMAT) 

#Try implanting a cluster instead
ix <- which(rr!=1)
rr_bin <- rep(0,1040)
rr_bin[ix] <-1

overlap <- matrix(rr_bin, nrow=1)%*%sparsematrix
overlap <- ifelse(overlap!=0,1,0)
#identify which PCs overlap cluster
#simulate cluster risk here
qgamma(c(0.01, 0.99),2000,1000)

lambdahat <- rep(1, 66870)
lambdahat[which(overlap!=0)] <- rgamma(length(which(overlap!=0)), 2000, 1000)

outObs <- rbinom(66870, size = 1, prob=0.05)#rpois(n = 66870, lambda=1)
outExp <- outObs/lambdahat
#test <- outObs/lambdahat
#outExp <- rpois(n = 66870, lambda=outObs)

#E1 <- rpois(1040,2)
#outExp <- (sparsemat%*%E1)
#outObs<-  rpois(length(outExp), lambda = (outExp@x*lambda)) #outExp*lambdahat


#Create other things needed
#outExp <- sparsemat%*%Ex_mod
#outObs <- sparsemat%*%Yx_mod
#lambdahat <- outObs/outExp
Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
Lambda_dense <- as.matrix(Lambda)
Lambda_dense[Lambda_dense == 0] <- 1

#back find these
Ex_mod <- t(outExp)%*%t(sparsematrix)
Yx_mod <- t(outObs)%*%t(sparsematrix)
#####################################################################################
#####################################################################################
#####################################################################################
#CLUSTER DETECTION BY Potential Cluster
#####################################################################################
#####################################################################################
#####################################################################################

sim_superclust_pc<- detectclusters(sparsematrix, Ex_mod@x, Yx_mod@x,
                                   numCenters, Time, maxclust,
                                   bylocation = FALSE, model="poisson",
                                   overdisp.est = NULL)
thetas <- seq(0.1,3, length=1000)
null <- dpoisson_theta(Yx_mod, 1, Ex_mod, ix=1:dim(sparsematrix)[1])

#BIC
clusterlocs.bic <- clusterlocs_ident(sim_superclust_pc, selection="selection.bic", sparseMAT)
modelthetas.bic <- extract_thetai(clusterlocs.bic$param_ix, clusterlocs.bic$w_q, Lambda_dense)
#calculate null and model curves
out.bic <- dpoisson_theta(Yx_mod@x, thetas, Ex_mod@x, ix = which(clusterlocs.bic$pclocs!=0))
bounds.bic <- calcbounds(null, out.bic, Yx_mod@x, Ex_mod@x,
                         modelthetas.bic$thetai, modelthetas.bic$thetaa, clusterlocs.bic$param_ix, 
                         clusterlocs.bic$w_q,sparsematrix, adjustedbounds = TRUE) 
if(sim_superclust_pc$selection.bic!=sim_superclust_pc$selection.aic){
    #AIC
    clusterlocs.aic <- clusterlocs_ident(sim_superclust_pc, selection="selection.aic", sparseMAT)
    modelthetas.aic <- extract_thetai(clusterlocs.aic$param_ix, clusterlocs.aic$w_q, Lambda_dense)
    #calculate null and model curves
    out.aic <- dpoisson_theta(Yx_mod, thetas, Ex_mod, ix = which(clusterlocs.aic$pclocs!=0))
    bounds.aic <- calcbounds(null, out.aic, Yx_mod@x, Ex_mod@x,
                             modelthetas.aic$thetai, 
                             modelthetas.aic$thetaa, clusterlocs.aic$param_ix, 
                             clusterlocs.aic$w_q,sparsematrix, adjustedbounds = TRUE) 
} else {
    print("BIC and AIC select same.")
    modelthetas.aic <- vector(mode = "list", length = 1)
    modelthetas.aic[[1]] <- -999999
    names(modelthetas.aic) <- c("thetaa")  
    bounds.aic <- vector(mode = "list", length = 2)
    bounds.aic[[1]] <- c(-999999, -999999)
    bounds.aic[[2]] <- c(-999999, -999999)
    names(bounds.aic) <- c("profiled", "ma_adjusted")
}

row <- cbind(simID = c(1,1),thetaa.bic= modelthetas.bic$thetaa,
             thetaa.aic = modelthetas.aic$thetaa,
             profiled.bic = bounds.bic$profiled,
             adjusted.bic = bounds.bic$ma_adjusted,
             profiled.aic = bounds.aic$profiled,
             adjusted.aic = bounds.aic$ma_adjusted)

master <- rbind(master, row)


pdfname <- paste0("simID_", s, ".pdf")
plotmapAllIC(res.bic = sim_superclust_pc$wLambda[sim_superclust_pc$selection.bic,],
             res.aic = sim_superclust_pc$wLambda[sim_superclust_pc$selection.aic,],
             oracle = Yx_mod@x/Ex_mod@x,
             pdfname = pdfname,
             genpdf = FALSE,
             maxrr = 2,
             minrr = 0.5,
             obs = FALSE)

