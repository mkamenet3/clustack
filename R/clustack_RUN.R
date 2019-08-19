set.seed(20190326)
source("R/clustack.R")


#0) Setup
#set data
cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath
centroids <- japanbreastcancer
periods <- japanbreastcancer
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]
#set global
rMax <- 20 
Time <- 5
maxclust <- 10
locLambdas <- vector("list", maxclust)
#create set of potential clusters based on distances
potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
n <- length(x)
init <- clusso::setVectors(periods$period, expected, cases,covars=NULL, Time)
E1 <- init$E0
Ex <- clusso::scale(init, Time)
Yx <- init$Y.vec
vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)    
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparseMAT <- spacetimeMat(potentialclusters, numCenters, Time) 

#res_loc <- detectclusters(sparseMAT, Ex, Yx,numCenters,Time, maxclust=4, bylocation=TRUE)
sparsemat <- Matrix::t(sparseMAT) #66870x1040
out <- poisLik(Ex, Yx, sparsemat)
Lik <- out$Lik
Lambda_dense <- out$Lambda_dense

res <- bylocation(Lik, sparsemat, locLambdas, Lambda_dense, maxclust)

###################################################
#PLOTTING
lapply(res, function(x){plotmap((x), genpdf = FALSE)})
lapply(res, function(x){summary((x@x))})
#lapply(res_loc, function(x){plotmap(x)})

# #Saving plots for each maxclust
# vec <- 1:maxclust
# #lapply(res, function(x) summary(x@x))
# lapply(1:maxclust, function(i){
#     name <- paste0("maxclust15_byloc_2019_13_8_nclust", vec[i], ".pdf")
#     plotmap(res[[i]], pdfname = name)
# })


# #################################################################
# ############################################
# #08/14/19
# ############################################
#  maxclust <- 10
#  locLambdas <- vector("list", maxclust)
# 
# outExp <- sparsemat%*%Ex #XE
# outObs <- sparsemat%*%Yx #XY
# lambdahat <- outObs/outExp
# Lambda <- as.vector(lambdahat)*sparsemat
# Lik <- ((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs
# outlogLik <- log(Lik)
# 
# outlogLik_scaled <- outlogLik - max(outlogLik)
# Lik <- exp(outlogLik_scaled)
# 
# Lambda_dense <- as.matrix(Lambda)
# Lambda_dense[Lambda_dense == 0] <- 1
# wi <- likweights(Lik) #tiny weights
# for (i in 1:maxclust){
#     #find location wiht largest weight
#     wi_loc <-t(wi)%*%sparsemat
#     maxloc <- which.max(as.vector(wi_loc))
#     message(paste0("Location identified: ",(maxloc)))
#     #find all potential clusters that overlap that location
#     locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
#     pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
#     #partition weights vector st for pclocmax=1, those weights sum to 1
#     Lik[which(pclocmax!=0)] <-1
#     #reweight whats in the cluster so that sums to 1
#     wi[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
#     
#     out <- t(wi)%*%Lambda_dense
#     locLambdas[[i]] <- out 
#     #e) Set Lik for overlapping locations to zero
#     Lik[which(pclocmax!=0)] <-0
#     #f) Recalculate scaled likelihoods
#     wi <- likweights(Lik)
#     
# }
# 
# lapply(locLambdas, function(x) summary(x@x))
# lapply(locLambdas, function(x) table(x@x))
# lapply(locLambdas, function(x){plotmap(x, genpdf = FALSE)})
# vec <- 1:maxclust
# #lapply(res, function(x) summary(x@x))
# lapply(1:maxclust, function(i){
#     name <- paste0("maxclust10_byloc_2019_14_8_nclust", vec[i], ".pdf")
#     plotmap(res[[i]], pdfname = name)
# })
