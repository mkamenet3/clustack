set.seed(20190326)
source("R/clustack.R")


#0) Setup
load("data/japanbreastcancer.RData")
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
#maxclust <- 10
locLambdas <- vector("list", maxclust)
#create set of potential clusters based on distances
potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
n <- length(x)
init <- clusso::setVectors(factor(jbc$period), expected, cases,covars=NULL, Time)
E1 <- init$E0
Ex <- clusso::scale(init, Time)
Yx <- init$Y.vec
vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)    
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparseMAT <- spacetimeMat(potentialclusters, numCenters, Time) 


###################################################
#RUN SL
###################################################
#By location
maxclust <- 11
test_loc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = TRUE)

plotmap(test_loc$wLambda[1,],genpdf = FALSE)
plotmap(test_loc$wLambda[5,],genpdf = FALSE)
#sapply(1:maxclust, function(i) summary(test_loc$wLambda[i,]))

#By PC
test_pc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = FALSE)
plotmap(test_loc$wLambda[3,],genpdf = FALSE)
sapply(1:maxclust, function(i) plotmap(test_pc$wLambda[i,],genpdf = FALSE))
sapply(1:maxclust, function(i) summary(test_pc$wLambda[i,]))

###################################################
#Debugging Start-Up
###################################################

# maxclust <- 1040
# outExp <- sparsemat%*%Ex
# outObs <- sparsemat%*%Yx
# #calc Lambda
# lambdahat <- outObs/outExp
# Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
# Lambda_dense <- as.matrix(Lambda)
# Lambda_dense[Lambda_dense == 0] <- 1
# #Get scaled likelihood
# Lik <- ((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs
# outlogLik <- log(Lik)
# outlogLik_scaled <- outlogLik-max(outlogLik)
# Lik <- Likorig <- exp(outlogLik_scaled)
# test <- bycluster(Lik, Lambda_dense, sparsemat,maxclust)
# # sapply(1:maxclust, function(i) plotmap(test[i,],genpdf = FALSE))
# # sapply(1:maxclust, function(i) summary(test[i,]))
