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
maxclust <- 1040
test_loc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = TRUE, model="poisson")
#BIC
plotmap(test_loc$wLambda[test_loc$selection.bic,],genpdf = FALSE)
summary(test_loc$wLambda[test_loc$selection.bic,])

#AIC/AICc
plotmap(test_loc$wLambda[test_loc$selection.aic,],genpdf = FALSE)
summary(test_loc$wLambda[test_loc$selection.aic,])



#AICc
plotmap(test_loc$wLambda[test_loc$selection.aicc,],genpdf = FALSE)
summary(test_loc$wLambda[test_loc$selection.aicc,])

###################
#By PC
test_pc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = FALSE, model="poisson")

#BIC and AIC/AICc all select same thing
plotmap(test_pc$wLambda[test_pc$selection.bic,],genpdf = FALSE)
plotmap(test_pc$wLambda[test_pc$selection.aic,],genpdf = FALSE)

summary(test_pc$wLambda[test_pc$selection.bic,])

