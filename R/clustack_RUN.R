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
maxclust <- 10
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
test_loc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = TRUE)
lapply(test_loc, function(x){plotmap((x), genpdf = FALSE)})
lapply(test_loc, function(x){plotmap((x), genpdf = FALSE, maxrr = 1.5, minrr = 0.5)})
lapply(test_loc, function(x){summary((x))})

#By PC
test_pc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = FALSE)
lapply(test_pc, function(x){plotmap((x), genpdf = FALSE, maxrr=2)})
lapply(test_pc, function(x){summary((x))})


######################################################################################################

###################################################
#Test out new funcs
###################################################
sparsemat <- Matrix::t(sparseMAT) #66870x1040
out <- poisLik(Ex, Yx, sparsemat)	# out <- poisLik(Ex, Yx, sparsemat)
Lik <- out$Lik	# Lik <- out$Lik
Lambda_dense <- out$Lambda_dense	# Lambda_dense <- out$Lambda_dense
#by loc
res <- bylocation(Lik, sparsemat, locLambdas, Lambda_dense, maxclust)



#lapply(res_pc, function(x){plotmap((x), genpdf = FALSE)})	
#by pc
#res_pc <- bycluster(Lik, sparsemat, locLambdas, Lambda_dense, maxclust)




