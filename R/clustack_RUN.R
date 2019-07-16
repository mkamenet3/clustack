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
maxclust <- 5
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


###################################################
#PLOTTING
lapply(res, function(x){plotmap(exp(x))})
lapply(res_loc, function(x){plotmap(x)})

#Saving plots for each maxclust
vec <- 1:maxclust
lapply(res, function(x) summary(x@x))
lapply(1:maxclust, function(i){
    name <- paste0("maxclust15_byloc_2019_10_7_nclust", vec[i], ".pdf")
    plotmap(res[[i]], pdfname = name)
})


