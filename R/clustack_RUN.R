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
sparsemat <- Matrix::t(sparseMAT) #66870x1040
out <- poisLik(Ex, Yx, sparsemat)
outLik <- out$outLik
lambdahat <- out$lambdahat
Lambda <- lambdahat*sparsemat #big Lambda matrix
locLambdas <- vector("list", maxclust)

res <- bylocation(outLik, sparsemat, locLambdas, Lambda, maxclust)

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


########################################
#LONG FORM Update 7-9-2019
########################################
outExp <- sparsemat%*%Ex
outObs <- sparsemat%*%Yx;outObs[outObs==0] <- 0.00001 #TOFIX
lambdahat<- outObs/outExp#log(outObs/outExp)

outLambda <- ifelse(outObs==0, 0, outObs/outExp)
#outlogLik <- log(outObs*outLambda + (sum(outObs)-outObs)*log((sum(outObs)-outObs)/(sum(outExp)-outExp)))
#outlogLik <- log(outExp*outLambda + (sum(outObs)-outObs)*log((sum(outObs)-outObs)/(sum(outExp)-outExp)))


outlogLik <- outObs * log(outLambda * outExp) - (outLambda * outExp)
#outlogLik <- outObs * ifelse(is.na(log(outLambda * outExp)),0,log(outLambda * outExp)) - (outLambda * outExp)

if(any(is.na(outlogLik))){
    ix <- which(is.na(outlogLik))
    outlogLik[ix] <- (outLambda[ix]*outExp[ix])
}

#outlogLik <- log(outObs*log(outObs/outExp) + (sum(outObs)-outObs)*log((sum(outObs)-outObs)/(sum(outExp)-outExp)))
#outlogLik <- log(outObs*outLambda + (sum(outObs)-outObs)*log((sum(outObs)-outObs)/(sum(outExp)-outExp)))
outlogLik_scaled <- outlogLik-max(outlogLik)
outLik <- exp(outlogLik_scaled)
Lambda <- outLambda*sparsemat
# Lambda <- as.vector(lambdahat)*sparsemat #t(lambdahat)%*%sparsemat #1x1040 (lambda for each location)
# outlogLik <- log(outObs*log(lambdahat) + (sum(outObs)-outObs)*log((sum(outObs)-outObs)/(sum(outExp)-outExp)))
# outlogLik_scaled <- outlogLik-max(outlogLik)
# outLik <- exp(outlogLik_scaled)
wi <- likweights2(outLik) #66870x1

lwi <- t(wi)%*%Lambda
maxlwi <- which.max(lwi@x)
message(paste0("Location identified: ",(maxlwi)))
locmax <- rep(0,numCenters*Time); locmax[maxlwi] <-1
findoverlap <- sparsemat%*%locmax
ix <- which(findoverlap!=0)
wxi <- likweights2(outLik[ix])
wi[ix] <- wxi
#tomap <- exp(t(wi)%*%Lambda)
tomap <- exp(t(wi)%*%Lambda)


#ixmap <- which(ifelse(findoverlap!=0,1,0)%*%sparsemat!=0)
#locLambdas[ixmap] <- sum(Xi*outLambda[ix])


#########
lwi <- t(wi)%*%sparsemat
maxlwi <- which.max(lwi@x) #552
message(paste0("Location identified: ",(maxlwi)))
#2) Find All PC's that Cover the max location
locmax <- rep(0,numCenters*Time); locmax[maxlwi] <-1
findoverlap <- sparsemat%*%locmax
ix <- which(findoverlap!=0)
wxi <- likweights(outLik[ix])
sum(wxi)

plotmap(log(tomap))


########################################
#LONG FORM Update 7-16-2019
########################################
maxclust <- 4
locLambdas <- vector("list", maxclust)

outExp <- sparsemat%*%Ex #XE
outObs <- sparsemat%*%Yx #XY
lambdahat <- outObs/outExp
Lambda <- as.vector(lambdahat)*sparsemat
#calc Lik
#outlogLik <- outObs * log(lambdahat*outExp) - (lambdahat*outExp)

outlogLik <- log(outObs*outLambda + (sum(outObs)-outObs)*log((sum(outObs)-outObs)/(sum(outExp)-outExp)))
if(any(is.na(outlogLik))){
    ix <- which(is.na(outlogLik))
    outlogLik[ix] <- (lambdahat[ix]*outExp[ix])
} #568 zeros after this in log lik
outlogLik_scaled <- outlogLik - max(outlogLik)
Lik <- exp(outlogLik_scaled)

for(i in 1:maxclust){
    wi <- likweights(Lik)
    lwi <- t(wi)%*%sparsemat
    maxlwi <- which.max(lwi@x) #552
    message(paste0("Location identified: ",(maxlwi)))
    
    locmax <- rep(0,numCenters*Time); locmax[maxlwi] <-1
    findoverlap <- sparsemat%*%locmax
    ix <- which(findoverlap!=0) #all PCs that overlap
    
    
    #upweight
    Lik[ix] <- 1
    wxi <- likweights(Lik[ix])
    print(sum(wxi))
    
    wi[ix] <- wxi
    #4) Find weighted Lambda by locations
    locLambdas[[i]] <- t(wi)%*%Lambda #locLambdas[[i]] <- lwi
    #5) Set Likelihood in cluster i to zero
    Lik[ix] <- 0
}







