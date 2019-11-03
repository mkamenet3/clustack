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
maxclust <- 11
test_loc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = TRUE)
lapply(test_loc, function(x){plotmap((x), genpdf = FALSE)})
lapply(test_loc, function(x){plotmap((x), genpdf = FALSE, maxrr = 1.5, minrr = 0.5)})
lapply(test_loc, function(x){summary((x))})
##plot overlapping clusters
plotmap(test_loc$selection, genpdf = FALSE)

#By PC
test_pc <-detectclusters(sparseMAT, Ex, Yx, numCenters, Time, maxclust, bylocation = FALSE)
lapply(test_pc, function(x){plotmap((x), genpdf = FALSE, maxrr=2)})
lapply(test_pc, function(x){summary((x))})
##plot overlapping clusters
plotmap(test_pc$selection, genpdf = FALSE)


######################################################################################################

###################################################
#Test out new funcs
###################################################
# sparsemat <- Matrix::t(sparseMAT) #66870x1040
# out <- poisLik(Ex, Yx, sparsemat)	# out <- poisLik(Ex, Yx, sparsemat)
# Lik <- out$Lik	# Lik <- out$Lik
# Lambda_dense <- out$Lambda_dense	# Lambda_dense <- out$Lambda_dense
# #by loc
# res <- bylocation(Lik, sparsemat, locLambdas, Lambda_dense, maxclust)



#lapply(res_pc, function(x){plotmap((x), genpdf = FALSE)})	
#by pc
#res_pc <- bycluster(Lik, sparsemat, locLambdas, Lambda_dense, maxclust)


#What's wrong with my likelihood

outExp <- sparsemat%*%Ex
outObs <- sparsemat%*%Yx
#calc Lambda
lambdahat <- outObs/outExp
Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
Lambda_dense <- as.matrix(Lambda)
Lambda_dense[Lambda_dense == 0] <- 1
#Get scaled likelihood
Lik <- ((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs
outlogLik <- log(Lik)
outlogLik_scaled <- outlogLik-max(outlogLik)
Lik <- Likorig <- exp(outlogLik_scaled)
log(sum(Likorig))
#[1] 8.179592e-06

wi <- likweights(Lik) #tiny weights
# 
# wiMAT <-  matrix(rep(NA, (ncol(sparsemat)+1)*nrow(sparsemat)),ncol=(ncol(sparsemat)+1))
# LikMAT <- matrix(rep(NA, (ncol(sparsemat)+1)*nrow(sparsemat)),ncol=(ncol(sparsemat)+1))
#store initial weights into wiMAT first columns
#LikMAT[,1] <- Lik@x
#wiMAT[,1] <- wi@x
a <- rep(NA, 11)
for (i in 1:maxclust){
    message(paste0("Searching for cluster ",i))
    #find location with largest weight
    wi_loc <- t(wi)%*%sparsemat
    maxloc <- which.max(as.vector(wi_loc))
    message(paste0("Location identified: ",(maxloc)))
    #find all potential clusters that overlap that location
    locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
    pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
    #partition weights vector st for pclocmax=1, those weights sum to 1
    Lik[which(pclocmax!=0)] <-1
    Likorig[which(pclocmax!=0)] <-1
    #print(sum(Lik@x))
    print(paste0("Unlogged orig: ",log(sum(Likorig@x)))) #7.254885
    a[i] <- sum(Likorig@x)
    print(paste0("Unlogged1: ",sum(Lik@x))) 
    #8.136226
    #7.842671
    # 7.853993
    #7.644919
    #7.352441
    #8.286269
    
    
    
    #LikMAT[,(i+1)] <- Lik@x
    #reweight whats in the cluster so that sums to 1
    wi[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
    out <- t(wi)%*%Lambda_dense
    print(paste0("something: ",sum(Lik)))
    
    print(paste0("withall out: ",dpoisson(Yx, out@x,Ex)))
    #idx <- which(round(out@x,5)!=1)
    
    #print(paste0("withall out: ",dpoisson(Yx, ifelse(round(out@x,5)==1,0, out@x),Ex)))
   # wiMAT[,(i+1)] <- wi@x
    #only keep elements inside cluster
    ix <- ifelse(t(matrix(pclocmax,ncol=1))%*%sparsemat!=0,1,0)
    outID <- ifelse(ix*out==0,1,ix*out)
    locLambdas[[i]] <- outID 
    #e) Set Lik for overlapping locations to zero
    Lik[which(pclocmax!=0)] <-0
    
    
    print(paste0("Un-logged2: ", sum(Lik)))
    #a[i] <- sum(Lik)
    #print(paste0("Logged2: ",log(sum(Likorig))))
    print(paste0("Onlythose inside cluster:",dpoisson(Yx, locLambdas[[i]],Ex)))
    #a[i] <- dpoisson(Yx, locLambdas[[i]],Ex)
    # [1] -52.8276
    #-110.4767
    #-152.6296
    # -207.8002
    #-216.5021
    #-232.7549
    #-254.6563
    
    #f) Recalculate scaled likelihoods
    wi <- likweights(Lik)
    
    #g) Store weights into wiMAT and Lik into LikMAT(i+1 because this wi corresponds to next iteration)
    #LikMAT[,(i+1)] <- Lik@x
    #wiMAT[,(i+1)] <- wi@x
}


#############################################################
#Debugging 2019-11-03
##Issue in what i'm reweighting
############################################################
sparsemat <- Matrix::t(sparseMAT) 
out <- poisLik(Ex, Yx, sparsemat)
#Lik0 <- out$Lik
outExp <- sparsemat%*%Ex
outObs <- sparsemat%*%Yx
#calc Lambda
lambdahat <- outObs/outExp
Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
Lambda_dense <- as.matrix(Lambda)
Lambda_dense[Lambda_dense == 0] <- 1
#1) Get scaled likelihood
Lik0 <- ((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs
outlogLik <- log(Lik0)
outlogLik_scaled <- outlogLik-max(outlogLik)
Lik <- Likorig <- exp(outlogLik_scaled)

#2) Get (scaled) Likelihood-based weights
wi <- wi_0 <- likweights(Lik)
sum(wi@x) #weights sum to 1

#3) Find max loc
#find location with largest weight
wi_loc <- t(wi)%*%sparsemat
maxloc <- which.max(as.vector(wi_loc))
message(paste0("Location identified: ",(maxloc)))
#find all potential clusters that overlap that location
locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
pclocmax <- as.vector(t(locmax)%*%t(sparsemat)) #all locations that overlap max location

#upweight weights
wi_0[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
test <- t(matrix(wi_0@x,ncol = 1))%*%Lambda_dense #Ok this is correct
