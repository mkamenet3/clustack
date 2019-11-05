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

sumex <- function(vec){
    return(sprintf("%.16f",sum(vec)))
}

#ITER 1
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
    ixpclocmax1 <- which(pclocmax!=0)
    
    #upweight weights
    Lik[which(pclocmax!=0)] <-1
    wi_0[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
    locLambdas <- t(matrix(wi_0@x,ncol = 1))%*%Lambda_dense #Ok this is correct
    
    #Likelihood
    #NULL: dpoisson(Yx, rep(1,1040), Ex) #8138.282
    dpoisson(Yx, as.vector(locLambdas), Ex) # 8152.666
    
    #Remove the Likelihood from being there
    takeout <- which(pclocmax!=0)
    Lik[takeout] <- 0
    #Lik[which(pclocmax!=0)] <-0
###############################    
#ITER 2
    wi2 <- wi_0
    
    wi2[which(pclocmax==0)] <- likweights(Lik[which(pclocmax==0)]) #these weights sum to 2 now which makes sense I think
    sprintf("%.16f",sum(wi2)) #this is all 1
    #3) Find max loc
    #find location with largest weight
    wtemp <- rep(0,66870)
    wtemp[which(pclocmax==0)] <- wi2[which(pclocmax==0)]
    sumex(wtemp)
    
    
    #
    #test <- t(wtemp)%*%sparsemat
    
    wi2_loc <-t(wtemp)%*%sparsemat
    maxloc2 <- which.max(as.vector(wi2_loc))
    # 
    # wi2_loc <- t(wi2)%*%sparsemat
    # maxloc2 <- which.max(as.vector(wi2_loc))
    message(paste0("Location identified: ",(maxloc2))) #this is correct 331
    #find all potential clusters that overlap that location
    locmax2 <- rep(0,numCenters*Time); locmax2[maxloc2] <-1; locmax2 <- matrix(locmax2,ncol=1)
    pclocmax2 <- as.vector(t(locmax2)%*%t(sparsemat)) #all locations that overlap max location
    
    Lik[which(pclocmax2!=0)] <-1
    wtemp[which(pclocmax2!=0)] <- likweights(Lik[which(pclocmax2!=0)])
    
    #wi_0[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
    locLambdas2 <- t(matrix(wtemp,ncol = 1))%*%Lambda_dense #Ok this is correct
    
    #ixpclocmax2 <- which(pclocmax2!=0)
    #ix <-ixpclocmax1[which(ixpclocmax1%in%ixpclocmax2==FALSE)]
    #upweight weights
    Lik[ix] <-1
    wi2[ix] <- likweights(Lik[ix])
    #wi_0[which(pclocmax2!=0)] <- likweights(Lik[which(pclocmax2!=0)])
    #sprintf("%.16f",sum(wi_0))
    sumex(wi2)
    locLambdas2 <- t(matrix(wi2@x,ncol = 1))%*%Lambda_dense #Ok this is correct
    #locLambdas2 <- t(matrix(wi_0,ncol = 1))%*%Lambda_dense #Ok this is correct
    dpoisson(Yx, as.vector(locLambdas2), Ex) #8150.041
    
    
#############################################################
#Debugging 2019-11-04
##It's the weight vector
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
Lik <- Likorig <- exp(outlogLik_scaled) #this Lik is Lik* (because its scaled)

sumex <- function(vec){
    return(sprintf("%.16f",sum(vec)))
}
#############
dpoisson(Yx, rep(1, 1040),Ex) # 8138.282
##ITER 1
#############
wfinal<- rep(0,66870)
w0 <- likweights(Lik)
#Find maxloc
wi_loc <- t(w0)%*%sparsemat
maxloc <- which.max(as.vector(wi_loc))
message(paste0("Location identified: ",(maxloc)))
#find all potential clusters that overlap that location
locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
pclocmax <- as.vector(t(locmax)%*%t(sparsemat)) #all locations that overlap max location [1x66870]
ix1 <- which(pclocmax!=0)
#upweight Lik* to 1 in these overlapping Pc's
Lik[which(pclocmax!=0)] <-1
#In all pclocmax !=0, reweight so that those sum to 1
w0_1 <- w0
w0_1[ix1] <- likweights(Lik[ix1])    
sumex(w0_1)
wfinal[ix1] <- likweights(Lik[ix1]); sumex(wfinal)
#Find weighted RR
locLambdas1 <- t(matrix(wfinal,ncol = 1))%*%Lambda_dense #Ok this is correct
sum(Lik0*wfinal)
#locLambdas1 <- t(matrix(w0_1@x,ncol = 1))%*%Lambda_dense #Ok this is correct
#Remove those pcs that overlap (pclocmax!=0) from being eligible again
Lik[ix1] <-0
#Now, reweight only Lik that were not in C1
#w1 <- w0_1
w0_1[-ix1] <- likweights(Lik[-ix1])
# w00[which(pclocmax==0)] <- likweights(Lik[which(pclocmax==0)])
# w1[which(pclocmax==0)] <- likweights(Lik[which(pclocmax==0)])
sumex(w1);sumex(w0_1)
dpoisson(Yx, as.vector(locLambdas1), Ex) #8152.666

#############
##ITER 2
#############
#Find second maxloc
#Find maxloc
wi_loc <- t(w0_1)%*%sparsemat
maxloc <- which.max(as.vector(wi_loc))
message(paste0("Location identified: ",(maxloc)))
#find all potential clusters that overlap that max location (but are not any of the PCsfrom before)
locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
pclocmax2 <- as.vector(t(locmax)%*%t(sparsemat)) #all locations that overlap max location [1x66870]
ix2 <- which(pclocmax2!=0)

pclocuniq <- ix2[which(!ix2 %in% ix1)]#intersect(which(pclocmax!=0), which(pclocmax2!=0))
#upweight Lik* to 1 in these overlapping Pc's
Lik[pclocuniq] <-1
#In all pclocmax !=0, reweight so that those sum to 1
w0_1[-pclocuniq] <- likweights(Lik[-pclocuniq])
wfinal[pclocuniq] <- likweights(Lik[pclocuniq]); sumex(wfinal)
sum(Lik0*wfinal)

locLambdas2 <- t(matrix(wfinal,ncol = 1))%*%Lambda_dense #Ok this is correct
dpoisson(Yx, as.vector(locLambdas2)/sum(wfinal), Ex) #8150.634

#############
##ITER 3
#############
wi_loc <- t(w0_1)%*%sparsemat
maxloc <- which.max(as.vector(wi_loc))
message(paste0("Location identified: ",(maxloc)))
#find all potential clusters that overlap that max location (but are not any of the PCsfrom before)
locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
pclocmax3 <- as.vector(t(locmax)%*%t(sparsemat)) #all locations that overlap max location [1x66870]
ix3 <- which(pclocmax3!=0)

pclocuniq3 <- ix3[which(!ix3 %in% unique(c(ix2,ix1)))]#intersect(which(pclocmax!=0), which(pclocmax2!=0))
#upweight Lik* to 1 in these overlapping Pc's
Lik[pclocuniq3] <-1
#In all pclocmax !=0, reweight so that those sum to 1
w0_1[-pclocuniq3] <- likweights(Lik[-pclocuniq3]) #I think this is wrong
wfinal[pclocuniq3] <- likweights(Lik[pclocuniq3]); sumex(wfinal)
sum(Likorig*wfinal)
locLambdas3 <- t(matrix(wfinal,ncol = 1))%*%Lambda_dense #Ok this is correct
dpoisson(Yx, as.vector(locLambdas3)/sum(wfinal), Ex) #8149.634

#############
##ITER 4
#############
wi_loc <- t(w0_1)%*%sparsemat
maxloc <- which.max(as.vector(wi_loc))
message(paste0("Location identified: ",(maxloc))) #353
#find all potential clusters that overlap that max location (but are not any of the PCsfrom before)
locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
pclocmax4 <- as.vector(t(locmax)%*%t(sparsemat)) #all locations that overlap max location [1x66870]
ix4 <- which(pclocmax4!=0)

pclocuniq4 <- ix4[which(!ix4 %in% unique(c(ix2,ix1,ix3)))]#intersect(which(pclocmax!=0), which(pclocmax2!=0))
#upweight Lik* to 1 in these overlapping Pc's
Lik[pclocuniq4] <-1
#In all pclocmax !=0, reweight so that those sum to 1
w0_1[-pclocuniq4] <- likweights(Lik[-pclocuniq4]) #I think this is wrong
wfinal[pclocuniq4] <- likweights(Lik[pclocuniq4]); sumex(wfinal)
#sum(Likorig*wfinal)
locLambdas4 <- t(matrix(wfinal,ncol = 1))%*%Lambda_dense #Ok this is correct
dpoisson(Yx, as.vector(locLambdas4)/sum(wfinal), Ex) #8148.047


#############
##ITER 5
#############
wi_loc <- t(w0_1)%*%sparsemat
maxloc <- which.max(as.vector(wi_loc))
message(paste0("Location identified: ",(maxloc))) #353
#find all potential clusters that overlap that max location (but are not any of the PCsfrom before)
locmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
pclocmax5 <- as.vector(t(locmax)%*%t(sparsemat)) #all locations that overlap max location [1x66870]
ix5 <- which(pclocmax5!=0)

pclocuniq5 <- ix5[which(!ix5 %in% unique(c(ix2,ix1,ix3,ix4)))]#intersect(which(pclocmax!=0), which(pclocmax2!=0))
#upweight Lik* to 1 in these overlapping Pc's
Lik[pclocuniq5] <-1
#In all pclocmax !=0, reweight so that those sum to 1
w0_1[-pclocuniq5] <- likweights(Lik[-pclocuniq5]) #I think this is wrong
wfinal[pclocuniq5] <- likweights(Lik[pclocuniq5]); sumex(wfinal)
#sum(Likorig*wfinal)
locLambdas5 <- t(matrix(wfinal,ncol = 1))%*%Lambda_dense #Ok this is correct
dpoisson(Yx, as.vector(locLambdas5)/sum(wfinal), Ex) #8148.047

# w1_1 <- w1; sumex(w1_1)
# w1_1[pclocuniq] <- likweights(Lik[pclocuniq]); sumex(w1_1)
# w00[pclocuniq] <-  likweights(Lik[pclocuniq]); sumex(w00)

# a <- rep(NA,10)
# b <- c(2,3,9)
# a[c(2,6,8)] <- b
# a <- c(2,3,5,6)
# b <- c(5,7,8,9)
# intersect(a,b)

#############################################################
#Debugging 2019-11-04
#CLEAN UP CLEAN UP CLEAN THIS SHIT UP AND ITERATE OVER
############################################################   
#init
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
Lik <- Likorig <- exp(outlogLik_scaled) #this Lik is Lik* (because its scaled)


foo <- function(maxclust,...){
    print(paste0("Init logLik: ",dpoisson(Yx, rep(1, 1040),Ex)))
    wfinal<- rep(0,66870)
    #takeout <- rep(0,66870)
    takeout <- NULL
    #begin iter
    for (i in 1:maxclust){
        wtmp <- likweights(Lik)    
        #Find maxloc
        wi_loc <- t(wtmp)%*%sparsemat
        maxloc <- which.max(as.vector(wi_loc))
        message(paste0("Location identified: ",(maxloc)))
        #find all potential clusters that overlap that location
        locmax <- rep(0,numCenters*Time); 
            locmax[maxloc] <-1;
            locmax <- matrix(locmax,ncol=1)
        pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
        #all locations that overlap max location [1x66870]
        ix1 <- which(pclocmax!=0)
        takeout <- unique(c(takeout, ix1))
        #upweight Lik* to 1 in these overlapping Pc's
        Lik[ix1] <-1
        #In all pclocmax !=0, reweight so that those sum to 1
        wfinal[ix1] <- likweights(Lik[ix1])
        wtmp[ix1] <- likweights(Lik[ix1])
        
        print(paste0("Sum wfinal: ",sumex(wfinal)))
        print(paste0("Sum wtmp: ",sumex(wtmp)))
        #Find weighted RR
        locLambdas1 <- t(matrix(wfinal,ncol = 1))%*%Lambda_dense 
        #poisson logLik
        print(paste0("logLike: ",dpoisson(Yx, as.vector(locLambdas1), Ex)))
        #Remove those pcs that overlap (pclocmax!=0) from being eligible again
        Lik[takeout] <-0 #check this
        #Now, reweight only Lik that were not in C1
        wtmp[-takeout] <- likweights(Lik[-takeout])
        print(paste0("sum wtmp: ", sumex(wtmp)))
    }
}
foo(3)


#############################################################
#Debugging 2019-11-05
#WHAT ATTEMPT AMI ON?
#I don't think sum of wi will give me total number of clusters anymore
#Nor do I think likelihood needs to be increasing
############################################################  

#init
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
Lik <- Likorig <- exp(outlogLik_scaled) #this Lik is Lik* (because its scaled)

wtMAT <- matrix(rep(NA, 1040*66870), ncol=1040)
###################
#Iter1
###################
foo <- function(Lik, wtMAT, maxclust){
    for(i in 1:maxclust){
        wtmp <- likweights(Lik)
        #Find maxloc
        wi_loc <- t(wtmp)%*%sparsemat
        maxloc <- which.max(as.vector(wi_loc))
        message(paste0("Location identified: ",(maxloc)))
        #find all potential clusters that overlap that location
        locmax <- rep(0,numCenters*Time); 
        locmax[maxloc] <-1;
        locmax <- matrix(locmax,ncol=1)
        pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
        ix <- which(pclocmax!=0) #indices of all PCs that overlap max location
        #upweight Lik* to 1 in all Pcs that overlap max cluster
        Lik[ix] <-1
        #reweight so that everything inside the PCs that overlap max cluster sum to 1
        wtmp[ix] <- likweights(Lik[ix]) 
        wtMAT[,i] <- wtmp@x
        #set Lik in everything inside the Pcs that overlap max cluster to 0
        Lik[ix] <-0    
    }
    return(wtMAT)
    
}
    
a <-foo(Lik, wtMAT, 1040)
test <- t(a)%*%Lambda_dense
likMat <- sapply(1:1040, function(i) dpoisson(Yx, test[i,], Ex))









