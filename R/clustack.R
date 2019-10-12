# clustack functions
#Load Libraries
library(clusso)
library(testthat)
library(Matrix)

########################################
#Functions
########################################


colorsgrey <- function (x) {
    y = colorRamp(RColorBrewer::brewer.pal(9, "Greys")[1:9])(x)
    rgb(y[, 1], y[, 2], y[, 3], maxColorValue = 255)
}


#'@title plotmap
#'@param res result 
#'@param pdfname String for name of pdf to be generated.
#'@param genpdf Boolean. If results should be generated to console, set to \code{FALSE}. Default is \code{TRUE} for pdf to be generated.
#'@param maxrr For the color ramp, what is the maximum relative risk color. Default is for the ramp to be between 0 and 2. 
#'@param minrr For the color ramp, what is the minimum relative risk color. Default is for the ramp to be between 0 and 2. 
#'@return Maps for central region of Japan for each time period.
plotmap <- function(res, pdfname=NULL, genpdf=TRUE, maxrr=2, minrr=0){
    if(!is.null(maxrr)){
        maxrr=maxrr
    }
    else{
        maxrr=2
    }
    if(!is.null(minrr)){
        minrr=minrr
    }
    else{
        minrr=0
    }
    #cluster_ix <- redblue(log(2 *  pmax(1/2, pmin(res, 2)))/log(4))
    #cluster_ix <- redblue(log(maxrr *  pmax(1/maxrr, pmin(res, maxrr)))/log(maxrr^2))
    cluster_ix <- redblue(log(maxrr *  pmax(minrr, pmin(res, maxrr)))/log(maxrr^2))
    colors_0 <- matrix(cluster_ix, ncol=5, byrow = FALSE)
    if(genpdf==TRUE){
        pdf(pdfname, height=11, width=10)    
    }
    
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,1] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P2
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,2] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P3
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,3] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P4
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,4] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P5
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,5] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #legend
    par(fig=c(.35,.75,0,0.2), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    #text(seq(.6,1.4,length=5),rep(.45,5),seq(0,maxrr,length.out=5),srt=330,adj=0)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(minrr,maxrr,length.out=5),srt=330,adj=0)
    
    if(genpdf==TRUE){
        dev.off()    
    }
    
}

#' @title likeweights
#' @description Generate likelihood-based initial weights.
#' @param liki Likelihood for each potential cluster.
#' @return A vector of weights.
likweights <- function(liki){
    wi <- liki/sum(liki)
    if(any(is.na(liki))){
        print("There are NA likelihoods")
    }
    return(wi)
}

#'@title poisLik
#'@description Poisson-based likelihood
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@export
#'@return Likelihood for each potential cluster.
poisLik <-function(Ex, Yx, sparsemat){
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
    Lik <- exp(outlogLik_scaled)
    return(list(Lik=Lik,
                Lambda_dense=Lambda_dense))
}

#'@title bylocation
#'@description Cluster detection based on maximum location.
#'@param Lik Likelihood for each potential cluster.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param locLambdas Empty list filled by each iteration estimated weighted relative risks inside identified cluster.
#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param maxclust Maximum number of clusters to be detected.
#'@return List. First element is a Weighted relative risks for identified locations. Second element is a large matrix of weights for each iteration.
bylocation <- function(Lik, sparsemat, locLambdas, Lambda_dense,maxclust){
    wi <- likweights(Lik) #tiny weights
    #LikMAT <- matrix(rep(NA, (ncol(sparsemat)+1)*nrow(sparsemat)),ncol=(ncol(sparsemat)+1))
    wiMAT <-  matrix(rep(NA, (ncol(sparsemat)+1)*nrow(sparsemat)),ncol=(ncol(sparsemat)+1))
    #store initial weights into wiMAT first columns
    #LikMAT[,1] <- Lik@x
    wiMAT[,1] <- wi@x
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
        #reweight whats in the cluster so that sums to 1
        wi[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
        out <- t(wi)%*%Lambda_dense
        #only keep elements inside cluster
        ix <- ifelse(t(matrix(pclocmax,ncol=1))%*%sparsemat!=0,1,0)
        outID <- ifelse(ix*out==0,1,ix*out)
        locLambdas[[i]] <- outID 
        #e) Set Lik for overlapping locations to zero
        Lik[which(pclocmax!=0)] <-0
        #f) Recalculate scaled likelihoods
        wi <- likweights(Lik)
        #g) Store weights into wiMAT and Lik into LikMAT(i+1 because this wi corresponds to next iteration)
        #LikMAT[,(i+1)] <- Lik@x
        wiMAT[,(i+1)] <- wi@x
    }
    return(list(locLambdas=locLambdas,
                #LikMAT = LikMAT,
                wiMAT = wiMAT))
}

#'@title bycluster
#'@description Cluster detection based on maximum potential cluster.
#'@param Lik Likelihood for each potential cluster.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param locLambdas Empty list filled by each iteration estimated weighted relative risks inside identified cluster.
#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param maxclust Maximum number of clusters to be detected.
#'@return Weighted relative risks for identified locations.
bycluster <-  function(Lik, sparsemat, locLambdas, Lambda_dense,maxclust){
    wi <- likweights(Lik) #tiny weights
    #LikMAT <- matrix(rep(NA, (ncol(sparsemat)+1)*nrow(sparsemat)),ncol=(ncol(sparsemat)+1))
    wiMAT <-  matrix(rep(NA, (ncol(sparsemat)+1)*nrow(sparsemat)),ncol=(ncol(sparsemat)+1))
    #store initial weights into wiMAT first columns
    LikMAT[,1] <- Lik@x
    wiMAT[,1] <- wi@x
    for (i in 1:maxclust){
        message(paste0("Searching for cluster ",i))
        #find potential cluster with largest weight
        maxpc <- which.max(as.vector(wi))
        #maxloc <- which.max(as.vector(wi_loc))
        message(paste0("Potential cluster identified: ",(maxpc)))
        #find all potential clusters that overlap that PC
        pcmax <- rep(0,length(wi)); pcmax[maxpc] <-1; pcmax <- matrix(pcmax,ncol=1)
        #partition weights vector st for pclocmax=1, those weights sum to 1
        Lik[which(pcmax!=0)] <-1
        #reweight whats in the cluster so that sums to 1
        wi[which(pcmax!=0)] <- likweights(Lik[which(pcmax!=0)])
        out <- t(wi)%*%Lambda_dense
        #take only things inside the cluster
        ix <- t(pcmax)%*%sparsemat
        outID <- ifelse(ix*out==0,1,ix*out)
        locLambdas[[i]] <- outID
        #e) Set Lik for overlapping locations to zero
        Lik[which(pcmax!=0)] <-0
        #f) Recalculate scaled likelihoods
        wi <- likweights(Lik)
        #g) Store weights into wiMAT and Lik into LikMAT(i+1 because this wi corresponds to next iteration)
        #LikMAT[,(i+1)] <- Lik@x
        wiMAT[,(i+1)] <- wi@x
    }
    return(list(locLambdas=locLambdas,
                #LikMAT = LikMAT,
                wiMAT = wiMAT))
}
    
#'@title detectclusters
#'@description Detect disease clusters either by location or by potential cluster.
#'@param sparseMAT Large sparsematrix TODO
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param numCenters Number of centroids.
#'@param Time Number of time periods.
#'@param maxclust Maximum number of clusters allowed. TODO - allow this to be unknown.
#'@param bylocation If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{TRUE}.
#'@param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#'@return Returns list for each iteration with weighted relative risks by location inside identified cluster.
#'@export
detectclusters <- function(sparseMAT, Ex, Yx,numCenters,Time, maxclust,bylocation=TRUE, cv=FALSE){
    sparsemat <- Matrix::t(sparseMAT) 
    out <- poisLik(Ex, Yx, sparsemat)
    Lik <- out$Lik
    Lambda_dense <- out$Lambda_dense
    locLambdas <- vector("list", maxclust)
    if(bylocation==FALSE){
        message("Cluster detection by potential cluster")
        res <- bycluster(Lik, sparsemat, locLambdas, Lambda_dense, maxclust)
        #perform selection by IC/CV
        #selection <- clusterselect(..., cv=FALSE)
        return(res)
    }
    else{
        #default
        message("Cluster detection by location")
        res <- bylocation(Lik, sparsemat, locLambdas, Lambda_dense, maxclust) 
        #perform selection by IC/CV
        return(res)
    }
}

#wiMAT <- res$wiMAT
#a <- res$locLambdas[[1]] 
locLambdas <- res$locLambdas
clusterselect <- function(locLambdas,maxclust, numCenters, Time,cv=FALSE){
    # if(!is.null(maxclust)){
    #     maxclust=maxclust#+1
    # }
    # else{
    #     maxclust = length(locLambdas)#+1
    # }
    
    loglik <- sapply(1:maxclust, function(i) sum(dpoisson(Yx, locLambdas[[i]], Ex)))#sum(dpoisson(Yx, a, Ex))
    
    
    
    #loglik <- log(apply(wiMAT[,1:maxclust1],2,sum))
    K <- seq(from=1, to=maxclust)#rep(0,maxclust1)#seq(from=0, to=(maxclust1-1))
    if(cv==TRUE){
        #TODO
    }
    else{
        #calc
        PLL.bic <- (-2*loglik) + K*log(numCenters*Time)
        PLL.aic <- 2*K - 2*loglik
        PLL.aicc <- 2*(K) - 2*(loglik) +
            ((2*K*(K + 1))/(n_uniq*Time - K - 1))
        #select
        select.bic <- which.min(PLL.aic)
        select.aic <- which.min(PLL.aic)
        select.aicc <- which.min(PLL.aicc)
        print(c(select.bic,select.aic, select.aicc))
    }
    
}




