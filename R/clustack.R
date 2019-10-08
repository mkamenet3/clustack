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
plotmap <- function(res, pdfname=NULL, genpdf=TRUE){
    cluster_ix <- redblue(log(2 *  pmax(1/2, pmin(res, 2)))/log(4))
    colors_0 <- matrix(cluster_ix, ncol=5, byrow = FALSE)
    if(genpdf==TRUE){
        pdf(pdfname, height=11, width=10)    
    }
    
    
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,1] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,paste0("Period 1"),cex=1.00)
    
    #P2
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,2] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,paste0("Period 2"),cex=1.00)
    
    #P3
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,3] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,paste0("Period 3"),cex=1.00)
    
    #P4
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,4] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,paste0("Period 4"),cex=1.00)
    
    #P5
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,5] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,paste0("Period 5"),cex=1.00)
    
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
    
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
#'@param locLambdas TODO
#'@param Lambda_dense TODO
#'@param maxclust Maximum number of clusters to be detected.
#'@return Weighted relative risks for identified locations.
bylocation <- function(Lik, sparsemat, locLambdas, Lambda_dense,maxclust){
    wi <- likweights(Lik) #tiny weights
    for (i in 1:maxclust){
        message(paste0("Searching for cluster ",i))
        #find location with largest weight
        wi_loc <-t(wi)%*%sparsemat
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
        locLambdas[[i]] <- out 
        #e) Set Lik for overlapping locations to zero
        Lik[which(pclocmax!=0)] <-0
        #f) Recalculate scaled likelihoods
        wi <- likweights(Lik)
    }
    return(locLambdas=locLambdas)
}

#'@title bycluster
#'@description Cluster detection based on maximum potential cluster.
#'@param Lik Likelihood for each potential cluster.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param locLambdas TODO
#'@param Lambda_dense TODO
#'@param maxclust Maximum number of clusters to be detected.
#'@return Weighted relative risks for identified locations.
bycluster <-  function(Lik, sparsemat, locLambdas, Lambda_dense,maxclust){
    wi <- likweights(Lik) #tiny weights
    for (i in 1:maxclust){
        message(paste0("Searching for cluster ",i))
        #find potential cluster with largest weight
        #wi_loc <- t(wi)%*%sparsemat
        maxpc <- which.max(as.vector(wi))
        #maxloc <- which.max(as.vector(wi_loc))
        message(paste0("Potential cluster identified: ",(maxpc)))
        #find all potential clusters that overlap that PC
        pcmax <- rep(0,length(wi)); pcmax[maxpc] <-1; pcmax <- matrix(pcmax,ncol=1)
        pclocmax_locs <- as.vector((as.vector(t(pcmax))%*%sparsemat)%*%t(sparsemat))
        pclocmax <- ifelse(pclocmax_locs!=0,1,0)
        
        
        #pcmax <- rep(0,numCenters*Time); locmax[maxloc] <-1; locmax <- matrix(locmax,ncol=1)
        #pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
        #partition weights vector st for pclocmax=1, those weights sum to 1
        Lik[which(pclocmax!=0)] <-1
        #reweight whats in the cluster so that sums to 1
        wi[which(pclocmax!=0)] <- likweights(Lik[which(pclocmax!=0)])
        
        out <- t(wi)%*%Lambda_dense
        locLambdas[[i]] <- out 
        #e) Set Lik for overlapping locations to zero
        Lik[which(pclocmax!=0)] <-0
        #f) Recalculate scaled likelihoods
        wi <- likweights(Lik)
    }
    return(locLambdas=locLambdas)
}
    
#'@title detectclusters
#'@description TODO
#'@param sparseMAT Large sparsematrix TODO
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param numCenters Number of centroids.
#'@param Time Number of time periods.
#'@param maxclust Maximum number of clusters allowed. TODO - allow this to be unknown.
#'@param bylocation If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{TRUE}.
#'@ return TODO
detectclusters <- function(sparseMAT, Ex, Yx,numCenters,Time, maxclust,bylocation=TRUE){
    #sparsemat large sparsematrix where rows are potential clusters and columns are space-time locations
    #Ex expected counts
    #Yx observed counts
    #numCenters number of centroids
    #Time number of time periods
    #maxclust maximum number of clusters to search for
    #bylocation if TRUE, clusters identified by location; if FALSE clusters identified by potential cluster
    
    #output empties
    sparsemat <- Matrix::t(sparseMAT) #66870x1040
    out <- poisLik(Ex, Yx, sparsemat)
    Lik <- out$Lik
    Lambda_dense <- out$Lambda_dense
    #outLik <- out$outLik
    #lambdahat <- out$lambdahat
    #Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
    locLambdas <- vector("list", maxclust)
    if(bylocation==FALSE){
        message("Cluster detection by potential cluster")
        res <- bycluster(outLik, sparsemat, locLambdas, Lambda_dense, maxclust)
        return(res)
    }
    else{
        #default
        #outLik, sparsemat, locLambdas, Lambda,maxclust
        message("Cluster detection by location")
        res <- bylocation(outLik, sparsemat, locLambdas, Lambda_dense, maxclust) #Lik, sparsemat, locLambdas, Lambda_dense, maxclust)
        return(res)
    }
}







