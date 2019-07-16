# clustack functions
#Load Libraries
library(clusso)
library(testthat)
library(Matrix)

########################################
#Functions
########################################

colorsmk <- function (x) {
    y = colorRamp(RColorBrewer::brewer.pal(9, "Greys")[1:9])(x)
    rgb(y[, 1], y[, 2], y[, 3], maxColorValue = 255)
}

plotmap <- function(res, pdfname=NULL){
    #cluster_ix <- redblue(log(2 *  pmax(1/2, pmin(res, 2)))/log(4))
    colmax <- 4;cluster_ix <- colorsmk(log(colmax *  pmax(1/colmax, pmin(res, colmax)))/log(colmax^2))
    
    colors_0 <- matrix(cluster_ix, ncol=5, byrow = FALSE)
    #pdf(pdfname, height=11, width=10)
    
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,1] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 1"),cex=1.00)
    
    #P2
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,2] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 2"),cex=1.00)
    
    #P3
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,3] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 3"),cex=1.00)
    
    #P4
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,4] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 4"),cex=1.00)
    
    #P5
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colors_0[,5] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 5"),cex=1.00)
    
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    #rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=pal(50),border=F)
    
    #rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=colorsmk(0:50/50),border=F)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(0,4,length.out=5),srt=330,adj=0)
    
    #dev.off()
}

likweights <- function(liki){
    # minmod <- min(liki, na.rm=TRUE)
    # print(minmod)
    # deltai <- liki-minmod
    wi <- liki/sum(liki)
    if(any(is.na(liki))){
        print("There are NA likelihoods")
    }
    return(wi)
}

poisLik <-function(Ex, Yx, sparsemat){
    #outLambda observed over expected;exponentiated rates
    #Ex expected counts
    #Yx observed counts
    #sparsemat large sparsematrix where rows are potential clusters and columns are space-time locations
    #return likelihood for each potential cluster
    outExp <- sparsemat%*%Ex
    outObs <- sparsemat%*%Yx
    lambdahat <- ifelse(outObs==0, 0, outObs/outExp)
    outlogLik <- outObs * log(lambdahat * outExp) - (lambdahat * outExp)
    if(any(is.na(outlogLik))){
        ix <- which(is.na(outlogLik))
        #print(ix)
        outlogLik[ix] <- (lambdahat[ix]*outExp[ix])
    }
    outlogLik_scaled <- outlogLik-max(outlogLik)
    outLik <- exp(outlogLik_scaled)
    return(list(outLik=outLik,
                lambdahat=lambdahat))
}

bylocation <- function(outLik, sparsemat, locLambdas, Lambda,maxclust){
    for(i in 1:maxclust){
        message(paste0("Searching for cluster ",i))
        #01) Get weight for each PC
        wi <- likweights(outLik)
        #1) Find LOCATION with highest weight
        lwi <- t(wi)%*%Lambda
        maxlwi <- which.max(lwi@x) #552
        message(paste0("Location identified: ",(maxlwi)))
        #2) Find All PC's that Cover the max location
        locmax <- rep(0,numCenters*Time); locmax[maxlwi] <-1
        findoverlap <- sparsemat%*%locmax
        ix <- which(findoverlap!=0)
        wxi <- likweights(outLik[ix])
        wi[ix] <- wxi
        #4) Find weighted Lambda by locations
        locLambdas[[i]] <- t(wi)%*%Lambda #locLambdas[[i]] <- lwi
        #5) Set Likelihood in cluster i to zero
        outLik[ix] <- 0
        
    }
    return(locLambdas = locLambdas)
}

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
    outLik <- out$outLik
    lambdahat <- out$lambdahat
    Lambda <- lambdahat*sparsemat #big Lambda matrix
    locLambdas <- vector("list", maxclust)
    if(bylocation==FALSE){
        message("Cluster detection by potential cluster")
        res <- bypotentialcluster(outLik, sparsemat, Lambda, locLambdas,outLambda, maxclust)
        return(res)
    }
    else{
        #default
        #outLik, sparsemat, locLambdas, Lambda,maxclust
        message("Cluster detection by location")
        res <- bylocation(outLik, sparsemat, locLambdas, Lambda, maxclust)
        return(res)
    }
}







