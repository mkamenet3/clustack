# clustack functions
#Load Libraries
library(clusso)
library(testthat)
library(Matrix)


########################################
#Functions
########################################
plotmeanrr_stack <- function(ric, Time, sim.i,ic, flav){
    #color.ic <- sapply(1:Time, function(i) redblue(log(2 *pmax(1/2, pmin(ric[, i], 2)))/log(4)))
    color.ic <- sapply(1:Time, function(i) redblue(log(1.1 *pmax(1/1.1, pmin(ric[, i], 1.1)))/log(1.1^2)))
    pdf(paste0(sim.i,"_meanrr_pc_",flav, "_",ic,".pdf"), height=11, width=10)
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=color.ic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,paste0(ic," - ", flav),cex=1.00, srt=90)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=color.ic[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=color.ic[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=color.ic[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=color.ic[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    #rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(0.9,1.1,length.out=5),srt=330,adj=0)
    
    dev.off()
}

clusso_prob_clusteroverlap <- function(sparseMAT,lassoresult,selected,rr, risk.ratio,nsim,Time, ncentroids, pow){
    #DEFINE TRUTH
    if(risk.ratio==1){
        warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk.ratio,1,0)    
    }
    #Take out the time vectors - only keep cluster part of matrix
    sparseMAT_clusteronly <- sparseMAT[,-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    #select out my betas for each sim
    betaselect <- lassoresult$coefs.lasso.all[,selected]
    #binarize
    betaselect_bin <- ifelse(abs(betaselect >= 10e-6),1,0)
    #only take the clusters betas
    betaselect_bin_clusteronly <- betaselect_bin[-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    ##INCLUSTER
    if(pow==TRUE){
        #print("pow")
        clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
        incluster_sim <- clusteroverlap_bin %*% betaselect_bin_clusteronly
        incluster_sim_bin <- ifelse(incluster_sim !=0,1,0)
        return(idin = incluster_sim_bin)
    }
    else{
        ##OUTCLUSTER
        #print("fp")
        clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
        outcluster_sim <- clusteroverlap_bin %*% betaselect_bin_clusteronly
        outcluster_sim_bin <- ifelse(outcluster_sim !=0,1,0)
        return(idout = outcluster_sim_bin)
    }
}




colorsgrey <- function (x) {
    y = colorRamp(RColorBrewer::brewer.pal(9, "Greys")[1:9])(x)
    rgb(y[, 1], y[, 2], y[, 3], maxColorValue = 255)
}

#'@title probplotmapAllIC
#'@param res.bic result bic
#'#'@param res.aic result aic
#'@param oracle
#'@param pdfname String for name of pdf to be generated.
#'@param genpdf Boolean. If results should be generated to console, set to \code{FALSE}. Default is \code{TRUE} for pdf to be generated.
#'@return Maps for central region of Japan for each time period.
probplotmapAllIC <- function(colprob, pdfname=NULL, genpdf=TRUE, obs=NULL){
    if(!is.null(obs)){
        firstrow = "Observed"
    }
    else{
        firstrow="Oracle"
    }
    if(genpdf==TRUE){
        pdf(pdfname, height=11, width=10)    
    }
    ####################################
    #oracle or Obs
    ####################################
    #P1
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[1]][,1] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,paste0(firstrow),cex=1.00, srt=90)
    
    #P2
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[1]][,2] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P3
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[1]][,3] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P4
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[1]][,4] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P5
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[1]][,5] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    ####################################
    #BIC
    ####################################
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[2]][,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,'BIC',cex=1.00, srt=90)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[2]][,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 2 - BIC',cex=1.00)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[2]][,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 3 - BIC',cex=1.00)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[2]][,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 4 - BIC',cex=1.00)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[2]][,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 5 - BIC',cex=1.00)
    
    ####################################
    #AIC
    ####################################
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[3]][,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,'AIC',cex=1.00, srt=90)
    
    par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[3]][,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 2 - AIC',cex=1.00)
    
    par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[3]][,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 3 - AIC',cex=1.00)
    
    par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[3]][,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 4 - AIC',cex=1.00)
    
    par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=colprob[[3]][,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 5 - AIC',cex=1.00)
    
    
    
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=greys(0:50/50),border=F)
    text(seq(.6,1.4,length=6),rep(.45,6),c('0.0','0.2','0.4','0.6','0.8','1.0'),srt=330,adj=0)
    #text(seq(0.6, 1.4, length = 5), rep(0.45, 5), seq(0.5, 1.5, length.out = 5), srt = 330, adj = 0)
    
    if(genpdf==TRUE){
        dev.off()    
    }
    
}

#'@title plotmapAllIC
#'@param res.bic result bic
#'#'@param res.aic result aic
#'@param oracle
#'@param pdfname String for name of pdf to be generated.
#'@param genpdf Boolean. If results should be generated to console, set to \code{FALSE}. Default is \code{TRUE} for pdf to be generated.
#'@param maxrr For the color ramp, what is the maximum relative risk color. Default is for the ramp to be between 0 and 2. 
#'@param minrr For the color ramp, what is the minimum relative risk color. Default is for the ramp to be between 0 and 2. 
#'@return Maps for central region of Japan for each time period.
plotmapAllIC <- function(res.bic, res.aic, oracle ,pdfname=NULL, genpdf=TRUE, maxrr=2, minrr=0, obs=NULL){
    if(!is.null(obs)){
        firstrow = "Observed"
    }
    else{
        firstrow="Oracle"
    }
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
    
    cluster_ix.bic <- matrix(redblue(log(maxrr *  pmax(minrr, pmin(res.bic, maxrr)))/log(maxrr^2)), ncol=5, byrow=FALSE)
    cluster_ix.aic <- matrix(redblue(log(maxrr *  pmax(minrr, pmin(res.aic, maxrr)))/log(maxrr^2)), ncol=5, byrow=FALSE)
    if(!is.null(obs)){
        oracle_ix <- matrix(redblue(log(2*pmax(1/2,pmin(oracle,2)))/log(4)), ncol=5, byrow=FALSE)
    } else {
        oracle_ix <- matrix(redblue(log(maxrr *  pmax(minrr, pmin(oracle, maxrr)))/log(maxrr^2)), ncol=5, byrow=FALSE)    
    }
    
    #colors_0 <- matrix(cluster_ix, ncol=5, byrow = FALSE)
    if(genpdf==TRUE){
        pdf(pdfname, height=11, width=10)    
    }
    ####################################
    #oracle or Obs
    ####################################
    #P1
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=oracle_ix [,1] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,paste0(firstrow),cex=1.00, srt=90)
    
    #P2
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=oracle_ix[,2] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P3
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=oracle_ix[,3] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P4
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=oracle_ix[,4] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    
    #P5
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=oracle_ix[,5] ,border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    ####################################
    #BIC
    ####################################
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.bic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,'BIC',cex=1.00, srt=90)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.bic[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 2 - BIC',cex=1.00)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.bic[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 3 - BIC',cex=1.00)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.bic[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 4 - BIC',cex=1.00)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.bic[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 5 - BIC',cex=1.00)
    
    ####################################
    #AIC
    ####################################
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.aic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(270,4120,'AIC',cex=1.00, srt=90)
    
    par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.aic[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 2 - AIC',cex=1.00)
    
    par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.aic[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 3 - AIC',cex=1.00)
    
    par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.aic[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 4 - AIC',cex=1.00)
    
    par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=cluster_ix.aic[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    #text(355,4120,'Period 5 - AIC',cex=1.00)
    
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

####NOTE: this overdisp function is part of clusso. Omit after testing
overdisp <- function(offset_reg, sim = TRUE, overdispfloor = TRUE) {
    if(sim==TRUE){
        stopifnot(inherits(offset_reg[[1]], c("glm", "lm")))
        phi <- max(unlist(lapply(1:length(offset_reg), function(i) deviance(offset_reg[[i]])/df.residual(offset_reg[[i]]))))
    }
    else{
        stopifnot(inherits(offset_reg, c("glm", "lm")))
        phi <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
    }
    if(overdispfloor == TRUE & phi < 1){
        message(paste("Underdispersion detected (", phi,"). Setting phi to 1"))
        phi <- 1
    }
    if(overdispfloor == FALSE & phi < 1){
        message(paste("Underdispersion detect (", phi,").\n 'Floor' argument was set to FALSE to underdispersed model will be run"))
    }
    return(phi)
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
    outlogLik <- outObs*log(outObs/outExp) + (outExp-outObs)
    outlogLik <- ifelse(is.na(outlogLik@x),0, outlogLik@x) #likelihood should be zero here
    outlogLik_scaled <- outlogLik-max(outlogLik)
    Lik <- exp(outlogLik_scaled)
    return(list(Lik=Lik,
                Lambda_dense=Lambda_dense))
}

#'@title binomLik
#'@description Poisson-based likelihood
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@export
#'@return Likelihood for each potential cluster.
binomLik <-function(nx, Yx, sparsemat){
    #outExp <- sparsemat%*%Ex
    outnx <- sparsemat%*%nx
    outObs <- sparsemat%*%Yx
    #calc Lambda
    lambdahat <- outObs/outExp
    Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
    Lambda_dense <- as.matrix(Lambda)
    Lambda_dense[Lambda_dense == 0] <- 1
    #Get scaled likelihood
    
    #Lik <- (((outObs/nx)/(sum(outObs)/sum(nx)))^(outobs))*(((sum(outObs)-outObs)/(sum(nx)-nx))/(sum(outObs)/sum(nx)))^(sum(outObs)-outObs)
    #((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs #TODO CHECK THIS
    outnt <- sum(outnx)
    outObst <- sum(outObs)
    Lik <- (((outObs/outnx)^outObs)*(1-(outObs/outnx))^(outnx-outObs))*(((outObst - outObs)/(outnt-outnx))^(outObst-outObs))*(1-((outObst-outObs)/(outnt - outnx)))^(outnt - outnx - outObst+outObs)
    
    outlogLik <- log(Lik)
    outlogLik_scaled <- outlogLik-max(outlogLik)
    Lik <- exp(outlogLik_scaled)
    return(list(Lik=Lik,
                Lambda_dense=Lambda_dense))
}

#'@title bylocation
#'@description Cluster detection based on maximum location.
#'@param Lik Likelihood for each potential cluster.
#'#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param maxclust Maximum number of clusters to be detected.
#'@return List. First element is a Weighted relative risks for identified locations. Second element is a large matrix of weights for each iteration.
#'

bylocation <- function(Lik, Lambda_dense,sparsemat, maxclust){
    wtMAT <- matrix(rep(NA, maxclust*nrow(sparsemat)), ncol=maxclust)
    wt0 <- likweights(Lik)
    ixall <- NULL
    ix <- NULL
    maxlocs <- rep(NA, maxclust)
    stop <- FALSE
    for(i in 1:maxclust){
        wtmp <-wt0
        wtmp[ixall] <-0
        #Find maxloc
        wi_loc <- t(wtmp)%*%sparsemat
        maxloc <- which.max(as.vector(wi_loc))
        message(paste0("Location identified: ",(maxloc)))
        maxlocs[i] <- maxloc 
        #find all potential clusters that overlap that location
        locmax <- rep(0,numCenters*Time); 
        locmax[maxloc] <-1;
        locmax <- matrix(locmax,ncol=1)
        pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
        ix <- which(pclocmax!=0) #indices of all PCs that overlap max location
        #upweight Lik* to 1 in all Pcs that overlap max cluster
        #reweight everything inside cluster to sum to 1
        wtmp[ix] <- likweights(Lik[ix])
        #save wtmp vector into master weight matrix
        wtMAT[,i] <- wtmp
        if(all(ix %in% ixall)){
            stop <- TRUE
            maxi = i
            break
        }
        if(stop){break}
        #save everything that's already been considered
        ixall <- c(ixall, ix)
        maxi = i 
        #Lik[ixall] <-0
    }
    wLambda <- crossprod(wtMAT[,1:(maxi)], Lambda_dense)
    return(list(wLambda = wLambda,
                wtMAT = wtMAT,
                #wtMAT0 = wtMAT0,
                maxlocs = maxlocs))
}




#'@title bycluster
#'@description Cluster detection based on maximum potential cluster.
#'@param Lik Likelihood for each potential cluster.
#'#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param maxclust Maximum number of clusters to be detected.
#'@return Weighted relative risks for identified locations.
bycluster <-  function(Lik, Lambda_dense, sparsemat,maxclust){
    wtMAT <- matrix(rep(NA, maxclust*nrow(sparsemat)), ncol=maxclust)
    wtMAT0 <- matrix(rep(NA, maxclust*nrow(sparsemat)), ncol=maxclust)
    maxpcs <- rep(NA, maxclust)
    stop = FALSE
    ixall <- NULL
        for (i in 1:maxclust){
            Lik[ixall] <-0
            wtmp <- likweights(Lik)#wt0
            wtmp[ixall] <-0
            #find potential cluster with largest weight
            maxpc <- which.max(as.vector(wtmp))
            #maxloc <- which.max(as.vector(wi_loc))
            message(paste0("Potential cluster identified: ",(maxpc)))
            maxpcs[i] <- maxpc 
            
            wtMAT0[,i] <- wtmp
            #find all potential clusters that overlap that PC
            pcmax <- rep(0,length(wtmp)); pcmax[maxpc] <-1; pcmax <- matrix(pcmax,ncol=1)
            pcpcmax <- t(t(sparsemat)%*%pcmax)%*%t(sparsemat); pcpcmax <- ifelse(pcpcmax!=0,1,0)
            ix <- which(pcpcmax!=0)
            #upweight Lik* to 1 in all PCs that overlap max PC
            #reweight so that everything inside the PCs that overlap max cluster sum to 1
            #aa <- likweights(Lik[ix]) 
            wtmp[ix] <- likweights(Lik[ix]) 
            wtMAT[,i] <- wtmp
            if(all(ix %in% ixall)){
                stop = TRUE
                maxi = i
                break
            }
            if (stop){break}
            #save everything that's already been considered
            ixall <- c(ixall,ix)
            maxi = i
        }
    wLambda <- crossprod(wtMAT[,1:(maxi)], Lambda_dense)
    return(list(wLambda = wLambda,
                wtMAT = wtMAT,
                #wtMAT0 = wtMAT0,
                maxpcs = maxpcs))
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
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#'@param overdisp.est Overdispersion parameter estimated across all simulations (max).
#'@return Returns list for each iteration with weighted relative risks by location inside identified cluster.
#'@export
detectclusters <- function(sparseMAT, Ex, Yx,numCenters,Time, maxclust,bylocation=TRUE, model=c("poisson", "binomial"),
                           cv=FALSE,  overdisp.est) {
    if(is.null(overdisp.est)){
        quasi <- FALSE
        message(paste0("Model specified: ", model))
    } else{ quasi <- TRUE
    #print("Quasi model")
    message(paste0("Model specified: ", "quasi-",model))
    }
    sparsemat <- Matrix::t(sparseMAT) 
    if (model=="poisson"){
        out <- poisLik(Ex, Yx, sparsemat)  
        #out <- poisLik(Ex[[1]], YSIM[[1]], sparsemat)  
    }
    else if (model=="binomial"){
        out <- binomLik(Ex, Yx, sparsemat) #in dev
    }
    else {
        stop("Model not specified. Please indicate `poisson` or `binomial`.")
    }
    #message(paste0("Model specified: ", model))
    Lik <- out$Lik
    Lambda_dense <- out$Lambda_dense
    #locLambdas <- vector("list", maxclust)
    if(bylocation==FALSE){
        print("Cluster detection by potential cluster")
        res <- bycluster(Lik, Lambda_dense, sparsemat, maxclust)
        #perform selection by IC/CV
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, quasi,cv=FALSE,overdisp.est)
        return(list(wLambda = res[[1]],
                    loglik = selection$loglik,
                    selection.bic_orig = selection$select.bic,
                    selection.aic_orig = selection$select.aic,
                    selection.aicc_orig = selection$select.aicc,
                    selection.bic = ifelse(selection$select.bic==0,1, selection$select.bic),
                    selection.aic = ifelse(selection$select.aic==0,1, selection$select.aic),
                    selection.aicc = ifelse(selection$select.aicc==0,1,selection$select.aicc),
                    #sparsemat = res[[2]],
                    wtMAT = res[[2]],
                    #wtMAT0 = res[[3]],
                    maxpcs = res[[3]],
                    Lambda_dense = Lambda_dense))#,
                    #Lambda_sparse = res[[3]]))
    }
    else{
        #default
        print("Cluster detection by location")
        res <- bylocation(Lik, Lambda_dense, sparsemat, maxclust)
        #perform selection by IC/CV
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, quasi,cv=FALSE,overdisp.est)
        return(list(wLambda = res[[1]],
                    loglik = selection$loglik,
                    selection.bic_orig = selection$select.bic,
                    selection.aic_orig = selection$select.aic,
                    selection.aicc_orig = selection$select.aicc,
                    selection.bic = ifelse(selection$select.bic==0,1,selection$select.bic),
                    selection.aic = ifelse(selection$select.aic==0,1, selection$select.aic),
                    selection.aicc = ifelse(selection$select.aic==0,1,selection$select.aicc),
                    #sparsemat = res[[2]],
                    wtMAT = res[[2]],
                    maxlocs = res[[3]],
                    Lambda_dense = Lambda_dense))#,
                    #wtMAT0 = res[[3]]))#,
                    #Lambda_sparse = res[[3]]))
    }
}

#'@title clusterselect
#'@description Select the optimal number of (overlapping) clusters using either information criteria or cross-validation
#'@param wLambda 
#'@param maxclust TODO
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param numCenters TODO
#'@param Time TODO
#'@param quasi Boolean. \code{TRUE} indicates a quasi-Poisson or quasi-binomial model that accounts for overdispersion. \code{FALSE} indicates a Poisson or binomial model without adjustment for overdispersion.
#'@param covars Dataframe of additional covariates to be included in the model that are un-penalized by the LASSO.
#'@param Yx Number of observed cases for each space-time location.
#'@param cv TODO
#'@param overdisp.est Overdispersion parameter estimated across all simulations (max).
#'@return TODO
clusterselect <- function(wLambda,Yx, Ex, model,maxclust, numCenters, Time,quasi,cv=FALSE,overdisp.est){
    #preprocess wLambda
    remtab <- which(sapply(1:nrow(wLambda), function(k) all(is.nan(wLambda[k,])))==FALSE)
    
    #test <- t(a)%*%Lambda_dense
    loglik <- sapply(1:length(remtab), function(i) dpoisson(Yx, wLambda[i,], Ex))
    #loglik <- sapply(1:nrow(wLambda), function(i) dpoisson(YSIM[[1]], wLambda[i,], Ex[[1]]))
    #add null model loglik
    if (model=="poisson"){
        loglik <- c(dpoisson(Yx, rep(1,length(Yx)), Ex), loglik)    
        #loglik <- c(dpoisson(YSIM[[1]], rep(1,length(YSIM[[1]])), Ex[[1]]), loglik)
        if (quasi== TRUE){
            message(paste("Overdispersion estimate:", round(overdisp.est,4)))
            if(quasi == TRUE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
        } else {
            print("no overdispersion being estimated")
        }
        
    }
    # else if (model=="binomial"){
    #     #TODO
    #     #loglik <- c(dbin(Yx, rep(1,length(Yx)), Ex), loglik)
    #     # if (quasi== TRUE){
    #     #     #TODO
    #     # } else {
    #     #     #TODO
    #     # }
    # }
    
    #find K clusters
    K <- seq(from=0, to=nrow(wLambda))#seq(from=0, to=maxclust)
    if(cv==TRUE){
        #TODO
    }
    else{
        #calc
        if (quasi== TRUE){
            PLL.bic <- (-2*loglik/overdisp.est) + (K+1)*log(sum(Yx)) #(-2*loglik) + K*log(numCenters*Time)
            PLL.aic <- 2*(K+1) - 2*(loglik/overdisp.est)
            PLL.aicc <- 2*(K+1) - 2*(loglik/overdisp.est) +
                ((2*(K+1)*(K+1 + 1))/(sum(Yx) - K+1 - 1))
        } else {
            PLL.bic <- (-2*loglik) + K*log(sum(Yx)) #(-2*loglik) + K*log(numCenters*Time)
            PLL.aic <- 2*K - 2*loglik
            PLL.aicc <- 2*(K) - 2*(loglik) +
                ((2*K*(K + 1))/(sum(Yx) - K - 1))
        }
        
        #make sure that aicc is not overparameterized
        if(any(is.infinite(PLL.aicc))) {
            idx <- which(is.infinite(PLL.aicc))
            #replace everything in penalized loglik at infinite value and onward to be NA - these are nonsense solutions
            PLL.aicc[idx:length(PLL.aicc)] <- NA
        }
        #select
        select.bic <- which.min(PLL.bic)-1 #-1 because the first element corresponds to 0 clusters (null)
        select.aic <- which.min(PLL.aic)-1
        select.aicc <- which.min(PLL.aicc)-1 
        cat(paste0("\t","Num. selected clusters by BIC: ", (select.bic),"\n",
                   "\t","Num. selected clusters by AIC: ",(select.aic),"\n",
                   "\t","Num. selected clusters by AICc: ",(select.aicc)))
        # if(all.equal(select.bic, select.aic, select.aicc)){
        #     #TODO
        # }
        # else{
        #     #Do it for each criterion
        #    #TODO
        #     #BIC
        #    #TODO
        #     #AIC
        #    #TODO    
        #     #AICc
        #    #TODO
        #     # return(list())
        # }
        return(list(loglik = loglik,
                    select.bic = select.bic,
                    select.aic = select.aic,
                    select.aicc = select.aicc))
    }
    
}

####################################################################
#CLUSTACKBOUNDS FUNCTIONS
####################################################################

bucklandbounds <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est) {
    NT <- rowSums(sparsematrix)
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/NT[k])
    } else {
        varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/NT[k])
    }
    withintheta <- (thetai - thetaa)^2
    varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))
    var_thetaa <- as.vector(varthetas_w)
    UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
    LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
    return(list(buckland.LB = LBa,
                clusterMA = thetaa,
                buckland.UB = UBa))
}

maw1 <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est) {
    NT <- rowSums(sparsematrix)
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/NT[k])
    } else {
        varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/NT[k])
    }
    withintheta <- (thetai - thetaa)^2
    se_thetaa <- sum(w_q*sqrt(varthetai+withintheta))
    var_thetaa <- as.vector(se_thetaa)
    UBa = as.vector(thetaa) + 1.96*(se_thetaa)
    LBa = as.vector(thetaa) - 1.96*(se_thetaa)
    return(list(maw1.LB = LBa,
                clusterMA = thetaa,
                maw1.UB = UBa))
}

maw2 <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est) {
    NT <- rowSums(sparsematrix)
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/NT[k])
    } else {
        varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/NT[k])
    }
    withintheta <- (thetai - thetaa)^2
    var_thetaa <- sum(w_q*(varthetai + withintheta))
    var_thetaa <- as.vector(varthetaa)
    UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
    LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
    return(list(maw2.LB = LBa,
                clusterMA = thetaa,
                maw2.UB = UBa))
}


mata_tailareazscore <- function(thetaii, thetaaa, se.thetaii, w_q, alpha){
    thetaii <- as.vector(thetaii)
    zval <- (thetaaa - thetaii)/se.thetaii
    zpnorm <- pnorm(zval)
    w_zpnorm <- sum((w_q*zpnorm))-alpha
    
}

matabounds <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, transform=c("none","log", "sqrt")) {
    NT <- rowSums(sparsematrix)
    switch(transform,
           none = matabounds_none(thetai,thetaa, w_q,sparsematrix, overdisp.est, NT),
           log = matabounds_log(thetai,thetaa, w_q,sparsematrix, overdisp.est, NT),
           sqrt = matabounds_sqrt(thetai,thetaa, w_q,sparsematrix, overdisp.est, NT))
}

matabounds_none <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, NT) {
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/NT[k])
    } else {
        varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/NT[k])
    }
    mataLB <- uniroot(f=mata_tailareazscore, interval=c(-1, 3),
                      thetaii= thetai,
                      se.thetaii=sqrt(varthetai),
                      w_q=w_q, alpha=0.025, tol=1e-8)$root    
    
    mataUB <- uniroot(f=mata_tailareazscore, interval=c(-1, 3),
                      thetaii= thetai,
                      se.thetaii=sqrt(varthetai),
                      w_q=w_q, alpha=1-0.025, tol=1e-8)$root 
    return(list(mata.LB = mataLB,
                clusterMA = thetaa,
                mata.UB = mataUB))
}

matabounds_sqrt <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, NT) {
    Tvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(4*NT[k]))
    mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                      thetaii= sqrt(thetai),
                      se.thetaii=sqrt(Tvarthetai),
                      w_q=w_q, alpha=0.025, tol=1e-8)$root

    mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                      thetaii= sqrt(thetai),
                      se.thetaii=sqrt(Tvarthetai),
                      w_q=w_q, alpha=1-0.025, tol=1e-8)$root
    return(list(matasqrt.LB = (mataLB)^2,
                clusterMA = thetaa,
                matasqrt.UB = (mataUB)^2))
}

matabounds_log <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, NT) {
    logTvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*NT[k]))
    mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                      thetaii= log(thetai),
                      se.thetaii=sqrt(logTvarthetai),
                      w_q=w_q, alpha=0.025, tol=1e-8)$root

    mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                      thetaii= log(thetai),
                      se.thetaii=sqrt(logTvarthetai),
                      w_q=w_q, alpha=1-0.025, tol=1e-8)$root
    return(list(matalog.LB = exp(mataLB),
                clusterMA = thetaa,
                matalog.UB = exp(mataUB)))
}

selectuniqRR <- function(uniqRRs){
    clusterRR_i <- rep(NA, nrow(uniqRRs))
    flag1 <- sapply(1:nrow(uniqRRs), function(k) length(unique(uniqRRs[k,])))
    non1s <- uniqRRs[which(flag1==2),]
    clusterRR_i[which(flag1==2)] <- sapply(1:nrow(non1s), function(k) non1s[k,which(non1s[k,]!=1)])
    clusterRR_i[which(flag1==1)] <- 1    
    return(unlist(clusterRR_i))
    
}







