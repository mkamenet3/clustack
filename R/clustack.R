# clustack functions
#Load Libraries
library(clusso)
library(testthat)
library(Matrix)


########################################
#Functions
########################################
create_plotmeanrr_stack <- function(sim_superclust,IC, flav,Time, nsim, sim.i){
    if(IC=="aic"){
        selects <- sapply(1:nsim, function(i) sim_superclust[[i]]$selection.aic)    
    } else{
        selects <- sapply(1:nsim, function(i) sim_superclust[[i]]$selection.bic)
    }
    #browser()
    if(all(selects==0)){
        ric <- matrix(rep(1, ncol(sim_superclust[[1]]$Lambda_dense)), ncol=Time)
    } else {
        meanrr <- lapply(1:nsim, function(i) sim_superclust[[i]]$wLambda[selects[i],])
        names(meanrr) <- paste0("l",1:nsim)
        meanrr[which(selects==0)] <- NULL
        meanrr.df <- do.call(rbind, meanrr)
        meanrrs <- apply(meanrr.df,2, mean)
        ric <- matrix(meanrrs, ncol = Time)
    }
    plotmeanrr_stack(ric, Time, sim.i,ic=IC, flav=flav)
}

create_plotFPR_stack <- function(sim_superclust,IC, flav,Time, nsim, sim.i){
    if(IC=="aic"){
        selects <- sapply(1:nsim, function(i) sim_superclust[[i]]$selection.aic)    
    } else{
        selects <- sapply(1:nsim, function(i) sim_superclust[[i]]$selection.bic)
    }
    #browser()
    if(all(selects==0)){
        #ric <- matrix(rep(1, ncol(sim_superclust[[1]]$Lambda_dense)), ncol=Time)
        probs <- rep(0, ncol(sim_superclust[[1]]$Lambda_dense))
    } else {
        #find pc's that overlap maxloc
        if(flav=="loc"){
            maxid <- sapply(1:nsim, function(i) sim_superclust_loc[[i]]$maxlocs[selects[i]])    
            vec <- rep(0, 208 * Time)
            position <- list(vec)[rep(1, nsim)]
            simindicator <- mapply(reval, position, maxid)
            probs <- Matrix::rowSums(simindicator)/nsim
        } else{
            maxid <- sapply(1:nsim, function(i) sim_superclust_loc[[i]]$maxpcs[selects[i]])   
            ixout <- vector(mode = "list", length = nsim)
            #find pcs
            for(i in 1:length(maxid)){
                pcmax <- rep(0,66870)
                pcmax[maxid[[i]]] <-1; pcmax <- matrix(pcmax,ncol=1)    
                identlocs<- sparsematrix%*%pcmax
                ixout[[i]] <- which(identlocs@x!=0)
                
            }
            vec <- rep(0, 208 * Time)
            position <- list(vec)[rep(1, nsim)]
            simindicator <- mapply(reval, position, ixout)
            probs <- Matrix::rowSums(simindicator)/nsim
        }
        
    }
    plotmeanrr_stack(ric=matrix(probs, ncol=Time), Time, sim.i,ic=IC, flav=flav, greys=TRUE)
}

plotmeanrr_stack <- function(ric, Time, sim.i,ic, flav, greys){
    if(greys==TRUE){
        color.ic <- sapply(1:Time, function(i) greys(ric[,i]))
        pdf(paste0(sim.i,"_fpr_",flav, "_",ic,".pdf"), height=11, width=10)
        
    } else {
        color.ic <- sapply(1:Time, function(i) redblue(log(2 *pmax(1/2, pmin(ric[, i], 2)))/log(4)))
        pdf(paste0(sim.i,"_meanrr_",flav, "_",ic,".pdf"), height=11, width=10)
    }
   
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
    if(greys==TRUE){
        rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=greys(0:50/50),border=F)
        text(seq(.6,1.4,length=5),rep(.45,5),seq(0,1,length.out=5),srt=330,adj=0)
    } else {
        rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
        text(seq(.6,1.4,length=5),rep(.45,5),seq(0.5,2,length.out=5),srt=330,adj=0)
    }
    
    
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
                    selection.bic = selection$select.bic,
                    selection.aic = selection$select.aic,
                    selection.aicc = selection$select.aicc,
                    selection.bic_forceid = ifelse(selection$select.bic==0,1, selection$select.bic),
                    selection.aic_forceid = ifelse(selection$select.aic==0,1, selection$select.aic),
                    selection.aicc_forceid = ifelse(selection$select.aicc==0,1,selection$select.aicc),
                    # selection.bic_orig = selection$select.bic,
                    # selection.aic_orig = selection$select.aic,
                    # selection.aicc_orig = selection$select.aicc,
                    # selection.bic = ifelse(selection$select.bic==0,1, selection$select.bic),
                    # selection.aic = ifelse(selection$select.aic==0,1, selection$select.aic),
                    # selection.aicc = ifelse(selection$select.aicc==0,1,selection$select.aicc),
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
                    selection.bic = selection$select.bic,
                    selection.aic = selection$select.aic,
                    selection.aicc = selection$select.aicc,
                    selection.bic_forceid = ifelse(selection$select.bic==0,1,selection$select.bic),
                    selection.aic_forceid = ifelse(selection$select.aic==0,1, selection$select.aic),
                    selection.aicc_forceid = ifelse(selection$select.aic==0,1,selection$select.aicc),
                    # selection.bic_orig = selection$select.bic,
                    # selection.aic_orig = selection$select.aic,
                    # selection.aicc_orig = selection$select.aicc,
                    # selection.bic = ifelse(selection$select.bic==0,1,selection$select.bic),
                    # selection.aic = ifelse(selection$select.aic==0,1, selection$select.aic),
                    # selection.aicc = ifelse(selection$select.aic==0,1,selection$select.aicc),
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

# 


####################################################################
#CLUSTACKBOUNDS FUNCTIONS
####################################################################

bucklandbounds <- function(thetai,thetaa, w_q,sparsematrix, outExp,overdisp.est, transform=NULL, cellrates=FALSE) {
    if(cellrates == TRUE){
        if(!is.null(transform)){
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/thetai[k,]*outExp)
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) sqrt(varthetai[,k] + withintheta[k,]))
            varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
            var_thetaa <- as.vector(varthetas_w)
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
            
        } else {
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
            }
            withintheta <- (thetai - thetaa)^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) sqrt(varthetai[,k] + withintheta[k,]))
            varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
            var_thetaa <- as.vector(varthetas_w)
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    } else {
        if(!is.null(transform)){
            #print("log-scale")
            #log transform
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k]*outExp[k])))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))^2
            var_thetaa <- as.vector(varthetas_w)
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
            
        } else {
            #print("No transform")
            #no transform
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
            }
            withintheta <- (thetai - thetaa)^2
            varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))
            var_thetaa <- as.vector(varthetas_w)
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    }
    
    return(list(buckland.LB = LBa,
                clusterMA = thetaa,
                buckland.UB = UBa))
}


maw2 <- function(thetai,thetaa, w_q,sparsematrix, outExp, overdisp.est, transform=NULL, cellrates=FALSE) {
    if(cellrates==TRUE){
        if(!is.null(transform)){
            #print("log-scale")
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k,]*outExp))
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) (varthetai[,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        } else{
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
            }
            withintheta <- (thetai - thetaa)^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) (varthetai[,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    } else {
        if(!is.null(transform)){
            #print("log-scale")
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k]*outExp[k])))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(var_thetaa)
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        } else{
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
            }
            withintheta <- (thetai - thetaa)^2
            var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(var_thetaa)
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    }
    return(list(maw2.LB = LBa,
                clusterMA = thetaa,
                maw2.UB = UBa))
}


mata_tailareazscore <- function(thetaii, thetaaa, se.thetaii, w_q, alpha){
    thetaii <- as.vector(thetaii)
    zval <- (thetaaa - thetaii)/se.thetaii
    zpnorm <- pnorm(zval)
   # print(str(zpnorm))
#    print(str(w_q))
    w_zpnorm <- sum((w_q*zpnorm))-alpha
    
}

matabounds <- function(thetai,thetaa, w_q,sparsematrix, outExp, overdisp.est,transform=c("none","log", "sqrt"), cellrates=FALSE) {
    #print(cellrates)
    #NT <- rowSums(sparsematrix)
    if(is.null(transform)){
        none = matabounds_none(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates)
    }
    switch(transform,
           none = matabounds_none(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates),
           log = matabounds_log(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp,cellrates),
           sqrt = matabounds_sqrt(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp,cellrates))
}

matabounds_none <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
    if(cellrates==TRUE){
        #browser()
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
        } else {
            varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
        }
        mataLB <- sapply(1:ncol(thetai), 
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai[,k],
                          se.thetaii=sqrt(varthetai[k,]),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root)    
        
        mataUB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai[,k],
                          se.thetaii=sqrt(varthetai[k,]),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root) 
    } else{
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
        } else {
            varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
        }
        mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai,
                          se.thetaii=sqrt(varthetai),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root    
        
        mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai,
                          se.thetaii=sqrt(varthetai),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root 
    }
    
    return(list(mata.LB = mataLB,
                clusterMA = thetaa,
                mata.UB = mataUB))
}

matabounds_sqrt <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
   # print("sqrt version")
    if(cellrates==TRUE){
    #    print("cellrates true")
        #print(str(outExp))
        #browser()
        Tvarthetai <-  1/(4*outExp)
        # test <- sapply(1:5, 
        #                function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
        #                                    thetaii= sqrt(thetai[,k]),
        #                                    se.thetaii=sqrt(Tvarthetai[k]),
        #                                    w_q=w_q, alpha=0.025, tol=1e-8)$root)
        
        #TODO: fix herel it should be thetai[,k] where k corresponds to the kth cell of cellsix, not 1:5
        mataLB <- sapply(1:ncol(thetai), 
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai[,k]),
                          se.thetaii=sqrt(Tvarthetai[k]),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root)
        
        mataUB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai[,k]),
                          se.thetaii=sqrt(Tvarthetai[k]),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root)
    } else{
        Tvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(4*outExp[k]))
        mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai),
                          se.thetaii=sqrt(Tvarthetai),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root
        
        mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai),
                          se.thetaii=sqrt(Tvarthetai),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root
    }
    return(list(matasqrt.LB = (mataLB)^2,
                clusterMA = thetaa,
                matasqrt.UB = (mataUB)^2))
}

matabounds_log <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
    if(cellrates==TRUE){
        logTvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k,]*outExp))
        mataLB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai[,k]),
                          se.thetaii=sqrt(logTvarthetai[k,]),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root)
        
        mataUB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai[,k]),
                          se.thetaii=sqrt(logTvarthetai[k,]),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root)
    } else {
        logTvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))   
        mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai),
                          se.thetaii=sqrt(logTvarthetai),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root
        
        mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai),
                          se.thetaii=sqrt(logTvarthetai),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root
    }
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

##############################################################################
#SIMULATION HELPER FUNCTIONS
##############################################################################

#nonma
nonma <- function(cluster_thetaa,sim_superclust_loc,clusterRR_ilarge,wslarge,idix, IC){
    #MA by maxloc
    #cluster_thetaa_locs <- lapply(1:nsim, function(i) 1)
    #cluster_thetaa <- lapply(1:length(idix), function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
    if(IC=="aic") {
        clusterRRlarge <- lapply(1:length(idix), function(j) unique(sim_superclust_loc[[idix[[j]]]]$Lambda_dense[,sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.aic]])[2])
        se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.aic]]))
        
    } else {
        clusterRRlarge <- lapply(1:length(idix), function(j) unique(sim_superclust_loc[[idix[[j]]]]$Lambda_dense[,sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.bic]])[2])
        se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.bic]]))
        
    }
    nonma.theta.time <- system.time(nonma.theta <- lapply(1:length(idix), function(i) cbind(lb=cluster_thetaa[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                                                    clusterMA = cluster_thetaa[[i]],
                                                                                    ub=cluster_thetaa[[i]]+1.96*se_clusterRRlarge[[i]])))
    return(list(nonma.theta.time = nonma.theta.time[[3]],
                nonma.theta = nonma.theta))
}


#nonma_asymp
nonma_asymp <- function(cluster_thetaa,sim_superclust_loc,clusterRR_ilarge,wslarge,idix, IC){
    #MA by maxloc
    #cluster_thetaa_locs <- lapply(1:nsim, function(i) 1)
    #cluster_thetaa <- lapply(1:length(idix), function(i) sum(clusterRR_ilarge[[i]]*wslarge[[i]]))
    if(IC=="aic") {
        clusterRRlarge <- lapply(1:length(idix), 
                                 function(j) unique(sim_superclust_loc[[idix[[j]]]]$Lambda_dense[,sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.aic]])[2])
        se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(clusterRRlarge[[j]]/outObs[[j]]@x[sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.aic]]))
        
    } else {
        clusterRRlarge <- lapply(1:length(idix), 
                                 function(j) unique(sim_superclust_loc[[idix[[j]]]]$Lambda_dense[,sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.bic_forceid]])[2])
        se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(clusterRRlarge[[j]]/outObs[[j]]@x[sim_superclust_loc[[idix[[j]]]]$maxlocs[sim_superclust_loc[[idix[[j]]]]$selection.bic_forceid]]))
        
    }
    nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- lapply(1:length(idix), function(i) cbind(lbasymp=cluster_thetaa[[i]]-1.96*se_clusterRRlarge_asymp[[i]], clusterMA = cluster_thetaa[[i]],ubasymp=cluster_thetaa[[i]]+1.96*se_clusterRRlarge_asymp[[i]])))
    return(list(nonma_asymp.theta.time = nonma_asymp.theta.time[[3]],
                nonma_asymp.theta = nonma_asymp.theta))
}


calcbounds <- function(id, IC, sim_superclust){
    #do all diagnostics
    idix <- which(id!=0)
    print(idix)
    print(paste0("calculating bounds for ", IC))
    #prep
    if(IC=="aic"){
        wslarge <- lapply(1:length(idix), function(j) sim_superclust[[idix[[j]]]]$wtMAT[,sim_superclust[[idix[[j]]]]$selection.aic])
    } else {
        wslarge <- lapply(1:length(idix), function(j) sim_superclust[[idix[[j]]]]$wtMAT[,sim_superclust[[idix[[j]]]]$selection.bic])
    }
    clusterRR_uniqlarge <- lapply(1:length(idix), function(j) sapply(1:nrow(sim_superclust[[idix[[j]]]]$Lambda_dense), 
                                                                     function(k) unique(sim_superclust[[idix[[j]]]]$Lambda_dense[k,]))) 
    
    clusterRR_ilarge <- lapply(1:length(idix), function(i) rep(NA, 66870))
    #clusterRR_uniq_ilarge <- lapply(1:length(idix), function(j) as.matrix(do.call(rbind, clusterRR_uniqlarge[[idix[[j]]]]), ncol=2))
    clusterRR_uniq_ilarge <- lapply(1:length(idix), function(j) as.matrix(do.call(rbind, clusterRR_uniqlarge[[j]]), ncol=2))
    #clusterRR_ilarge <- lapply(1:length(idix), function(j) selectuniqRR(clusterRR_uniq_ilarge[[idix[[j]]]]))
    clusterRR_ilarge <- lapply(1:length(idix), function(j) selectuniqRR(clusterRR_uniq_ilarge[[j]]))
    cluster_thetaa <- lapply(1:length(idix), function(j) sum(clusterRR_ilarge[[j]]*wslarge[[j]]))
    
    
    
    #Perform
    outnonma.time <- system.time(outnonma <- nonma(cluster_thetaa, sim_superclust, clusterRR_ilarge, wslarge, idix, IC=IC))
    outnonma_asymp.time <- system.time(outnonma_asymp <- nonma_asymp(cluster_thetaa ,sim_superclust, clusterRR_ilarge, wslarge, idix, IC=IC))
    print("nonma finished")
    outbuck.theta.time <- system.time(outbuck.theta <- lapply(1:length(idix), function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                 thetaa = cluster_thetaa[[i]], 
                                                                                                 w_q=wslarge[[i]], 
                                                                                                 sparsematrix=t(sparsematrix), 
                                                                                                 outExp[[i]],overdisp.est = NULL)))
    outbuck.theta.time <- system.time(outbuck.theta <- lapply(1:length(idix), function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]],
                                                                                                 thetaa = cluster_thetaa[[i]], 
                                                                                                 w_q=wslarge[[i]], 
                                                                                                 sparsematrix=t(sparsematrix), 
                                                                                                 outExp[[i]],overdisp.est = NULL)))
    outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- lapply(1:length(idix), function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]],
                                                                                                         thetaa =cluster_thetaa[[i]], 
                                                                                                         w_q=wslarge[[i]], 
                                                                                                         sparsematrix=t(sparsematrix),
                                                                                                         outExp[[i]],
                                                                                                         overdisp.est = NULL, 
                                                                                                         transform=TRUE)))
    print("buckland finished")
    outmaw2.theta.time <- system.time(outmaw2.theta <- lapply(1:length(idix), function(i) maw2(thetai=clusterRR_ilarge[[i]], 
                                                                                       thetaa = cluster_thetaa[[i]], 
                                                                                       w_q=wslarge[[i]],
                                                                                       sparsematrix=t(sparsematrix ), 
                                                                                       outExp[[i]],overdisp.est = NULL)))
    
    outmaw2Tlog.theta.time <- system.time(outmaw2Tlog.theta  <- lapply(1:length(idix), function(i) maw2(thetai=clusterRR_ilarge[[i]], 
                                                                                                thetaa = cluster_thetaa[[i]], 
                                                                                                w_q=wslarge[[i]], 
                                                                                                sparsematrix=t(sparsematrix), 
                                                                                                outExp[[i]], 
                                                                                                overdisp.est = NULL,
                                                                                                transform=TRUE)))
    print("maw2 finished")
    outmata.theta.time <- system.time(outmata.theta <- lapply(1:length(idix), function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                             thetaa = cluster_thetaa[[i]], 
                                                                                             w_q=wslarge[[i]], 
                                                                                             sparsematrix=t(sparsematrix ), 
                                                                                             outExp = outExp[[i]],
                                                                                             overdisp.est = NULL,
                                                                                             transform="none")))
    outmataT.theta.time <- system.time(outmataT.theta <- lapply(1:length(idix), function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                               thetaa = cluster_thetaa[[i]], 
                                                                                               w_q=wslarge[[i]], 
                                                                                               sparsematrix=t(sparsematrix ), 
                                                                                               outExp = outExp[[i]],
                                                                                               overdisp.est = NULL, 
                                                                                               transform="sqrt")))
    outmataTlog.theta.time <- system.time(outmataTlog.theta <- lapply(1:length(idix), function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                     thetaa = cluster_thetaa[[i]], 
                                                                                                     w_q=wslarge[[i]], 
                                                                                                     sparsematrix=t(sparsematrix ),
                                                                                                     outExp = outExp[[i]], 
                                                                                                     overdisp.est = NULL,
                                                                                                     transform="log")))
    print("mata finished")
    #return
    return(list(
        outnonma = outnonma,
        outnonma_asymp = outnonma_asymp,
        outbuck.theta = outbuck.theta,
        outbuckTlog.theta = outbuckTlog.theta,
        outmaw2.theta = outmaw2.theta,
        outmaw2Tlog.theta = outmaw2Tlog.theta,
        outmata.theta = outmata.theta,
        outmataT.theta = outmataT.theta,
        outmataTlog.theta = outmataTlog.theta,
        
        outnonma.time = outnonma.time[[3]],
        outnonma_asymp.time = outnonma_asymp.time[[3]],
        outbuck.theta.time = outbuck.theta.time[[3]],
        outbuckTlog.theta.time = outbuckTlog.theta.time[[3]],
        outmaw2.theta.time = outmaw2.theta.time[[3]],
        outmaw2Tlog.theta.time = outmaw2Tlog.theta.time[[3]],
        outmata.theta.time = outmata.theta.time[[3]],
        outmataT.theta.time = outmataT.theta.time[[3]],
        outmataTlog.theta.time = outmataTlog.theta.time[[3]] ))
    
}


##################
#Comparisons
##################

make.name.tree <- function(x, recursive, what.names)
{
    if (!is.character(what.names) || length(what.names) != 1)
        stop("'what.names' must be a single string")
    what.names <- match.arg(what.names, c("inherited" , "full"))
    .make.name.tree.rec <- function(x, parent_name, depth)
    {
        if (length(x) == 0)
            return(character(0))
        x_names <- names(x)
        if (is.null(x_names))
            x_names <- rep.int(parent_name, length(x))
        else if (what.names == "full")
            x_names <- paste0(parent_name, x_names)
        else
            x_names[x_names == ""] <- parent_name
        if (!is.list(x) || (!recursive && depth >= 1L))
            return(x_names)
        if (what.names == "full")
            x_names <- paste0(x_names, ".")
        lapply(seq_len(length(x)),
               function(i) .make.name.tree.rec(x[[i]], x_names[i], depth + 1L))
    }
    .make.name.tree.rec(x, "", 0L)
}

unlist2 <- function(x, recursive=TRUE, use.names=TRUE, what.names="inherited")
{
    ans <- unlist(x, recursive, FALSE)
    if (!use.names)
        return(ans)
    if (!is.character(what.names) || length(what.names) != 1)
        stop("'what.names' must be a single string")
    what.names <- match.arg(what.names, c("inherited" , "full"))
    names(ans) <- unlist(make.name.tree(x, recursive, what.names), recursive, FALSE)
    ans
}
dpoisson_fstage <- function(x, lambda, log = FALSE) {
    if(log == FALSE) 
        return(lambda^x * exp(-lambda)/factorial(x))
    else
        return(x*ifelse(lambda==0,1,log(lambda))-lambda-log(factorial(x)))
}
#this function returns poisson distribution if not log; otherwise, return the log of the likelihood function#
###not sure why this is written like this


f=function(x,y) x*log(y+(x==0))    ## auxiliary function ; ##what does this do????

scales=function(E.MAT, Y.MAT, Time) as.vector(sapply(1:Time, function(i) E.MAT[,i]*sum(Y.MAT[,i])/sum(E.MAT[,i])))


st_mat <- function(X) {
    U <- sapply(1:ncol(X),function(i) clusso::prod_YxCpp(X[,i], clusters$last, clusters$center))   #minimizes rows of columns from 20800 -->8960; #some sort of cluster identification?
    W=NULL
    for(i in 1:5) {
        v=numeric(nrow(U)) #creates empty vector of 0's that is 1:8960
        for(j in 0:(5-i))
        {v=v+U[,i+j] #running total sum of all the rows from U
        W=c(W,v)
        }
    }
    return(W)
}
#sums U/X matrix and puts it into one long vector; according to this algorithm, there are 15 iteractions of i and j (from length(W)/8960)
#8960 is somehow an optimized number of counts from 104000 --> widdled down

cluster_model <- function(delta,Y,E0,sd, Time) {
    Ys=matrix(Y,ncol=Time) #this is 20800 x 5 matrix (5 time periods)
    n=nrow(Ys) #num. rows of Ys = 20800
    
    j.path = numeric(max.ndelta) #1:8000 of 0's vector
    
    E.path = matrix(nrow=n*Time,ncol=max.ndelta+1) #logi [1:104000, 1:8001] NA NA NA NA NA NA
    #nrow=20800*5, ncol=8001
    
    yx <- st_mat(Ys) #yx = W defined above in the function
    
    E.cur=scales(matrix(E0,ncol=Time),Ys, Time)
    E.path[,1]=E.cur #replace first column with the E.cur data
    
    for(ndelta in 1:max.ndelta)
    {      
        Ex <- st_mat(matrix(E.cur,ncol=Time)) #num[1:134400]
        
        deri=(yx-Ex)/sd #standardized; num[1:134400]
        jmax=which.max(abs(deri)) #this selects the max of 'deri' which is the standardarized difference of real-expected values of 'x', which i think is deaths?
        #max value of our deri vector
        sign=sign(deri[jmax]) #extract original sign of this
        
        j.path[ndelta]=jmax #finds the max value at each increment from 1-8000
        
        
        j=jmax %% nrow(clusters) #find remainder from jmax (max of difference) and number of rows in clusters
        #nrow(clusters)=8960;how many aren't max
        l=ceiling(jmax/nrow(clusters)) #sets ceiling for (number of clusters?)
        #col=maxCol(n, j, clusters$n, clusters$last)
        col=clusso::max_colCpp(n, j, clusters$n, clusters$last)
        cols=rep(col,end[l]-start[l]+1) #not sure what this does; sets max col for cols?
        xj=numeric(Time*n) #Time=5 time periods; n=208 from clusteres dataset
        xj[(1+(start[l]-1)*n) : (end[l]*n)] <- cols
        
        E.cur=E.cur*exp(sign*delta*xj) #E.cur*e^sign*delta*xj
        E.cur=scales(matrix(E.cur,ncol=Time),Ys, Time) #rescales down 
        
        E.path[,ndelta+1]=E.cur
    }
    K <- sapply(1:max.ndelta, function(i) length(unique(j.path[1:i]))) 
    K <- c(0,K)
    
    loglike <- sapply(1:(max.ndelta+1),function(i) sum(dpoisson_fstage(Y,E.path[,i],log=T)))    #this sets the log-likelihood as a function of (i) across max.ndelta+1, and uses the poisson distrib specified above
    nsize = sum(Y)
    PLL.bic  <- -2*(loglike) + ((K)*log(nsize))
    #PLL.bic=loglike-log(n*Time)/2*K #calculautes best BIC
    #best.bic=which.max(PLL.bic)
    best.bic=which.min(PLL.bic)
    E.bic=E.path[,best.bic]
    
    PLL.aic <-  2*(K) - 2*(loglike)
    best.aic=which.min(PLL.aic)
    # PLL.aic=loglike-K #AIC
    # best.aic=which.max(PLL.aic)
    E.aic=E.path[,best.aic]
    
    #PLL.aicc=loglike-K*n*Time/(n*Time-K-1) #AICc
    PLL.aicc <- 2*(K) - 2*(loglike) +
        ((2*K*(K + 1))/(nsize - K - 1))
    best.aicc=which.min(PLL.aicc)
    #best.aicc=which.max(PLL.aicc)
    E.aicc=E.path[,best.aicc]
    
    #convert to Relative Risks
    RRobs <- matrix(as.vector(Y)/as.vector(E0),ncol=Time)
    RRbic <- matrix(E.bic/as.vector(E0),ncol=Time)
    RRaic <- matrix(E.aic/as.vector(E0),ncol=Time)
    RRaicc <- matrix(E.aicc/as.vector(E0),ncol=Time)
    
    
    return(list(RRobs = RRobs , RRbic = RRbic, RRaic = RRaic, RRaicc = RRaicc,
                E.bic = E.bic,E.aic= E.aic,E.aicc = E.aicc,
                n.bic=best.bic-1,n.aic=best.aic-1,n.aicc=best.aicc-1,
                PLL.bic,PLL.aic,PLL.aicc, 
                K=K))
    #returns bic,aic,aicc
}




stepscan <- function(Yx, Ex, Time,sparsematrix, nsim, maxclust){
    #stop <- FALSE
    #browser()
    mlc <- matrix(rep(NA,1040*maxclust), ncol = maxclust)
    mlc_pvals <- rep(NA, maxclust)
    maxLiks <- rep(NA, maxclust)
    #pval <- 0
    #while(pval <= 0.05){
    for(m in 1:maxclust){
        #identify first MLC
        yy  <- matrix(Yx, nrow=1) %*% sparsematrix[,1:66870]
        ee  <- matrix(Ex, nrow=1) %*% sparsematrix[,1:66870]
        LLRa <- yy@x*log(yy@x/ee@x) + (sum(Yx) - yy@x)*log((sum(Yx)- yy@x)/(sum(Ex)-ee@x))
        maxlL <- which.max(LLRa)
        maxLiks[m] <- maxlL
        LLRa_stat <- LLRa[maxlL]
        ixmaxL <- sparsematrix[,maxlL]
        RRmaxL <- yy@x[maxlL]/ee@x[maxlL]
        #MC test for cluster
        YSIM= rmultinom(nsim, sum(Yx), prob=Ex)
        #under the null
        yy0 <- sapply(1:nsim, function(i) matrix(YSIM[,i], nrow=1)%*% sparsematrix[,1:66870])
        ee0 <- ee
        LLR0 <- sapply(1:nsim, function(i) yy0[[i]]@x*log(yy0[[i]]@x/ee0@x) + (sum(YSIM[,i]) - yy0[[i]]@x)*log((sum(YSIM[,i])- yy0[[i]]@x)/(sum(Ex)-ee0@x)))
        #calculate pvalue
        pval <- sum(apply(LLR0,2, function(x) ifelse(any(x>= LLRa_stat , na.rm = TRUE),1,0)))/nsim
        #print(pval)
        mlc_pvals[m] <- pval
        #remove cluster
        rr <- rep(1,1040)
        rr[which(ixmaxL!=0)] <- RRmaxL
        mlc[,m] <- rr
        E.cur <- Ex*rr 
        Ex <- scales(matrix(E.cur,ncol=Time),matrix(Yx, ncol=Time), Time)    
        #str(Ex)
        if(pval > 0.05){
            break
        }
        next
    }
    
    return(list(clusters=mlc,
                maxLiks = maxLiks,
                pvals = mlc_pvals))
}


clusso_prob_clusteroverlap <- function(sparsematrix,lassoresult,selected,rr, risk,nsim,Time, ncentroids, pow){
    #DEFINE TRUTH
    if(risk==1){
        warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk,1,0)    
    }
    #Take out the time vectors - only keep cluster part of matrix
    sparsematrix_clusteronly <- sparsematrix[,-c(ncol(sparsematrix)-Time+1:ncol(sparsematrix))]
    #select out my betas for each sim
    betaselect <- lassoresult$coefs.lasso.all[,selected]
    #binarize
    betaselect_bin <- ifelse(abs(round(betaselect,4) >= 10e-3),1,0)
    #only take the clusters betas
    betaselect_bin_clusteronly <- betaselect_bin[-c(ncol(sparsematrix)-Time+1:ncol(sparsematrix))]
    ##INCLUSTER
    if(pow==TRUE){
        #print("pow")
        clusteroverlap <- t(rrmatvec) %*% sparsematrix_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
        incluster_sim <- clusteroverlap_bin %*% betaselect_bin_clusteronly
        incluster_sim_bin <- ifelse(incluster_sim !=0,1,0)
        return(idin = incluster_sim_bin)
    }
    else{
        ##OUTCLUSTER
        #print("fp")
        clusteroverlap <- t(rrmatvec) %*% sparsematrix_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
        outcluster_sim <- clusteroverlap_bin %*% betaselect_bin_clusteronly
        outcluster_sim_bin <- ifelse(outcluster_sim !=0,1,0)
        return(idout = outcluster_sim_bin)
    }
}


step_clusterix <- function(sparsematrix, stepscan, numclustersid){
    ixids <- NULL
    if(numclustersid!=0){
        for(i in 1:numclustersid){
            ixid_i <- which(sparsematrix[,stepscan$maxLiks[which(stepscan$pvals>0.05)-i]]==1)
            ixids <- c(ixids, ixid_i)
        }
    } else{
        ixids <-0
    }

    return(unique(ixids))
    
}

spatscanfs_prob_clusteroverlap <- function(res_stepsscan, ixids,numclustersid ,sparsematrix,rr,risk,pow,nsim){
    #DEFINE TRUTH
    if(risk==1){
        warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk,1,0)    
    }
    # numclustersid <- lapply(1:nsim, function(i) which(res_stepsscan[[i]]$pvals>0.05)-1)
    #ixids <- lapply(1:nsim, function(i) step_clusterix(sparsematrix, res_stepsscan[[i]], numclustersid=numclustersid[[i]]))
    
    rrmatvec <- ifelse(as.vector(rr)==risk,1,0)  
    sparsematrix_clusteronly <- sparsematrix[,-c(ncol(sparsematrix)-Time+1:ncol(sparsematrix))]
    clusteroverlap <- t(rrmatvec) %*% sparsematrix_clusteronly
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) rep(0,1040))
    for(i in 1:nsim){
        betaselect_bin_clusteronly[[i]][ixids[[i]]] <- 1
    }
    overlap <- ifelse(as.vector(sparsematrix_clusteronly %*% clusteroverlap_bin)!=0,1,0)
    if(pow==TRUE){
        #power
        power <- sum(ifelse(unlist(lapply(1:nsim, function(i) betaselect_bin_clusteronly[[i]]%*%overlap))!=0,1,0))/nsim
        return(idin=power)
    } else{
        overlapfp <- ifelse(overlap==1,0,1)
        fp <- sum(ifelse(unlist(lapply(1:nsim, function(i) betaselect_bin_clusteronly[[i]]%*%overlapfp))!=0,1,0))/nsim
        return(idout=fp)
    }
}

reval <- function(probs, ix){
    probs[ix] <-1
    return(probs)
}

forwardstage_prob_clusteroverlap <- function(sim_stage,sparsematrix, rr, risk,pow, nsim){
    if(risk==1){
        warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk,1,0)    
    }
    bkg.bic <- lapply(1:nsim, function(i) sapply(1:Time, function(t) round(as.numeric(attributes(sort(table(sim_stage[[i]]$RRbic[,t]), decreasing = TRUE)[1])[[1]]),4)))
    bkg.aic <- lapply(1:nsim, function(i) sapply(1:Time, function(t) round(as.numeric(attributes(sort(table(sim_stage[[i]]$RRaic[,t]), decreasing = TRUE)[1])[[1]]),4)))
    sparsematrix_clusteronly <- sparsematrix[,-c(ncol(sparsematrix)-Time+1:ncol(sparsematrix))]
    
    clusteroverlap <- t(rrmatvec) %*% sparsematrix_clusteronly
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
    overlap <- ifelse(as.vector(sparsematrix_clusteronly %*% clusteroverlap_bin)!=0,1,0)
    #####################
    #Power
    #####################
    betaselect_bin.bic <- lapply(1:nsim, function(i) as.vector(sapply(1:Time, function(t) ifelse(abs(round(sim_stage[[i]]$RRbic[,t],4) > bkg.bic[[i]][t]),1,0))))
    betaselect_bin.aic <- lapply(1:nsim, function(i) as.vector(sapply(1:Time, function(t) ifelse(abs(round(sim_stage[[i]]$RRaic[,t],4) > bkg.aic[[i]][t]),1,0))))
    if(pow==TRUE){
        #power
        power.bic <- sum(ifelse(unlist(lapply(1:nsim, function(i) betaselect_bin.bic[[i]]%*%overlap))!=0,1,0))/nsim
        power.aic <- sum(ifelse(unlist(lapply(1:nsim, function(i) betaselect_bin.aic[[i]]%*%overlap))!=0,1,0))/nsim
        return(list(power.bic = power.bic, power.aic = power.aic,  
                    betaselect_bin.bic =  betaselect_bin.bic,
                    betaselect_bin.aic =  betaselect_bin.aic))
    }  else{
        #fp
        overlapfp <- ifelse(overlap==1,0,1)
        fp.bic <- sum(ifelse(unlist(lapply(1:nsim, function(i) betaselect_bin.bic[[i]]%*%overlapfp))!=0,1,0))/nsim
        fp.aic <- sum(ifelse(unlist(lapply(1:nsim, function(i) betaselect_bin.aic[[i]]%*%overlapfp))!=0,1,0))/nsim
        return(list(fp.bic = fp.bic, fp.aic = fp.aic,
                    betaselect_bin.bic =  betaselect_bin.bic,
                    betaselect_bin.aic =  betaselect_bin.aic))
    }
}

