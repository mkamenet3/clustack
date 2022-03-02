colorsgrey <- function (x) {
    y = colorRamp(RColorBrewer::brewer.pal(9, "Greys")[1:9])(x)
    rgb(y[, 1], y[, 2], y[, 3], maxColorValue = 255)
}




create_plotmeanrr_stack <- function(res,IC, flav,Time, nsim, sim.i, greys){
    if(IC=="aic"){
        selects <- sapply(1:nsim, function(i) res[[i]]$selection.aic)    
    } else{
        selects <- sapply(1:nsim, function(i) res[[i]]$selection.bic)
    }
    #browser()
    if(all(selects==0)){
        ric <- matrix(rep(1, ncol(res[[1]]$Lambda_dense)), ncol=Time)
    } else {
        meanrr <- lapply(1:nsim, function(i) res[[i]]$wLambda[selects[i],])
        names(meanrr) <- paste0("l",1:nsim)
        meanrr[which(selects==0)] <- NULL
        meanrr.df <- do.call(rbind, meanrr)
        meanrrs <- apply(meanrr.df,2, mean)
        ric <- matrix(meanrrs, ncol = Time)
    }
    plotmeanrr_stack(ric, Time, sim.i,ic=IC, flav=flav, greys)
}

create_plotFPR_stack <- function(res,IC, flav,Time, nsim, sim.i){
    if(IC=="aic"){
        selects <- sapply(1:nsim, function(i) res[[i]]$selection.aic)    
    } else{
        selects <- sapply(1:nsim, function(i) res[[i]]$selection.bic)
    }
    #browser()
    if(all(selects==0)){
        #ric <- matrix(rep(1, ncol(res[[1]]$Lambda_dense)), ncol=Time)
        probs <- rep(0, ncol(res[[1]]$Lambda_dense))
    } else {
        #find pc's that overlap maxloc
        if(flav=="loc"){
            maxid <- sapply(1:nsim, function(i) res[[i]]$maxid[selects[i]])    
            vec <- rep(0, 208 * Time)
            position <- list(vec)[rep(1, nsim)]
            simindicator <- mapply(reval, position, maxid)
            probs <- Matrix::rowSums(simindicator)/nsim
        } else{
            maxid <- sapply(1:nsim, function(i) res[[i]]$maxid[selects[i]])   
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
        g.ramp <- gray.colors(26, start=1, end =0)
        #color.ic <- sapply(1:Time, function(i) greys(ric[,i]))
        color.ic <- sapply(1:Time, function(i) g.ramp[100*ric[,i]])
        pdf(paste0(sim.i,"_fpr_",flav, "_",ic,".pdf"), height=11, width=10)
        
    } else {
        color.ic <- sapply(1:Time, function(i) redblue(log(1.5 *pmax(1/1.5, pmin(ric[, i], 1.5)))/log(1.5^2)))
        print(paste0(sim.i,"_meanrr_",flav, "_",ic,".pdf"))
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
        text(seq(.6,1.4,length=5),rep(.45,5),seq(0.5,1.5,length.out=5),srt=330,adj=0)
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
plotmap <- function(res, pdfname=NULL, genpdf=TRUE, maxrr=2){
    if(!is.null(maxrr)){
        maxrr=maxrr
    }
    else{
        maxrr=2
    }
    # if(!is.null(minrr)){
    #     minrr=minrr
    # }
    # else{
    #     minrr=0
    # }
    #cluster_ix <- redblue(log(2 *  pmax(1/2, pmin(res, 2)))/log(4))
    #cluster_ix <- redblue(log(maxrr *  pmax(1/maxrr, pmin(res, maxrr)))/log(maxrr^2))
    cluster_ix <- redblue(log(maxrr *  pmax(1/maxrr, pmin(res, maxrr)))/log(maxrr^2))
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
    len <- length(redblue(15:50/50))
    rect(seq(.6,1.4,length=len)[-len],.5,seq(.65,1.4,length=len)[-1],.62,col=redblue(15:50/50),border=F)
    text(seq(.6,1.4,length=7),rep(.45,5),seq(1/maxrr,maxrr,length.out=7),srt=330,adj=0)   
    
    # plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    # rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    # #text(seq(.6,1.4,length=5),rep(.45,5),seq(0,maxrr,length.out=5),srt=330,adj=0)
    # text(seq(.6,1.4,length=5),rep(.45,5),seq(minrr,maxrr,length.out=5),srt=330,adj=0)
    
    if(genpdf==TRUE){
        dev.off()    
    }
    
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
    for(c in 1:maxclust){ 
        #identify first MLC
        yy  <- matrix(Yx, nrow=1) %*% sparsematrix[,1:66870]
        ee  <- matrix(Ex, nrow=1) %*% sparsematrix[,1:66870]
        LLRa <- yy@x*log(yy@x/ee@x) + (ee@x - yy@x) #cluster only likelihood
        #LLRa <- yy@x*log(yy@x/ee@x) + (sum(Yx) - yy@x)*log((sum(Yx)- yy@x)/(sum(Ex)-ee@x))
        maxlL <- which.max(LLRa)
        maxLiks[c] <- maxlL
        LLRa_stat <- LLRa[maxlL]
        ixmaxL <- sparsematrix[,maxlL]
        RRmaxL <- yy@x[maxlL]/ee@x[maxlL]
        #MC test for cluster
        YSIM= rmultinom(nsim, sum(Yx), prob=Ex)
        #YSIM= rmultinom(nsim, sum(Ex), prob=Ex)
        #under the null
        yy0 <- sapply(1:nsim, function(i) matrix(YSIM[,i], nrow=1)%*% sparsematrix[,1:66870])
        ee0 <- ee
        #ee0 <- ee
        #LLR0 <- sapply(1:nsim, function(i) yy0[[i]]@x*log(yy0[[i]]@x/ee0@x) + (sum(YSIM[,i]) - yy0[[i]]@x)*log((sum(YSIM[,i])- yy0[[i]]@x)/(sum(Ex)-ee0@x)))
        
        #LLR0 <- sapply(1:nsim, function(i) yy0[[i]]@x*log(yy0[[i]]@x/ee0@x) + (sum(YSIM[,i]) - yy0[[i]]@x)*log((sum(YSIM[,i])- yy0[[i]]@x)/(sum(Ex)-ee0@x)))
        #calculate pvalue
        
        #test <- apply(LLR0,2, function(x) which(x > LLRa_stat))
        LLR0 <-  sapply(1:nsim, function(i) yy0[[i]]@x*log(yy0[[i]]@x/ee0@x) + (ee0@x - yy0[[i]]@x))
        
        
        pval <- sum(apply(LLR0,2, function(x) ifelse(any(x>= LLRa_stat , na.rm = TRUE),1,0)))/nsim
        #print(pval)
        mlc_pvals[c] <- pval
        #remove cluster
        rr <- rep(1,1040)
        rr[which(ixmaxL!=0)] <- RRmaxL
        mlc[,c] <- rr
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

step_clusterix <- function(sparsematrix, stepscan, numclustersid, thresh){
    ixids <- NULL
    #if(numclustersid!=0){
    if(length(numclustersid)!=0){
        #browser()
        for(i in 1:length(numclustersid)){
           # print(i)
            #browser()
            #for(j in 0:length(which(!is.na(stepscan$pvals)))){
            #for(j in 0:(length(which(!is.na(stepscan$pvals)))-1)){
            for(j in 0:length(which(stepscan$pvals < thresh))){
                #browser()
                #print(j)
                #ixid_i <- which(sparsematrix[,stepscan$maxLiks[which(stepscan$pvals>thresh)-(i+j)]]==1)
                #ixid_i <- which(sparsematrix[,stepscan$maxLiks[which(stepscan$pvals>thresh)-(j+1)]]==1)
                ixid_i <- which(sparsematrix[,stepscan$maxLiks[max(which(stepscan$pvals < thresh))-j]]==1)
                ixids <- c(ixids, ixid_i)
            }
        }
    } else{
        ixids <-0
    }
    
    return(unique(ixids))
    
}


spatscanfs_prob_clusteroverlap <- function(res_stepsscan, ixids,numclustersid ,sparsematrix,rr,risk,pow){
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
    betaselect_bin_clusteronly <- rep(0,1040)
    betaselect_bin_clusteronly[ixids] <- 1
    # for(i in 1:nsim){
    #     betaselect_bin_clusteronly[[i]][ixids[[i]]] <- 1
    # }
    overlap <- ifelse(as.vector(sparsematrix_clusteronly %*% clusteroverlap_bin)!=0,1,0)
    if(pow==TRUE){
        #power
        power <- ifelse(unlist(betaselect_bin_clusteronly%*%overlap)!=0,1,0)
        return(idin=power)
    } else{
        overlapfp <- ifelse(overlap==1,0,1)
        fp <- ifelse(unlist(betaselect_bin_clusteronly%*%overlapfp)!=0,1,0)
        return(idout=fp)
    }
}

reval <- function(probs, ix){
    probs[ix] <-1
    return(probs)
}

forwardstage_prob_clusteroverlap <- function(sim_stage,sparsematrix, rr, risk,pow){
    if(risk==1){
        warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk,1,0)    
    }
    bkg.bic <-sapply(1:Time, function(t) round(as.numeric(attributes(sort(table(sim_stage$RRbic[,t]), decreasing = TRUE)[1])[[1]]),4))
    bkg.aic <- sapply(1:Time, function(t) round(as.numeric(attributes(sort(table(sim_stage$RRaic[,t]), decreasing = TRUE)[1])[[1]]),4))
    sparsematrix_clusteronly <- sparsematrix[,-c(ncol(sparsematrix)-Time+1:ncol(sparsematrix))]
    
    clusteroverlap <- t(rrmatvec) %*% sparsematrix_clusteronly
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
    overlap <- ifelse(as.vector(sparsematrix_clusteronly %*% clusteroverlap_bin)!=0,1,0)
    #####################
    #Power
    #####################
    betaselect_bin.bic <- as.vector(sapply(1:Time, function(t) ifelse(abs(round(sim_stage$RRbic[,t],4) > bkg.bic[t]),1,0)))
    betaselect_bin.aic <- as.vector(sapply(1:Time, function(t) ifelse(abs(round(sim_stage$RRaic[,t],4) > bkg.aic[t]),1,0)))
    if(pow==TRUE){
        #power
        power.bic <- ifelse(unlist(betaselect_bin.bic%*%overlap)!=0,1,0)
        power.aic <- ifelse(unlist(betaselect_bin.aic%*%overlap)!=0,1,0)
        return(list(power.bic = power.bic, power.aic = power.aic,  
                    betaselect_bin.bic =  betaselect_bin.bic,
                    betaselect_bin.aic =  betaselect_bin.aic))
    }  else{
        #fp
        overlapfp <- ifelse(overlap==1,0,1)
        fp.bic <- ifelse(unlist(betaselect_bin.bic%*%overlapfp)!=0,1,0)
        fp.aic <- ifelse(unlist(betaselect_bin.aic%*%overlapfp)!=0,1,0)
        return(list(fp.bic = fp.bic, fp.aic = fp.aic,
                    betaselect_bin.bic =  betaselect_bin.bic,
                    betaselect_bin.aic =  betaselect_bin.aic))
    }
}

