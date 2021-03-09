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
#'@param byloc If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{TRUE}.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#'@param overdisp.est Overdispersion parameter estimated across all simulations (max).
#'@return Returns list for each iteration with weighted relative risks by location inside identified cluster.
#'@export
detectclusters <- function(sparseMAT, Ex, Yx,numCenters,Time, maxclust,byloc=TRUE, model=c("poisson", "binomial"),
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
    if(byloc==FALSE){
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
                    #maxpcs = res[[3]],
                    maxid = res[[3]],
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
                    #maxlocs = res[[3]],
                    maxid = res[[3]],
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

bucklandbounds.cluster <- function(thetaa, thetai,res, w_q, outExp_out ,IC,transform="none",overdisp.est) {
    #browser()
    if(transform=="log"){
        if(!is.null(overdisp.est)){
            varthetai <- overdisp.est*(1/(thetai*outExp_out ))
        } else {
            varthetai <- 1/(thetai*outExp_out )
        }
        withintheta <- (log(thetai) - log(thetaa))^2
        varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))^2
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
        LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        
    } else {
        if(!is.null(overdisp.est)){
            varthetai <- overdisp.est*(thetai/outExp_out )
        } else {
            varthetai <- (thetai/outExp_out)
        }
        withintheta <- (thetai - thetaa)^2
        varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
        LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
    }
    
    
    return(list(buckland.LB = LBa,
                clusterMA = thetaa,
                buckland.UB = UBa))
}

bucklandbounds.cells <- function(thetaa, res, w_q, outExp_out ,IC,transform="none",tsparsematrix,overdisp.est, cellsix_out) {
    #browser()
    thetai <- res$Lambda_dense
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    print(str(thetaa))
    if(transform=="log"){
        print("transform")
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(tsparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp_out )))
        } else {
            varthetai <- sapply(1:nrow(tsparsematrix), function(k) 1/(thetai[k,]*outExp_out))
        }
        withintheta <- sapply(1:length(cellsix_out), function(j) (log(thetai[,cellsix_out][,j]) - log(thetaa[j])))^2
        wtnbtn <- sapply(1:nrow(tsparsematrix), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
        varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
        LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        
    } else {
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(tsparsematrix), function(k) overdisp.est*thetai[k,]/outExp_out)
        } else {
            varthetai <- sapply(1:nrow(tsparsematrix), function(k) thetai[k,]/outExp_out)
        }
        withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
        wtnbtn <- sapply(1:nrow(tsparsematrix), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
        varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
        LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
    }
    return(list(buckland.LB = LBa,
                clusterMA = thetaa,
                buckland.UB = UBa))
}

maw2.cluster <- function(thetaa,thetai, w_q, outExp, transform,overdisp.est) {
    if(transform=="log"){
        #print("log-scale")
        if(!is.null(overdisp.est)){
            varthetai <-  overdisp.est*(1/(thetai*outExp))
        } else {
            varthetai <-1/(thetai*outExp)
        }
        withintheta <- (log(thetai) - log(thetaa))^2
        var_thetaa <- sum(w_q*(varthetai + withintheta))
        var_thetaa <- as.vector(var_thetaa)
        UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
        LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
    } else{
        if(!is.null(overdisp.est)){
            varthetai <- overdisp.est*(thetai/outExp)
        } else {
            varthetai <- thetai/outExp
        }
        withintheta <- (thetai - thetaa)^2
        var_thetaa <- sum(w_q*(varthetai + withintheta))
        var_thetaa <- as.vector(var_thetaa)
        UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
        LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
    }
    return(list(maw2.LB = LBa,
                clusterMA = thetaa,
                maw2.UB = UBa))
}
#cluster_thetaa,
# clusterRR_ilarge,
# w_q=wslarge,
# outExp,
# transform="none",
# tsparsematrix=t(sparsematrix),
# overdisp.est,
# cellsix_out
maw2.cells <- function(thetaa,res, w_q, outExp, IC, transform,tsparsematrix, overdisp.est, cellsix_out){
    thetai <- res$Lambda_dense
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
        if(transform=="log"){
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(tsparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
            } else {
                varthetai <- sapply(1:nrow(tsparsematrix), function(k) 1/(thetai[k,]*outExp))
            }
            withintheta <- sapply(1:length(cellsix_out), function(j) (log(thetai[, cellsix_out][,j]) - log(thetaa[j])))^2
            wtnbtn <- sapply(1:nrow(tsparsematrix), function(k) (varthetai[cellsix_out,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        } else{
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(tsparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
            } else {
                varthetai <- sapply(1:nrow(tsparsematrix), function(k) thetai[k,]/outExp)
            }
            withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
            wtnbtn <- sapply(1:nrow(tsparsematrix), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
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

matabounds.cluster <- function(thetaa,thetai, w_q, outExp, transform=c("none","log", "sqrt"), overdisp.est) {
    if(is.null(transform)){
        none = matabounds_none(thetaa,thetai, w_q, outExp, overdisp.est)
    }
    switch(transform,
           none = matabounds_none.cluster(thetaa,thetai, w_q,  outExp, overdisp.est),
           log = matabounds_log.cluster(thetaa, thetai,w_q, outExp, overdisp.est),
           sqrt = matabounds_sqrt.cluster(thetaa,thetai, w_q,  outExp, overdisp.est))
}


matabounds.cells <- function(thetaa,res, w_q, outExp, IC, transform=c("none","log", "sqrt"),tsparsematrix,overdisp.est, cellsix_out) {
    if(is.null(transform)){
        none = matabounds_none.cells(thetaa, res, w_q, outExp, IC, tsparsematrix, overdisp.est, cellsix_out)
    }
    switch(transform,
           none = matabounds_none.cells(thetaa, res, w_q, outExp, IC, tsparsematrix, overdisp.est, cellsix_out),
           log = matabounds_log.cells(thetaa, res, w_q, outExp, IC, tsparsematrix, overdisp.est, cellsix_out),
           sqrt = matabounds_sqrt.cells(thetaa, res, w_q, outExp, IC, tsparsematrix, overdisp.est, cellsix_out))
}

matabounds_none.cluster <- function(thetaa,thetai, w_q, outExp, overdisp.est) {
   # browser()
    if(!is.null(overdisp.est)){
        varthetai <- overdisp.est*(thetai/outExp)
    } else {
        varthetai <- thetai/outExp
    }
    mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                      thetaii= thetai,
                      se.thetaii=sqrt(varthetai@x),
                      w_q=w_q, alpha=0.025, tol=1e-8)$root    
    
    mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                      thetaii= thetai,
                      se.thetaii=sqrt(varthetai@x),
                      w_q=w_q, alpha=1-0.025, tol=1e-8)$root 
    return(list(mata.LB = mataLB,
                clusterMA = thetaa,
                mata.UB = mataUB))
}

matabounds_none.cells <- function(thetaa, res, w_q,outExp, IC, tsparsematrix, overdisp.est,cellsix_out) {
    #browser()
    thetai <- res$Lambda_dens
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(tsparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
    } else {
        varthetai <- sapply(1:nrow(tsparsematrix), function(k) thetai[k,]/outExp)
    }
    mataLB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetaii= thetai[,cellsix_out[k]],
                                         se.thetaii=sqrt(varthetai[cellsix_out[k],]),
                                         w_q=w_q, alpha=0.025, tol=1e-8)$root)
    
    mataUB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetaii= thetai[,cellsix_out[k]],
                                         se.thetaii=sqrt(varthetai[cellsix_out[k],]),
                                         w_q=w_q, alpha=1-0.025, tol=1e-8)$root)
    
    return(list(mata.LB = mataLB,
                clusterMA = thetaa,
                mata.UB = mataUB))
}

# matabounds_none <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
#     if(cellrates==TRUE){
#         #browser()
#         if(!is.null(overdisp.est)){
#             varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
#         } else {
#             varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
#         }
#         mataLB <- sapply(1:ncol(thetai), 
#                          function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
#                           thetaii= thetai[,k],
#                           se.thetaii=sqrt(varthetai[k,]),
#                           w_q=w_q, alpha=0.025, tol=1e-8)$root)    
#         
#         mataUB <- sapply(1:ncol(thetai),
#                          function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
#                           thetaii= thetai[,k],
#                           se.thetaii=sqrt(varthetai[k,]),
#                           w_q=w_q, alpha=1-0.025, tol=1e-8)$root) 
#     } else{
#         if(!is.null(overdisp.est)){
#             varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
#         } else {
#             varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
#         }
#         mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
#                           thetaii= thetai,
#                           se.thetaii=sqrt(varthetai),
#                           w_q=w_q, alpha=0.025, tol=1e-8)$root    
#         
#         mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
#                           thetaii= thetai,
#                           se.thetaii=sqrt(varthetai),
#                           w_q=w_q, alpha=1-0.025, tol=1e-8)$root 
#     }
#     
#     return(list(mata.LB = mataLB,
#                 clusterMA = thetaa,
#                 mata.UB = mataUB))
# }

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
#cellsix if TRUE, then perform bounds for cluster RR. Otherwise cellsix is the index of each cell and bounds calcuated for each cell
#nonma

nonma.cluster <- function(cluster_thetaa,res, clusterRR_ilarge,wslarge,idix, outExp_out,IC, transform, byloc=TRUE){
    if(byloc==FALSE){
        #browser()
        print("by pc")
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- unique(res$Lambda_dense[res$maxid[res$selection.aic],])[2]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outExp_out[res$maxid[res$selection.aic]]))
                
            } else {
                clusterRRlarge <- unique(res$Lambda_dense[res$maxid[res$selection.bic],])[2]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out[res$maxid[res$selection.bic]])
                
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=exp(log(cluster_thetaa)-1.96*se_clusterRRlarge), 
                                                                 clusterMA = cluster_thetaa,
                                                                 ub=exp(log(cluster_thetaa)+1.96*se_clusterRRlarge)))
        } else {
            if(IC=="aic") {
                clusterRRlarge <- unique(res$Lambda_dense[res$maxid[res$selection.aic],])[2]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out[res$maxid[res$selection.aic]])
                
            } else {
                clusterRRlarge <- unique(res$Lambda_dense[res$maxid[res$selection.bic],])[2]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out[res$maxid[res$selection.bic]])
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=cluster_thetaa-1.96*se_clusterRRlarge, 
                                                                 clusterMA = cluster_thetaa,
                                                                 ub=cluster_thetaa+1.96*se_clusterRRlarge))
        }
    } else {
        print("by loc")
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- unique(res$Lambda_dense[,res$maxid[res$selection.aic]])[2]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outExp_out[res$maxid[res$selection.aic]]))
                
            } else {
                clusterRRlarge <- unique(res$Lambda_dense[,res$maxid[res$selection.bic]])[2]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outExp_out[res$maxid[res$selection.bic]]))
                
            }
            nonma.theta.time <-system.time(nonma.theta <-  cbind(lb=exp(log(cluster_thetaa)-1.96*se_clusterRRlarge), 
                                                                 clusterMA = cluster_thetaa,
                                                                 ub=exp(log(cluster_thetaa)+1.96*se_clusterRRlarge)))
        } else {
            if(IC=="aic") {
                clusterRRlarge <-unique(res$Lambda_dense[,res$maxid[res$selection.aic]])[2]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out[res$maxid[res$selection.aic]])
                
            } else {
                clusterRRlarge <- unique(res$Lambda_dense[,res$maxid[res$selection.bic]])[2]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out[res$maxid[res$selection.bic]])
                
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=cluster_thetaa-1.96*se_clusterRRlarge, 
                                                                 clusterMA = cluster_thetaa,
                                                                 ub=cluster_thetaa+1.96*se_clusterRRlarge))
            
        }
    }
    return(list(nonma.theta.time = nonma.theta.time[[3]],
                nonma.theta = nonma.theta))
}

nonma.cells <- function(cluster_thetaa,res, clusterRR_ilarge,wslarge,idix, outExp_out,IC, transform, byloc=TRUE, cellrisk_wt_out=NULL, cellsix_out){
    if(byloc==FALSE){
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outExp_out))
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
                
            } else {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outExp_out))
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=exp(log(clusterMA)-1.96*se_clusterRRlarge), 
                                                                 clusterMA = clusterMA,
                                                                 ub=exp(log(clusterMA)+1.96*se_clusterRRlarge)))
        }else {
            if(IC=="aic") {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out)
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
            } else {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outExp_out)
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=clusterMA-1.96*se_clusterRRlarge, 
                                                                 clusterMA = clusterMA,
                                                                 ub=clusterMA+1.96*se_clusterRRlarge))
        } 
        

    }else {
        if(transform=="log"){
            if(IC=="aic") {
                se_clusterRRlarge <- sqrt(1/(cellrisk_wt_out*outExp_out))
            } else {
                se_clusterRRlarge <- sqrt(1/(cellrisk_wt_out*outExp_out))
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=exp(log(cellrisk_wt_out)-1.96*se_clusterRRlarge), 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=exp(log(cellrisk_wt_out) +1.96*se_clusterRRlarge)))
        } else {
            if(IC=="aic") {
                se_clusterRRlarge <- sqrt(cellrisk_wt_out/outExp_out)
            } else {
                se_clusterRRlarge <- sqrt(cellrisk_wt_out/outExp_out)
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=cellrisk_wt_out-1.96*se_clusterRRlarge, 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=cellrisk_wt_out +1.96*se_clusterRRlarge))
        }
    }
    return(list(nonma.theta.time = nonma.theta.time[[3]],
                nonma.theta = nonma.theta))
}




#nonma_asymp
nonma_asymp.cluster <- function(cluster_thetaa,res, clusterRR_ilarge,wslarge,idix, outObs_out,IC, transform, byloc=TRUE){
    if(byloc==FALSE){
        print("by pc")
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.aic],][2]
                se_clusterRRlarge_asymp <- sqrt(1/(cluster_thetaa*outObs_out[res$maxid[res$selection.aic]]))
            } else {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.bic],][2]
                se_clusterRRlarge_asymp <- sqrt(1/(clusterRRlarge*outObs_out[res$maxid[res$selection.bic]]))
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lbasymp=exp(log(cluster_thetaa)-1.96*se_clusterRRlarge_asymp), 
                                                                             clusterMA = cluster_thetaa,
                                                                             ubasymp=exp(log(cluster_thetaa)+1.96*se_clusterRRlarge_asymp)))
        } else {
            if(IC=="aic") {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.aic],][2]
                se_clusterRRlarge_asymp <- sqrt(clusterRRlarge/outObs_out[res$maxid[res$selection.aic]])
            } else {
                clusterRRlarge <-  res$Lambda_dense[res$maxid[res$selection.bic],][2]
                se_clusterRRlarge_asymp <- sqrt(clusterRRlarge/outObs_out[res$maxid[res$selection.bic]])
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lbasymp=cluster_thetaa-1.96*se_clusterRRlarge_asymp, 
                                                                             clusterMA = cluster_thetaa,
                                                                             ubasymp=cluster_thetaa+1.96*se_clusterRRlarge_asymp))
        }
    } else {
        print("by LOC")
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- unique(res$Lambda_dense[,res$maxid[res$selection.aic]])[2]
                se_clusterRRlarge_asymp <- sqrt(1/(clusterRRlarge*outObs_out[res$maxid[res$selection.aic]]))    
            } else {
                clusterRRlarge <- unique(res$Lambda_dense[,res$maxid[res$selection.bic]])[2]
                se_clusterRRlarge_asymp <- sqrt(1/(clusterRRlarge*outObs_out[res$maxid[res$selection.bic]]))    
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lbasymp=exp(log(cluster_thetaa)-1.96*se_clusterRRlarge_asymp), 
                                                                             clusterMA = cluster_thetaa,
                                                                             ubasymp=exp(log(cluster_thetaa)+1.96*se_clusterRRlarge_asymp)))
            
        } else {
            if(IC=="aic") {
                se_clusterRRlarge_asymp <- sqrt(cluster_thetaa/outObs_out[res$maxid[res$selection.aic]]) 
            } else{
                se_clusterRRlarge_asymp <- sqrt(cluster_thetaa/outObs_out[res$maxid[res$selection.bic]]) 
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lbasymp=cluster_thetaa-1.96*se_clusterRRlarge_asymp, 
                                                                             clusterMA = cluster_thetaa,
                                                                             ubasymp=cluster_thetaa+1.96*se_clusterRRlarge_asymp))
            
        }
        
    }
    return(list(nonma_asymp.theta.time = nonma_asymp.theta.time[[3]],
                nonma_asymp.theta = nonma_asymp.theta))
}


nonma_asymp.cells <-function(cluster_thetaa,res, clusterRR_ilarge,wslarge,idix, outObs_out,IC, transform, byloc=TRUE, cellrisk_wt_out, cellsix_out){
    if(byloc==FALSE){
        print("by PC")
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outObs_out))
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
                
            } else {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_clusterRRlarge <- sqrt(1/(clusterRRlarge*outObs_out))
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=exp(log(clusterMA)-1.96*se_clusterRRlarge), 
                                                                 clusterMA = clusterMA,
                                                                 ub=exp(log(clusterMA)+1.96*se_clusterRRlarge)))
        }else {
            if(IC=="aic") {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outObs_out)
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
            } else {
                clusterRRlarge <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_clusterRRlarge <- sqrt(clusterRRlarge/outObs_out)
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=clusterMA-1.96*se_clusterRRlarge, 
                                                                 clusterMA = clusterMA,
                                                                 ub=clusterMA+1.96*se_clusterRRlarge))
        } 
        
        
    }else {
        print("by loc")
        print(cellrisk_wt_out)
        if(transform=="log"){
            if(IC=="aic") {
                se_clusterRRlarge <- sqrt(1/(cellrisk_wt_out*outObs_out))
            } else {
                se_clusterRRlarge <- sqrt(1/(cellrisk_wt_out*outObs_out))
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=exp(log(cellrisk_wt_out)-1.96*se_clusterRRlarge), 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=exp(log(cellrisk_wt_out) +1.96*se_clusterRRlarge)))
        } else {
            if(IC=="aic") {
                se_clusterRRlarge <- sqrt(cellrisk_wt_out/outObs_out)
            } else {
                se_clusterRRlarge <- sqrt(cellrisk_wt_out/outObs_out)
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=cellrisk_wt_out-1.96*se_clusterRRlarge, 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=cellrisk_wt_out +1.96*se_clusterRRlarge))
        }
        
    }
    
    return(list(nonma_asymp.theta.time = nonma_asymp.theta.time[[3]],
                nonma_asymp.theta = nonma_asymp.theta))
}


calcbounds <- function(idix, IC, res, byloc, Ex, Obs,target=c("cluster", "cells"), cellsix=NULL, sparsematrix){
    if(is.null(cellsix) | target=="cluster"){
        print("no cell rates")
    }
    if(!is.null(cellsix) & is.null(sparsematrix)){
        stop("You must provide the sparsematrix when calculating rates for each cell.")
    }
    if(IC=="aic"){
        wslarge <- res$wtMAT[,res$selection.aic]
    } else {
        wslarge <- res$wtMAT[,res$selection.bic]
    }
    clusterRR_uniqlarge <- sapply(1:nrow(res$Lambda_dense), function(k) unique(res$Lambda_dense[k,]))
    #clusterRR_ilarge <- rep(NA, dim(res$Lambda_dense)[1])
    clusterRR_uniq_ilarge <- as.matrix(do.call(rbind, clusterRR_uniqlarge), ncol=2)
    clusterRR_ilarge <- selectuniqRR(clusterRR_uniq_ilarge)
    cluster_thetaa <- sum(clusterRR_ilarge*wslarge)
    if(target=="cells" & is.null(cellsix)){
        stop("For cell-wise estimates, you must provide the index of the cells of interest")
    }
    switch(target,
           cluster = calcbounds.cluster(idix, IC, res, byloc, Ex, Obs,wslarge, cluster_thetaa,clusterRR_ilarge, sparsematrix),
           cells = calcbounds.cells(idix, IC, res, byloc, Ex, Obs,wslarge, cluster_thetaa,clusterRR_ilarge, sparsematrix, cellsix))
}


calcbounds.cluster <- function(idix, IC, res, byloc, Ex, Obs,wslarge, cluster_thetaa,clusterRR_ilarge, sparsematrix){
    print("calcbounds.cluster")
    print(str(sparsematrix))
    if(byloc==TRUE){
        print(paste0("byloc", byloc))
        outExp_out <- Ex
        outObs_out <- Obs
        outExp <- t(sparsematrix)%*%Ex
        outObs <- t(sparsematrix)%*%Obs
    }
    else if(byloc==FALSE){
        print(paste0("byloc", byloc))
        outExp <- t(sparsematrix)%*%Ex
        outObs <- t(sparsematrix)%*%Obs
        outExp_out <- outExp@x
        outObs_out <- outObs@x
    }
    outnonma.time <- system.time(outnonma <- nonma.cluster(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                   outExp_out, IC=IC, transform="none", byloc))
    outnonmaTlog.time <- system.time(outnonmaTlog <- nonma.cluster(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                           outExp_out, IC=IC, transform="log", byloc))
    
    
    outnonma_asymp.time <- system.time(outnonma_asymp <- nonma_asymp.cluster(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                                             outObs_out, IC=IC, transform="none", byloc))
    outnonma_asympTlog.time <- system.time(outnonma_asympTlog <- nonma_asymp.cluster(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                                                     outObs_out, IC=IC, transform="log", byloc))
    print("nonma finished")
    outbuck.theta.time <- system.time(outbuck.theta <- bucklandbounds.cluster(cluster_thetaa,
                                                                              clusterRR_ilarge,
                                                                            res,
                                                                            w_q=wslarge,
                                                                            outExp,
                                                                            IC=IC,
                                                                            transform="none",
                                                                            overdisp.est))
    outbuckTlog.theta.time  <- system.time(outbuckTlog.theta <- bucklandbounds.cluster(cluster_thetaa,
                                                                              clusterRR_ilarge,
                                                                              res,
                                                                              w_q=wslarge,
                                                                              outExp,
                                                                              IC=IC,
                                                                              transform="log",
                                                                              overdisp.est))
    print("buckland finished")
    
    outmaw2.theta.time <- system.time(outmaw2.theta <- maw2.cluster(cluster_thetaa,
                                                                    clusterRR_ilarge,
                                                                    w_q=wslarge,
                                                                    outExp,
                                                                    transform="none",
                                                                    overdisp.est))
    
    outmaw2Tlog.theta.time <- system.time(outmaw2Tlog.theta <-maw2.cluster(cluster_thetaa,
                                                                    clusterRR_ilarge,
                                                                    w_q=wslarge,
                                                                    outExp,
                                                                    transform="log",
                                                                    overdisp.est))
   
     print("maw2 finished")
    outmata.theta.time <- system.time(outmata.theta <- matabounds.cluster(cluster_thetaa,
                                                                          clusterRR_ilarge,
                                                                  w_q=wslarge,
                                                                  outExp,
                                                                  transform="none",
                                                                  overdisp.est))
    
    #thetaa,thetai, w_q, outExp, overdisp.est,transform=c("none","log", "sqrt")
    
    # outmataTlog.theta.time <- system.time(outmataTlog.theta <- matabounds.cluster(thetai=clusterRR_ilarge,
    #                                                                       thetaa = cluster_thetaa,
    #                                                                       w_q=wslarge,
    #                                                                       sparsematrix=t(sparsematrix ),
    #                                                                       outExp,
    #                                                                       overdisp.est = NULL,
    #                                                                       transform="log"))
    # print("mata finished")
    # 
    # 
    return(list(
        outnonma = outnonma,
        outnonmaTlog = outnonmaTlog,
        outnonma_asymp = outnonma_asymp,
        outnonma_asympTlog = outnonma_asympTlog,
        outbuck.theta = outbuck.theta,
        outbuckTlog.theta = outbuckTlog.theta,
        outmaw2.theta = outmaw2.theta,
        outmaw2Tlog.theta = outmaw2Tlog.theta,
        outmata.theta = outmata.theta,
        # outmataTlog.theta = outmataTlog.theta,
        
        outnonma.time = outnonma.time[[3]],
        outnonmaTlog.time = outnonmaTlog.time[[3]],
        outnonma_asymp.time = outnonma_asymp.time[[3]],
        outnonma_asympTlog.time = outnonma_asympTlog.time[[3]],
        outbuck.theta.time = outbuck.theta.time[[3]],
        outbuckTlog.theta.time = outbuckTlog.theta.time[[3]],
        outmaw2.theta.time = outmaw2.theta.time[[3]],
        outmaw2Tlog.theta.time = outmaw2Tlog.theta.time[[3]],
        outmata.theta.time = outmata.theta.time[[3]]#,
        # outmataTlog.theta.time = outmataTlog.theta.time[[3]]
    ))
}
    
calcbounds.cells <- function(idix, IC, res, byloc, Ex, Obs,wslarge, cluster_thetaa,clusterRR_ilarge, sparsematrix, cellsix){
    print("calcbounds.cells")
    cellrisk_wt_out <- rep(NA, length(cellsix))
    cellsix_out <- cellsix
    if(byloc==TRUE){
        print(paste0("byloc", byloc))
        outExp <- t(sparsematrix)%*%Ex
        outObs <- t(sparsematrix)%*%Obs
        outExp_out <- Ex[cellsix]
        outObs_out <- Obs[cellsix]
        for(i in 1:length(cellsix)){
            cellsixvec <- rep(0,dim(res$wLambda)[2])
            cellsixvec[cellsix[i]] <-1
            overlapid <- matrix(cellsixvec, nrow=1)%*%t(res$Lambda_dense)
            cellrisk_wt <- overlapid%*%wslarge 
            cellrisk_wt_out[[i]] <- cellrisk_wt
        }
    }
    else if (byloc==FALSE){
        print(paste0("byloc: ", byloc))
        outExp <- t(sparsematrix)%*%Ex
        outObs <- t(sparsematrix)%*%Obs
        print(paste0("cellsix: ", cellsix))
        outExp_out <- outExp@x[cellsix]
        outObs_out <- outObs@x[cellsix]
        cellrisk_wt_out <- NULL
    }    
        outnonma.time <- system.time(outnonma <- nonma.cells(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                       outExp_out, IC=IC, transform="none", byloc, cellrisk_wt_out, cellsix_out))
        outnonmaTlog.time <- system.time(outnonmaTlog <- nonma.cells(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                                outExp_out, IC=IC, transform="log", byloc, cellrisk_wt_out, cellsix_out))
        
        outnonma_asymp.time <- system.time(outnonma_asymp <- nonma_asymp.cells(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                                               outObs_out, IC=IC, transform="none", byloc, cellrisk_wt_out, cellsix_out))
        outnonma_asympTlog.time <- system.time(outnonma_asympTlog <- nonma_asymp.cells(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix,
                                                                                       outObs_out, IC=IC, transform="log", byloc, 
                                                                                       cellrisk_wt_out, cellsix_out))
        print("nonma finished")
        
        outbuck.theta.time <- system.time(outbuck.theta <- bucklandbounds.cells(cluster_thetaa,
                                                                                res,
                                                                                w_q=wslarge,
                                                                                Ex,
                                                                                IC=IC,
                                                                                transform="none",
                                                                                tsparsematrix=t(sparsematrix),
                                                                                overdisp.est, cellsix_out))
        outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- bucklandbounds.cells(cluster_thetaa,
                                                                                        res,
                                                                                        w_q=wslarge,
                                                                                        Ex,
                                                                                        IC=IC,
                                                                                        transform="log",
                                                                                        tsparsematrix=t(sparsematrix),
                                                                                        overdisp.est, cellsix_out))
        print("buckland finished")
        outmaw2.theta.time <- system.time(outmaw2.theta <- maw2.cells(cluster_thetaa,
                                                                      res,
                                                                w_q=wslarge,
                                                                Ex,
                                                                IC,
                                                                transform="none",
                                                                tsparsematrix=t(sparsematrix),
                                                                overdisp.est,
                                                                cellsix_out))
        outmaw2Tlog.theta.time <- system.time(outmaw2Tlog.theta <- maw2.cells(cluster_thetaa,
                                                                      res,
                                                                      w_q=wslarge,
                                                                      Ex,
                                                                      IC,
                                                                      transform="log",
                                                                      tsparsematrix=t(sparsematrix),
                                                                      overdisp.est,
                                                                      cellsix_out))

        # print("maw2 finished")
        outmata.theta.time <- system.time(outmata.theta <- matabounds.cells(cluster_thetaa,
                                                                            res,
                                                                      w_q=wslarge,
                                                                      Ex,
                                                                      IC,
                                                                      transform="none",
                                                                      tsparsematrix=t(sparsematrix),
                                                                      overdisp.est,
                                                                      cellsix_out))
        #thetaa, res, w_q,outExp, IC, tsparsematrix, overdisp.est,cellsix_out
        # outmataTlog.theta.time <- system.time(outmataTlog.theta <- matabounds.cells(thetai=clusterRR_ilarge,
        #                                                                       thetaa = cluster_thetaa,
        #                                                                       w_q=wslarge,
        #                                                                       sparsematrix=t(sparsematrix ),
        #                                                                       outExp,
        #                                                                       overdisp.est = NULL,
        #                                                                       transform="log"))
        # print("mata finished")
        
        
        return(list(
            outnonma = outnonma,
            outnonmaTlog = outnonmaTlog,
            outnonma_asymp = outnonma_asymp,
            outnonma_asympTlog = outnonma_asympTlog,
            outbuck.theta = outbuck.theta,
            outbuckTlog.theta = outbuckTlog.theta,
            outmaw2.theta = outmaw2.theta,
            outmaw2Tlog.theta = outmaw2Tlog.theta,
            outmata.theta = outmata.theta,
            # outmataTlog.theta = outmataTlog.theta,
            
            outnonma.time = outnonma.time[[3]],
            outnonmaTlog.time = outnonmaTlog.time[[3]],
            outnonma_asymp.time = outnonma_asymp.time[[3]],
            outnonma_asympTlog.time = outnonma_asympTlog.time[[3]],
            outbuck.theta.time = outbuck.theta.time[[3]],
            outbuckTlog.theta.time = outbuckTlog.theta.time[[3]],
            outmaw2.theta.time = outmaw2.theta.time[[3]],
            outmaw2Tlog.theta.time = outmaw2Tlog.theta.time[[3]],
            outmata.theta.time = outmata.theta.time[[3]]#,
            # outmataTlog.theta.time = outmataTlog.theta.time[[3]]
        ))
}

    
