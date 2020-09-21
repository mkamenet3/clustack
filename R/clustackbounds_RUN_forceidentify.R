#Calculating Clustack Bounds
#Simulation Study
#M.Kamenetsky
#v1: 2020-06-30
#V2: 20200-7-18

rm(list=ls())
library(dplyr)
library(tidyr)
library(clusso)
library(ggplot2)
library(forcats)



####################################################################
set.seed(20200630)
####################################################################
#library(clusso)
library(MASS)
library(clusso)
source("clustack.R")
load("../data/japanbreastcancer.RData")
# source("R/clustack.R")
# load("../clustack/data/japanbreastcancer.RData")
cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath


#0) Setup
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]


#set global
rMax <- 20 
Time <- 5
#create set of potential clusters based on distances
potentialclusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparsematrix <- spacetimeMat(potentialclusters, numCenters, Time) 
#set maxclust
maxclust <- 12

##############################
#ARGS
##############################
#Arguments passed will only control theta parameter
##R CMD BATCH "--args arg1 arg2" myscript.R -> 
#args<- commandArgs(TRUE)
theta <- 60#as.numeric(args[1])
nsim <- 50#as.numeric(args[2])
#Fixed Params
cent <- 150
tim <- c(1:5)

#Sim-Over Params
risks <- c(1.1, 1.5, 2.0)
radii <- c(9, 11, 18)
# 
# # #tester
 risk = 1.5
# # cent = 150
 rad = 9
# # tim <- c(1:5)

########################################################################################################################
#20200811: Separating out model uncertainty from detection uncertainty
#Population Size
#RR
########################################################################################################################
clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
n <- length(x)
init <- setVectors(japanbreastcancer$period,
                   japanbreastcancer$expdeath, japanbreastcancer$death,
                   covars=NULL, Time=Time)
E1 <- init$E0
#Ex <- scale(init, Time)
Yx <- init$Y.vec
vectors <- list(Period = init$Year, Ex=init$E0, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)
n_uniq <- length(unique(clusters$center))
numCenters <- n_uniq
tmp <- clusters[clusters$center==cent,]
rr = matrix(1, nrow=n, ncol=Time)
###############################
#9 km cluster
###############################
risk <-1
cluster9 <- tmp[(tmp$r <= 9),]
rr9<-rr
rr9[cluster9$last, tim[1]:tail(tim, n=1)] <- NA
E1 <- as.vector(rr)*init$E0
incluster9 <- which(is.na(rr9))
YSIM <- rpois(length(E1), lambda = E1)
Ex <- unlist(scale_sim(list(YSIM), init, 1, Time))

poisson.test(x=sum(YSIM[incluster9])*2, T=sum(Ex[incluster9]), r=1)

ncell_out <- 0
for (i in 1:length(incluster9)){
    a <- poisson.test(x=sum(YSIM[incluster9[i]])*8, T=sum(Ex[incluster9[i]]), r=1)
    if(a$conf.int[1] < 1 & a$conf.int[2] > 1){
        print(a)
        ncell_out <- ncell_out +1
    }
}
poisson.test(x=sum(YSIM[incluster9])*8, T=sum(Ex[incluster9]), r=1)



###########################
#population
###########################
E1 <- init$E0*10
#Ex <- scale(init, Time)
Yx <- init$Y.vec
vectors <- list(Period = init$Year, Ex=init$E0, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)
n_uniq <- length(unique(clusters$center))
numCenters <- n_uniq
tmp <- clusters[clusters$center==cent,]
rr = matrix(1, nrow=n, ncol=Time)
YSIM <- rpois(length(E1), lambda = E1)
Ex <- unlist(scale_sim(list(YSIM), init, 1, Time))


poisson.test(x=sum(YSIM[incluster9])*2, T=sum(Ex[incluster9]), r=1)

ncell_out <- 0
for (i in 1:length(incluster9)){
    a <- poisson.test(x=sum(round(YSIM[incluster9[i]]*2)), T=sum(Ex[incluster9[i]]), r=1,
                      alternative = "greater")
    if(a$conf.int[1] < 1 & a$conf.int[2] > 1){
        print(a)
        ncell_out <- ncell_out +1
    }
    #print(a)
}
poisson.test(x=sum(round(YSIM[incluster9]*2)), T=sum(Ex[incluster9]), r=1, alternative = "greater")



###########################
#POWER
###########################
# count data, x= population mean
nPois<-function(x0,x1){        
    n<-4/(sqrt(x0)-sqrt(x1))^2    
    n    
}
nPois(1,2)


#takes background rate x into account
nbrPois<-function(x,x0,x1){        
    n<-4/(sqrt(x+x0)-sqrt(x+x1))^2    
    n    
}
nbrPois(1.5,1,2)


# n beyond background rate needed
nbbrPois<-function(x){             
    n<-4*sqrt(x)
    n
}
#if there are usually 20 cases in a pop, would need to see 18 more cases in the population

nbbrPois(20)
nbbrPois(sum(init$Y.vec[incluster9]))
#if there are usually 226 cases in the cluster, then I would need to see 60 more cases in the cluster to detect beyond the background rate

nRR<-function(rr, p0){
    n<-4/(p0*((sqrt(rr)-1)^2))
    n   
}

nRR(3, .01)



nRRoutcomes<-function(rr){
    m<-4/((sqrt(rr)-1)^2)
    m   
}
#So to detect a RR of 3, need (4/sqrt(3)-1)^2 ~=8 outcomes in the unexposed and 3*8 in exposed group
nRRoutcomes(3)

nRRoutcomes(2)
#I need 2*23 in exposed group

aa <- nbbrPois(init$E0[incluster9])
#prevalence:
nRRoutcomes(2)
nRRoutcomes(4)

#####################
#Power Analysis
#rate above background needed for a statistically significant result
set.seed(1)
N <- 1040
relrisk <-2
tmp <- clusters[clusters$center==cent,]
cluster9 <- tmp[(tmp$r <= 9),]
rr = matrix(1, nrow=n, ncol=Time)
rr9<-rr
rr9[cluster9$last, tim[1]:tail(tim, n=1)] <- relrisk
E1 <- as.vector(rr9)*init$E0
incluster9 <- which(rr9==relrisk)

#rr9[cluster9$last, tim[1]:tail(tim, n=1)] <- relrisk
reps <-100
ee1 <- as.vector(rr9)*init$E0
yy <- rpois(length(ee1), lambda = ee1)
ee <- unlist(scale_sim(list(yy), init, 1, Time))
#poisson.test(x=sum(yy[incluster9]),T=sum(ee[incluster9]), r=1, alternative = "greater")
yy1 <- replicate(reps,rpois(length(ee1), lambda = ee1))


x <-sum(yy[incluster9])
T <- sum(ee[incluster9])
alpha <- 0.05
#m <- r * T
p.L <- function(x, alpha) {
    if (x == 0) 
        0
    else qgamma(alpha, x)
}
pl <- p.L(x, 0.05)
p.U <- function(x, alpha) qgamma(1 - alpha, x + 1)
pu <- p.U(x, 0.05)

CINT <- c(p.L(x, alpha), p.U(x, alpha))/T
##############################################################
#asymptotic Poisson CI
t <- sum(ee[incluster9])
lambda <- sum(yy[incluster9])/t
round(lambda + c(-1, 1) * qnorm(0.975) * sqrt(lambda/t), 3)
#[1] 0.007 0.099

#lambdavals <- seq(0.005, 0.1, by = 0.01)
lambdavals <- seq(1,10, by=0.1)
nosim <- 1000
#t <- 100
coverage <- sapply(lambdavals, function(lambda) {
    lhats <- replicate(100,rpois(1040, lambda = lambda*ee)/ee)
    #lhats <- rpois(nosim, lambda = lambda * t)/t
    ll <- lhats - qnorm(0.975) * sqrt(lhats/ee)
    ul <- lhats + qnorm(0.975) * sqrt(lhats/ee)
    pow <- sapply(1:1040, function(i) mean(ll[i,] < lambda & ul[i,] > lambda))
    #apply(mean(ll < lambda & ul > lambda)
})

coverage2 <- sapply(lambdavals, function(lambda) {
    lhats <- rpois(100, lambda = lambda*sum(ee))/sum(ee)
    #lhats <- rpois(nosim, lambda = lambda * t)/t
    ll <- lhats - qnorm(0.975) * sqrt(lhats/sum(ee))
    ul <- lhats + qnorm(0.975) * sqrt(lhats/sum(ee))
    #pow <- sapply(1:1040, function(i) mean(ll[i,] < lambda & ul[i,] > lambda))
    (mean(ll < lambda & ul > lambda))
})

coverage_q <- coverage[incluster9,]
tcoverage <- t(coverage_q)
dat <- data.frame(tcoverage) %>%
    #slice(incluster9) %>%
    mutate(lambdavals = lambdavals,
           cluster = coverage2) %>%
    gather(key="variable", value="coverage", -lambdavals)
    
    
ggplot() +
    geom_line(data=subset(dat, variable!="cluster"), aes(x=lambdavals, y=coverage, group=variable), alpha=0.5) +
    geom_line(data=subset(dat, variable=="cluster"), aes(x=lambdavals, y=coverage, group=variable), col="red", size=1, alpha=0.7) +
    theme_minimal() +
    ggtitle("Proposed lambda vals")

###################################################################################################################
popdouble <- seq(1, 10, by=0.1)
#lambdavals <- seq(1,10, by=0.1)
nosim <- 1000
#t <- 100
lam <- 2
coverage <- sapply(popdouble, function(pops) {
    lhats <- replicate(100,rpois(1040, lambda = lam*(ee*pops))/(ee*pops))
    #lhats <- rpois(nosim, lambda = lambda * t)/t
    ll <- lhats - qnorm(0.975) * sqrt(lhats/(ee*pops))
    ul <- lhats + qnorm(0.975) * sqrt(lhats/(ee*pops))
    pow <- sapply(1:1040, function(i) mean(ll[i,] < lam & ul[i,] > lam))
    #apply(mean(ll < lambda & ul > lambda)
})

coverage2 <- sapply(lambdavals, function(lambda) {
    lhats <- rpois(100, lambda = lambda*sum(ee))/sum(ee)
    #lhats <- rpois(nosim, lambda = lambda * t)/t
    ll <- lhats - qnorm(0.975) * sqrt(lhats/sum(ee))
    ul <- lhats + qnorm(0.975) * sqrt(lhats/sum(ee))
    #pow <- sapply(1:1040, function(i) mean(ll[i,] < lambda & ul[i,] > lambda))
    (mean(ll < lambda & ul > lambda))
})

coverage_q <- coverage[incluster9,]
tcoverage <- t(coverage_q)
dat <- data.frame(tcoverage) %>%
    #slice(incluster9) %>%
    mutate(popdouble = popdouble,
           cluster = coverage2) %>%
    gather(key="variable", value="coverage", -popdouble)


ggplot() +
    geom_line(data=subset(dat, variable!="cluster"), aes(x=popdouble, y=coverage, group=variable), alpha=0.5) +
    geom_line(data=subset(dat, variable=="cluster"), aes(x=popdouble, y=coverage, group=variable), col="red", size=1, alpha=0.7) +
    theme_minimal() +
    ggtitle("Pop Doubling at lambda 2")



##############################
master <- NULL

for(rad in radii){
    for(risk in risks){
        for(s in 1:nsim){
            print(paste0("SimID: ",s,
                         "Risk: ", risk,
                         "Radius: ", rad))
            print(paste0("Theta: ", theta))
            print(str(theta))
            simid <- s
            #put the cluster in
            clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
            n <- length(x)
            init <- setVectors(japanbreastcancer$period,
                               japanbreastcancer$expdeath, japanbreastcancer$death,
                               covars=NULL, Time=Time)
            E1 <- init$E0
            #Ex <- scale(init, Time)
            Yx <- init$Y.vec
            vectors <- list(Period = init$Year, Ex=init$E0, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)
            n_uniq <- length(unique(clusters$center))
            numCenters <- n_uniq
            tmp <- clusters[clusters$center==cent,]
            cluster <- tmp[(tmp$r <= rad),]
            rr = matrix(1, nrow=n, ncol=Time)
            rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
            E1 <- as.vector(rr)*init$E0
            message(paste("Running model for periods",tim[1],"through", tail(tim, n=1)))
            #simulate data here
            #if (theta=="Inf"){
            if (is.infinite(theta)){
                print("Inf theta: pure poisson")
                #YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
                YSIM <- rpois(length(E1), lambda = E1)
            } else {
                print("OverD: NB")
                YSIM <- MASS::rnegbin(E1, theta = theta)
            }
            Ex <- unlist(scale_sim(list(YSIM), init, 1, Time))
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #SET-UP
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #run superclust by location for each sim (apply)
            if (is.infinite(theta)){
                overdisp.est <- NULL
            } else {
                offset_reg <- glm(YSIM ~ as.factor(rep(c("1","2","3","4","5"),each=length(Ex)/Time)) + offset(log(Ex)),
                                  family=quasipoisson)
                overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = TRUE)
                
            }
            outEx <- matrix(Ex, nrow=1)%*%sparsematrix
            outYx <- matrix(YSIM, nrow=1)%*%sparsematrix
            lambdahat <- YSIM/Ex
            
            Lambda <- as.vector(lambdahat)*sparsematrix #big Lambda matrix
            Lambda_dense <- as.matrix(Lambda)
            Lambda_dense[Lambda_dense == 0] <- 1
            
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #CLUSTER DETECTION BY Potential Cluster
            #####################################################################################
            #####################################################################################
            #####################################################################################
            sim_superclust_pc<- detectclusters(sparsematrix, Ex, YSIM,
                                               numCenters, Time, maxclust,
                                               bylocation = FALSE, model="poisson",
                                               overdisp.est = overdisp.est)
            #print(str(sim_superclust_pc))
            #BIC
            if(sim_superclust_pc$selection.bic!=0){
                clusterlocs.bic <- clusterlocs_ident(sim_superclust_pc, selection="selection.bic", sparsematrix)
                modelthetas.bic <- extract_thetai(clusterlocs.bic$param_ix, clusterlocs.bic$w_q, Lambda_dense)
                #calculate null and model curves
                bounds.bic.log <- calcbounds(Yx, Ex,
                                             modelthetas.bic$thetai, 
                                             modelthetas.bic$thetaa, 
                                             clusterlocs.bic$param_ix, 
                                             clusterlocs.bic$w_q,
                                             sparsematrix, 
                                             overdisp.est=overdisp.est)
            } else {
                modelthetas.bic <- list(thetai = rep(1, 1040), thetaa = rep(1, 1040))
                bounds.bic.log <- vector(mode = "list", length = 2)
                bounds.bic.log[[1]] <- rep(NA, 1040)
                bounds.bic.log[[2]] <- rep(NA, 1040)
                names(bounds.bic.log) <- c("ma_adjusted.LB", "ma_adjusted.UB")
            }
            print("BIC finished")
            
            if(sim_superclust_pc$selection.aic!=0) {
                if(sim_superclust_pc$selection.bic!=sim_superclust_pc$selection.aic){
                    #AIC
                    clusterlocs.aic <- clusterlocs_ident(sim_superclust_pc, selection="selection.aic", sparsematrix)
                    modelthetas.aic <- extract_thetai(clusterlocs.aic$param_ix, clusterlocs.aic$w_q, Lambda_dense)
                    #calculate null and model curves
                    bounds.aic.log <- calcbounds(Yx, Ex,
                                                 modelthetas.aic$thetai, 
                                                 modelthetas.aic$thetaa, 
                                                 clusterlocs.aic$param_ix, 
                                                 clusterlocs.aic$w_q,
                                                 sparsematrix, 
                                                 overdisp.est=overdisp.est) 
                } else {
                    print("BIC and AIC select the same path.")
                    bounds.aic.log <- bounds.bic.log
                    modelthetas.aic <- modelthetas.bic
                }
            } else {
                modelthetas.aic <- list(thetai = rep(1, 1040), thetaa = rep(1, 1040))
                bounds.aic.log <- vector(mode = "list", length = 2)
                #as.vector(bounds.bic.log$ma_adjusted.LB)
                bounds.aic.log[[1]] <- rep(NA, 1040)
                bounds.aic.log[[2]] <- rep(NA, 1040)
                names(bounds.aic.log) <- c("ma_adjusted.LB", "ma_adjusted.UB")
                
            }
            print("AIC finished")
            ########################
            #Calculate MATA-Bounds
            ########################
            #BIC
            if(sim_superclust_pc$selection.bic!=0){
                if(!is.null(overdisp.est)) {
                    se.thetai.bic <- sqrt(apply(modelthetas.bic$thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))*sqrt(overdisp.est)
                } else{
                    se.thetai.bic <- sqrt(apply(modelthetas.bic$thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))   
                }
                modelthetas_thetai.bic.log <- log(modelthetas.bic$thetai)
                mata_st.bic <- sapply(1:1040, function(x) mata_solvebounds(modelthetas_thetai.bic.log[,x],
                                                                           se.thetaii = se.thetai.bic[,x],
                                                                           w_q =  clusterlocs.bic$w_q,
                                                                           alpha = 0.025,
                                                                           tol=1e-8))
                mata_st.bic <- t(mata_st.bic)
            } else {
                mata_st.bic <- matrix(rep(NA, 1040*2),ncol=2)
            }
            print("BIC MATA finished")

            #AIC
            if(sim_superclust_pc$selection.aic!=0) {
                if(sim_superclust_pc$selection.bic!=sim_superclust_pc$selection.aic){
                    
                    if(!is.null(overdisp.est)) {
                        se.thetai.aic <- sqrt(apply(modelthetas.aic$thetai, 2, function(x) 1/outEx@x[clusterlocs.aic$param_ix]))*sqrt(overdisp.est)   
                    } else {
                        se.thetai.aic <- sqrt(apply(modelthetas.aic$thetai, 2, function(x) 1/outEx@x[clusterlocs.aic$param_ix]))    
                    }
                    modelthetas_thetai.aic.log <- log(modelthetas.aic$thetai)
                    mata_st.aic <- sapply(1:1040, function(x) mata_solvebounds(modelthetas_thetai.aic.log[,x],
                                                                               se.thetaii = se.thetai.aic[,x],
                                                                               w_q =  clusterlocs.aic$w_q,
                                                                               alpha = 0.025,
                                                                               tol=1e-8))
                    mata_st.aic <- t(mata_st.aic)
                    
                } else {
                    print("BIC and AIC select the same path.")
                    mata_st.aic <- mata_st.bic
                }
                
            } else {
                mata_st.aic <- matrix(rep(NA, 1040*2),ncol=2)
            }
            print("AIC MATA finished")
            #Create simrow
            row <- cbind(simID = rep(s,1040),
                         risk = rep(risk, 1040),
                         radius = rep(rad, 1040),
                         theta = rep(theta, 1040),
                         YSIM = YSIM,
                         Ex = Ex,
                         select.bic = rep(sim_superclust_pc$selection.bic,1040),
                         select.aic = rep(sim_superclust_pc$selection.aic,1040),
                         select.bic_orig = rep(sim_superclust_pc$selection.bic_orig,1040),
                         select.aic_orig = rep(sim_superclust_pc$selection.aic_orig,1040),
                         thetaa.bic= as.vector(modelthetas.bic$thetaa),
                         thetaa.aic = as.vector(modelthetas.aic$thetaa),
                         adjusted.bic.log.LB = as.vector(bounds.bic.log$ma_adjusted.LB),
                         adjusted.bic.log.UB = as.vector(bounds.bic.log$ma_adjusted.UB),
                         adjusted.aic.log.LB = as.vector(bounds.aic.log$ma_adjusted.LB),
                         adjusted.aic.log.UB = as.vector(bounds.aic.log$ma_adjusted.UB),
                         mata.bic.LB = mata_st.bic[,1],
                         mata.bic.UB = mata_st.bic[,2],
                         mata.aic.LB = mata_st.aic[,1],
                         mata.aic.UB = mata_st.aic[,2])
            #print(row)
            master <- rbind(master, row)
        }
    }
}



write.csv(master, file=paste0("clustackbounds_sim_bypc", "_theta_forceidentify", theta,".csv"))


