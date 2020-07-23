#Analyze clustackbounds simulation
#M.Kamenetsky
#2020-07-18

master <- read.csv(master, file="clustackbounds_sim_bypc.csv")


##########################################################################################
#POST ANALYSIS
##########################################################################################
master2 <- as.data.frame(master)
master2$locid <- rep(1:1040, times=nsim)
master2$locidF <- as.factor(rep(1:1040, times=2))
clusterix <- which(rr==risk)
master2 <- master2 %>%
    mutate(clusterix = ifelse(locid %in% clusterix,1,0)) 
master2 <- master2 %>%
    mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB & clusterix==1), 1, 0),
           coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB & clusterix==1), 1, 0),
           coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB & clusterix==1), 1, 0),
           coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB & clusterix==1), 1, 0))

################################
#Coverage probability
################################



ggplot(data=master2,aes(x=locid, y=thetaa.bic, group=factor(simID), color=factor(simID))) +
    theme_bw() +
    geom_ribbon(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB, fill=factor(simID)), alpha=0.1) +
    geom_line()

master2 %>%
    arrange(thetaa.bic) %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic, group=factor(simID), color=factor(simID))) +
    theme_bw() +
    geom_ribbon(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB, fill=factor(simID)), alpha=0.1) +
    geom_line()

master2 %>%
    dplyr::filter(coverage.bic==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=locidF, y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB,color=factor(simID)), alpha=0.9, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("nsim=10, estimates & bounds where RR=2 within bounds (adequate coverage)")

master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.bic.log.LB, ymax = adjusted.bic.log.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (Buckland, BIC)")
master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.aic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.aic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.aic.log.LB, ymax = adjusted.aic.log.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (Buckland, AIC)")


master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.bic.LB, ymax = mata.bic.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (MATA, BIC)")
master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.aic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.aic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.aic.LB, ymax = mata.aic.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (MATA, AIC)")

##########################################
#coverage probability
##########################################
#make a quick and dirty plot
dat <- cbind.data.frame(thetaa=as.vector(modelthetas.bic$thetaa), UBa = as.vector(bounds.bic$ma_adjusted[[2]]), 
                        LBa=as.vector(bounds.bic$ma_adjusted[[1]]))

dat <- cbind.data.frame(UBa= UBa, LBa=LBa, thetaa=as.vector(thetaa))

dat <- dat %>%arrange(thetaa) %>%mutate(id = 1:nrow(.))
ggplot(data=dat) +
    geom_line(aes(x=id, y=thetaa), size=2) +
    theme_bw() +
    geom_line(aes(x=id, y=LBa), col="blue", size=1.5, alpha=0.5) +
    geom_line(aes(x=id, y=UBa), col="blue", size=1.5, alpha=0.5) +
    geom_hline(yintercept = 1, col="red", linetype="dashed", size=1) +
    geom_vline(xintercept = 823, col="red", linetype="dashed", size=1)


# #https://github.com/jlaake/RMark/blob/master/RMark/R/mata.wald.r
# thetai <- log(modelthetas.bic$thetai)
# thetaa <- log(modelthetas.bic$thetaa)
# se.thetai <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))
# w_q <- clusterlocs.bic$w_q
# thetaaa <- log(as.vector(thetaa)[827])
# thetaii <- log(thetai[,827])
# se.thetaii <- se.thetai[,827]
# 
# ########################################################
# #MATA-Intervals
# ########################################################
# thetai <- as.matrix(thetai)
# diff1 = apply(thetai,1, function(x) (as.vector(thetaa)-x))
# zquant <- sapply(1:ncol(diff1), function(x) diff1[,x]/se.thetai[x,])
# test1 <- apply(zquant, 2, function(x) pnorm(zquant))

######################################################
mataL <- uniroot(f=mata_tailareazscore, interval=c(-3, 3),
                 thetaii= thetaii,
                 se.thetaii=se.thetaii,
                 w_q=w_q, alpha=0.025, tol=1e-8)$root

mataU <- uniroot(f=mata_tailareazscore, interval=c(-3, 3),
                 thetaii= thetaii,
                 se.thetaii=se.thetaii,
                 w_q=w_q, alpha=1-0.025, tol=1e-8)$root
exp(cbind(mataL, mataU))



mata_solvebounds(thetaii = thetaii, se.thetaii = se.thetaii, w_q, alpha = 0.025, tol=1e-8)
mata_solvebounds(thetaii = thetaii, se.thetaii = se.thetaii, w_q, alpha = 0.025, tol=1e-8)
#do this for all ST locations

mata_st <- sapply(1:1040, function(x) mata_solvebounds(thetai[,x],
                                                       se.thetaii = se.thetai[,x],
                                                       w_q = w_q,
                                                       alpha = 0.025,
                                                       tol=1e-8))
mata_test <- sapply(1:10, function(x) mata_solvebounds(thetai[,x],
                                                       se.thetaii = se.thetai[,x],
                                                       w_q = w_q,
                                                       alpha = 0.025,
                                                       tol=1e-8))
mata_stt <- t(mata_st)
mata_df <- as.data.frame(mata_stt)
names(mata_df) <- c("mataLB", "mataUB")
mata_df$locid <- 1:1040

mata_df <- mata_df %>%
    mutate(clusterix = ifelse(locid %in% clusterix,1,0)) 
mata_df <- mata_df %>%
    mutate(coverage.bic = ifelse((risk < mataUB & risk > mataLB & clusterix==1), 1, 0))





