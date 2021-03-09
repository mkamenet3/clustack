##########################################################################
##########################################################################
##########################################################################
#SCRIPT TO COMBINE CONDOR RESULTS
#M.Kamenetsky
#2021-02-28

#The goal of the script is to combine the condor results for both
#the main simulation and clustackbounds.
#The clustack condor sim results are stored on my external D drive (the big one)
##########################################################################
##########################################################################
##########################################################################
library(tidyverse)
##########################################################################
#CLUSTACK SIMULATION
##########################################################################
csvs <- list.files(path = "D:/RESEARCH/clustack/simulations/outcondor/singlelcuster/",
                 pattern="*.csv", full.names = TRUE)

# csvsInf <- list.files(path = "D:/RESEARCH/clustack/simulations/outcondor/",
#                    pattern="thetaInf_singlecluster_iter\\d*_startsim\\d*.csv", full.names = TRUE)
# csvs60 <- list.files(path = "D:/RESEARCH/clustack/simulations/outcondor/",
#                    pattern="theta60_singlecluster_iter\\d*_startsim\\d*.csv", full.names = TRUE)
clustacksim <- do.call(rbind, lapply(csvs, read.csv))
write.csv(clustacksim, file="../../data/clustacksim_raw.csv")
write.csv(clustacksim, file="D:/RESEARCH/clustack/simulations/outcondor/clustacksim_raw.csv")


clustackout <-   clustacksim %>%
    dplyr::select(-X) %>%
    rename(radius =X.1,
           RR = X.2,
           pop = X.3,
           theta = X.4) %>%
    mutate(incluster = if_else((method=="loc" | method=="pc") & incluster > 1e-1, 1,
                               if_else((method=="loc" | method=="pc") & incluster <= 1e-1, 0,
                                       if_else((method!="loc" | method!="pc"), incluster, NA_real_))),
           outcluster = if_else((method=="loc" | method=="pc") & outcluster > 1e-1, 1,
                                if_else((method=="loc" | method=="pc") & outcluster <= 1e-1,0,
                                        if_else((method!="loc" | method!="pc"), outcluster, NA_real_)))) %>%
    relocate(c(incluster, outcluster), .before=outcluster) %>% 
    #dplyr::select(-c(incluster.bic, outcluster.bic)) %>%
    dplyr::select(stsmodel:iter)
    
    
clustack2 <- clustackout %>%
    group_by(stsmodel, IC, radius, RR, pop, theta, timeperiod, method) %>%
    summarize(inclust = mean(incluster),
           outclust = mean(outcluster))
clustack2$method <- factor(clustack2$method, levels=c("loc", "pc", "clusso", "fstagewise", "stepscan"))


#PLOTS

powfp_plots <- function(mytheta){
    if(mytheta==Inf){
        mod <- "Poisson"
    } else{
        mod <- "Quasi-Poisson"
    }
    large_st <- ggplot(data=subset(clustack2, pop==150 & stsmodel=="ST" &  theta==mytheta)) +
         geom_line(col="blue",aes(x=RR, y=inclust, group=radius), size=1) +
         geom_line(col="red",aes(x=RR, y=outclust, group=radius), size=1) +
         facet_grid(radius~method+IC) + theme_bw() +
        ggtitle(paste0("Large Pop., Spatio-Temporal, ", mod)) +
        ylab("% Detection")
    
    small_st <- ggplot(data=subset(clustack2, pop==35 & stsmodel=="ST"&  theta==mytheta)) +
            geom_line(col="blue",aes(x=RR, y=inclust, group=radius), size=1) +
            geom_line(col="red",aes(x=RR, y=outclust, group=radius), size=1) +
            facet_grid(radius~method+IC) + theme_bw()+
        ggtitle(paste0("Small Pop., Spatio-Temporal, ", mod)) +
        ylab("% Detection")
    
    
    large_space <- ggplot(data=subset(clustack2, pop==150 & stsmodel=="space"&  theta==mytheta)) +
            geom_line(col="blue",aes(x=RR, y=inclust, group=radius), size=1) +
            geom_line(col="red",aes(x=RR, y=outclust, group=radius), size=1) +
            facet_grid(radius~method+IC) + theme_bw() +
        ggtitle(paste0("Large Pop., Spatial, ", mod)) +
        ylab("% Detection")
    
    small_space <- ggplot(data=subset(clustack2, pop==35 & stsmodel=="space"&  theta==mytheta)) +
            geom_line(col="blue",aes(x=RR, y=inclust, group=radius), size=1) +
            geom_line(col="red",aes(x=RR, y=outclust, group=radius), size=1) +
            facet_grid(radius~method+IC) + theme_bw()+
        ggtitle(paste0("Small Pop., Spatial, ", mod)) +
        ylab("% Detection")
    return(list(large_st,
                small_st,
                large_space,
                small_space))
    
}

##INF
powfp_plots(mytheta = Inf)
#QP
powfp_plots(mytheta = 60)
# test <- clustack2 %>%
#     filter(RR == 2, stsmodel=="ST",
#            radius==11, pop==150, method=="loc")



# ggplot(data=mastercoverage_select, aes(x=risk, y=prob)) +
#     geom_line(col="blue", size=2) +
#     facet_grid(radius~IC) +
#     xlab("Relative risk (RR)") +
#     ylab("MATA-Interval Coverage Probability") +
#     theme_bw()+
#     ggtitle("MATA-Interval Coverage") +
#     scale_x_continuous(breaks=c(1.1, 1.5,2)) +
#     ylim(0,1)+
#     theme(plot.title = element_text(hjust = 0.5),
#           plot.caption = element_text(face="italic")) +
#     labs(caption="Coverage probability: at least 1 ST cell was within the bounds. \nAll simulations considered (including those with no detection)")


#test <- read.csv("tocondor/fromcondor/thetaInf_singlecluster_iter4_startsim1.csv")
##########################################################################
#CLUSTACK BOUNDS
##########################################################################
csvs_bounds <- list.files(path = "D:/RESEARCH/clustack/simulations/outcondor/bounds/",
                   pattern="*.csv", full.names = TRUE)
#clustacksim_bounds <- do.call(rbind, lapply(csvs_bounds, read.csv))
#write.csv(clustacksim, file="D:/RESEARCH/clustack/simulations/outcondor/clustacksim_bounds_raw.csv")

clustackbounds2 <- clustacksim_bounds %>%
    dplyr::select(-X)%>%
    group_by(stsmodel, rad, risk, cent, theta, method) %>%
    summarize(across(nonma.bic.LB:matalog.aic.time, ~sum(.x, na.rm = TRUE)/100)) 
    # summarise_each(.vars = letters[1:4],
    #                .funs = c(mean="mean"))

clean <- clustackbounds2 %>%
    mutate(clusterMA_aic = clusterMA.aic.1,
           clusterMA_bic = clusterMA.bic.1,
           pcloc = method) %>%
    dplyr::select(-c(clusterMA.bic.1,clusterMA.bic.2,clusterMA.bic.3,
                     clusterMA.bic.4, clusterMA.bic.5, clusterMA.bic.6,
                     clusterMA.bic.7, clusterMA.bic.8, clusterMA.bic.9, clusterMA.bic.10,
                     clusterMA.aic.1,clusterMA.aic.2,clusterMA.aic.3,
                     clusterMA.aic.4, clusterMA.aic.5, clusterMA.aic.6,
                     clusterMA.aic.7, clusterMA.aic.8, clusterMA.aic.9, clusterMA.aic.10)) %>%
    pivot_longer(cols = ends_with("B") | ends_with("time"),
                 names_to="mtd",
                 values_to="Bound") %>%
    separate(mtd, sep="([.])", into=c("method", "IC", "LBUB"), remove=FALSE)  %>%
    mutate(method = recode_factor(method, 
                                  'outbuck'= "buck",
                                  "outbucklog" = "bucklog",
                                  "outmaw2" = "maw2",
                                  "outmaw2log" = "maw2log",
                                  "outmata" = "mata",
                                  "outmataTlog" = "matalog")) %>%
    dplyr::select(-mtd) %>%
    pivot_wider(names_from="LBUB",
                values_from="Bound") %>%
    mutate(method = as.factor(method)) %>%
    mutate(method= fct_relevel(method, "nonma", "nonma_asymp","buck","bucklog", "maw2","maw2log", "mata", "matalog"))


clean_all <- clustacksim_bounds %>%
    mutate(clusterMA_aic = clusterMA.aic.1,
           clusterMA_bic = clusterMA.bic.1,
           pcloc = method) %>%
    dplyr::select(-c(clusterMA.bic.1,clusterMA.bic.2,clusterMA.bic.3,
                     clusterMA.bic.4, clusterMA.bic.5, clusterMA.bic.6,
                     clusterMA.bic.7, clusterMA.bic.8, clusterMA.bic.9, clusterMA.bic.10,
                     clusterMA.aic.1,clusterMA.aic.2,clusterMA.aic.3,
                     clusterMA.aic.4, clusterMA.aic.5, clusterMA.aic.6,
                     clusterMA.aic.7, clusterMA.aic.8, clusterMA.aic.9, clusterMA.aic.10)) %>%
    pivot_longer(cols = ends_with("B") | ends_with("time"),
                 names_to="mtd",
                 values_to="Bound") %>%
    separate(mtd, sep="([.])", into=c("method", "IC", "LBUB"), remove=FALSE)  %>%
    mutate(method = recode_factor(method, 
                                  'outbuck'= "buck",
                                  "outbucklog" = "bucklog",
                                  "outmaw2" = "maw2",
                                  "outmaw2log" = "maw2log",
                                  "outmata" = "mata",
                                  "outmataTlog" = "matalog")) %>%
    dplyr::select(-mtd) %>%
    pivot_wider(names_from="LBUB",
                values_from="Bound") %>%
    mutate(method = as.factor(method)) %>%
    mutate(method= fct_relevel(method, "nonma", "nonma_asymp","buck","bucklog", "maw2","maw2log", "mata", "matalog")) 



#test <-  clean_all %>% group_by(stsmodel, risk, cent, theta, method, pcloc, IC) %>% mutate(ID = group_indices())
    

ggplot(data=subset(clean, (stsmodel=="space" & pcloc=="LOC" & cent==150 & theta==Inf &IC == "bic"))) +
    geom_line(aes(x=rad, y=LB), col="goldenrod", size=1) +
    geom_line(aes(x=rad, y=UB), col="goldenrod", size=1) +
    geom_line(aes(x=rad, y=clusterMA_bic), col="black", size=1) +
    geom_line(aes(x=rad, y=risk), col="black", size=1, linetype="dashed") +
    #geom_hline(aes(yintercept=clusterMA_aic, col="black"), size=1) +
    #geom_line(aes(x=ecount, y=maw2.LB, group=simID, color=as.factor(cellid))) +
    #geom_line(aes(x=ecount, y=cellrisk, group=simID), color="black") +
    facet_grid(as.factor(risk)~ method) + theme_bw() +
    ggtitle("Large Pop, Spatial, Quasi-Poisson, by PC+BIC") +
    ylab("Detection") + ylim(0,3)


# 
# ggplot(data=subset(test, (stsmodel=="ST" & pcloc=="PC" & cent==150 & theta==60))) +
#     geom_point(aes(x=rad, y=LB), col="goldenrod", size=1) +
#     #geom_line(aes(x=rad, y=UB), col="goldenrod", size=1) +
#     #geom_line(aes(x=rad, y=risk), col="black", size=1) +
#     #geom_hline(aes(yintercept=clusterMA_aic, col="black"), size=1) +
#     #geom_line(aes(x=ecount, y=maw2.LB, group=simID, color=as.factor(cellid))) +
#     #geom_line(aes(x=ecount, y=cellrisk, group=simID), color="black") +
#     facet_grid(as.factor(risk)~ method +IC) + theme_bw() +
#     ggtitle("Large Pop, Spatio-Temporal, Quasi-Poisson, by PC")
# 
