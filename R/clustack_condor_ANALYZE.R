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
csvs <- list.files(path = "D:/RESEARCH/clustack/simulations/outcondor/",
                 pattern="thetaInf_singlecluster_iter\\d*_startsim\\d*.csv", full.names = TRUE)
clustacksim <- do.call(rbind, lapply(csvs, read.csv))

clustackout <-   clustacksim %>%
    dplyr::select(-X) %>%
    rename(radius =X.1,
           RR = X.2,
           pop = X.3,
           theta = X.4) %>%
    mutate(incluster = if_else((method=="loc" | method=="pc") & incluster.bic > 1e-1, 1,
                               if_else((method=="loc" | method=="pc") & incluster.bic <= 1e-1, 0,
                                       if_else((method!="loc" | method!="pc"), incluster.bic, NA_real_))),
           outcluster = if_else((method=="loc" | method=="pc") & outcluster.bic > 1e-1, 1,
                                if_else((method=="loc" | method=="pc") & outcluster.bic <= 1e-1,0,
                                        if_else((method!="loc" | method!="pc"), outcluster.bic, NA_real_)))) %>%
    relocate(c(incluster, outcluster), .before=outcluster.bic) %>% 
    #dplyr::select(-c(incluster.bic, outcluster.bic)) %>%
    dplyr::select(stsmodel:iter)
    
    
clustack2 <- clustackout %>%
    group_by(stsmodel, IC, radius, RR, pop, theta, timeperiod, method) %>%
    summarize(inclust = mean(incluster),
           outclust = mean(outcluster))
#PLOTS

(large_st <- ggplot(data=subset(clustack2, pop==150 & stsmodel=="ST")) +
    geom_line(col="blue",aes(x=RR, y=inclust, group=radius)) +
    geom_line(col="red",aes(x=RR, y=outclust, group=radius)) +
    facet_grid(radius~IC+method) + theme_bw())

(small_st <- ggplot(data=subset(clustack2, pop==35 & stsmodel=="ST")) +
    geom_line(col="blue",aes(x=RR, y=inclust, group=radius)) +
    geom_line(col="red",aes(x=RR, y=outclust, group=radius)) +
    facet_grid(radius~IC+method) + theme_bw())


(large_space <- ggplot(data=subset(clustack2, pop==150 & stsmodel=="space")) +
    geom_line(col="blue",aes(x=RR, y=inclust, group=radius)) +
    geom_line(col="red",aes(x=RR, y=outclust, group=radius)) +
    facet_grid(radius~IC+method) + theme_bw())

(small_space <- ggplot(data=subset(clustack2, pop==35 & stsmodel=="space")) +
    geom_line(col="blue",aes(x=RR, y=inclust, group=radius)) +
    geom_line(col="red",aes(x=RR, y=outclust, group=radius)) +
    facet_grid(radius~IC+method) + theme_bw())

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