### 1) Install packages
#install.packages("survival")
#install.packages("tidyverse")
#install.packages("ggfortify")
#install.packages("survminer")


### 2) Load packages
library(survival)
library(tidyverse)
library(ggfortify)
library(survminer)

### 3) Activate survivalFun
survivalFun <- function(data=dat, title){
  
  noofgroups <- length(unique(data$trt))
  
  dat2 <- data[rep(1:nrow(data),data$no),]%>%select(-no) # expand the raw data into individual level
  rownames(dat2)<-1:nrow(dat2) # rename rows
  
  
  ### log-rank test
  LR <- survdiff(Surv(time, delta) ~ trt, data = dat2) 
  # estimate the survival function of the two groups and do the log-rank (LR) test
  
  LRp <- formatC((1 - pchisq(LR$chisq,df=1)), format = "e", digits = 3) 
  # calculate the LR test p-value and convert it to scientific notation
  
  LRpairwise <- pairwise_survdiff(Surv(time, delta) ~ trt, data = dat2)
  # calculate pairwise p-value using Benjamini-Hochberg (BH) method for multiple groups
  
  ### cox regression
  cox <- coxph(Surv(time, delta) ~ trt, data = dat2, ties = "breslow")
  coxsummary <- summary(cox)
  coxsummary
  
  LRTp <- formatC(coxsummary$logtest[3], format = "e", digits = 3) # likelihood ratio test p-value
  Waldp <- formatC(coxsummary$waldtest[3], format = "e", digits = 3) # wald test p-value
  Scorep <- formatC(coxsummary$sctest[3], format = "e", digits = 3) # score test p-value, which should be identical to log-rank when the each individual has unique event time
  
  tab <- data.frame("Likelihood-ratio"=LRTp, "Wald"=Waldp, "Score"=Scorep, "Log-rank"=LRp)
  rownames(tab) <- ""
  tab <- as.table(as.matrix(tab))
  
  KM <- survfit(Surv(time, delta) ~ trt, data = dat2)
  # estimate Kaplan-Meier (KM) survival function for plot
  KMsum <- summary(KM)


  plotcolor <- c("blue","orange","grey","purple","green","yellow","red")
  
  phcheck <- plot(KM,fun = "cloglog", col = plotcolor[1:noofgroups], main="Kaplan-Meier survival cloglog plot")
  
  plab <- ifelse(LRp<1e-50, paste("LR p-value:\n< 1e-50"), paste("LR p-value:\n",LRp))
  
  if(noofgroups <= 2){
    plot <- autoplot(KM, conf.int = F, surv.size = 1.5, censor.size = 5,
                      xlab = "Days", ylab = "Survival Probability", main = paste(title),
                      ylim = c(0,1)) + 
      theme_bw() + 
      scale_colour_manual(values = plotcolor[1:noofgroups])+
      guides(col=guide_legend(title = "Genetic Background")) +    # change the legend title
      theme(legend.position = "bottom",     # legend position
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16),
            legend.title = element_text(size=14),
            legend.text = element_text(size=12)) +
      annotate(geom="text", x=10, y=0.18, label=plab, size=4) 
  } else 
    {
    plot <- autoplot(KM, conf.int = F, surv.size = 1.5, censor.size = 5,
                      xlab = "Days", ylab = "Survival Probability", main = paste(title),
                      ylim = c(0,1)) + 
      theme_bw() + 
      scale_colour_manual(values = plotcolor[1:noofgroups]) + 
      guides(col=guide_legend(title = "Genetic Background",nrow=2, byrow=T)) +
      theme(legend.position = "bottom",     # legend position
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16),
            legend.title = element_text(size=14),
            legend.text = element_text(size=12),
            legend.key.size = unit(0.5,"line")) +
      annotate(geom="text", x=10, y=0.18, label=plab, size=4)
  }
  
  return(list(KMplot=plot, coxtable=tab, pairwise=LRpairwise, coxsummary=coxsummary, 
              KMsummary=KMsum, median=KM))
}
# survivalFun is a function that can be used to analyze the data in the right format


### 4) Choose working directory
#setwd("C:/Users/yuany/Desktop/Survival") # find the directory where the data located


### 5) Load data
#dat <- read.csv("Yolk-gal4-hAAT.csv", header = T) # read data
#dat <- read.csv("Lsp2-gal4-hAAT.csv", header = T)
#dat <- read.csv("Cg-gal4-hAAT.csv", header = T)
#dat <- read.csv("Lsp2-gal4-D6.csv", header = T)
dat <- read.csv("yolk.csv", header = T)


dat$trt <- relevel(dat$trt, ref = as.vector(unique(dat$trt)[grep("BackgroundxYolk-Gal4",unique(dat$trt))])) 
dat$trt <- relevel(dat$trt, ref ="Control") 
# relevel the trt column so that background with "CTR" is the reference level

### 6) Run analysis using data
output <- survivalFun(title = "Project 2.1 Trial 2 Female")
# after loaded the data, we can directly use survivalFun to do the analysis
# change the title for the plot
# the plot is to check whether the data is proportional hazard. This is the assumption for cox model
# the two lines should be parallel without any crossing if proportional hazard 
# if yes, we can interpret the coefficient as following; if not, we should only focus on log-rank test


output$coxtable
# this table contain the p-value for cox regression model (likelihood ratio test, wald test, and score test) and log-rank test

output$KMsummary
# estimated survival over time

output$KMplot
# this is the plot for survival function

output$pairwise
# this is the pairwise comparison for multiple groups with log-rank test, Benjamini-Hochberg (BH) method

output$coxsummary
# this is the summary of cox regression, where we find the coefficient
# interpretation of this parameter is: compared to strain1, strain 2 will increase/decrease the hazard of dying at any time by x times

output$median
# estimated median survival time with 95% CI

### 7) Output plot with publication quality
ggsave("Yolk.png", width=5, height=4, dpi=600)
# save the plot in png format, don't forget to change the png name




