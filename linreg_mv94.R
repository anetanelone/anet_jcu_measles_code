# Author: Amet J. N. Anelone
# Code for the interview for the position of Infectious Disease Modeler, vacancy reference  number 18889, James Cook University.
# Analysis of the dynamics of cell-associated viremia  following changes in measles infection doses
#Corresponding Journal article: Anelone, A.J. and Clapham, H.E., 2024. Measles Infection Dose Responses: Insights from Mathematical Modeling. Bulletin of Mathematical Biology, 86(7), p.85. DOI : 10.1007/s11538-024-01305-0
#cat(  "set-up") ##Set up
cat(" Set workingdirectory to source file location ")
setwd("~/Anet_MEV_NUS/mv_rcode")
library(ggplot2)
# Assign contrasting colors to improve data visualization, and individuals with color blindness.
# Each color corresponds to an infection dose TCID
custom_colors_TCID50 <- c( "red",   "blue", "orange","magenta"  , "green")
# cat("Get data") ##Get  data
data_MV1994 <- read.csv("datamv94.csv", header = TRUE, sep = ",") #  Read the file from current working directory
str(data_MV1994)
cat("Plot data") ##plot data
fig_data_MV1994 <-  ggplot( data = data_MV1994 ,  aes(x = day, y =  viremia,  color= as.factor(TCID50), fill=as.factor(TCID50),   group= as.factor(TCID50)) ) +
  geom_point( aes(  shape= as.factor(TCID50) ),  size=2 )+   
  scale_shape_manual( values=c( 18, 8,17,19,  15)) +
  scale_fill_manual(values= custom_colors_TCID50 ) +
  scale_color_manual(values= custom_colors_TCID50 ) +
  labs( shape = bquote(  TCID[50] ~ "="), fill = bquote(  TCID[50] ~ "="), color = bquote(  TCID[50] ~ "=") ) +
  #geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.35) +
  geom_line(  size =0.75) + 
  #geom_point(data = data94,  aes(x = t, y = mes,   color=as.factor(TCID)),   size=2) +
  #geom_line(mapping = aes( y = si,m_ivl),linetype = "dashed",size=2, color = "blue") + 
  labs(x = "Days post MV infection", y =bquote(  "Log" ~  TCID[50] ~  "/" ~  10^6 ~  "PBMC"))+
  ggtitle("Infectious viremia ")  +
  #annotate(geom="text", x=25, y=4, label=mcqii, color="black")+
  theme(legend.position="bottom")   +
  #theme(legend.position = "none") 
  theme(legend.text = element_text(colour="black",  face="bold")) +
  geom_hline(yintercept = 0.3, linetype = "dotdash", color = "darkgray", size=1) +
  scale_x_continuous(  breaks=c(0,  3,  5,6,  7,  9, 11, 13, 14, 17, 18 ))+
  scale_y_continuous(limits = c(0,4),breaks=c(0.3,1,2,  3,3.4,  4 ))
#Assess code from here ##Analysis  ############################################
cat("Data analysis of measles viral growth and decay to study measles infection dose response")
data_MV1994_uncensored <- subset( data_MV1994,  data_MV1994$viremia >0.3   )# Consider detectable viremia, by removing  data at the limit of detection 0.3
day_peak_10e4 <- 7;day_peak_10e2 <- 9;day_peak_10e0 <- 13 # day of the peak viremia
data_MV1994_uncensored$aligned_peakdays <-      data_MV1994_uncensored$day
data_MV1994_uncensored$aligned_peakdays[which( data_MV1994_uncensored$TCID50 ==10^4)] <-  (day_peak_10e2 -day_peak_10e4)+    data_MV1994_uncensored$day[which( data_MV1994_uncensored$TCID50 ==10^4)] 
data_MV1994_uncensored$aligned_peakdays[which( data_MV1994_uncensored$TCID50 ==10^0)] <-  (day_peak_10e2 -day_peak_10e0)+    data_MV1994_uncensored$day[which( data_MV1994_uncensored$TCID50 ==10^0)] 
group_viremia <- data.frame(   viremia10e4 = data_MV1994_uncensored$viremia[which( data_MV1994_uncensored$TCID50 ==10^4)],    viremia10e3 = data_MV1994_uncensored$viremia[which( data_MV1994_uncensored$TCID50 ==10^3)],   viremia10e2 = data_MV1994_uncensored$viremia[which( data_MV1994_uncensored$TCID50 ==10^2)], viremia10e1 = data_MV1994_uncensored$viremia[which( data_MV1994_uncensored$TCID50 ==10^1)] )  
 cor( group_viremia) # Correlation matrix for the viremia of each infection dose
lm_MV_growth  <-lm( data_MV1994_uncensored$viremia[ which(data_MV1994_uncensored$aligned_peakdays < day_peak_10e2 )]  ~ data_MV1994_uncensored$aligned_peakdays[ which(data_MV1994_uncensored$aligned_peakdays < day_peak_10e2 )]  )  
MV_growth_slope <- summary(lm_MV_growth  )$coefficients[2, "Estimate"]
MV_growth_slope_standard_errors <- summary(lm_MV_growth  )$coefficients[2, "Std. Error"]
MV_growth_slope_pvalue <-   summary(lm_MV_growth  )$coefficients[2, 4]
lm_MV_decline <-lm( data_MV1994_uncensored$viremia[ which(data_MV1994_uncensored$aligned_peakdays >= day_peak_10e2 )]  ~ data_MV1994_uncensored$aligned_peakdays[ which(data_MV1994_uncensored$aligned_peakdays >= day_peak_10e2 )]  )  
MV_decline_slope <- summary(lm_MV_decline  )$coefficients[2, "Estimate"]
MV_decline_slope_standard_errors <- summary(lm_MV_decline  )$coefficients[2, "Std. Error"]
MV_decline_slope_pvalue <-   summary(lm_MV_decline  )$coefficients[2, 4]
data_MV1994_uncensored$dynamics  <- ifelse( data_MV1994_uncensored$aligned_peakdays < day_peak_10e2,"Growth","Decline")
tau <- 1 # day according to estimations in the literature
MV_reproductiveratio <- (1+ MV_growth_slope / abs(MV_decline_slope) )* exp (MV_growth_slope *tau) #Estimation of the within-host reproductive ratio for measles using the derivation in Ribeiro, R.M., Qin, L., Chavez, L.L., Li, D., Self, S.G. and Perelson, A.S., 2010. Estimation of the initial viral growth rate and basic reproductive number during acute HIV-1 infection. Journal of virology, 84(12), pp.6096-6102.
cat("Plot results") ##plot results
fig_results <-       ggplot( data = data_MV1994_uncensored,  aes(x = aligned_peakdays, y = viremia,  color= as.factor(dynamics), fill=as.factor(dynamics),   group=as.factor(dynamics) ) ) +
                     geom_point( aes(  shape=as.factor(dynamics) ),  size=2 )+           scale_shape_manual( values=c( 16,   8)) +
                     labs( shape = bquote( "Dynamics"), fill = bquote( "Dynamics"), color = bquote("Dynamics") ) +
                     labs(x = " Day", y =bquote(  "Log" ~  TCID[50] ~  "/" ~  10^6 ~  "PBMC"))+     ggtitle("Aligning the peaks of viremia for MV infection dose response ")  +
                     theme(legend.position="bottom")   +       theme(legend.text = element_text(colour="black",  face="bold")) +
                     scale_y_continuous(limits = c(0,6), breaks=c(0.3,1,2,  3,3.4,  4 ))+      scale_x_continuous(  breaks=c(  5,6,  7,  9, 11, 13, 14 ))+
                     geom_smooth( method = "lm", se = TRUE) +
                     annotate("text", x = 9, y =5,  label = paste(" Growth rate =", round(MV_growth_slope, 2)),)+  annotate("text", x = 9, y =4,  label = paste(" Decline rate =", round(abs(MV_decline_slope), 2)),)+   annotate("text", x = 9, y =6,  label = paste("In vivo reproductive ratio =", round( MV_reproductiveratio, 2)),)
save_plot("f_data_MV1994_JCU.pdf",  fig_results  )
