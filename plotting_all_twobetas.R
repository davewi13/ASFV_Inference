# Set herd number
herd_no <- "two_betas_freq"
nchains <- 50

source('//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Run_on_azog/guinatPriors.R')

process_chains <- function(herd){
  if(herd == 1){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd1/")
    list2env(herd1_guinat_results,globalenv())
  } else if (herd == 2){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd2/")
    list2env(herd2_guinat_results,globalenv())
  } else if (herd == 3){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd3/")
    list2env(herd3_guinat_results,globalenv())
  } else if (herd == 4){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd4/")
    list2env(herd4_guinat_results,globalenv())
  } else if (herd == 5){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd5/")
    list2env(herd5_guinat_results,globalenv())
  } else if (herd == 6){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd6/")
    list2env(herd6_guinat_results,globalenv())
  } else if (herd == 7){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd7/")
    list2env(herd7_guinat_results,globalenv())
  } else if (herd == 8){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd8/")
    list2env(herd8_guinat_results,globalenv())
  } else if (herd == 9){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd9/")
    list2env(herd9_guinat_results,globalenv())
  } else if (herd == 10){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Herd10/")
  } else if (herd == "exp_mean"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_exp_mean/")
  } else if (herd == "exp_shape"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_exp_shape/")
  } else if (herd == "inf_mean"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_inf_mean/")
  } else if (herd == "inf_shape"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_inf_shape/")
  } else if (herd == "beta"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_beta/")
  } else if (herd == "all_uninf"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_all_uninf/")
  } else if (herd == "two_betas"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_twobetas/")
  } else if (herd == "two_betas_freq"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_twobetas_freqdep/")
  } else if (herd == "simulated_data"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_simulated_data4/")
  } else if (herd == "fd_noext"){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_fd_noext/")
  } else {
    stop("Invalid herd number")
  }
  filelist <- list.files()
  for(i in 1:nchains){
    load(filelist[i])
    samps <- results$samps#[1:401,]
    #samps$no_initial <- results$no_initial
    if(i != 1){
      results_df <- rbind(results_df, samps)
    } else {
      results_df <- samps
    }
  }
  results_df$Chain <- rep(1:nchains,each=dim(samps)[1])
  results_df$Iteration <- rep(1:dim(samps)[1],times=nchains)
  return(results_df)
}

results_df <- process_chains(herd=herd_no)

# Priors used
x <- seq(0,100,0.01)
prior_beta_int <- dgamma(x,2,1)
prior_beta_ext <- dgamma(x,1,1/1500)
prior_exp_mean <- dgamma(x,10,10/6.25)
prior_exp_shape <- dgamma(x,5,5/19.39)
prior_inf_mean <- dgamma(x,10,10/9.12)
prior_inf_shape <- dgamma(x,5,5/22)
library(LaplacesDemon)

#################################################
#################################################
##                                             ##
##      Code chunk to plot simulated data      ##
##                                             ##
#################################################
#################################################

library(LaplacesDemon)
data_summary <- function(x) {
  m <- median(x)
  ymin <- p.interval(x)[1]
  ymax <- p.interval(x)[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

sim_true_beta_int <- data.frame(beta_int = c(0.4,0.3,0.25,0.25,0.35,0.45,0.3,0.35,0.4), Herd=factor(1:9))
sim_true_beta_ext <- data.frame(beta_ext = c(0.0002,0.0002,0.0001,0.0002,0.0001,0.0002,0.0001,0.0002,0.0001), Herd=factor(1:9))
gamma_means <- data.frame(means = c(6.25, 9.12), stage = c("Latent", "Infectious"))
gamma_shapes <- data.frame(shapes = c(19.39, 22), stage = c("Latent", "Infectious"))

results_df <- subset(results_df, Iteration >= 1500)

beta_int_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Prior = factor(rep(c(F,T),each=dim(results_df)[1]*9)),
                          Value = c(results_df$beta1_int,
                                    results_df$beta2_int,
                                    results_df$beta3_int,
                                    results_df$beta4_int,
                                    results_df$beta5_int,
                                    results_df$beta6_int,
                                    results_df$beta7_int,
                                    results_df$beta8_int,
                                    results_df$beta9_int,
                                    rgamma(50200*9,shape=2,rate=1)))

beta_ext_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_ext,
                                    results_df$beta2_ext,
                                    results_df$beta3_ext,
                                    results_df$beta4_ext,
                                    results_df$beta5_ext,
                                    results_df$beta6_ext,
                                    results_df$beta7_ext,
                                    results_df$beta8_ext,
                                    results_df$beta9_ext))

means_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                        Iteration = rep(1:dim(results_df)[1],4),
                        Prior = factor(rep(c(F,T),each=dim(results_df)[1]*2)),
                        Value = c(results_df$exp_mean,
                                  results_df$inf_mean,
                                  rgamma(50200,shape=10,rate=10/6.25),
                                  rgamma(50200,shape=10,rate=10/9.12)))

shapes_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                       Iteration = rep(1:dim(results_df)[1],4),
                       Prior = factor(rep(c(F,T),each=dim(results_df)[1]*2)),
                       Value = c(results_df$exp_shape,
                                 results_df$inf_shape,
                                 rgamma(50200,shape=5,rate=5/19.39),
                                 rgamma(50200,shape=5,rate=5/22)))

p1 <- ggplot(subset(beta_int_df, Prior == "FALSE"), aes(x=Herd, y=Value)) +
  geom_violin(position="identity", alpha=0.5, col="#56B4E9", fill="#56B4E9") +
  #scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  #scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Within-herd transmission rate '*beta['2'])) +
  theme_bw() +
  ylim(0,1) +
  theme(legend.position = "none") +
  geom_point(data = sim_true_beta_int, mapping = aes(x=Herd, y=beta_int), col="black", shape=3, size=3)

p2 <- ggplot(beta_ext_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('External transmission rate '*beta['1'])) +
  theme_bw() +
  geom_point(data = sim_true_beta_ext, mapping = aes(x=Herd, y=beta_ext), col="black", shape=3, size=3)

p3 <- ggplot(means_df, aes(x=Stage, y=Value, col=Prior)) +
  geom_violin(aes(fill=Prior), position="identity", alpha=0.5) +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(data=subset(means_df, Prior == "FALSE"), fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Mean stage duration') + 
  theme_bw() +
  ylim(0,20) +
  theme(legend.position = "none") +
  geom_point(data = gamma_means, mapping = aes(x=stage, y=means), col="black", shape=3, size=3)

p4 <- ggplot(shapes_df, aes(x=Stage, y=Value, col=Prior)) +
  geom_violin(aes(fill=Prior), position="identity", alpha=0.5) +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(data = subset(shapes_df, Prior == "FALSE"), fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Shape parameter') + 
  ylim(0,50) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = gamma_shapes, mapping = aes(x=stage, y=shapes), col="black", shape=3, size=3)

library(gridExtra)
png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/Sim_study_results.png", height=6, width=6, units="in", res=300)
grid.arrange(grobs=list(p1,p3,p2,p4), widths=c(2,1), layout.matrix=matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T))
dev.off()



#################################################
#################################################
##                                             ##
##         Code chunk to plot ASF data         ##
##                                             ##
#################################################
#################################################


results_df <- subset(results_df, Iteration > 3000)
R0_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_int*results_df$inf_mean,
                                    results_df$beta2_int*results_df$inf_mean,
                                    results_df$beta3_int*results_df$inf_mean,
                                    results_df$beta4_int*results_df$inf_mean,
                                    results_df$beta5_int*results_df$inf_mean,
                                    results_df$beta6_int*results_df$inf_mean,
                                    results_df$beta7_int*results_df$inf_mean,
                                    results_df$beta8_int*results_df$inf_mean,
                                    results_df$beta9_int*results_df$inf_mean))

beta_int_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_int,
                                    results_df$beta2_int,
                                    results_df$beta3_int,
                                    results_df$beta4_int,
                                    results_df$beta5_int,
                                    results_df$beta6_int,
                                    results_df$beta7_int,
                                    results_df$beta8_int,
                                    results_df$beta9_int))

beta_ext_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_ext,
                                    results_df$beta2_ext,
                                    results_df$beta3_ext,
                                    results_df$beta4_ext,
                                    results_df$beta5_ext,
                                    results_df$beta6_ext,
                                    results_df$beta7_ext,
                                    results_df$beta8_ext,
                                    results_df$beta9_ext))

means_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                       Iteration = rep(1:dim(results_df)[1],2),
                       Value = c(results_df$exp_mean,
                                 results_df$inf_mean))

shapes_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                       Iteration = rep(1:dim(results_df)[1],2),
                       Value = c(results_df$exp_shape,
                                 results_df$inf_shape))

p1 <- ggplot(R0_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('R'[0])) +
  ylim(0,40) +
  geom_hline(yintercept=9.12) +
  theme_bw()
  
p2 <- ggplot(beta_int_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Within-herd transmission '*beta['2'])) +
  ylim(0,6) +
  geom_hline(yintercept=1) +
  geom_hline(yintercept=0.24, lty=2) +
  geom_hline(yintercept=5.57, lty=2) +
  theme_bw()

p3 <- ggplot(beta_ext_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('External transmission '*beta['1'])) +
  theme_bw() 

p4 <- ggplot(means_df, aes(x=Stage, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Mean stage duration') +
  geom_segment(aes(x = 1.5, y = 6.25, xend = 2.5, yend = 6.25),col="black") +
  geom_segment(aes(x = 1.5, y = 3, xend = 2.5, yend = 3),col="black", lty=2) +
  geom_segment(aes(x = 1.5, y = 10.7, xend = 2.5, yend = 10.7),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 9.12, xend = 1.5, yend = 9.12),col="black") +
  geom_segment(aes(x = 0.5, y = 4.37, xend = 1.5, yend = 4.37),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 15.58, xend = 1.5, yend = 15.58),col="black", lty=2) +
  theme_bw() 

p5 <- ggplot(shapes_df, aes(x=Stage, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Shape parameter') +
  geom_segment(aes(x = 1.5, y = 22, xend = 2.5, yend = 22),col="black") +
  geom_segment(aes(x = 1.5, y = 7.14, xend = 2.5, yend = 7.14),col="black", lty=2) +
  geom_segment(aes(x = 1.5, y = 45.06, xend = 2.5, yend = 45.06),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 19.39, xend = 1.5, yend = 19.39),col="black") +
  geom_segment(aes(x = 0.5, y = 6.30, xend = 1.5, yend = 6.30),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 39.72, xend = 1.5, yend = 39.72),col="black", lty=2) +
  theme_bw() 


png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/ASF_results.png", height=6, width=6, units="in", res=300)
grid.arrange(grobs=list(p1,p4,p2,p5,p3), widths=c(2,1), layout.matrix=matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T))
dev.off()

#################################################
#################################################
##                                             ##
##       Code chunk to plot uninf priors       ##
##                                             ##
#################################################
#################################################


results_df <- process_chains("exp_mean")
results_df$exp_shape <- process_chains("exp_shape")$exp_shape
results_df$inf_mean <- process_chains("inf_mean")$inf_mean
results_df$inf_shape <- process_chains("inf_shape")$inf_shape
results_df[,10:18] <- process_chains("beta")[,10:18]
results_df <- subset(results_df, Iteration > 2000)

R0_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                    Iteration = rep(1:dim(results_df)[1], 9),
                    Value = c(results_df$beta1_int*results_df$inf_mean,
                              results_df$beta2_int*results_df$inf_mean,
                              results_df$beta3_int*results_df$inf_mean,
                              results_df$beta4_int*results_df$inf_mean,
                              results_df$beta5_int*results_df$inf_mean,
                              results_df$beta6_int*results_df$inf_mean,
                              results_df$beta7_int*results_df$inf_mean,
                              results_df$beta8_int*results_df$inf_mean,
                              results_df$beta9_int*results_df$inf_mean))

beta_int_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_int,
                                    results_df$beta2_int,
                                    results_df$beta3_int,
                                    results_df$beta4_int,
                                    results_df$beta5_int,
                                    results_df$beta6_int,
                                    results_df$beta7_int,
                                    results_df$beta8_int,
                                    results_df$beta9_int))

beta_ext_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_ext,
                                    results_df$beta2_ext,
                                    results_df$beta3_ext,
                                    results_df$beta4_ext,
                                    results_df$beta5_ext,
                                    results_df$beta6_ext,
                                    results_df$beta7_ext,
                                    results_df$beta8_ext,
                                    results_df$beta9_ext))

means_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                       Iteration = rep(1:dim(results_df)[1],2),
                       Value = c(results_df$exp_mean,
                                 results_df$inf_mean))

shapes_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                        Iteration = rep(1:dim(results_df)[1],2),
                        Value = c(results_df$exp_shape,
                                  results_df$inf_shape))

p1 <- ggplot(R0_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('R'[0])) +
  ylim(0,40) +
  geom_hline(yintercept=9.12) +
  theme_bw()

p2 <- ggplot(beta_int_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Within-herd transmission '*beta['2'])) +
  ylim(0,6) +
  geom_hline(yintercept=10) +
  geom_hline(yintercept=0.5, lty=2) +
  geom_hline(yintercept=19.5, lty=2) +
  theme_bw()

p3 <- ggplot(beta_ext_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('External transmission '*beta['1'])) +
  theme_bw() 

p4 <- ggplot(means_df, aes(x=Stage, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Mean stage duration') +
  geom_segment(aes(x = 1.5, y = 6.25, xend = 2.5, yend = 6.25),col="black") +
  geom_segment(aes(x = 1.5, y = 0.75, xend = 2.5, yend = 0.75),col="black", lty=2) +
  geom_segment(aes(x = 1.5, y = 17.4, xend = 2.5, yend = 17.4),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 9.12, xend = 1.5, yend = 9.12),col="black") +
  geom_segment(aes(x = 0.5, y = 1.1, xend = 1.5, yend = 1.1),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 25.4, xend = 1.5, yend = 25.4),col="black", lty=2) +
  theme_bw() 

p5 <- ggplot(shapes_df, aes(x=Stage, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Shape parameter') +
  geom_segment(aes(x = 1.5, y = 22, xend = 2.5, yend = 22),col="black") +
  geom_segment(aes(x = 1.5, y = 2.7, xend = 2.5, yend = 2.7),col="black", lty=2) +
  geom_segment(aes(x = 1.5, y = 61.3, xend = 2.5, yend = 61.3),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 19.39, xend = 1.5, yend = 19.39),col="black") +
  geom_segment(aes(x = 0.5, y = 2.4, xend = 1.5, yend = 2.4),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 54, xend = 1.5, yend = 54),col="black", lty=2) +
  theme_bw() 


png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/ASF_results_uninf_priors.png", height=6, width=6, units="in", res=300)
grid.arrange(grobs=list(p1,p4,p2,p5,p3), widths=c(2,1), layout.matrix=matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T))
dev.off()


#################################################
#################################################
##                                             ##
##       Code chunk to plot no ext trans       ##
##                                             ##
#################################################
#################################################

results_df <- subset(results_df, Iteration > 3000)

R0_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                    Iteration = rep(1:dim(results_df)[1], 9),
                    Value = c(results_df$beta1_int*results_df$inf_mean,
                              results_df$beta2_int*results_df$inf_mean,
                              results_df$beta3_int*results_df$inf_mean,
                              results_df$beta4_int*results_df$inf_mean,
                              results_df$beta5_int*results_df$inf_mean,
                              results_df$beta6_int*results_df$inf_mean,
                              results_df$beta7_int*results_df$inf_mean,
                              results_df$beta8_int*results_df$inf_mean,
                              results_df$beta9_int*results_df$inf_mean))

beta_int_df <- data.frame(Herd = factor(rep(1:9,each=dim(results_df)[1])),
                          Iteration = rep(1:dim(results_df)[1], 9),
                          Value = c(results_df$beta1_int,
                                    results_df$beta2_int,
                                    results_df$beta3_int,
                                    results_df$beta4_int,
                                    results_df$beta5_int,
                                    results_df$beta6_int,
                                    results_df$beta7_int,
                                    results_df$beta8_int,
                                    results_df$beta9_int))
beta_int_df <- beta_int_df[beta_int_df$Value != 0,]

means_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                       Iteration = rep(1:dim(results_df)[1],2),
                       Value = c(results_df$exp_mean,
                                 results_df$inf_mean))

shapes_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                        Iteration = rep(1:dim(results_df)[1],2),
                        Value = c(results_df$exp_shape,
                                  results_df$inf_shape))

p1 <- ggplot(R0_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('R'[0])) +
  ylim(0,40) +
  geom_hline(yintercept=9.12) +
  theme_bw()

p2 <- ggplot(beta_int_df, aes(x=Herd, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Within-herd transmission '*beta['2'])) +
  ylim(0,6) +
  geom_hline(yintercept=1) +
  geom_hline(yintercept=5.57, lty=2) +
  geom_hline(yintercept=0.24, lty=2) +
  theme_bw()

p3 <- ggplot(means_df, aes(x=Stage, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Mean stage duration') +
  geom_segment(aes(x = 1.5, y = 6.25, xend = 2.5, yend = 6.25),col="black") +
  geom_segment(aes(x = 1.5, y = 0.75, xend = 2.5, yend = 0.75),col="black", lty=2) +
  geom_segment(aes(x = 1.5, y = 17.4, xend = 2.5, yend = 17.4),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 9.12, xend = 1.5, yend = 9.12),col="black") +
  geom_segment(aes(x = 0.5, y = 1.1, xend = 1.5, yend = 1.1),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 25.4, xend = 1.5, yend = 25.4),col="black", lty=2) +
  theme_bw() 

p4 <- ggplot(shapes_df, aes(x=Stage, y=Value)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Shape parameter') +
  geom_segment(aes(x = 1.5, y = 22, xend = 2.5, yend = 22),col="black") +
  geom_segment(aes(x = 1.5, y = 2.7, xend = 2.5, yend = 2.7),col="black", lty=2) +
  geom_segment(aes(x = 1.5, y = 61.3, xend = 2.5, yend = 61.3),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 19.39, xend = 1.5, yend = 19.39),col="black") +
  geom_segment(aes(x = 0.5, y = 2.4, xend = 1.5, yend = 2.4),col="black", lty=2) +
  geom_segment(aes(x = 0.5, y = 54, xend = 1.5, yend = 54),col="black", lty=2) +
  theme_bw() 

png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/ASF_noext.png", height=6, width=6, units="in", res=300)
grid.arrange(grobs=list(p1,p3,p2,p4), widths=c(2,1), layout.matrix=matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T))
dev.off()

#################################################
#################################################
##                                             ##
##  Code chunk to plot fitted distributions    ##
##                                             ##
#################################################
#################################################

x <- seq(1,16,0.01)
y.latent <- dgamma(x, shape=18.4, rate=18.4/8.1)
y.infectious <- dgamma(x, shape=17.4, rate=17.4/8.8)
df <- data.frame(x=x, y.latent=y.latent, y=y.infectious)
library(gridExtra)

xmin <- qgamma(0.025, shape=18.4, rate=18.4/8.1)
xmax <- qgamma(0.975, shape=18.4, rate=18.4/8.1)
shade <- rbind(c(xmin,0),
               subset(df, (x > xmin) & (x < xmax)), 
               c(xmax, 0))
ymin <- 0
ymax <- max(y.latent)
p1 <- ggplot(df, aes(x=x, y=y.latent)) + 
  geom_line(col="#56B4E9") +
  ylab("Density") + 
  xlab("Days") + 
  ggtitle("Latent period") +
  theme_bw()#+ 
  #geom_polygon(data = shade, aes(x, y.latent), col="#56B4E9", fill="#56B4E9", alpha=0.5)

df <- data.frame(x=x, y.infectious=y.infectious)
xmin <- qgamma(0.025, shape=17.4, rate=17.4/8.8)
xmax <- qgamma(0.975, shape=17.4, rate=17.4/8.8)
shade <- rbind(c(xmin,0),
               subset(df, (x > xmin) & (x < xmax)), 
               c(xmax, 0))
ymin <- 0
ymax <- max(y.infectious)
p2 <- ggplot(df, aes(x=x, y=y.infectious)) + 
  geom_line(col="#56B4E9") +
  ylab("") + 
  xlab("Days") + 
  ggtitle("Infectious period") +
  theme_bw()# + 
  #geom_polygon(data = shade, aes(x, y.infectious), col="#56B4E9", fill="#56B4E9", alpha=0.5)

png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/Fitted_distributions.png", height=4.5, width=9, units="in", res=300)
grid.arrange(grobs=list(p1,p2),ncol=2)
dev.off()

x <- seq(1,16,0.01)
y.latent <- dgamma(x, shape=18.4, rate=18.4/8.1)
y.infectious <- dgamma(x, shape=17.4, rate=17.4/8.8)
df <- data.frame(x=rep(x,2), y=c(y.latent, y.infectious), type=rep(c("Latent","Infectious"), each=length(x)))
p1 <- ggplot(df, aes(x=x, y=y, col=type)) + 
  geom_line() +
  ylab("Density") + 
  xlab("Days") + 
  theme_bw() +
  scale_color_discrete(name = "Stage")

png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/Fitted_distributions.png", height=3, width=4.5, units="in", res=300)
p1
dev.off()

#################################################
#################################################
##                                             ##
##  Find proportion infected at positive ASFV  ##
##                                             ##
#################################################
#################################################


process_chains <- function(){
    setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_twobetas_freqdep/")

    df <- data.frame(Herd = rep(1:9,50), Proportion=NA)
    filelist <- list.files()
      for(i in 1:nchains){
        load(filelist[i])
        times <- results$times
        for(j in 1:length(times)){
          temp <- times[[j]]
            if(j == 1){
              infection <- sum(temp$Infection < 1)/1614#22)#/1614
            } else if(j == 2){
              infection <- sum(temp$Infection < 2)/1949#11)#/1949
            } else if(j == 3){
              infection <- sum(temp$Infection < 5)/1753#11)#/1753
            } else if(j == 4){
              infection <- sum(temp$Infection < 1)/1833#8)#/1833
            } else if(j == 5){
              infection <- sum(temp$Infection < 4)/1320#8)#/1320
            } else if(j == 6){
              infection <- sum(temp$Infection < 5)/600#9)#/600
            } else if(j == 7){
              infection <- sum(temp$Infection < 1)/600#9)#/600
            } else if(j == 8){
              infection <- sum(temp$Infection < 6)/600#9)#/600
            } else if(j == 9){
              infection <- sum(temp$Infection < 2)/2145#17)#/2145
            }
          df$Proportion[(i-1)*9+j] <- infection
        }
      }
      return(df)
}

results_df <- process_chains()
results_df$Herd <- factor(results_df$Herd)
library(dplyr)
results_df %>%
  group_by(Herd) %>%
  summarise(Proportion = median(Proportion))

# Output data for table
apply(results_df[,10:18], 2, data_summary)







# Set herd number
nchains <- 50
source('//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Run_on_azog/guinatPriors.R')

process_chains <- function(){

  results_df <- data.frame(Herd=factor(rep(1:9,50)), Intro_times=NA, Intro_n=NA)
  setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_fd_noext/")

  filelist <- list.files()
  for(i in 1:nchains){
    load(filelist[i])
    results_df$Intro_times[((i-1)*9+1):(i*9)] <- lapply(results$times, min)
    for(j in 1:9){
      temp <- results$times[[j]]
      results_df$Intro_n[(9*(i-1)+j)] <- sum(temp$Exposure < min(temp$Infection))
    }
  }
  return(results_df)
}

results_df <- process_chains()
results_df$Intro_times <- as.numeric(results_df$Intro_times)
results_df$Intro_n <- as.numeric(results_df$Intro_n)

p1 <- ggplot(results_df, aes(x=Herd, y=Intro_times)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Time of introduction, T'[0])) +
  ylim(-35,0) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("")

p2 <- ggplot(results_df, aes(x=Herd, y=Intro_n)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Initial external infections, n'[0])) +
  theme_bw()

library(gridExtra)
png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/Introduction_inferences_noext.png", height=6, width=9, units="in", res=600)
p1
dev.off()

# Redo relative to positive diagnosis
# Set herd number
nchains <- 50
source('//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Run_on_azog/guinatPriors.R')

process_chains <- function(){
  
  results_df <- data.frame(Herd=factor(rep(1:9,50)), Timepoint="Diagnosis", Intro_times=NA, Intro_n=NA)
  #results_df <- rbind(results_df, data.frame(Herd=factor(rep(1:9,50)), Timepoint="First Mortality", Intro_times=NA, Intro_n=NA))
  setwd("C:/Users/dewing/HPC_results/Fixed_multiple_inits/Combined_herds_twobetas_freqdep/")
  
  filelist <- list.files()
  for(i in 1:nchains){
    load(filelist[i])
    results_df$Intro_times[((i-1)*9+1):(i*9)] <- as.numeric(lapply(results$times, min))
    for(j in 1:9){
      temp <- results$times[[j]]
      results_df$Intro_n[(9*(i-1)+j)] <- sum(temp$Exposure < min(temp$Infection))
    }
  }
  #for(i in (nchains+1):(2*nchains)){
  #  load(filelist[i-nchains])
  #  results_df$Intro_times[((i-1)*9+1):(i*9)] <- as.numeric(lapply(results$times, min)) - c(0,1,4,0,2,4,0,5,1)
  #  for(j in 1:9){
  #    temp <- results$times[[j]]
  #    results_df$Intro_n[(9*(i-1)+j)] <- sum(temp$Exposure < min(temp$Infection))
  #  }
  #}
  return(results_df)
}

results_df <- process_chains()
results_df$Intro_times <- as.numeric(results_df$Intro_times)
results_df$Intro_n <- as.numeric(results_df$Intro_n)

pd = position_dodge(0.9)
#p1 <- ggplot(results_df, aes(x=Herd, y=Intro_times, col=Timepoint, fill=Timepoint)) +
#  geom_violin(alpha=0.5) + 
#  stat_summary(fun.data=data_summary, geom="pointrange", col="black", position=pd) +
#  ylab(expression('Time of introduction, T'[0])) +
#  ylim(-45,0) +
#  theme_bw() +
#  theme(axis.text.x = element_blank(),
#        axis.ticks.x = element_blank()) +
#  xlab("") +
#  scale_fill_manual(values=c("#56B4E9", "#FF99CC")) +
#  scale_color_manual(values=c("#56B4E9", "#FF99CC"))

p1 <- ggplot(results_df, aes(x=Herd, y=Intro_times)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Time of introduction, T'[0])) +
  ylim(-25,0) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("")

time.before.diagnosis <- data.frame(Herd=factor(1:9), Time=c(36,23,20,22,21,19,20,19,32))
p1 <- p1 +
  geom_text(data=time.before.diagnosis, aes(x=Herd, label=Time, y=-2))

p2 <- ggplot(results_df, aes(x=Herd, y=Intro_n)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Initial external infections, n'[0])) +
  theme_bw()

library(gridExtra)
png("//212.219.57.223/HomeDirVol/dewing/My_files/Research/Stochastic_SIR_Modelling/Plots_for_paper/Introduction_inferences.png", height=6, width=9, units="in", res=600)
#p1
#grid.arrange(p1,p2, nrow=2, widths=c(0.86,0.14), layout_matrix = rbind(c(1,1),c(2,NA)))
grid.arrange(p1, p2, nrow=2)
dev.off()
