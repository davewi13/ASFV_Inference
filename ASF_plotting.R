# Set herd number
nchains <- 5

# This function assumes results files for each of n chains herds are in their own directories as specified below
process_chains <- function(herd){
  results_df <- NULL
  if(herd == 1){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h1")
  } else if (herd == 2){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h2")
  } else if (herd == 3){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h3")
  } else if (herd == 4){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h4")
  } else if (herd == 5){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h5")
  } else if (herd == 6){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h6")
  } else if (herd == 7){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h7")
  } else if (herd == 8){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h8")
  } else if (herd == 9){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_R0_h9")
  } else if (herd == "other"){}
  filelist <- list.files()
  for(i in 1:nchains){
    data <- read.table(filelist[i], header=T, sep="\t")
    data$Chain <- i
    results_df <- rbind(results_df, data)
  }
  results_df$Herd <- herd
  return(results_df)
}

# Combine data from all herds together into one results file
results_df <- NULL
for(i in 1:9){
  temp <- process_chains(herd=i)
  results_df <- rbind(results_df,temp)
}


#################################################
#################################################
##                                             ##
##      Code chunk to plot Figures 3, 4, 6     ##
##                                             ##
#################################################
#################################################

# Function to get data summary
library(LaplacesDemon)
data_summary <- function(x) {
  m <- median(x)
  ymin <- p.interval(x)[1]
  ymax <- p.interval(x)[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Store true values in the case of simulated data
sim_true_R0_int <- data.frame(R0 = c(8,5,11,4,3,14,5,7,6), Herd=factor(1:9))
sim_true_beta_ext <- data.frame(beta_ext = c(0.01,0.5,3,0.1,4,0.05,6,0.01,1), Herd=factor(1:9))
sim_true_T0 <- data.frame(T0 = rep(0,9), Herd=factor(1:9))
gamma_means <- data.frame(means = c(rep(7,9),rep(9,9)), stage = c(rep("Latent",9), rep("Infectious",9)), herd = c(1:9, 1:9))
gamma_shapes <- data.frame(shapes = c(rep(19.39,9),rep(22,9)), stage = c(rep("Latent",9), rep("Infectious",9)), herd = c(1:9, 1:9))

# Remove burn-in period
results_df <- subset(results_df, State >= max(results_df$State)*0.1)

# Plotting indicator to show these are results and not simulations from the prior
results_df$Prior.log <- F

# Herd sizes and times of initial removals (commented out the simulated data case)
herd_sizes <- c(1614,1949,1753,1833,1320,600,600,600,2145)
#herd_sizes <- c(500,500,500,1000,1000,1000,2000,2000,2000)
init_times <- c(116,117,120,116,118,120,116,121,117)
#init_times <- rep(100,9)

# Get everything scaled and located in the right place
for(i in 1:9){
  results_df[results_df$Chain == i,]$beta1_0 <- exp(results_df[results_df$Chain == i,]$beta1_0)*herd_sizes[i]
  results_df[results_df$Chain == i,]$Tinit.Herd0 <- results_df[results_df$Chain == i,]$Tinit.Herd0-init_times[i]
}

# Simulate loads of values from the priors for the shaded areas in the plots
prior_df <- data.frame(State=1:90000,
                       mu_L=rgamma(90000, 10, 10/6.25),
                       sh_L=rgamma(90000, 5, 5/19.39),
                       mu_I=rgamma(90000, 10, 10/9.12),
                       sh_I=rgamma(90000, 5, 5/22.00),
                       beta1_0=runif(90000, 0, 1),
                       R0_0=runif(90000, 0, 30),
                       Ninf.Herd0=NA,
                       Tinit.Herd0=NA,
                       Init.infected0=NA,
                       Ext.infected0=NA,
                       Li.Herd0=NA,
                       Pr.index.Herd0=NA,
                       Prior=NA,
                       Chain=NA,
                       Herd=rep(1:9, each=10000),
                       Prior.log=T)

# Combine together
results_df <- rbind(results_df, prior_df)

results_df$Herd <- factor(results_df$Herd)

# Switch to T/F for simulated/observed data
sim_data <- F

# Create the plots
library(ggplot2)
p1 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=R0_0)) +
  geom_violin(position="identity", alpha=0.5, col="#56B4E9", fill="#56B4E9", scale="width") +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab("R0") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = sim_true_R0_int, mapping = aes(x=Herd, y=R0), col="black", shape=3, size=3*sim_data)
p1 <- p1 + geom_violin(data=subset(results_df, (Prior.log == T) & (R0_0 <= 30)), aes(x=Herd, y=R0_0), alpha=0.2, col="#808080", fill="#808080", scale="width")

p2 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=beta1_0)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5, scale="width") + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Herd size times '*beta['1i'])) +
  theme_bw() +
  ylim(0,20) +
  geom_point(data = sim_true_beta_ext, mapping = aes(x=Herd, y=beta_ext), col="black", shape=3, size=3*sim_data)

p3 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=mu_L)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", position="identity", alpha=0.5, scale="width") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Mean latent period '*mu['L'])) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = subset(gamma_means, stage == "Latent"), mapping = aes(x=factor(herd), y=means), col="black", shape=3, size=3*sim_data)
p3 <- p3 + geom_violin(data=subset(results_df, (Prior.log == T) & (mu_L <= 15)), aes(x=Herd, y=mu_L), alpha=0.2, col="#808080", fill="#808080", scale="width")

p4 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=mu_I)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", position="identity", alpha=0.5, scale="width") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Mean infectious period '*mu['I'])) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = subset(gamma_means, stage == "Infectious"), mapping = aes(x=factor(herd), y=means), col="black", shape=3, size=3*sim_data)
p4 <- p4 + geom_violin(data=subset(results_df, (Prior.log == T) & (mu_I <= 15)), aes(x=Herd, y=mu_I), alpha=0.2, col="#808080", fill="#808080", scale="width")

p5 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=sh_L)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", position="identity", alpha=0.5, scale="width") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Latent shape '*k['L'])) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = subset(gamma_shapes, stage == "Latent"), mapping = aes(x=factor(herd), y=shapes), col="black", shape=3, size=3*sim_data)
p5 <- p5 + geom_violin(data=subset(results_df, (Prior.log == T) & (sh_L >= 5)), aes(x=Herd, y=sh_L), alpha=0.2, col="#808080", fill="#808080", scale="width")

p6 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=sh_I)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", position="identity", alpha=0.5, scale="width") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Infectious shape '*k['I'])) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = subset(gamma_shapes, stage == "Infectious"), mapping = aes(x=factor(herd), y=shapes), col="black", shape=3, size=3*sim_data)
p6 <- p6 + geom_violin(data=subset(results_df, (Prior.log == T) & (sh_I >= 5)), aes(x=Herd, y=sh_I), alpha=0.2, col="#808080", fill="#808080", scale="width")

p7 <- ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=Tinit.Herd0)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", position="identity", alpha=0.5, scale="width") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Time of initial infection') + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = sim_true_T0, mapping = aes(x=Herd, y=T0), col="black", shape=3, size=3*sim_data)

# Store the plots - some editing of grid.arrange will be required depending on exactly which plot is needed.
library(gridExtra)
png("C:/Users/dewing/Exact_Bayesian_Inference/Submission to Interface/Plots_for_paper_resubmission/ASF_R0.png", height=7.5, width=6, units="in", res=600)
grid.arrange(grobs=list(p1,p2,p3,p4,p5,p6,p7), nrow=4, ncol=2, layout.matrix=rbind(c(1,2),c(3,4),c(5,6),c(7,7)))
dev.off()

#################################################
#################################################
##                                             ##
##      Code chunk to plot Figure 7            ##
##                                             ##
#################################################
#################################################

png("C:/Users/dewing/Exact_Bayesian_Inference/Submission to Interface/Plots_for_paper_resubmission/Introduction_inferences.png", height=4.5, width=4.5, units="in", res=600)
ggplot(subset(results_df, Prior.log == F), aes(x=Herd, y=Ext.infected0)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", position="identity", alpha=0.5, scale="width") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  scale_color_manual(values=c("#56B4E9", "#D3D3D3")) +
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab('Initial external infections') + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

#################################################
#################################################
##                                             ##
##      Code chunk to plot Figure 5            ##
##                                             ##
#################################################
#################################################

# Assumes the PPC files are stored in separate directories as above
process_ppc <- function(herd,nchains){
  results_df <- NULL
  if(herd == 1){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h1_ppc")
  } else if (herd == 2){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h2_ppc")
  } else if (herd == 3){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h3_ppc")
  } else if (herd == 4){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h4_ppc")
  } else if (herd == 5){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h5_ppc")
  } else if (herd == 6){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h6_ppc")
  } else if (herd == 7){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h7_ppc")
  } else if (herd == 8){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h8_ppc")
  } else if (herd == 9){
    setwd("C:/Users/dewing/Exact_Bayesian_Inference/Final_model/Final_mod_h9_ppc")
  } else if (herd == "other"){}
  filelist <- list.files()
  for(i in 1:nchains){
    data <- read.table(filelist[i], header=T, sep="\t")
    data$Chain <- i
    data <- data[order(data$State, data$Removal),]
    data <- subset(data, Observed == 1)
    data$State <- factor(data$State)
    data$CumRem <- 1
    for(i in 2:length(data$State)){
      if(data$State[i] == data$State[i-1]){
        data$CumRem[i] <- data$CumRem[i-1]+1
      }
    }
    results_df <- rbind(results_df, data)
  }
  index <- results_df$Removal[1]
  for(i in 2:length(results_df$Removal)){
    if((results_df$Herd[i] != results_df$Herd[i-1]) || (results_df$State[i] != results_df$State[i-1])){
      index <- results_df$Removal[i]
    }
    results_df$Removal[i] <- results_df$Removal[i]-index
  }
  results_df$Removal[1] <- 0
  results_df$State <- as.numeric(as.character(results_df$State))
  results_df <- subset(results_df, State > 200000)
  results_df$Herd <- herd
  results_df$Herd <- factor(results_df$Herd)
  return(results_df)
}

# Combine all the results
results_df <- NULL
for(i in 1:9){
  temp <- process_ppc(i, 5)
  results_df <- rbind(results_df,temp)
}

results_df$Chain <- factor(results_df$Chain)
results_df$Chain.State <- interaction(results_df$State,results_df$Chain)

# Input the observations
Times <- c(0:23)
Herd <- 1:9
asf_dat <- expand.grid(Times, Herd)
asf_dat$CumRem <- c(c(1,1,1,2,2,2,3,4,4,4,4,5,6,6,10,16,21,27,28,36,36,37,48,61),
                    c(1,2,3,5,8,9,18,24,28,32,43,46,51,58,60,63,rep(NA,8)),
                    c(1,2,3,5,7,12,16,25,26,29,32,37,43,rep(NA,11)),
                    c(1,3,8,15,20,27,34,42,46,51,55,65,67,81,99,102,rep(NA,8)),
                    c(4,10,14,18,24,32,38,39,42,47,56,64,83,rep(NA,11)),
                    c(1,2,4,6,7,7,7,8,9,17,rep(NA,14)),
                    c(2,2,3,3,5,6,10,12,20,32,40,46,51,51,rep(NA,10)),
                    c(2,2,5,8,12,15,15,16,16,rep(NA,15)),
                    c(7,14,17,19,19,24,25,25,31,31,32,36,36,46,82,106,123,141,148,188,230,261,273,NA))
asf_dat$CumRem <- as.numeric(asf_dat$CumRem)
names(asf_dat)[c(1,2)] <- c("Times", "Herd")

p1 <- ggplot(asf_dat, aes(x=Times, y=CumRem)) +
  facet_wrap(~Herd, scales = "free") +
  geom_line(lwd=3, col="red") +
  geom_line(data=results_df, aes(x=Removal, y=CumRem, grouping=Chain.State),alpha=0.1) +
  scale_color_discrete("viridis") +
  xlim(0,24) + 
  theme(legend.position = "none") +
  xlab("Day") +
  ylab("Cumulative Removals") +
  ggtitle("Including external transmission") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# The code above will produce the PPC plots for a given model run - to exactly reproduce Figure 5 run this twice for the cases with and without beta1 and store the results as p1 and p2 then plot using grid.arrange.
png("C:/Users/dewing/Exact_Bayesian_Inference/Submission to Interface/Plots_for_paper_resubmission/ASF_ppc.png", height=9, width=6, units="in", res=600)
grid.arrange(p1,nrow=1,ncol=1)
dev.off()

#################################################
#################################################
##                                             ##
##      Code chunk to plot Figure 1            ##
##                                             ##
#################################################
#################################################

asf_dat$Herd <- factor(asf_dat$Herd)
head(asf_dat)
asf_dat_sum <- data.frame(Herd = factor(1:9), Obs = c(22,10,7,8,6,5,9,4,16))

png("C:/Users/dewing/Exact_Bayesian_Inference/Submission to Interface/Plots_for_paper_resubmission/Data_plot.png", height=4.5, width=4.5, units="in", res=600)
ggplot(asf_dat, aes(x=Times, y=CumRem)) +
  geom_line(lwd=1.5) +
  facet_wrap(~Herd, scales="free") +
  geom_vline(data=asf_dat_sum, aes(xintercept=Obs),col="red") +
  ylab("Cumulative Removals") +
  xlab("Day") + 
  theme_bw()
dev.off()


