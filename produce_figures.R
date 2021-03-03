# Set number of chains you have run
nchains <- 50

# Set the working directory
setwd("./results.R")

process_chains <- function(nchains=50){
# This is a short function to combine the results across the chains into a single data frame for plotting
  filelist <- list.files()
  for(i in 1:nchains){
    load(filelist[i])
    samps <- results$samps
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

# Run the processing function above
results_df <- process_chains(herd=herd_no)

#########################################################
#########################################################
##                                                     ##
##      Code chunk to plot simulated data - Fig 3      ##
##                                                     ##
#########################################################
#########################################################

# Note that this assumes you have run the code above on the correct set of chains

library(LaplacesDemon)
data_summary <- function(x) {
# Short function to add data summary
  m <- median(x)
  ymin <- p.interval(x)[1]
  ymax <- p.interval(x)[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# True values for annotating the plots
sim_true_beta_int <- data.frame(beta_int = c(0.4,0.3,0.25,0.25,0.35,0.45,0.3,0.35,0.4), Herd=factor(1:9))
sim_true_beta_ext <- data.frame(beta_ext = c(0.0002,0.0002,0.0001,0.0002,0.0001,0.0002,0.0001,0.0002,0.0001), Herd=factor(1:9))
gamma_means <- data.frame(means = c(6.25, 9.12), stage = c("Latent", "Infectious"))
gamma_shapes <- data.frame(shapes = c(19.39, 22), stage = c("Latent", "Infectious"))

# Remove first 75% of observations
results_df <- subset(results_df, Iteration >= 1500)

# Create a data frame for plotting within-herd transmission (with shaded prior in background)
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

# Create a data frame for plotting external transmission
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

# Create a data frame for plotting means (with shaded prior in background)
means_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                        Iteration = rep(1:dim(results_df)[1],4),
                        Prior = factor(rep(c(F,T),each=dim(results_df)[1]*2)),
                        Value = c(results_df$exp_mean,
                                  results_df$inf_mean,
                                  rgamma(50200,shape=10,rate=10/6.25),
                                  rgamma(50200,shape=10,rate=10/9.12)))

# Create a data frame for plotting shape parameters (with shaded prior in background)
shapes_df <- data.frame(Stage = rep(rep(c("Latent", "Infectious"), each=dim(results_df)[1]),2),
                       Iteration = rep(1:dim(results_df)[1],4),
                       Prior = factor(rep(c(F,T),each=dim(results_df)[1]*2)),
                       Value = c(results_df$exp_shape,
                                 results_df$inf_shape,
                                 rgamma(50200,shape=5,rate=5/19.39),
                                 rgamma(50200,shape=5,rate=5/22)))

# Make the plots
p1 <- ggplot(subset(beta_int_df, Prior == "FALSE"), aes(x=Herd, y=Value)) +
  geom_violin(position="identity", alpha=0.5, col="#56B4E9", fill="#56B4E9") +
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

# Arrange into a grid and save
library(gridExtra)
png("./Sim_study_results.png", height=6, width=6, units="in", res=300)
grid.arrange(grobs=list(p1,p3,p2,p4), widths=c(2,1), layout.matrix=matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T))
dev.off()



########################################################
########################################################
##                                                    ##
##         Code chunk to plot ASF data - Fig 4        ##
##                                                    ##
########################################################
########################################################

# Note that this assumes you have run the code in the top section of this file on the correct set of chains

# Just like above, this removes the burn-in observations, makes data frames of results and then produces and arranges the plots
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


png("./ASF_results.png", height=6, width=6, units="in", res=300)
grid.arrange(grobs=list(p1,p4,p2,p5,p3), widths=c(2,1), layout.matrix=matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T))
dev.off()

#########################################################
#########################################################
##                                                     ##
##  Code chunk to plot fitted distributions - Fig 6    ##
##                                                     ##
#########################################################
#########################################################

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

png("./Fitted_distributions.png", height=3, width=4.5, units="in", res=300)
p1
dev.off()

#########################################################
#########################################################
##                                                     ##
##   Code chunk to plot introduction times - Fig 5     ##
##                                                     ##
#########################################################
#########################################################


process_chains2 <- function(nchains=50){
# Slightly modified processing function to give out introduction times and numbers of externally infected

  results_df <- data.frame(Herd=factor(rep(1:9,50)), Intro_times=NA, Intro_n=NA)
  setwd("./Results")

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

# Process results and plot
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

time.before.diagnosis <- data.frame(Herd=factor(1:9), Time=c(36,23,20,22,21,19,20,19,32))
p1 <- p1 +
  geom_text(data=time.before.diagnosis, aes(x=Herd, label=Time, y=-2))

p2 <- ggplot(results_df, aes(x=Herd, y=Intro_n)) +
  geom_violin(col="#56B4E9", fill="#56B4E9", alpha=0.5) + 
  stat_summary(fun.data=data_summary, geom="pointrange", col="black") +
  ylab(expression('Initial external infections, n'[0])) +
  theme_bw()

library(gridExtra)
png("./Introduction_inferences.png", height=6, width=9, units="in", res=600)
p1
dev.off()


