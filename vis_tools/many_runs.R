
{
  library(ggplot2)
  demography <- read.csv('demography.csv')
  dynamics <- read.csv('dynamics.csv')
  metadata <- read.csv('metadata.csv')
  
  mean_eff_mig_total <- aggregate(eff_mig_mean ~ run_id, data=dynamics, FUN = mean)
  mean_eff_mig_total$eff_mig_avg <- mean_eff_mig_total$eff_mig_mean
  mean_eff_mig_total$eff_mig_mean <- NULL
  dynamics <- merge(dynamics, mean_eff_mig_total, by="run_id")
  #demography <- merge(demography, metadata, by="run_id")
  #dynamics <- merge(dynamics, metadata, by="run_id")
  
  data <- merge(demography, dynamics, by=c("run_id","generation"))
  data <- merge(data, metadata, by="run_id")
}

ggplot(subset(data, fitness == 0), aes(generation, global_ld_mean, group=(run_id), color=(MEAN_LOCUS_WEIGHT))) + geom_point() + geom_smooth(method="lm")

ggplot(subset(data), aes(generation, fst_mean, group=(run_id), color=(ENV_FACTOR_H_VALUE))) + geom_point() + geom_line()

ld_model <- (lm(global_ld_mean ~ generation*MEAN_LOCUS_WEIGHT, data = data))



ggplot(subset(data), aes(generation, prop_of_k_mean, group=(run_id), color=(ENV_FACTOR_H_VALUE))) + geom_point() + geom_line()



#ggplot(data, aes(generation, eff_mig_mean, color=as.factor(run_id), group=as.factor(run_id))) + geom_smooth() + geom_point()


