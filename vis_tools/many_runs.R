
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

ggplot(subset(data, fitness==1), aes(generation, global_ld_mean, group=(run_id), color=(NUM_PATCHES))) + geom_point() + geom_line()

ggplot(subset(data), aes(generation, fst_mean, group=(run_id), color=(ENV_FACTOR_H_VALUE))) + geom_point() + geom_line()

ggplot(subset(data, fitness==1), aes(generation, n_loci_fixed, group=(run_id), color=(PATCH_DECAY))) + geom_point() + geom_line()


ggplot(data, aes(generation, eff_mig_mean, group=(run_id), color=(INCIDENCE_FUNCTION_DECAY))) + geom_point() + geom_line()

ggplot(subset(data), aes(generation, prop_of_k_mean, group=(run_id), color=(ENV_FACTOR_H_VALUE))) + geom_point() + geom_line()

aggregate_by <- function(arg){
  d <- subset(data, fitness == 1)
  col = d[[arg]]
  agg_data <- do.call(data.frame, aggregate(global_ld_mean ~ generation*col, data = d, FUN = function(x) c(mn = mean(x), sd = sd(x))))
  agg_data$ymax <- agg_data$global_ld_mean.mn + agg_data$global_ld_mean.sd
  agg_data$ymin <- agg_data$global_ld_mean.mn - agg_data$global_ld_mean.sd
  ggplot(agg_data, aes(generation, global_ld_mean.mn, color = col, group = col )) + geom_ribbon(aes(ymin=ymin, ymax=ymax + global_ld_mean.sd)) + geom_line()
}

aggregate_by("PATCH_DECAY")

