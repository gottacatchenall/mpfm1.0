
{
  library(ggplot2)
  library(gmodels)
  library(dplyr)
  demography <- read.csv('demography.csv')
  #dynamics <- read.csv('dynamics.csv')
  metadata <- read.csv('metadata.csv')
  
  #mean_eff_mig_total <- aggregate(eff_mig_mean ~ run_id, data=dynamics, FUN = mean)
  #mean_eff_mig_total$eff_mig_avg <- mean_eff_mig_total$eff_mig_mean
  #mean_eff_mig_total$eff_mig_mean <- NULL
  #dynamics <- merge(dynamics, mean_eff_mig_total, by="run_id")
  #demography <- merge(demography, metadata, by="run_id")
  #dynamics <- merge(dynamics, metadata, by="run_id")
  
  #data <- merge(demography, dynamics, by=c("run_id","generation"))
  data <- merge(metadata, demography, by="run_id")
}
ggplot(subset(data), aes(generation, global_ld_mean, group=(run_id), color=(PATCH_K_MEAN))) + geom_point() + geom_line()

ggplot(subset(data), aes(generation, fst_mean, group=(run_id), color=(ENV_FACTOR_H_VALUE))) + geom_point() + geom_line()

ggplot(subset(data, fitness==1), aes(generation, n_loci_fixed, group=(run_id), color=(PATCH_DECAY))) + geom_point() + geom_line()


ggplot(data, aes(generation, eff_mig_mean, group=(run_id), color=(INCIDENCE_FUNCTION_DECAY))) + geom_point() + geom_line()

ggplot(subset(data), aes(generation, prop_of_k_mean, group=(run_id), color=(ENV_FACTOR_H_VALUE))) + geom_point() + geom_line()

aggregate_by_ld <- function(arg, xlabel="Generation", ylabel="Mean LD", col_name="col_name"){
    d <- subset(data)
    col = d[[arg]]
    agg_data <- do.call(data.frame, aggregate(global_ld_mean ~ generation*col, data = d, FUN = ci))
    agg_data$ymax <- agg_data$global_ld_mean.CI.upper
    agg_data$ymin <- agg_data$global_ld_mean.CI.lower 
    ggplot(agg_data, aes(generation, global_ld_mean.Estimate, color =  as.factor(col), group =  as.factor(col) )) + geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.5) + geom_line() + xlab(xlabel) + ylab(ylabel) + labs(color=col_name)
}

aggregate_by_local_ld <- function(arg, xlabel="Generation", ylabel="Mean LD", col_name="col_name"){
  d <- subset(data)
  col = d[[arg]]
  agg_data <- do.call(data.frame, aggregate(local_ld_mean ~ generation*col, data = d, FUN = ci))
  agg_data$ymax <- agg_data$local_ld_mean.CI.upper
  agg_data$ymin <- agg_data$local_ld_mean.CI.lower 
  ggplot(agg_data, aes(generation, local_ld_mean.Estimate, color =  as.factor(col), group =  as.factor(col) )) + geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.5) + geom_line() + xlab(xlabel) + ylab(ylabel) + labs(color=col_name)
}


aggregate_by_n_fixed <- function(arg, xlabel="Generation", ylabel="Number of Loci Fixed", col_name="col_name"){
  d <- subset(data)
  col = d[[arg]]
  agg_data <- do.call(data.frame, aggregate(n_loci_fixed ~ generation*col, data = d, FUN = ci))
    agg_data$ymax <- agg_data$n_loci_fixed.CI.upper
    agg_data$ymin <- agg_data$n_loci_fixed.CI.lower
    ggplot(agg_data, aes(generation, n_loci_fixed.Estimate, color =  as.factor(col), group = as.factor(col) )) + geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.5) + geom_line() + xlab(xlabel) + ylab(ylabel) + labs(color=col_name)
}



mult <- c("INCIDENCE_FUNCTION_DECAY", "NUM_PATCHES")
plot_df <-data %>% group_by(NUM_PATCHES) %>%
  do(
    plots = ggplot(data = .) + aes(x = generation, y = get_ld_mean, group=run_id) +
      geom_point() + ggtitle(.$NUM_PATCHES)
  )


aggregate_by_n_fixed("N_INDIVIDUALS", col_name="Total K")
aggregate_by_ld("N_INDIVIDUALS", col_name="Total K")


aggregate_by_n_fixed("NUM_PATCHES", col_name = "Number of Patches")
aggregate_by_ld("NUM_PATCHES", col_name = "Number of Patches")

aggregate_by_n_fixed("INCIDENCE_FUNCTION_DECAY", col_name = "Dist. Decay Strength")
aggregate_by_ld("INCIDENCE_FUNCTION_DECAY", col_name = "Dist. Decay Strength")

aggregate_by_local_ld("INCIDENCE_FUNCTION_DECAY", col_name = "Dist. Decay Strength")



# 2deme
aggregate_by_ld("PATCH_K_MEAN", col_name = "Dist. Decay Strength")
aggregate_by_n_fixed("PATCH_K_MEAN", col_name = "Dist. Decay Strength")


aggregate_by_n_fixed("MIGRATION_RATE", col_name = "Dist. Decay Strength")
aggregate_by_ld("MIGRATION_RATE", col_name = "Dist. Decay Strength")

aggregate_by_ld("EF_DIFFERENCE", col_name = "Dist. Decay Strength")


##
ggplot(subset(data, fitness==1), aes(generation, global_ld_mean, color=as.factor(NUM_PATCHES))) + geom_point(aes(group=run_id, size=0.5),alpha=0.05) + geom_line(aes(group=run_id),alpha=0.10)+ geom_smooth(aes(group=as.factor(interaction(N_INDIVIDUALS,INCIDENCE_FUNCTION_DECAY,NUM_PATCHES))), alpha=1) + facet_wrap(N_INDIVIDUALS ~ INCIDENCE_FUNCTION_DECAY) + labs(title = 'LD') + theme(legend.position="none")


