

ggplot(mpg, aes(cty, hwy, color = class)) + geom_density2d()

#Read data
{
library(vegan)
library(ggplot2)
patch_data <- read.csv('patches.csv')
att_migration <- read.csv('attempted_migration.csv')
success_migration <- read.csv('successful_migration.csv')
eff_migration <- read.csv('eff_migration.csv')
linkage <- read.csv('linkage.csv')
gbl_linkage <- read.csv('global_linkage.csv')
alleles <- read.csv('allele_freq.csv')
f_st <- read.csv('f_st.csv')
#colonization <- read.csv('colonization.csv')
#extinction <- read.csv('extinction.csv')
patch_data$patch_num <- factor(patch_data$patch_num)
linkage$D = abs(linkage$D)
gbl_linkage$D = abs(gbl_linkage$D)
linkage$patch_num = factor(linkage$patch_num)
alleles$patch_num <- factor(alleles$patch_num)
eff_migration$patch_num <- factor(eff_migration$patch_num)
success_migration$patch_from_num <- factor(success_migration$patch_from_num)
success_migration$patch_to_num <- factor(success_migration$patch_to_num)
att_migration$patch_from_num <- factor(att_migration$patch_from_num)
att_migration$patch_to_num <- factor(att_migration$patch_to_num)
}
run_name <- 'N=20 K=300, H=0.2, dynamic ef,'


# Mean Fitness
mean_w_by_gen <- aggregate(mean_w ~ generation*patch_num , data = patch_data, FUN=mean)

ggplot(mean_w_by_gen, aes(generation, mean_w, color = patch_num, group=patch_num))  + geom_point(size=0.5) + guides(col = guide_legend(nrow = 10)) + labs(y='Mean w', title=paste('Mean Fitness -- ',run_name))

# Effective Migration Total
eff_mig_by_gen <- aggregate(eff_migration ~ generation , data = eff_migration, FUN=mean)
ggplot(eff_mig_by_gen, aes(generation, eff_migration)) + geom_point() + geom_smooth() + guides(col = guide_legend(nrow = 10))+ labs(y='Effective Migration', title=paste('Effective Migration -- ',run_name))+ scale_colour_brewer(palette="Set1")

# Effective Migration by patch
eff_mig_by_gen <- aggregate(eff_migration ~ generation*patch_num , data = eff_migration, FUN=mean)
ggplot(eff_mig_by_gen, aes(generation, eff_migration,color = patch_num, group=patch_num)) + geom_line() + guides(col = guide_legend(nrow = 10)) + labs(y='Effective Migration', title=paste('Effective Migration -- ',run_name))

eff_mig_by_patch <- aggregate(eff_migration ~ generation*patch_num , data = eff_migration, FUN=mean)

ggplot(eff_mig_by_patch, aes(eff_migration, color = patch_num)) + geom_density()

## F_st
f_st_by_gen <- aggregate(F_st ~ (generation) , data = subset(f_st), FUN=mean)
ggplot(f_st_by_gen, aes(generation, F_st)) + geom_point() + geom_smooth(span=0.1) + guides(col = guide_legend(nrow = 10)) + labs(y='F_ST', title=paste('F_st -- ',run_name))


f_st_by_locus <- aggregate(F_st ~ locus*generation, data = subset(f_st), FUN=mean)
ggplot(f_st_by_locus, aes(F_st, color=locus, group=locus)) + geom_density() + xlim(-.05,1)


# Avg LD
mean_linkage_by_gen <- aggregate(D ~ generation*type , data = linkage, FUN=mean)
ggplot(mean_linkage_by_gen, aes(generation, D, color = type, group=type)) + geom_point() + geom_smooth(span=0.4) + guides(col = guide_legend(nrow = 10))+ labs(y='D', title=paste('Within-deme LD -- ',run_name))+ scale_colour_brewer(palette="Set1")

# global LD
mean_gbl_linkage_by_gen <- aggregate(D ~ generation*type , data = gbl_linkage, FUN=mean)
ggplot(mean_gbl_linkage_by_gen, aes(generation, D, color = type, group=type)) + geom_point()+ geom_smooth(span=0.4) + guides(col = guide_legend(nrow = 10))+ labs(y='D', title=paste('Global LD -- ',run_name)) + scale_colour_brewer(palette="Set1")
#+ geom_vline(aes(xintercept=generation), extinction)

## F_st
f_st_by_gen <- aggregate(F_st ~ generation , data = subset(f_st), FUN=mean)
ggplot(f_st_by_gen, aes(generation, F_st)) + geom_point(color='dodgerblue') + guides(col = guide_legend(nrow = 10)) + labs(y='F_ST', title=paste('F_st -- ',run_name))  

## Allele freq
## 
allele_freq <- aggregate(frequency ~ generation*patch_num , data = alleles, FUN=mean)
# look at linkage within groups of loci for the same ef
ggplot(allele_freq, aes(generation, frequency, color = patch_num, group=patch_num)) + geom_point() + geom_line() + guides(col = guide_legend(nrow = 10)) + labs(y='Mean Allele Frequency', title=paste('Mean Allele Frequency -- ',run_name)) 





## One locus one patch
{
one_patch <- subset(alleles, locus == 27 & patch_num=="12")
one_patch$allele_val <- factor(one_patch$allele_val)
ggplot(one_patch, aes(generation, frequency, color = allele_val, group=allele_val)) + geom_point() + geom_line()+ guides(col = guide_legend(nrow = 10)) + labs(y='Allele Frequency', title=paste('Fitness Allele -- locus=39, patch=3 --',run_name))
}


## allele dissimilarity

# subset by gen prior

avg_num_shared <- function(gen){
  s <- 0
  ct <- 0
  for (p1 in levels(alleles$patch_num)){
    p1_al <- subset(alleles, patch_num == p1 & generation == gen)
    p1_als = unique(p1_al$allele_val)
    for (p2 in levels(alleles$patch_num)){
      p2_al <- subset(alleles, patch_num == p2 & generation == gen)
      p2_als = unique(p2_al$allele_val)
      s <- s + length(intersect(p2_als, p1_als))/min(length(p2_als), length(p1_als))
      ct <- ct + 1
    }
  }
  
  if (ct > 0){
    return(s/ct)
  }
  else{
    return(0)
  }
}


gens <- c()
shared <- c()
i <- 1
for (gen in levels(factor(alleles$generation))){
  print(gen)
  if (as.numeric(gen) %% 50 == 0){
    gens[i] = as.numeric(gen)
    shared[i] <- avg_num_shared(as.numeric(gen))
    i <- i + 1
  }
}
shared_df <- data.frame(gens,shared)

ggplot(shared_df, aes(gens, shared)) + geom_point() + guides(col = guide_legend(nrow = 10)) + labs(y='Prop of Alleles Shared', title=paste('Proportion of Alleles Shared --',run_name))

### Prop of K
### 
### % 18 , 15, 0
mean_k_by_gen <- aggregate(prop_of_k ~ generation*patch_num, data=patch_data,FUN=mean)
ggplot(mean_k_by_gen, aes(generation, prop_of_k, color=patch_num, group=patch_num)) + geom_point() + geom_line()+ guides(col = guide_legend(nrow = 10))


summary(lm(mean_k_by_gen$prop_of_k ~ mean_k_by_gen$generation))

## Migration
## 
{
migration_by_gen <- aggregate(prop_of_new_patch_migrants ~ generation*patch_from_num*patch_to_num , data = att_migration, FUN=mean)
pairs = c()
j <- 1
for (i in c(1:nrow(migration_by_gen))){
 
  from <- as.numeric(migration_by_gen$patch_from_num[i])
  to <- as.numeric(migration_by_gen$patch_to_num[i])
  if ( from != to){
    this_pair_coord <- c(from, to)
    new = TRUE
    for (pair in pairs){
      if (pair == this_pair_coord){
          new = FALSE
      }
    }
    
    if (new == TRUE){
      pairs[[j]] <- this_pair_coord
      j <- j+1
    }
  }
}

migration_pair <- c()
generation <- c()
prop_new_migrants <- c()

i <- 1
for (pair in pairs){
  from <- pair[1]
  to <- pair[2]
  this_pair <- subset(migration_by_gen, as.numeric(migration_by_gen$patch_from_num) == from & as.numeric(migration_by_gen$patch_to_num) == to)
  
  this_pair_by_gen <- aggregate(prop_of_new_patch_migrants ~ generation, data = this_pair, FUN=mean)
    
  for (j in c(1:nrow(this_pair_by_gen))){
    migration_pair[i] <- paste(from, "->",to)
    generation[i] <- this_pair_by_gen$generation[j]
    prop_new_migrants[i] <- this_pair_by_gen$prop_of_new_patch_migrants[j]
    i <- i+1
  }
}

migr_df <- data.frame(factor(migration_pair), generation, prop_new_migrants)

ggplot(migr_df, aes(generation, prop_new_migrants, color=migration_pair, group=migration_pair)) + geom_point() + geom_line()+ guides(col = guide_legend(nrow = 10))+ labs(y='attempted migration as prop of new patch', title=paste('Attempted Migration -- ',run_name))
}



### LD By Patch
+ geom_density2d()

