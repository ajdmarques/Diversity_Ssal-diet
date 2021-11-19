#  #########

# Packages
# if you want to learn more about dplyr and tidyverse, check here https://www.tidyverse.org/
library(dplyr)
library(tidyverse)
library(iNEXT)

# Setting directories ####
main.dir <- "./" # set the path to your working dir
setwd(main.dir)

# Reload saved image #########
session::restore.session(paste0(main.dir,"02-diet.R"))

#####
for ( j in seq_along(group.list)){
  # Select which group to use
  group.X <- group.list[[j]]
  # Identify columns in source table
  match(as.factor(names(group.list)),colnames(merged.table)) -> col.indx
  group.names <- as.factor(colnames(merged.table[,col.indx]))
  # and which filter to apply
  filter.X <- filter.1
  
## Generate Richness input data
# First we want to detect species presence (0,1) in each sample and by group
richness.raw <- merged.table %>%
  # Apply filters
  filter(diet.status == "ok" & 
             amphibian.species == filter.X) %>%
  select(prey.id, sample, 
         #fragment, region, 
         #as.factor(colnames(merged.table[,col.indx])),
         group.names,
         total.diet.reads) %>%
  group_by(prey.id, sample) %>%
  #unique.data.frame() %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  pivot_wider(names_from = prey.id, values_from = total.diet.reads, values_fn = mean) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  as.data.frame()

# List prey names
as.factor(levels(merged.table$prey.id)) -> prey.name

# Calculate Richness as the sum of prey species
cbind.data.frame(richness.raw, 
                 rowSums(richness.raw[,colnames(richness.raw)%in%prey.name])) ->
  richness.in
colnames(richness.in) <- c(as.character(colnames(richness.raw)),"Richness")
write.csv(richness.in, paste0(output.dir,group.names[j],"_rich_table.csv"))
  
##Build GLM model##
glm(Richness ~ richness.in[,as.character(group.names[j])], 
    data=richness.in, family=poisson) -> rich.glm
summary(rich.glm)

# Prey composition ####
sp_composition <- merged.table %>%
  # Apply filters
  filter(diet.status == "ok" &
           amphibian.species == filter.X) %>%
   dplyr::select(sample, prey.id, blankfiltered.reads) %>%
   group_by(sample, prey.id) %>%
   #summarise(reads = sum(blankfiltered.reads)) %>%
   pivot_wider(names_from = prey.id, values_from = blankfiltered.reads, values_fn = sum) %>%
   mutate_if(is.numeric, replace_na, 0)
 
## Proportion of species reads
sp_composition_reads <- data.frame(Sample = sp_composition$sample, 
                                   sp_composition[-1]/rowSums(sp_composition[-1]))
 
## Presence/absence ####
sp_composition_pa <- sp_composition %>% 
   mutate_if(is.numeric, ~1 * (. > 0)) 
 
## Proportion of species richness ####
sp_composition_prop <- data.frame(sample = sp_composition_pa$sample, 
                                   sp_composition_pa[-1]/rowSums(sp_composition_pa[-1]))
 
# Join all species composition metrics in a list
sp_composition <- list(composition_tot_reads = sp_composition,
                       composition_pa = sp_composition_pa,
                        composition_prop_reads = sp_composition_reads,
                        composition_prop_sp = sp_composition_prop)

# write prey composition to file
write.csv(sp_composition, paste0(output.dir,group.names[j],"sp_composition.csv"))


## Identify presence of each prey by taxonomic hierarchy
prey.data <- merged.table %>%
  # Apply filters
  filter(diet.status == "ok" & 
           amphibian.species == filter.X) %>%
  select(prey.id,
         Phylum, Class, Order, Family,
         blankfiltered.reads) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  group_by(Phylum, Class, Order, Family, prey.id) %>%
  summarize(count=sum(blankfiltered.reads)) %>%
  arrange(Phylum, Class, Order, Family) %>%
  as.data.frame()
## Save as table
write.csv(prey.data, paste0(output.dir,"prey.count.csv"))

}

# Save image
save.image(paste0(main.dir,"03-diet.R"), compress = "gzip")
