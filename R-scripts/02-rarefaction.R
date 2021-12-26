# Rarefaction and extrapolation of species diversity (Hill numbers) #########
# From this point onwards, needs to be reviewed
# Here we will use package 'iNEXT' and we will need to convert our data to 'iNEXT' format (check package vignette for more information)

# Packages
# if you want to learn more about dplyr and tidyverse, check here https://www.tidyverse.org/
library(dplyr)
library(tidyverse)
library(iNEXT)

# Setting directories ####
main.dir <- "./" # set the path to your working dir
setwd(main.dir)

# set input and output directories
input.dir <- file.path(paste0(main.dir,"input/"))
output.dir <- file.path(paste0(main.dir,"output/"))

## Identify all fragments used.
fragments <- c("LXT","fwh1")

## Select the merged file for each fragment
path <- list() 
for (i in seq_along(fragments)){
  path[[i]] <-  list.files(paste0(main.dir,"output/",fragments[i]),
                           pattern="merged.table",
                           full.names = TRUE)
}
# Read files
files <- list()
for (i in seq_along(fragments)){
  files[[i]] <-  read.table(path[[i]],
                            header = TRUE,
                            fill = TRUE,
                            na.strings = c(""," ", "na"),
                            sep = "\t")
}
# Designate the otu fragment
for (i in seq_along(fragments)){
  files[[i]]$fragment <- fragments[i]
}
# Merge the two fragments into a common table
files[[1]] -> merged.table
for (i in seq_along(fragments)[-1]){
  bind_rows(merged.table, files[[i]]) -> merged.table
}
dim(merged.table)

#####
## Identify which groups you will be comparing
#####
# Consider making this a list of treatments in a  parameter file.
# Example. Select only one species
levels(merged.table$amphibian.species)[2] -> filter.1

# Example: Select several habitats
# Decide whether to generate new variable
region.setting = TRUE
if (region.setting == TRUE){
  # Identify habitat types
  levels(merged.table$habitat) -> var.habitat
  # Set region names
  as.factor(c("Morrazo","Morrazo","Ons","Porto")) -> var.region #Must be same length as var.habitat
  # Name new variable based on another conditions
  merged.table$habitat -> merged.table$region
  library(plyr)
  mapvalues(merged.table$region,
          from = as.character(var.habitat), 
          to = as.character(var.region)) -> merged.table$region
  as.factor(na.exclude(unique(merged.table$region))) -> var.region
  # Set region and variable groups
  levels(as.factor(merged.table$region)) -> group.1
}
# group by fragment
levels(as.factor(merged.table$fragment)) -> group.2

# List groups
group.list <- list(group.1,group.2) 
names(group.list) <- c('region','fragment')

## Create a variable for the final prey ID
merged.table %>%
  na_if('NA') %>%
  ## Identify Family when Genus is unknown
  mutate(prey.id = coalesce(Genus,Family,Order,Class,Phylum))  -> merged.table

# Write table
write.table(merged.table, file = paste0(output.dir, 
                                        "merged.table.txt"),
            sep = '\t',
            row.names = FALSE)
############################################################
## iNEXT used to estimate rarifaction curves for each group.
############################################################
for ( j in seq_along(group.list)){
# Select which group to use
group.X <- group.list[[j]]
# and which filter to apply
filter.X <- filter.1

# First we want to detect species presence (0,1) in each sample and by group
incidence_raw <- list()
# Index of variable
match(as.factor(names(group.list)[j]),colnames(merged.table)) -> col.indx
for (i in seq_along(group.X)){
  incidence_raw[[i]] <- merged.table %>%
    # Apply filters
    filter(diet.status == "ok" & 
             amphibian.species == filter.X &
             merged.table[,col.indx] == group.X[i]
             ) %>%
    # Select the level of prey
    select(prey.id, sample, total.diet.reads) %>%
    group_by(prey.id, sample) %>%
    unique.data.frame() %>%
    #summarise(reads = sum(total.diet.reads)) %>%
    mutate_if(is.numeric, ~1 * (. > 0)) %>%
    pivot_wider(names_from = sample, values_from = total.diet.reads, values_fn = mean) %>%
    mutate_if(is.numeric, replace_na, 0) %>%
    as.data.frame()
}

## iNetx requires that species names are st as row names
for (i in seq_along(group.X)){
  row.names(incidence_raw[[i]]) <- incidence_raw[[i]]$prey.id
  incidence_raw[[i]]$prey.id <- NULL
}

## create iNEXT incidence_raw
names(incidence_raw) <- group.X

## create iNEXT incidence_frequency object
incidence_freq <- lapply(incidence_raw, as.incfreq)
 
# we recommend saving your iNEXT objects into your output 
save(incidence_raw, incidence_freq, file = paste0(output.dir, names(group.list)[j],'.iNEXT.data.R'))

## specific filter for fragments
if (identical(group.X,group.2)==T){
  ## Include only samples with both fragments
  frag_comp_sample <- colnames(incidence_raw[[1]])[colnames(incidence_raw[[1]]) %in% 
                                                     colnames(incidence_raw[[2]])]
  incidence_raw[[1]] <- incidence_raw[[1]][,frag_comp_sample]
  incidence_raw[[2]] <- incidence_raw[[2]][,frag_comp_sample]
}

## run iNEXT
out1 <- list()
out1[[j]] <- iNEXT(incidence_raw, q = 0, datatype = 'incidence_raw')
write.csv(out1[[j]]$AsyEst, paste0(output.dir,names(group.list)[j],"_richness.csv"))

out2 <- list()
out2[[j]] <- iNEXT(incidence_freq, q = 0, datatype = 'incidence_freq')

##Run iNEXT
inext.plot.unit.div <- list()
inext.plot.unit.cov <- list()
inext.plot.cov.div <- list()
ggiNEXT(out1[[j]], type = 1) -> inext.plot.unit.div[[j]]
ggiNEXT(out1[[j]], type = 2) -> inext.plot.unit.cov[[j]]
ggiNEXT(out1[[j]], type = 3) -> inext.plot.cov.div[[j]]
jpeg(paste0(output.dir,names(group.list)[j],".inext.plot.jpeg"))
ggpubr::ggarrange(inext.plot.unit.div[[j]],
                  inext.plot.unit.cov[[j]],
                  inext.plot.cov.div[[j]],
                  common.legend = T,
                  nrow = 3)
dev.off()
jpeg(paste0(output.dir,names(group.list)[j],".coverage.diversity.plot.jpeg"))
inext.plot.cov.div[[j]]
dev.off()

}
# Save image
session::save.session(paste0(main.dir,"02-diet.R"), compress = "gzip")
 
