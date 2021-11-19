####PerMANOVA - Prey composition####

##Load packages##
library(vegan)
library(dplyr)
library(magrittr)
library(ggplot2)

# Setting directories ####
main.dir <- "./" # set the path to your working dir
setwd(main.dir)

# set input and output directories
in.dir <- file.path(paste0(main.dir,"input/"))
out.dir <- file.path(paste0(main.dir,"output/"))

# Reload saved image #########
session::restore.session(paste0(main.dir,"03-diet.R"))

############################################################
for ( j in seq_along(group.list)){
  # Select which group to use
  group.X <- group.list[[j]]
  # and which filter to apply
  filter.X <- filter.1
  # Names of group
  group.names
############################################################

# List prey names
as.factor(levels(merged.table$prey.id)) -> prey.name

##Convert read to presence
pres.data <- merged.table %>%
  # Apply filters
  filter(diet.status == "ok" & 
           amphibian.species == filter.X) %>%
  select(prey.id, sample, 
         #region, fragment, # groups being considered 
         group.names,
         blankfiltered.reads) %>%
  group_by(prey.id, sample) %>%
  #unique.data.frame() %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  pivot_wider(names_from = prey.id, values_from = blankfiltered.reads, values_fn = mean) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  as.data.frame()

## Number of reads per prey
rra.data <- merged.table %>%
  # Apply filters
  filter(diet.status == "ok" & 
           amphibian.species == filter.X) %>%
  select(prey.id, sample, 
         #fragment, region, # groups being considered
         group.names,
         blankfiltered.reads) %>%
  group_by(prey.id, sample) %>%
  #unique.data.frame() %>%
  #mutate_if(is.numeric, ~1 * (. > 0)) %>%
  pivot_wider(names_from = prey.id, values_from = blankfiltered.reads, values_fn = sum) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  as.data.frame()

## Proportion of species reads
rra.data <- data.frame(rra.data[!colnames(rra.data) %in% prey.name], 
                                   rra.data[colnames(rra.data) %in% prey.name]/
                                     rowSums(rra.data[colnames(rra.data) %in% prey.name]))

##Identify the most abundant prey among all samples
sort(colSums(pres.data[,colnames(pres.data) %in% prey.name]))

if (identical(group.X,group.2)==T){
## Include only samples with both fragments
  frag_comp_samples <- pres.data[
    pres.data$fragment == fragments[[1]],]$sample[
      pres.data[pres.data$fragment == fragments[[1]],]$sample %in%
        pres.data[pres.data$fragment == fragments[[2]],]$sample
      ]
# Remove samples from presence/absence dataframe
  pres.data <- pres.data[pres.data$sample %in% frag_comp_samples,]
# Remove samples from RRA dataframe
  rra.data <- rra.data[rra.data$sample %in% frag_comp_samples,]
}

##Create distance matrix##
Y.data <- pres.data[,colnames(pres.data) %in% prey.name] #select only columns with prey

#Presence/absence#
dist.Y<-vegdist(Y.data, method="jaccard")

#Abundance (RRA or weighted occurrences)#
dist.Y<-vegdist(Y.data, method="bray")

######
##Build PerMANOVA model##
######
adonis(dist.Y ~ pres.data[,as.character(group.names[j])], # Verify the group
       data=pres.data, method="binomial", permutations=1000)->perm1 
perm1

##Abundance (RRA or weighted occurrences)#
adonis(dist.Y ~ rra.data[,as.character(group.names[j])], # Verify the group
       data=pres.data, method="bray", permutations=1000)->perm2 
perm2

##Test dispersion among groups##
betadisper(dist.Y, pres.data[,as.character(group.names[j])], # Verify the group
           type="centroid")->bdisp.var1
bdisp.var1

anova(bdisp.var1)
jpeg(paste0(out.dir,"permova-dispersion.",group.names[j],".jpeg"))
#par(mfrow=c(1,2))
plot(bdisp.var1, 
     #col = c("red","blue","green"), 
     main = NULL, 
     label = T, label.cex = 0.5,
     sub = NULL, 
     hull = F, ellipse = T, 
     segments = T)
dev.off()

##Assess which prey are most different##
simper(Y.data, 
       pres.data[,as.character(group.names[j])],
       permutations=1000)->sim1
summary(sim1)

##Identify the names of the comparisons
names(sim1) -> sim.names
##Colour palette with black:
cbbPalette <- c("#000000", "#FFFFFF")
##Create lists
sim.df <- list()
sim.df.plot <- list()
for (i in seq_along(sim.names)){
  ##Isolate single comparison
  as.data.frame(summary(sim1)[i]) -> sim.df[[i]]
  ##Identify variables being compared 
  as.factor(unlist(strsplit(sim.names[i], split='_', fixed=TRUE))) -> ab.title
  ##Retain only significant prey p < 0.05
  sim.df[[i]][sim.df[[i]][,7] < 0.05,] -> sim.df[[i]]
  ##Prepare data into ggplot compatible table with variables "A" and "B"
  cbind.data.frame(rep("A", nrow(sim.df[[i]])),
                   rownames(sim.df[[i]]),
                   sim.df[[i]][,4]) -> a
  c("condition","prey","avg") -> colnames(a)
  cbind.data.frame(rep("B", nrow(sim.df[[i]])),
                   rownames(sim.df[[i]]),
                   -1*sim.df[[i]][5]) -> b
  c("condition","prey","avg") -> colnames(b)
  rbind.data.frame(a,b) -> ab
  ##Pass to ggplot
  ab %>%
    ggplot(aes(x = prey
               , y = avg
               , fill = condition)) +
    geom_col(stat="identity"
             , colour="black") +                    #Use bar-plot
    coord_flip() +                                  #Flip x and y axes
    expand_limits(y=c(-1,1)) +                      #Standardize y-axis
    ggtitle(paste(ab.title[2],"and",ab.title[1])) + #Set title
    scale_fill_manual(values=cbbPalette,            #Set legend colour 
                      name="Region",               #Set legend title
                      labels=ab.title) +             #Set legend names
    scale_x_discrete(name="") +                     #Set x label
    scale_y_continuous(name="") +     #Set y label
    theme(panel.grid.minor=element_blank(),         #Remove gridlines 
          panel.grid.major=element_blank(),
          plot.title = element_text(size = 10)) -> sim.df.plot[[i]]
}

if (identical(group.X,group.1)==T){
jpeg(paste0(out.dir,"simper_",group.names[j],".jpeg"))
# Group plots together
ggpubr::ggarrange(sim.df.plot[[1]],
                  sim.df.plot[[2]],
                  sim.df.plot[[3]],
                  common.legend = F,
                  nrow = 3)
dev.off()
}
}
