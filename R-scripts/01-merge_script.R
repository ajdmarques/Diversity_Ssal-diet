#####################################
### Metabarcoding training School ###

# Packages
# if you want to learn more about dplyr and tidyverse, check here https://www.tidyverse.org/
library(dplyr)
library(tidyverse)
library(iNEXT)

# Table Preparation #########

# Setting directories ####
main.dir <- "./" # set the path to your working dir

setwd(main.dir)

# name of fragment
fragment = c("LXT","fwh1")

# Repeat for each fragment
for (i in seq_along(fragment)){

# set input and output directories
input.dir <- file.path(paste0(main.dir,"input/",fragment[i]))

output.dir <- file.path(paste0(main.dir,"output/",fragment[i]))

# create directory to save results
dir.create(output.dir, showWarnings = FALSE,
           recursive = TRUE)

## Load txt files to a list ####

# path to folders containing txt files
path <- list.files(path = input.dir, 
                   pattern="\\.csv$", 
                   full.names = TRUE)

files <- list()

# Load OUT table
files[[1]] <- read.table(path[grep("read-count",path)],
               header = TRUE,
               fill = TRUE,
               na.strings = c(""," ", "na"),
               sep = ",")

# Load OTU id
files[[2]] <- read.table(path[grep("boldigger",path)],
                         header = TRUE,
                         fill = TRUE,
                         na.strings = c("", " ", "na"),
                         sep = ",")
colnames(files[[2]])[1] <- "otu" 
files[[2]][,1] <- as.factor(stringi::stri_replace_all(files[[2]][,1], 
                                                      replacement="",fixed = ">" ))

# Load Sample info
files[[3]] <- read.csv(path[grep("sample-info",path)],
                         header = TRUE,
                         fill = TRUE,
                         na.strings = c("", " ", "na"),
                         sep = ",",
                         quote = '')
colnames(files[[3]])[1] <- "sample" 

length(files) # your list should contain three objects (files)

# set files list names
file.names <- c('otu.table.txt','otu.id.txt','sample.info.txt')

names(files) <- file.names

# NOTE: Be careful to check your file names as they might need to be updated accordingly 'files$NAME'
### Step 1 ####
# Pivot OTU table So each line is the reads for each sample
merged.table <- files$otu.table.txt %>% 
  pivot_longer(c(-(1:5)), 
               names_to = 'sample', 
               values_to = 'reads')%>% 
  # Remove lines where the OTU is not found in a sample
  filter(reads != 0)
 
### Step 2 ####
# Join the sample info by header sample
merged.table <- left_join(merged.table, 
                        files$sample.info.txt, 
                        by = 'sample')

### Step 3 ####
# Join the OTU id by OTU name
merged.table <- left_join(merged.table, 
                        files$otu.id.txt, 
                        by = 'otu')

### Step 4 ####
# In original script, foudn to not be necessary
# Make new columns to identify extraction and PCR blanks
#merged.table$type == 'extblank' -> merged.table$extblank
#merged.table$type == 'pcrblank' -> merged.table$pcrblank
# Create a new column for Extraction and PCR blanks.
#merged.table <- merged.table %>% 
#  unite(extblankotu, c('extblank', 'otu'), remove = FALSE, sep = '') %>%
#  unite(pcrblankotu, c('pcrblank', 'otu'), remove = FALSE, sep = '') %>%
#  relocate(extblankotu, .after = extblank) %>%
#  relocate(pcrblankotu, .after = pcrblank)

### Step 5 & 6 ####
# Identify blanks/negatives as a new column
blanks <- merged.table %>% filter(type == 'extblank' | type == 'pcrblank') %>%
  select(c('sample', 'otu', 'reads', 'type')) %>%
  unite(blank, c('sample', 'otu'), remove = FALSE, sep = '') %>%
  relocate(blank, .after = reads)

head(blanks)

### Step 7 ####
# Identify the number of reads for each extraction blank 
tmp <- blanks %>% 
  filter(type == 'extblank')%>%
  select(otu, reads) %>%
  rename(extblankreads = reads) %>%
  group_by(otu) %>%
  summarise(sum(extblankreads)) 
colnames(tmp) <- c('otu','extblankreads')
 
merged.table <- merge(merged.table, tmp , 
            by = 'otu', all = T)

# Identify the number of reads for each PCR blank
tmp <- blanks %>% 
  filter(type == 'pcrblank') %>%
  select(otu, reads) %>%
  rename(pcrblankreads = reads) %>%
  group_by(otu) %>%
  summarise(sum(pcrblankreads)) 
colnames(tmp) <- c('otu','pcrblankreads')

merged.table <- merge(merged.table, tmp , 
                      by = 'otu', all = T)

rm(tmp)

# replace NA in 'extblankreads' and 'pcrblankreads' by 0
merged.table <- merged.table %>%
  mutate_at(vars(extblankreads, pcrblankreads), ~replace_na(.,0))

### Step8 ####
# Remove a number of reads from the sample 
# equal to the higher value of blank reads between Extraction and PCR Blanks
merged.table <- merged.table %>% 
  mutate(blankfiltered.reads =  case_when(extblankreads >= pcrblankreads ~
                                            reads-extblankreads,
                                          extblankreads < pcrblankreads ~ 
                                            reads-pcrblankreads,
                               TRUE ~ as.numeric(reads)))

# Remove lines with no reads after filtering
merged.table <- merged.table %>%         
           filter(blankfiltered.reads >  0)

### Step9 ####
# Designate which Phylum are dietary
levels(merged.table$Phylum)
c("Annelida","Arthropoda","Mollusca") -> diet.phylum
merged.table$Phylum %in% diet.phylum -> merged.table$diet

# Summarize the number of reads for dietary taxa per sample
diet_reads <- merged.table %>% 
  filter(diet == TRUE) %>% 
  drop_na(blankfiltered.reads) %>%
  group_by(sample) %>% 
  summarize(diet.reads = sum(blankfiltered.reads))

### Step10 ####
# Add the number of dietary reads per sample to each line
merged.table <- merged.table %>% left_join(diet_reads,
                        by = 'sample') %>%
  rename(total.diet.reads = diet.reads)

### Step11 ####
# Identify distribution of read lengths
bp.length <- ggplot(merged.table, aes(x=type, y=seq_length)) + 
  geom_boxplot()
# Make a function to find the MODE
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
bp.length + geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(title="Length of reads",x="OTUs", y = "Length (bp)") +
  stat_summary(fun=getmode, geom="point", shape=23, size=4) +
  theme_minimal() -> bp1

# Filter reads that are blanks
merged.table <- merged.table %>% mutate(diet.status =  case_when(type == 'extblank' | 
                                                              type == 'pcrblank' |
                                                              # Non-diet OTUs
                                                              diet == FALSE |
                                                              # not target length
                                                              #!seq_length == getmode(merged.table$seq_length)  |
                                                              # low match percentage
                                                              as.numeric(as.character(merged.table$Similarity))  < 90 |
                                                              # Fewer than 100 reads
                                                              total.diet.reads < 100 |
                                                              # Less than 1% of the total reads
                                                              blankfiltered.reads < 0.01*total.diet.reads ~ 
                                                              'delete',
                                                            TRUE ~ 'ok'
                                                            ))
## Step12 ###
# Filtered reads as proportion of total dietary reads
merged.table$'%reads' <-  merged.table$blankfiltered.reads/merged.table$total.diet.reads
# Set errors as NA
merged.table$'%reads'[merged.table$'%reads' <= 0] <- NA

# Filter reads by a given condition
merged.table <- merged.table %>% 
  filter(diet.status == 'ok')

write.table(merged.table, file = paste0(output.dir, 
                                      "/merged.table.",fragment[i],".txt"),
            sep = '\t',
            row.names = FALSE)
# Write a summary of the filtering 
{
  sink(paste0(output.dir,"/merge-log.txt"))
  cat(paste("Target Fragment:",fragment[i],"\n"))
  cat(paste("Initial OTUs:",dim(files[[1]])[1]),"\n")
  cat("Sample Submissions:","\n")
  cat(levels(files[[3]]$sample),sep = ",","\n")
  cat("Contamination OTUs:","\n")
  print(files[[2]][files[[2]]$otu %in% blanks$otu,], row.names=F)
  cat("Deitary groups:","\n")
  cat(diet.phylum, sep = ",","\n")
  cat("Average % Similarity of OTU assignments:","\n")
  cat(mean(as.numeric(as.character(
    files[[2]][files[[2]]$otu %in% unique(merged.table$otu),]$Similarity))), "+-",
    sd(as.numeric(as.character(
      files[[2]][files[[2]]$otu %in% unique(merged.table$otu),]$Similarity))),"\n"
  )
  cat("Total dietary reads per sample:","\n")
  print(as.data.frame(diet_reads), rown.names=F)
  cat(paste("Median read length:",getmode(merged.table$seq_length)),"\n")
  cat(paste("Final OTUs:",length(unique(merged.table$otu))),"\n")
sink()
}
jpeg(paste0(output.dir,"/otu-length-boxplot.jpeg"))
bp1
dev.off()
}

