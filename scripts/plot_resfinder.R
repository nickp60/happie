library(ggplot2)
library(tidyverse)

classes <- read.csv("https://bitbucket.org/genomicepidemiology/resfinder_db/raw/fc516c662c9e535f8cc3022c43891fa9f80dc537/notes.txt", header=F, stringsAsFactors = F)
classes_df <- data.frame("class"=character(), "name"=character(), stringsAsFactors = F)
class = NA
for (i in 1:nrow(classes)){
  t <- classes[i, 1]
  if (startsWith(t, "#")){
    class = gsub("#", "", t)
  } else {
    tmp <- data.frame("class"=class, "name"=t)
    classes_df <- rbind(classes_df, tmp)
  }
} 
classes_df <- classes_df %>% separate(name, into = c("GENE", "res", "note"), sep = ":")
# get rid of annoying leading white space
classes_df$GENE <- as.character(classes_df$GENE)
classes_df$GENE <-  gsub("^[ \\(\\)]*", "", classes_df$GENE)

dat <- read.csv2("~/Desktop/2018-12-02_happie_resfinder.tab", sep="\t", header = F, comment.char = "#", stringsAsFactors = F)
headers <- c("FILE", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS", "PERCCOVERAGE", "PERCIDENTITY", "DATABASE", "ACCESSION", "PRODUCT")
colnames(dat) <- headers

# make a simplified genome column, making everything lowercase, removing everything after the first bit
(classes_df$gene <- gsub("(.*?)[^a-z0-9](.*)", "\\1", tolower(classes_df$GENE)))
(dat$gene <- gsub("(.*?)[^a-z0-9](.*)", "\\1", tolower(dat$GENE)))

str(classes_df$gene)
str(dat$gene)
table(dat$gene %in% classes_df$gene)
dat$gene[!dat$gene %in% classes_df$gene]

dat2 <- left_join(dat, classes_df, by="gene")

ggplot(dat2, aes(x=class)) + geom_histogram(stat="count") + coord_flip() + scale_y_log10()

ggplot(dat_sum, aes(x=n)) + geom_histogram(stat="count") + coord_flip()
