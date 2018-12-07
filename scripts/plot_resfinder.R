library(ggplot2)
library(tidyverse)
if(!require(devtools)) install.packages("devtools")
if(!require(ggpubr)) devtools::install_github("kassambara/ggpubr")
library(ggpubr)

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


ggplot(dat2, aes(x=class)) + geom_histogram(stat="count") + coord_flip() + scale_y_log10()

ggplot(dat_sum, aes(x=n)) + geom_histogram(stat="count") + coord_flip()


#####################
mobile <- read.csv2("~/GitHub/happie/2018-12-04-happie-mobile-coords", header = F, sep = "\t", stringsAsFactors = F)
mobilenames <- c("path", "tool", "feature", "contig", "start", "end")
colnames(mobile) <- mobilenames
mobile$source <- ifelse(grepl("Lys", mobile$path), "soil", "clinical")

headers <- c("path", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS", "PERCCOVERAGE", "PERCIDENTITY", "DATABASE", "ACCESSION", "PRODUCT")
abricate <- read.csv2("~/GitHub/happie/2018-12-04-happie-abricate", header = F, sep = "\t", comment.char = "#", stringsAsFactors = F)
colnames(abricate) <- headers
abricate$source <- ifelse(grepl("Lys", abricate$path), "soil", "clinical")

dat <- abricate[abricate$DATABASE=="resfinder", ]


su <- mobile %>% 
  group_by(source, path, feature) %>%
  summarise(n=n()) %>%
  as.data.frame()
  
ggplot(su, aes(x=feature, y=n, fill=source)) +
  geom_boxplot() + stat_compare_means(n ~ feature+source)

ggboxplot(su, x = "source", y = "n",
          color = "source", palette = "jco", facet.by = "feature")+ 
#  stat_compare_means(label.y = c(29, 35, 40))+
  stat_compare_means(label.y = 45)
for (i in c("prophages")){
  print(t.test(
    su[su$source=="soil" & su$feature==i, "n"], 
    su[su$source=="clinical" & su$feature==i, "n"], 
    alternative = "two.sided", 
    var.equal = FALSE)
  )
}
classes <- read.csv("https://bitbucket.org/genomicepidemiology/resfinder_db/raw/fc516c662c9e535f8cc3022c43891fa9f80dc537/notes.txt", header=F, stringsAsFactors = F)
classes_df <- data.frame("class"=character(), "name"=character(), stringsAsFactors = F)
class = NA
for (i in 1:nrow(classes)){
  t <- classes[i, 1]
  if (startsWith(t, "#")){
    class = gsub(":", "", gsub("#", "", t))
  } else {
    tmp <- data.frame("class"=class, "name"=t)
    classes_df <- rbind(classes_df, tmp)
  }
} 
classes_df <- classes_df %>% separate(name, into = c("GENE", "res", "note"), sep = ":")

# make a simplified genome column, making everything lowercase, removing everything after the first bit
classes_df$gene <- sapply(strsplit(gsub(" ", "", gsub("\\(", "", gsub("\\)", "", classes_df$GENE))), "\\, |\\,| |\\-|1|2|3|4|5|6|7|8|9|0|_"), function(x){x[[1]]})
dat$gene <- sapply(strsplit(gsub(" ", "", gsub("\\(", "", gsub("\\)", "", dat$GENE))), "\\, |\\,| |\\-|1|2|3|4|5|6|7|8|9|0|_"), function(x){x[[1]]})

str(classes_df$gene)
str(dat$gene)
table(dat$gene %in% classes_df$gene)
dat$gene[!dat$gene %in% classes_df$gene]

dat2 <- left_join(dat, classes_df, by="gene")

dat2 <- dat2 %>%
  group_by(source) %>%
  mutate(nsource=n()) %>%
  group_by(source, class) %>%
  mutate(n=n(), perc= n/nsource * 100) %>% as.data.frame()

ggplot(dat2, aes(x=gene, y=path)) + geom_tile()
ggplot(dat2, aes(x=class, y=perc, group=source) )+
  #geom_boxplot()+
  geom_point() + 
  coord_flip() + scale_y_log10()

dat3 <- dat2[!duplicated(dat2[, c("source", "class")]), c("source", "class", "perc")]
t.test(dat)


#install.packages("caret")  

library(caret)
set.seed(12345)
test_i <- createDataPartition(y=dat2$source, times = 1, p = .33)$Resample1
test <-  dat2[train_i,]
train_hold <- dat2[-test_i,]
train_i <- createDataPartition(y=train_hold$source, times = 1, p = .5)$Resample1
train <- train_hold[train_i, ]
hold <- train_hold[-train_i, ]
  

control <- trainControl(method="repeatedcv", number=10, repeats=3)
metric <- "Accuracy"
mtry <- sqrt(ncol(train))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(source~., data=train, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)



