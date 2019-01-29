library(tidyverse)
if(!require(devtools)) install.packages("devtools")
if(!require(ggpubr)) devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(caret)
########################  Putting the "fun" in "functions" ############
remove_linear_combos <- function(df, ignor){
  outcomes <- df[, ignor]
  dat = df %>%
    select(- ignor) %>%
    select(-findLinearCombos(.)$remove)
  dat[, ignor] = outcomes
  return(dat)
}
make_partitioned_df <- function(df, interest, p1=.33, p2=.5){
  train_i <- createDataPartition(y=df[, interest], times = 1, p = p1)$Resample1
  traindf <- df[train_i, ]
  traindf$partition <- "train"
  test_hold <- df[-train_i, ]
  test_i <- createDataPartition(y=test_hold[, interest], times = 1, p = p2)$Resample1
  testdf <-  test_hold[test_i, ]
  testdf$partition <- "test"
  holddf <- test_hold[-test_i, , ]
  holddf$partition <- "holdout"
  return(rbind(traindf, testdf, holddf))
}



# classes <- read.csv("https://bitbucket.org/genomicepidemiology/resfinder_db/raw/fc516c662c9e535f8cc3022c43891fa9f80dc537/notes.txt", header=F, stringsAsFactors = F)
# classes_df <- data.frame("class"=character(), "name"=character(), stringsAsFactors = F)
# class = NA
# for (i in 1:nrow(classes)){
#   t <- classes[i, 1]
#   if (startsWith(t, "#")){
#     class = gsub("#", "", t)
#   } else {
#     tmp <- data.frame("class"=class, "name"=t)
#     classes_df <- rbind(classes_df, tmp)
#   }
# } 
# classes_df <- classes_df %>% separate(name, into = c("GENE", "res", "note"), sep = ":")
# # get rid of annoying leading white space
# classes_df$GENE <- as.character(classes_df$GENE)
# classes_df$GENE <-  gsub("^[ \\(\\)]*", "", classes_df$GENE)
# 
# dat <- read.csv2("~/GitHub/FB/Ecoli_comparative_genomics/results/2018-12-02_happie_resfinder.tab", sep="\t", header = F, comment.char = "#", stringsAsFactors = F)
# headers <- c("FILE", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS", "PERCCOVERAGE", "PERCIDENTITY", "DATABASE", "ACCESSION", "PRODUCT")
# colnames(dat) <- headers
# 
# 
# ggplot(dat2, aes(x=class)) + geom_histogram(stat="count") + coord_flip() + scale_y_log10()
# 
# ggplot(dat_sum, aes(x=n)) + geom_histogram(stat="count") + coord_flip()
# 

#####################
mobile <- read.csv2("~/GitHub/happie/2018-12-04-happie-mobile-coords", header = F, 
                    sep = "\t", stringsAsFactors = F, 
                    col.names =  c("path", "tool", "feature", "contig", "start", "end"))
mobile$source <- ifelse(grepl("Lys", mobile$path), "soil", "clinical")

abheaders <- c("path", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS", "PERCCOVERAGE", "PERCIDENTITY", "DATABASE", "ACCESSION", "PRODUCT")
abricate <- read.csv2("~/GitHub/happie/2018-12-04-happie-abricate", header = F, 
                      sep = "\t", comment.char = "#", stringsAsFactors = F,
                      col.names = abheaders)
abricate$short_name <- gsub(".fasta", "", 
                            gsub("\\/(.*?)\\/(.*?)\\/(.*?)\\/(.*?)\\/(.*)", "\\4", abricate$path))

table(abricate$DATABASE)

abricate$source <- ifelse(grepl("Lys", abricate$path), "soil", "clinical")

dat <- abricate %>% filter(DATABASE=="resfinder")


su <- mobile %>% 
  group_by(source, path, feature) %>%
  summarise(n=n()) %>%
  as.data.frame()
  
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
# make a simplified gene column, making everything lowercase, removing everything after the first bit
classes_df$gene <- sapply(strsplit(gsub(" ", "", gsub("\\(", "", gsub("\\)", "", classes_df$GENE))), "\\, |\\,| |\\-|1|2|3|4|5|6|7|8|9|0|_"), function(x){x[[1]]})
# drop complex gene names
classes_df <- classes_df %>% select(-GENE, -note) %>% distinct()

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

ggplot(dat2, aes(x=source)) + geom_histogram(stat="count") + labs(title="total AMR hits")
ggplot(dat2 %>% distinct(path, source), aes(x=source)) + geom_histogram(stat="count")  + labs(title="Organisms with AMR hits")
ggplot(dat2, aes(x=gene, y=short_name)) + geom_tile() + labs(title="AMR heatmap")
# ggplot(dat2, aes(x=class, y=perc, group=source))+
#   geom_point() + 
#   coord_flip() + scale_y_log10()
# ggplot(dat2, aes(x=class, y=perc, color=source, group=interaction(source, class)))+
#   geom_boxplot()+
#   geom_point() + 
#   coord_flip() + scale_y_log10()



#install.packages("caret")  
dat3 <-dat2 %>% select(short_name, source, class) %>% distinct()
wd3 <- dat3 %>% group_by(short_name, class) %>% mutate(classn=n()) %>% ungroup() %>% spread_(key_col = "class", value_col="classn", fill=0) %>% select(-short_name)
wd3$source <- ifelse(grepl("clinical", wd3$source), 0, 1)
str(wd3)
set.seed(12345)
test_i <- createDataPartition(y=wd3$source, times = 1, p = .33)$Resample1
test <-  wd3[test_i,]
train_hold <- wd3[-test_i,]
train_i <- createDataPartition(y=train_hold$source, times = 1, p = .5)$Resample1
traindf <- train_hold[train_i, ]
holddf <- train_hold[-train_i, ]
  

control <- trainControl(method="repeatedcv", number=10, repeats=3)
metric <- "Accuracy"
mtry <- sqrt(ncol(traindf))
tunegrid <- expand.grid(.mtry=mtry)
# rf_default <- train(
#   source~., 
#   data=traindf,
#   method="rf")
# print(rf_default)

ctrl <- 
  trainControl(
    method = "repeatedcv", 
    number = 10,
    repeats = 5,
    classProbs = F,
    summaryFunction = twoClassSummary)
# this fails -- we dont have enough positive hits for out soil saples
train(source ~ .,
      data=traindf,
      method = "fr",
      family = "binomial",
      metric = "ROC",
      trControl = ctrl)
glm_1 <- train(source~., method="glmboost", data=traindf)

summary(glm_1$finalModel)


########################  VF
raw_vf <- abricate %>% filter(DATABASE =="vfdb") %>% select_("short_name", "source", "GENE", "PRODUCT")
# get keywords
#download.file("http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz", "./VFs.xls.gz")
#system("gunzip VFs.xls.gz")
# cant seem to read in the xl file directly; open, export as csv
# desc_file <- readxl::read_excel(path = "./VFs.xls", sheet = 1, skip = 1)
desc_file <- read.csv("~/GitHub/happie/validation/VFs.csv", skip = 1, stringsAsFactors = F) %>% 
  select_("Keyword", "VFID")

# not all these have matches, cause the E. coli ones dont have VF ids
raw_vf$id <- gsub(".*(VF....).*", "\\1", raw_vf$PRODUCT)
vf_with_match <- left_join(raw_vf[grepl("VF...", raw_vf$id), ], desc_file, by = c("id" = "VFID"))
vf_no_match <- raw_vf[!grepl("VF...", raw_vf$id), ]
vf_no_match$id <- NA
vf_no_match<- vf_no_match %>% 
  mutate(
    Keyword = 
      case_when(
        grepl("fimbria", vf_no_match$PRODUCT, ignore.case = T) ~ "Adherence",
        grepl("Type III secretion", vf_no_match$PRODUCT, ignore.case = T) ~ "Secretion; Type III secretion",
        grepl("T3SS", vf_no_match$PRODUCT, ignore.case = T) ~ "Secretion; Type III secretion",
        grepl("esterase", vf_no_match$PRODUCT, ignore.case = T) ~ "Metabolic adaptation",
        grepl("ATP binding cassette", vf_no_match$PRODUCT, ignore.case = T) ~ "ATP transporter",
        grepl("glucosyltransferase", vf_no_match$PRODUCT, ignore.case = T) ~ "Toxin; Intracellular toxin; Glucosyltransferase
",
        grepl("colinisation", vf_no_match$PRODUCT, ignore.case = T) ~ "Adherence; Type IV pili",
        grepl("Chaperone", vf_no_match$PRODUCT, ignore.case = T) ~ "Intracellular growth; Protein folding",
        grepl("negative regulator", vf_no_match$PRODUCT, ignore.case = T) ~ "Regulation",
        grepl("chemotaxis", vf_no_match$PRODUCT, ignore.case = T) ~ "Adherence",
        grepl("coli surface antigen", vf_no_match$PRODUCT, ignore.case = T) ~ "Adherence",
        TRUE ~ "unknown"
      )
  )
  
vf <- rbind(vf_with_match, vf_no_match)
vf$gene<-vf$Keyword
vf <- vf %>% select(-PRODUCT, -id, -GENE, -Keyword)

# simplify gene
#vf$gene <- gsub("(.*?)[\\/-].*", "\\1", vf$GENE)
#vf <-vf[!duplicated(vf), colnames(vf)!="GENE"]

vf <- vf %>% group_by(short_name, gene)%>% mutate(genen=n()) %>% distinct()
vfs <- vf %>% transform(gene = strsplit(gene, "; ")) %>% unnest(gene) %>% 
  group_by(short_name, gene) %>% mutate(genet=sum(genen)) %>% 
  distinct(short_name, source, gene, genet)  %>%
  transform(source=as.factor(source)) %>%
  spread(key = gene, value = genet, fill=0 ) %>% as.data.frame()
vfs <- vf %>% transform(gene = strsplit(gene, "; ")) %>% unnest(gene) %>% 
  group_by(short_name, gene) %>% mutate(genet=1) %>% 
  distinct(short_name, source, gene, genet)  %>%
  transform(source=as.factor(source)) %>%
  spread(key = gene, value = genet, fill=0 ) %>% as.data.frame()

rownames(vfs) <- vfs$short_name
vfs <- vfs %>% select(-short_name)
#vfs$source<- as.factor(ifelse(grepl("Lys", rownames(vfs)), "soil", "clinical"))
#vfs$source <- ifelse(grepl("clinical", vfs$source), 0, 1)
#vfs$source <- ifelse(grepl("clinical", vfs$source), 0, 1)
str(vfs)

set.seed(12345)
vfs_trim <-vfs %>% do(remove_linear_combos(., "source"))
mlvf <- make_partitioned_df(df=vfs_trim, interest="source", p1=.33, p2=.5)
table(mlvf$partition)
metric <- "Accuracy"
mtry <- sqrt(ncol(traindf))
tunegrid <- expand.grid(.mtry=mtry)
ctrl <- 
  trainControl(
    method = "repeatedcv", 
    #number = 10,
    repeats = 5,
    classProbs = T,
    summaryFunction = twoClassSummary)

fit <- train(source~.,
             data=mlvf %>% 
               filter(partition %in% c("train", "holdout")) %>% 
               select(-partition) %>%
               do(remove_linear_combos(., "source")),
             method = "rf",
             family = "binomial",
             metric = "ROC",
             trControl = ctrl)
fit <- train(source~.,
             data=mlvf,
             method = "rf",
             family = "binomial",
             metric = "ROC",
             trControl = ctrl)
fit$finalModel
summary(fit)
varImp(fit)
confusionMatrix(predict(fit, 
                        mlvf %>% filter(partition %in% c("train", "holdout")) %>% 
                          select(-partition)), 
                mlvf[mlvf$partition %in% c("train", "holdout"), "source"])

confusionMatrix(predict(fit, 
                        mlvf %>% filter(partition %in% c("test")) %>% 
                          select(-partition)), 
                mlvf[mlvf$partition %in% c("test"), "source"])

fit

#################  plotting size of pangenomes
mob <- mobile %>% group_by(path, feature) %>% mutate(size=sum(end-start)) %>% select(path, size, source, feature) %>% distinct()
ggplot(mob, aes(x=source, y=size, fill=feature)) + geom_boxplot()


#### ML on presence absense


pg_raw <- read.table("~/GitHub/FB/Ecoli_comparative_genomics/results/2018-12-05-happie-roary/2018-12-05-roary/gene_presence_absence.Rtab", header = T)
str(pg_raw)
#pg_raw <- read.csv("~/GitHub/FB/Ecoli_comparative_genomics/results/2018-12-05-happie-roary/2018-12-05-roary/gene_presence_absence.csv")
#pg <-pg_raw[,c(1, 15:ncol(pg_raw))]

#  Fitting right away leeps crashing, so lets thin the data by removing singletons and dupletons
pg <-  pg_raw %>% mutate(n=rowSums(pg_raw %>% select(-Gene))) %>%filter(n > 2) %>% select(-n)

pg <- pg %>%
  #head() %>%
  gather("strain", "present", 2:ncol(.)) %>%
  spread(Gene, present) %>%
  mutate(source = ifelse(grepl("Lys", strain), 1, 0)) #%>% View()
mlpg <- make_partitioned_df(df=pg, interest="source", p1=.33, p2=.5)
#findLinearCombos(pg %>% select(-strain, -source))

ctrl <- 
  trainControl(
    method = "repeatedcv", 
    number = 10,
    repeats = 5,
    classProbs = T,
    summaryFunction = twoClassSummary)
thisdata=mlpg %>% 
  filter(partition=="train") %>% 
  select(-partition) %>%
  select(-strain) %>% do(remove_linear_combos(., ignor = "source"))
pgfit <- train(source~.,
             data=thisdata%>% transform(source=factor(source, labels=c("clinical", "soil"))),
             method = "rpart",
             #family = "binomial",
             metric = "ROC",
             trControl = ctrl)
summary(pgfit)
pgfit$finalModel
confusionMatrix(predict(pgfit, 
                        mlpg %>% filter(partition=="test") %>% 
                          select(-partition)), 
                mlvf[mlvf$partition=="test", "source"])

rcfit <- rpart::rpart(source~., data = thisdata, method = "class")
summary(rcfit)
plot(rcfit )
text(rcfit, use.n = TRUE)


