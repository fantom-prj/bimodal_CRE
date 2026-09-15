library(dplyr) 
library(magrittr)
library(edgeR)
library(knitr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggrepel)
library(ggalt)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(grid)
library(data.table)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

##########################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")
CRE_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/")
path_fig2_data=paste0(primary_folder,"Fig2/data/")

#===============================================================================
setwd(CRE_folder)

#===============================================================================
#define TATA box from non-transcribed CRE

# bedtools intersect with TATA
setwd(paste0(CRE_folder,"CRE_bed/"))
system("bedtools intersect -wb -a untranscribed_aCRE.midpoint_summit250.bed.gz -b hg38.TATAbox.main.bed.gz | gzip > untranscribed_aCRE.midpoint_summit250.TATA.bed.gz")

#==================================
setwd(paste0(CRE_folder,"CRE_bed/"))
ump=read.delim("untranscribed_aCRE.midpoint_summit250.bed.gz", header=F, stringsAsFactors = F)
all=ump

TATAump=read.delim("untranscribed_aCRE.midpoint_summit250.TATA.bed.gz", header=F, stringsAsFactors = F)
length(unique(TATAump$V4)) #413633

TATA=TATAump[which(TATAump$V11 >=3),] #score >=3

TATA=left_join(TATA, all[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
TATA$locS=TATA$V2-TATA$V2_ori
TATA$locE=TATA$V3-TATA$V2_ori

#==================================
#bedtools
setwd(paste0(CRE_folder,"CRE_bed/INR/"))
system("sh bed_intersect250.sh")

#==================================
path1=paste0(CRE_folder,"CRE_bed/INR/")
files=list.files(path=path1, pattern="bed.gz")
files=files[grep("250",files)]
files=files[grep("untranscribed",files)]

INR=read.delim(paste0(path1,files[1]), header=F, stringsAsFactors = F)
for (i in 2:length(files)){
  INR1=read.delim(paste0(path1,files[i]), header=F, stringsAsFactors = F)
  INR=rbind(INR,INR1)}
INR=left_join(INR, all[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))

INR=INR[which(INR$V11 >= 48.5),]
INR$locS=INR$V2-INR$V2_ori
INR$locE=INR$V3-INR$V2_ori

all=inner_join(TATA[,c(4,10,11:15)], INR[,c(4,10,11:15)], by="V4", copy=F)
all=all[which(all$V12.x==all$V12.y),]
all$distance=all$locS.y-all$locS.x #TATA start to INR start distance, add 2 refer to TSS
all$distance[which(all$V12.x=="-")]=all$locE.x[which(all$V12.x=="-")]-all$locE.y[which(all$V12.x=="-")]
all1=all[which(all$V12.x == "+" & all$locS.x < all$locS.y),]
all2=all[which(all$V12.x == "-" & all$locS.x > all$locS.y),]
all1=all1[which(all1$distance > 0 & all1$distance < 60),]
all2=all2[which(all2$distance > 0 & all2$distance < 60),]

#width of TATA (avoid edge TATA and edge INR)
all1=all1[which(all1$locE.x-all1$locS.x == 7),]
all1=all1[which(all1$locE.y-all1$locS.y == 7),]
all2=all2[which(all2$locE.x-all2$locS.x == 7),]
all2=all2[which(all2$locE.y-all2$locS.y == 7),]

all=rbind(all1,all2)
colnames(all)=c("CREID","TATA_ID","TATA_score","TATA_strand","TATA_region501_start","TATA_locS","TATA_locE","INR_ID","INR_score","INR_strand","INR_region501_start","INR_locS","INR_locE","distance")
all0 = all |> dplyr::filter(distance >= 10, distance <= 45) 
write.table(all0, gzfile(paste0(CRE_folder,"CRE_bed/TATA_INR_prediction_untranscribed_CRE_midpoint250.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#do the same for transcribed CREs, to predict score cutoff
setwd(paste0(CRE_folder,"CRE_bed/"))
system("bedtools intersect -wb -a transcribed_aCRE.midpoint_summit250.bed.gz -b hg38.TATAbox.main.bed.gz | gzip > transcribed_aCRE.midpoint_summit250.TATA.bed.gz")

#=====
ump=read.delim(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.midpoint_summit250.bed.gz"), header=F, stringsAsFactors = F)
all=ump
TATAump=read.delim(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.midpoint_summit250.TATA.bed.gz"), header=F, stringsAsFactors = F)
length(unique(TATAump$V4)) #37203

TATA=TATAump[which(TATAump$V11 >=3),]

TATA=left_join(TATA, all[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
TATA$locS=TATA$V2-TATA$V2_ori
TATA$locE=TATA$V3-TATA$V2_ori

path1=paste0(CRE_folder,"CRE_bed/INR/")
files=list.files(path=path1, pattern="bed.gz")
files=files[grep("250",files)]
files=files[1:9] #take transcribed alone

INR=read.delim(paste0(path1,files[1]), header=F, stringsAsFactors = F)
for (i in 2:length(files)){
  INR1=read.delim(paste0(path1,files[i]), header=F, stringsAsFactors = F)
  INR=rbind(INR,INR1)}
INR=left_join(INR, all[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))

INR=INR[which(INR$V11 >= 48.5),]
INR$locS=INR$V2-INR$V2_ori
INR$locE=INR$V3-INR$V2_ori

all=inner_join(TATA[,c(4,10,11:15)], INR[,c(4,10,11:15)], by="V4", copy=F)
all=all[which(all$V12.x==all$V12.y),]
all$distance=all$locS.y-all$locS.x
all$distance[which(all$V12.x=="-")]=all$locE.x[which(all$V12.x=="-")]-all$locE.y[which(all$V12.x=="-")]
all1=all[which(all$V12.x == "+" & all$locS.x < all$locS.y),]
all2=all[which(all$V12.x == "-" & all$locS.x > all$locS.y),]
all1=all1[which(all1$distance > 0 & all1$distance < 60),]
all2=all2[which(all2$distance > 0 & all2$distance < 60),]

#width of TATA (avoid edge TATA and edge INR)
all1=all1[which(all1$locE.x-all1$locS.x == 7),]
all1=all1[which(all1$locE.y-all1$locS.y == 7),]
all2=all2[which(all2$locE.x-all2$locS.x == 7),]
all2=all2[which(all2$locE.y-all2$locS.y == 7),]

all=rbind(all1,all2)
colnames(all)=c("CREID","TATA_ID","TATA_score","TATA_strand","TATA_region501_start","TATA_locS","TATA_locE","INR_ID","INR_score","INR_strand","INR_region501_start","INR_locS","INR_locE","distance")
#bed files look for TSS from SCAFE
all$chr=sapply(strsplit(all$CREID,"-"),"[",1)
all$start=all$INR_region501_start+all$INR_locS
all$end=all$INR_region501_start+all$INR_locE
all$ID=paste0(all$chr,"_",all$start,"_",all$end,"_",all$INR_strand)
all_bed=unique(all[,c("chr","start","end","ID","start","INR_strand")])
all_bed$start.1="."
write.table(all_bed[order(all_bed$chr, all_bed$start),], gzfile(paste0(CRE_folder,"CRE_bed/TATA_INR_bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(all, gzfile(paste0(CRE_folder,"CRE_bed/TATA_INR_prediction_transcribed_CRE_midpoint250.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============
#identify INR region with transcription
setwd(paste0(primary_folder,"Data_and_code/SCAFE/sc5end/CTSS_in_merged_cluster"))
system(paste0("for file in *.bed.gz; do bedtools intersect -wa -wb -a ",CRE_folder,"CRE_bed/TATA_INR_bed.gz -b \"$file\" -s | gzip > \"TATA_INR_${file%.CTSS.bed.gz}.bed.gz\"; done")) 
files=list.files(pattern="TATA_INR")
data=data.frame()
for (i in 1: length(files)){
  data1=read.delim(files[i], header=F, stringsAsFactors = F)
  data1=data1%>%group_by(V4)%>%dplyr::summarise(CTSS=sum(V11),unGCTSS=sum(V10))
  data=rbind(data,data1)}
data=data%>%group_by(V4)%>%dplyr::summarise(CTSS=sum(CTSS),unGCTSS=sum(unGCTSS))
need=data$V4[which(data$unGCTSS>=3)] # at least 3 reads from iPSC/NSC/Neuron as positive transcription
all=read.delim(paste0(CRE_folder,"CRE_bed/TATA_INR_prediction_transcribed_CRE_midpoint250.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all$true=0
all$true[which(all$ID %in% need)]=1

# test for best distance
#remove redundancy
all0 = all |> dplyr::filter(distance >= 10, distance <= 45) |> dplyr::group_by(ID) |> dplyr::slice_max(TATA_score, n = 3, with_ties = FALSE) |> dplyr::ungroup()
library(splines)
fit_splines <- glm(true ~ TATA_score + INR_score + splines::ns(distance, df = 3), data = all0, family = binomial())
summary(fit_splines)
dd = data.frame(distance = seq(min(all0$distance), max(all0$distance), by = 1),
                TATA_score = median(all0$TATA_score),
                INR_score = median(all0$INR_score))
dd$prob = predict(fit_splines, newdata = dd, type = "response")

write.table(dd,gzfile(paste0(CRE_folder,"CRE_bed/TATA_INR_distance.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
dd95 = dd[dd$prob >= 0.95 * max_prob, ]
range(dd95$distance) #25-33

ggplot(dd, aes(distance, prob)) +
  geom_line(linewidth=0.2)+
  geom_vline(xintercept=c(25,33), linetype="dashed", linewidth=0.2)+
  labs(x = "TATA–Inr distance (nt)",y = "Predicted probability") +
  theme1

write.table(dd,gzfile(paste0(path_fig2_data,"TATA_INR_distance.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#-> for fig ex2j

#===============
# define TATA and INR motif score cutoff
#remove redundancy
all1 = all |> dplyr::filter(distance >= 25, distance <= 33) |> dplyr::group_by(ID) |> dplyr::slice_max(TATA_score, n = 1, with_ties = FALSE) |> dplyr::ungroup()

#separate train and test
set.seed(123)
id <- sample(seq_len(nrow(all1)), size = 0.7 * nrow(all1))
train <- all1[id, ]
test  <- all1[-id, ]

fit <- glm(true ~ TATA_score + INR_score, data = train, family = binomial())
summary(fit)

library(pROC)
test$score <- predict(fit, newdata = test, type = "link")
test$score1=test$TATA_score+test$INR_score
roc1 <- roc(test$true, test$score)
auc1=auc(roc1) #0.5809
roc_df <- data.frame(
  specificity = roc1$specificities,
  sensitivity = roc1$sensitivities,
  threshold = roc1$thresholds)
write.table(roc_df,gzfile(paste0(CRE_folder,"CRE_bed/TATA_ROC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

tata1=ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "blue", size=0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.2)+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", signif(auc1, 3)), size = 2, color = "black")+
  labs(title = "ROC Curve for detecting TATA-box", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme1
pdf(paste0(CRE_folder,"CRE_bed/TATA_ROC.pdf"), width = 2, height = 1.7)
print(tata1)
dev.off() 

cutoff1=coords(roc1, "best", ret = "threshold", best.method = "youden")$threshold #-2.363024
cutoff95=roc_df$threshold[which.min(abs(roc_df$specificity - 0.95))] #-1.96714
length(which(test$score>cutoff95))/nrow(test) #0.05777

roc2 <- roc(test$true, test$score1)
auc2=auc(roc2) #0.5803
cutoff2=coords(roc2, "best", ret = "threshold", best.method = "youden")$threshold #54.8540995

#======================================
#use cutoff95
unt_all=read.delim(paste0(CRE_folder,"CRE_bed/TATA_INR_prediction_untranscribed_CRE_midpoint250.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
unt_all$score <- predict(fit, newdata = unt_all, type = "link")
unt_all1 = unt_all |> dplyr::filter(distance >= 25, distance <= 33) |> dplyr::group_by(INR_region501_start,INR_locS) |> dplyr::slice_max(TATA_score, n = 1, with_ties = FALSE) |> dplyr::ungroup()
unt_all2=unt_all1[which(unt_all1$score > cutoff95),]
length(unique(unt_all2$CREID))#24526
write.table(unt_all2, gzfile(paste0(CRE_folder,"CRE_bed/TATA_INR_prediction_untranscribed_CRE_midpoint250_about_cutoff95.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
