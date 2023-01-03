rm(list=ls())
library(magrittr)
library(stringr)
library(lance)
library(GEOquery)
library(dplyr)
library(data.table)
library(IRanges)
require(rtracklayer)

Sys.setenv("VROOM_CONNECTION_SIZE"=131072 * 50)
gs = getGEO("GSE58294", destdir = ".", getGPL = T, AnnotGPL = T)
df.exp = exprs(gs[[1]]) %>% as.data.frame %>% lc.tableToNum()
df.pheno = pData(gs[[1]]) %>% as.data.frame

df.annot = getGEO("GPL570", destdir = ".") %>% Table
#df.annot$gene = sapply(df.annot$gene_assignment, function(x){
#  x %>% str_split(" /// ") %>% unlist %>% (dplyr::first) %>% str_split(" // ") %>% unlist %>% (dplyr::last) %>% as.numeric()
#}) %>% as.character()
df.annot = subset(df.annot, !str_detect(ENTREZ_GENE_ID, "///"))
df.exp$gene = df.annot$ENTREZ_GENE_ID[match(rownames(df.exp), df.annot$ID)] %>% as.character
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
df.exp = subset(df.exp, gene != "")

write.table(df.exp, "01.Data/03.Validation.Expression.xls", sep = "\t", col.names = T, row.names = F)
write.table(df.pheno, "01.Data/04.Validation.MetaInfo.xls", sep = "\t", col.names = T, row.names = F)

rm(list=ls())
library(reshape2)
models = readRDS("model.RDS")
df.exp = read.delim2("01.Data/03.Validation.Expression.xls", row.names = 1) %>% lc.tableToNum()
hub.gene = read.delim2("05.HubGene/03.HubGene.xls")
df.exp = df.exp[hub.gene$GeneID %>% as.character, ]
rownames(df.exp) = hub.gene$Gene

df.exp.plot = cbind(gene = rownames(df.exp), df.exp)
df.exp.plot = reshape2::melt(df.exp.plot, id.vars = "gene", variable.name = "sample", value.name = "value")
df.exp = t(df.exp) %>% as.data.frame()

df.pheno = read.delim2("01.Data/04.Validation.MetaInfo.xls")
df.pheno$group = ifelse(df.pheno$group.ch1 == "Control", "Control", "Stroke")
df.pheno$sample = df.pheno$geo_accession
df.pheno = df.pheno[c("sample","group")]
df.exp = df.exp[df.pheno$sample,]
colnames(df.exp) = make.names(colnames(df.exp))

df.exp.plot$group = df.pheno$group[match(df.exp.plot$sample, df.pheno$sample)]
ggplot() +
  geom_boxplot(data = df.exp.plot, aes(x = group, y = value, color = group)) + 
  facet_wrap(gene ~ ., scales = "free", ncol = 2) +
  theme_bw() + 
  xlab(NULL) + 
  ylab("Expression") + 
  scale_color_manual(values = c("darkgreen", "darkorange"), guide = guide_legend(ncol = 1), labels = c("Low","High")) +
  theme(text = element_text(size = 16), legend.title = element_blank(), legend.position = "none")

library(glmnet)
library(pROC)
library(ggplot2)
df.pheno$pred.en = predict(models, newx = as.matrix(df.exp), type = "response")[,1]
roc.en = roc.curve(df.pheno$pred.en, weights.class0 = df.pheno$group == "Stroke", curve = T)
png("08.Model/07.Validation.png", width = 4, height = 5, res = 300, units = "in", bg = "white")
plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
abline(0,1); dev.off()
pdf("08.Model/07.Validation.pdf", width = 4, height = 5)
plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
abline(0,1); dev.off()

