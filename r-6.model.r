rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(reshape2)
library(randomForest)
library(glmnet)
library(ggplot2)
library(corrplot)
library(reshape2)
library(PRROC)
library(dcurves)
library(caret)

df.exp = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% lc.tableToNum()
hub.gene = read.delim2("05.HubGene/03.HubGene.xls")
df.exp = df.exp[hub.gene$GeneID %>% as.character, ]
rownames(df.exp) = hub.gene$Gene
df.exp = t(df.exp) %>% as.data.frame()
df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
df.exp$group = df.pheno$group[match(rownames(df.exp), df.pheno$geo_accession)]
df.exp$group = factor(df.exp$group)
colnames(df.exp) = make.names(colnames(df.exp))

set.seed(12345)
res.en = cv.glmnet(as.matrix(df.exp[-ncol(df.exp)]), df.exp$group, family = "binomial", nfolds = 10)
ggsave("08.Model/01.LASSO.CV.png", plot(res.en), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("08.Model/01.LASSO.CV.pdf", plot(res.en), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("08.Model/02.LASSO.Coef.png", plot(res.en$glmnet.fit, xvar = 'lambda'), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("08.Model/02.LASSO.Coef.pdf", plot(res.en$glmnet.fit, xvar = 'lambda'), width = 6, height = 6, dpi = 300, units = "in", bg = "white")

res.en = glmnet(as.matrix(df.exp[-ncol(df.exp)]), df.exp$group, family = "binomial", lambda = res.en$lambda.min)
cs.en = confusion.glmnet(res.en, newx = as.matrix(df.exp[-ncol(df.exp)]), newy = df.exp$group) %>% as.matrix()
cs.en = cs.en[2:1,2:1]
df.en = predict.glmnet(res.en, newx = as.matrix(df.exp[-ncol(df.exp)]), type = "link")[,1]
df.en = data.frame(sample = rownames(df.exp), score = df.en)
df.en$group = df.pheno$group[match(df.en$sample, df.pheno$geo_accession)]

png("08.Model/03.Confusion.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
corrplot(cs.en, method = "color", is.corr = F, tl.col = "black", 
         tl.cex = 1.5, cl.pos = "n", 
         col = colorRampPalette(c("white","royalblue"))(10))
corrplot(cs.en, method = "number", is.corr = F, cl.pos = "n",
         col = "black", number.digits = 0,
         number.cex = 2, add = T, bg = "transparent", tl.pos = 'n')
dev.off()
pdf("08.Model/03.Confusion.pdf", width = 4, height = 4)
corrplot(cs.en, method = "color", is.corr = F, tl.col = "black", 
         tl.cex = 1.5, cl.pos = "n", 
         col = colorRampPalette(c("white","royalblue"))(10))
corrplot(cs.en, method = "number", is.corr = F, cl.pos = "n",
         col = "black", number.digits = 0,
         number.cex = 2, add = T, bg = "transparent", tl.pos = 'n')
dev.off()

roc.en = roc.curve(df.en$score, weights.class0 = df.en$group == "Stroke", curve = T)
png("08.Model/04.ROC.png", width = 4, height = 5, res = 300, units = "in", bg = "white")
plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
abline(0,1); dev.off()
pdf("08.Model/04.ROC.pdf", width = 4, height = 5)
plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
abline(0,1); dev.off()

pr.en = pr.curve(df.en$score, weights.class0 = df.en$group == "Stroke", curve = T)
png("08.Model/05.PR.png", width = 4, height = 5, res = 300, units = "in", bg = "white")
plot(pr.en, auc.main = T, legend = F, color = F, asp = 1);dev.off()
pdf("08.Model/05.PR.pdf", width = 4, height = 5)
plot(pr.en, auc.main = T, legend = F, color = F, asp = 1);dev.off()

df.en$prob = (df.en$score - min(df.en$score))/(max(df.en$score) - min(df.en$score))
dca(group ~ prob, data = df.en, 
    label = list(prob = "LASSO")) %>% plot(smooth = T) + 
  theme(aspect.ratio = 1, legend.position = "top", text = element_text(size = 12))
ggsave("08.Model/06.DCA.png", width = 4, height = 4, dpi = 300)
ggsave("08.Model/06.DCA.pdf", width = 4, height = 4, dpi = 300)

saveRDS(res.en, "model.RDS")










