rm(list = ls())
library(magrittr)
library(stringr)
library(org.Hs.eg.db)
library(GSVA)

df = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% apply(., c(1,2), as.numeric) %>% as.data.frame
gene.annot = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(df), keytype = "ENTREZID", columns = "SYMBOL")
gene.annot = na.omit(gene.annot)
df = df[gene.annot$ENTREZID,]
rownames(df) = gene.annot$SYMBOL
gene.set = read.delim2("ssgsea.genesets.tsv")
gene.set.list = split(as.matrix(gene.set)[,1], gene.set[,2])
df.im = gsva(as.matrix(df), gene.set.list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
df.im = as.data.frame(df.im)
df.im = cbind(Immune = rownames(df.im), df.im)  
write.table(df.im, "03.Immune/01.ImmuneProfile.xls", sep = "\t", quote = F, col.names = T, row.names = F)

rm(list = ls())
df.im = read.delim2("03.Immune/01.ImmuneProfile.xls", row.names = 1) %>% apply(.,c(1,2),as.numeric) %>% as.data.frame
df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
df.pheno$group = factor(df.pheno$group, levels = c("Control", "Stroke"))

df.im = df.im[df.pheno$geo_accession]
p.diff = apply(df.im, 1, function(x){
  wilcox.test(as.numeric(x)~df.pheno$group)$p.value
})
df.diff = data.frame(Immune = rownames(df.im), p.value = p.diff, p.adj = p.adjust(p.diff, method = "bonferroni"))
write.table(df.diff, "03.Immune/02.ImmuneDiff.xls", sep = "\t", quote = F, col.names = T, row.names = F)

library(pheatmap)
col.mat = df.pheno[c("geo_accession", "group")]
col.mat$group = factor(col.mat$group, levels = c("Control", "Stroke"))
col.mat = col.mat[order(col.mat$group, decreasing = F), ]
df.im = df.im[col.mat$geo_accession]
rownames(col.mat) = col.mat[[1]]
colnames(col.mat)[2] = "Group"
col.mat = col.mat[-1]
p = pheatmap(df.im, cluster_cols = F, show_colnames = F, annotation_col = col.mat,
             border_color = "white", scale = "row",
             annotation_colors = list(Group = c(Control = "darkgreen",Stroke = "darkorange")))
png("03.Immune/03.Heatmap.png", width = 10, height = 5, units = "in", res = 300)
print(p)
dev.off()
pdf("03.Immune/03.Heatmap.pdf", width = 10, height = 5)
print(p)
dev.off()

library(ggplot2)
library(reshape2)
rm(list = ls())
df.im = read.delim2("03.Immune/01.ImmuneProfile.xls", row.names = 1) %>% apply(.,c(1,2),as.numeric) %>% as.data.frame
df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
df.pheno$group = factor(df.pheno$group)

df.im = cbind(im = rownames(df.im), df.im)
df.im4plot = melt(df.im, id.vars = "im", variable.name = "sample", value.name = "value")

df.im4plot$group = df.pheno[match(df.im4plot$sample, df.pheno$geo_accession),]$group
df.im4plot$group = factor(df.im4plot$group, levels = c("Control", "Stroke"))

df.diff = read.delim2("03.Immune/02.ImmuneDiff.xls", row.names = 1) %>% apply(.,c(1,2),as.numeric) %>% as.data.frame
df.sig = subset(df.diff, p.adj < 0.05)
df.sig = cbind(im = rownames(df.sig), df.sig)
df.sig$p = signif(df.sig$p.value,2) %>% paste0("p = ", .)
df.im4plot = subset(df.im4plot, im %in% rownames(df.sig))

ggplot() + 
  geom_boxplot(data = df.im4plot, mapping = aes(x = im, y = value, color = group)) +
  geom_text(inherit.aes = F, data = df.sig, mapping = aes(x = im, y = 1, label = p), angle = 60, color = "red") +
  xlab("Immune Cells") + ylab("Score") +
  scale_color_manual(values = c("darkblue","darkred")) +
  theme_bw() + ylim(c(0,1.2)) +
  theme(axis.text.x = element_text(angle = 60, hjust=.95, vjust=.95, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16), legend.text = element_text(size = 14),
        legend.title = element_blank(), legend.position = "top")
ggsave("03.Immune/04.Boxplot.png", width = 12, height = 8, units = "in", dpi = 300)
ggsave("03.Immune/04.Boxplot.pdf", width = 12, height = 8, units = "in", dpi = 300)




