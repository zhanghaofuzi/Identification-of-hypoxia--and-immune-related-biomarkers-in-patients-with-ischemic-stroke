rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(limma)
library(org.Hs.eg.db)

df = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% lc.tableToNum
df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df = df[df.pheno$geo_accession]
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
design.mat = cbind(Control = ifelse(df.pheno$group == "Control", 1, 0), 
                   Patient = ifelse(df.pheno$group == "Control", 0, 1))
contrast.mat = makeContrasts(contrasts="Patient-Control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
fit = fit[c(1,3,4)]
fit.annot = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(fit), keytype = "ENTREZID", columns = "SYMBOL")
fit.annot = subset(fit.annot, !duplicated(fit.annot$ENTREZID))
fit = cbind(GeneID = fit.annot$ENTREZID, Gene = fit.annot$SYMBOL, fit)
fit = na.omit(fit)
fit$Direction = ifelse(fit$logFC > 0, "Up", "No")
fit$Direction = ifelse(fit$logFC < 0, "Down", fit$Direction)
fit$Direction = ifelse(fit$adj.P.Val < 0.05, fit$Direction, "No")
fit$Direction = factor(fit$Direction, levels = c("Down", "No", "Up"))
fit = fit[order(fit$adj.P.Val, decreasing = F), ]
fit.out = subset(fit, Direction != "No")
write.table(fit.out, "02.DEG/01.DEG.xls", sep = "\t", quote = F, col.names = T, row.names = F)

library(ggplot2)
library(ggrepel)
ggplot(fit, aes(x = logFC, y = -log10(adj.P.Val), color = Direction)) + 
  geom_point(size = 1, alpha = 0.7, shape = 20) + 
  scale_color_manual(values=c("darkblue", "grey30", "darkred"), labels = c("Down", "No Change", "Up")) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype = 2) +
  # geom_vline(xintercept = c(-1,1), col="grey60", linetype = 2) +
  ggtitle("Stroke vs Control") +
  ylab("-log10(adj.P.Val)") +
  xlab("log2FC") +
  xlim(c(-2,2)) +
  geom_label_repel(
    data = fit[which(fit$Gene %in% hub_gene),],
    aes(label = fit[which(fit$Gene %in% hub_gene),]$Gene),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "grey",
    segment.size =0.5,
    segment.alpha = 0.5,
    show.legend = F,
    direction = "y",
    # hjust = 0,
    max.overlaps = 30,
    xlim = 1.6,
    # ylim = c(2, max_padj * 0.9)
  ) +
  theme_bw(base_size = 16) + 
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        legend.text = element_text(size = 12), 
        axis.text = element_text(face = "bold"), legend.position = c(0.1,0.9),
        legend.background = element_blank())
ggsave("02.DEG/02.Volcano_20221231_label.png", width = 10, height = 5, dpi = 300, bg = "white")
ggsave("02.DEG/02.Volcano_20221231_label.pdf", width = 10, height = 5, dpi = 300, bg = "white")

library(pheatmap)
fit = fit[order(fit$Direction, fit$P.Value, decreasing = F),]
fit = na.omit(fit)
fit.sub = rbind(subset(fit, Direction == "Up")[1:15,],
                subset(fit, Direction == "Down")[1:15,])
df.sub = df[rownames(fit.sub),]
rownames(df.sub) = fit.sub$Gene

col.mat = data.frame(Sample = df.pheno$geo_accession,
                     Group = ifelse(df.pheno$group =="Control", "Control", "Stroke"))
col.mat$Group = factor(col.mat$Group, levels = c("Control","Stroke"))
col.mat = col.mat[order(col.mat$Group),]
rownames(col.mat) = col.mat$Sample

row.mat = data.frame(Direction = factor(fit.sub$Direction, levels = c("Up", "Down")))
rownames(row.mat) = rownames(df.sub)
df.sub = df.sub[rownames(col.mat)]

p = pheatmap(df.sub, border_color = "black", scale = "row",
             color = colorRampPalette(c("green", "black", "red"))(50),
             cluster_rows = T, cluster_cols = F, show_colnames = F, 
             annotation_colors = list(Group = c(Control = "darkgreen", Stroke = "darkorange"),
                                      Direction = c(Up = "red", Down = "blue")),
             annotation_names_col = F, annotation_names_row = F,
             annotation_col = col.mat[-1], annotation_row = row.mat)

png("02.DEG/03.Heatmap.png", width = 8, height = 4, units = "in", res = 300)
print(p)
dev.off()
pdf("02.DEG/03.Heatmap.pdf", width = 8, height = 4)
print(p)
dev.off()

rm(list = ls())
ferr.gene = read.delim2("ferr.genes")
deg = read.delim2("02.DEG/01.DEG.xls")

p.venn = lc.vennFromList(list(DEG = deg$Gene, Ferroptosis = ferr.gene$SYMBOL), label.dist = -20)
png("02.DEG/04.Venn.png", width = 5, height = 5, res = 300, units = "in", bg = "white")
grid.newpage();grid.draw(p.venn);dev.off()
pdf("02.DEG/04.Venn.pdf", width = 5, height = 5)
grid.newpage();grid.draw(p.venn);dev.off()

deg.sub = subset(deg, GeneID %in% ferr.gene$ENTREZID)
write.table(deg.sub, "02.DEG/05.Ferroptosis.DEG.xls", quote = F, sep = "\t", row.names = F, col.names = T)

#deg.sub = deg.sub[order(deg.sub$adj.P.Val, decreasing = F),]
#deg.plot = rbind(subset(deg.sub, Direction == "Up")[1:10,],
#                 subset(deg.sub, Direction == "Down")[1:10,])

df = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% lc.tableToNum
df = df[as.character(deg.sub$GeneID),]

df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
col.mat = data.frame(Sample = df.pheno$geo_accession,
                     Group = factor(df.pheno$group, levels = c("Control", "Stroke")))
col.mat = col.mat[order(col.mat$Group),]
df = df[col.mat$Sample]
rownames(col.mat) = col.mat$Sample
col.mat = col.mat[-1]
row.mat = data.frame(Direction = factor(deg.sub$Direction, levels = c("Up", "Down")))
rownames(row.mat) = deg.sub$Gene
rownames(df) = deg.sub$Gene
p = pheatmap(df, border_color = "black", scale = "row",
             color = colorRampPalette(c("green","green","black","red","red"))(100),
             cluster_rows = T, cluster_cols = F, show_colnames = F, 
             annotation_colors = list(Group = c(Control = "darkgreen", Stroke = "darkorange"),
                                      Direction = c(Up = "red", Down = "blue")),
             annotation_names_col = F, annotation_names_row = F,
             annotation_col = col.mat, annotation_row = row.mat)
png("02.DEG/06.Ferroptosis.Heatmap.png", width = 8, height = 8, units = "in", res = 300)
print(p)
dev.off()
pdf("02.DEG/06.Ferroptosis.Heatmap.pdf", width = 8, height = 8)
print(p)
dev.off()
