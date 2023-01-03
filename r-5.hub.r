rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(reshape2)
library(org.Hs.eg.db)

df.deg = read.delim2("02.DEG/01.DEG.xls", row.names = 1)
hub.gene = read.delim2("hub.gene", header = F)[[1]]
hub.gene = AnnotationDbi::select(org.Hs.eg.db, hub.gene, "ENTREZID", "SYMBOL")
df.deg = df.deg[hub.gene$ENTREZID,]
df.deg = cbind(GeneID = rownames(df.deg), df.deg)
write.table(df.deg, "05.HubGene/03.HubGene.xls", sep = "\t", col.names = T, row.names = F)

rm(list = ls())
df = read.delim2("degree")
hub.gene = read.delim2("hub.gene", header = F)[[1]]
df = subset(df, gene %in% hub.gene)
ggplot(df, aes(x = degree, y = reorder(gene, degree), fill = gene)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = degree + 0.2, label = degree)) +
  scale_fill_aaas() +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,3)) +
  ylab("Hub Gene") + xlab("Degree") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("05.HubGene/04.Degree.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("05.HubGene/04.Degree.pdf", width = 6, height = 4, units = "in", dpi = 300)

rm(list = ls())
df.exp = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% lc.tableToNum()
hub.genes = read.delim2("05.HubGene/03.HubGene.xls")
df.exp = df.exp[as.character(hub.genes$GeneID),]
df.exp$gene = hub.genes$Gene
rownames(df.exp) = hub.genes$Gene
df.exp.plot = melt(df.exp, id.vars = "gene", variable.name = "sample", value.name = "value")

df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
df.exp.plot$group = df.pheno$group[match(df.exp.plot$sample, df.pheno$geo_accession)]
df.exp.plot$group = factor(df.exp.plot$group)
df.diff = read.delim2("02.DEG/01.DEG.xls", row.names = 1)
df.diff = df.diff[as.character(hub.genes$GeneID),]
df.diff = lc.tableToNum(df.diff)
df.diff$p = df.diff$adj.P.Val %>% signif(3) %>% paste0("p = ", .)
df.diff$label = paste0(df.diff$Gene,"\n",df.diff$p)

df.exp.plot = lc.tableToNum(df.exp.plot)
df.exp.plot$label = df.diff$label[match(df.exp.plot$gene, df.diff$Gene)]
df.exp.plot$group = factor(df.exp.plot$group)

library(ggplot2)
ggplot() +
  geom_boxplot(data = df.exp.plot, aes(x = group, y = value, color = group)) + 
  facet_wrap(label ~ ., scales = "free", ncol = 2) +
  theme_bw() + 
  xlab(NULL) + 
  ylab("Expression") + 
  scale_color_manual(values = c("darkgreen", "darkorange"), guide = guide_legend(ncol = 1), labels = c("Low","High")) +
  theme(text = element_text(size = 16), legend.title = element_blank(), legend.position = "none")
ggsave("05.HubGene/05.Boxplot.png", width = 5, height = 5, units = "in", dpi = 300)
ggsave("05.HubGene/05.Boxplot.pdf", width = 5, height = 5, units = "in", dpi = 300)


library(pROC)
res.roc = list()
df.exp = df.exp[df.pheno$geo_accession] %>% t %>% as.data.frame() %>% lc.tableToNum()
df.pheno$group = factor(df.pheno$group)
df.exp = apply(df.exp, 2, as.numeric)

apply(df.exp, 2, function(x){
  res = roc(df.pheno$group, x)
  res.roc[[length(res.roc) + 1]] <<- res
})
names(res.roc) = colnames(df.exp)

res.auc = sapply(res.roc, auc) %>% data.frame
res.auc = cbind(Gene = rownames(res.auc), AUC = res.auc[[1]]) %>% as.data.frame %>% lc.tableToNum
res.auc$AUC = round(res.auc$AUC, 3)
res.auc.plot = res.auc
res.auc.plot = res.auc.plot[order(res.auc.plot$AUC, decreasing = T),]

library(patchwork)
library(gridExtra)
library(ggsci)
p1 = ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + scale_color_aaas() +
  annotation_custom(grob = tableGrob(res.auc.plot, rows = NULL), xmin = 0.7, xmax = 0.9, ymin = 0, ymax = 0.65) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"), 
        legend.background = element_rect(fill = "transparent"),
        legend.position = "bottom", legend.text = element_text(size = 14))
ggsave("05.HubGene/06.ROC.png", width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("05.HubGene/06.ROC.pdf", width = 6, height = 6, dpi = 300, units = "in", bg = "white")

library(corrplot)
res.cor = cor(df.exp)
png("05.HubGene/07.Correlation.png", width = 8, height = 6, units = "in", bg = "white", res = 300)
col3 <- colorRampPalette(c("blue", "white","red"))
corrplot(res.cor, tl.col = "darkblue", col = col3(100), type = "lower",tl.cex = 1.6, cl.cex = 1.2)
corrplot(res.cor, tl.col = "darkblue", col = "black", type = "lower", number.cex = 1.2, bg = "transparent",
         tl.cex = 1.6, method = "number",diag = FALSE, tl.pos = "n", cl.pos = "n", add = T)
dev.off()
pdf("05.HubGene/07.Correlation.pdf", width = 8, height = 6)
col3 <- colorRampPalette(c("blue", "white","red"))
corrplot(res.cor, tl.col = "darkblue", col = col3(100), type = "lower",tl.cex = 1.6, cl.cex = 1.2)
corrplot(res.cor, tl.col = "darkblue", col = "black", type = "lower", number.cex = 1.2, bg = "transparent",
         tl.cex = 1.6, method = "number",diag = FALSE, tl.pos = "n", cl.pos = "n", add = T)
dev.off()

rm(list = ls())
library(GOSemSim)
hub.genes = read.delim2("hub.gene", header = F)[[1]]
hub.genes = AnnotationDbi::select(org.Hs.eg.db, hub.genes, "ENTREZID", "SYMBOL")
d1 = godata('org.Hs.eg.db', ont="MF", computeIC = F)
d2 = godata('org.Hs.eg.db', ont="BP", computeIC = F)
d3 = godata('org.Hs.eg.db', ont="CC", computeIC = F)

gene.pair = combn(as.character(hub.genes$ENTREZID), 2)

res.sim1 = apply(gene.pair, 2, function(x) geneSim(x[1],x[2],d1)$geneSim)
res.sim1 = rbind(data.frame(gene = t(gene.pair)[,1], sim = res.sim1), 
                 data.frame(gene = t(gene.pair)[,2], sim = res.sim1))
res.sim1$name = hub.genes$SYMBOL[match(res.sim1$gene, hub.genes$ENTREZID)]

res.sim2 = apply(gene.pair, 2, function(x) geneSim(x[1],x[2],d2)$geneSim)
res.sim2 = rbind(data.frame(gene = t(gene.pair)[,1], sim = res.sim2), 
                 data.frame(gene = t(gene.pair)[,2], sim = res.sim2))
res.sim2$name = hub.genes$SYMBOL[match(res.sim2$gene, hub.genes$ENTREZID)]

res.sim3 = apply(gene.pair, 2, function(x) geneSim(x[1],x[2],d3)$geneSim)
res.sim3 = rbind(data.frame(gene = t(gene.pair)[,1], sim = res.sim3), 
                 data.frame(gene = t(gene.pair)[,2], sim = res.sim3))
res.sim3$name = hub.genes$SYMBOL[match(res.sim3$gene, hub.genes$ENTREZID)]

res.sim = rbind(res.sim1, res.sim2, res.sim3)

gene.order = aggregate(.~name, FUN = median, data = res.sim[-1])
gene.order = gene.order[order(gene.order$sim, decreasing = F),]
res.sim$name = factor(res.sim$name, levels = gene.order$name)

library(ggpubr)
ggplot(res.sim, aes(y = name, x = sim, color = name)) +
  geom_boxplot(outlier.color = "transparent") + 
  geom_vline(xintercept = 0.75, color = "red", linetype = 2, size = 1) +
  theme_pubclean() + 
  ylab("Gene") + xlab(NULL) +
  theme(axis.line.x.bottom = element_line(size = 1), legend.position = "none", 
        axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 18))
ggsave("05.HubGene/08.GOSemSim.png", width = 6, height = 3, units = "in", dpi = 300, bg = "white")
ggsave("05.HubGene/08.GOSemSim.pdf", width = 6, height = 3, units = "in", dpi = 300, bg = "white")

rm(list = ls())
df.im = read.delim2("03.Immune/01.ImmuneProfile.xls", row.names = 1)
top.genes = read.delim2("05.HubGene/03.HubGene.xls")

df.exp = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% lc.tableToNum()
df.exp = df.exp[as.character(top.genes$GeneID),]
rownames(df.exp) = top.genes$Gene
df.exp = df.exp[colnames(df.im)]

df.exp = t(df.exp) %>% as.data.frame()
df.im = t(df.im) %>% as.data.frame()

com = rep(1:ncol(df.im), ncol(df.exp)) %>% matrix(ncol = 1)
com = cbind(rep(1:ncol(df.exp), each = ncol(df.im)), com)

res.cor = apply(com, 1, function(x){
  v1 = df.exp[[x[1]]] %>% as.numeric
  v2 = df.im[[x[2]]] %>% as.numeric
  a = cor.test(v1, v2)
  a = a[c("estimate","p.value")] %>% unlist
  a = c(colnames(df.exp)[x[1]], colnames(df.im)[x[2]], a)
}) %>% t %>% as.data.frame
res.cor = lc.tableToNum(res.cor)
colnames(res.cor) = c("gene", "im", "cor", "pv")
res.cor$pl = signif(res.cor$pv, 2) %>% paste0("p = ",.)

ggplot(data = res.cor, mapping = aes(x = cor, y = im)) + 
  geom_segment(aes(y = im, yend = im, x = 0, xend = cor), color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(aes(size = abs(cor), fill = pv), pch = 21, color = "white") +
  geom_text(aes(y = im, x = 0.95, label = signif(pv,3), color = pv < 0.05), hjust = 0) +
  scale_fill_gradient(low = "red", high = "green", 
                      guide = guide_colourbar(title = "p.value",barheight = unit(5, "in")), breaks = c(0.05, 0.1, 0.3, 0.5, 0.7,0.9,1)) +
  scale_color_manual(values = c("darkblue", "darkred"), guide = guide_none()) +
  scale_size(guide = guide_none()) +
  xlim(-1,1.2) +
  xlab("Correlation Coefficient") + 
  ylab(NULL) +
  facet_wrap(.~gene, ncol = 2) +
  theme_pubr(base_size = 14) + 
  theme(panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        legend.position = "right")
ggsave("05.HubGene/09.ImmuneCorrelation.png", width = 13, height = 10, units = "in", dpi = 300)
ggsave("05.HubGene/09.ImmuneCorrelation.pdf", width = 13, height = 10, units = "in", dpi = 300)

rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(ggrepel)
library(enrichplot)
library(aplot)

gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = 'PuOr') +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank())
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_brewer(palette = 'PuOr')
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

df.m = msigdbr()
df.kegg = subset(df.m, gs_subcat == "CP:KEGG")[c(3,5)]
df.exp = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% t %>% as.data.frame %>% lc.tableToNum()
hub.genes = read.delim2("05.HubGene/03.HubGene.xls")

hub = hub.genes$GeneID[1] %>% as.character()
hub.exp = df.exp[[hub]]
hub.cor = cor(df.exp, hub.exp) %>% as.data.frame %>% na.omit
hub.coreff = hub.cor[[1]]
names(hub.coreff) = rownames(hub.cor)
hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.05, seed = 1, pAdjustMethod = "BH", eps = 0)
top.kegg = res.kegg@result
top.kegg = c(top.kegg[order(top.kegg$NES, decreasing = T),]$Description[1:5],
             top.kegg[order(top.kegg$NES, decreasing = F),]$Description[5:1])
p = gsea.plot(res.kegg, top.kegg, hub.genes$Gene[1])
p
ggsave("06.GSEA/01.DDIT3.png", p, width = 12, height = 6, units = "in", limitsize = 300)
ggsave("06.GSEA/01.DDIT3.pdf", p, width = 12, height = 6, units = "in", limitsize = 300)

hub = hub.genes$GeneID[2] %>% as.character()
hub.exp = df.exp[[hub]]
hub.cor = cor(df.exp, hub.exp) %>% as.data.frame %>% na.omit
hub.coreff = hub.cor[[1]]
names(hub.coreff) = rownames(hub.cor)
hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.05, seed = 1, pAdjustMethod = "BH", eps = 0)
top.kegg = res.kegg@result
top.kegg = c(top.kegg[order(top.kegg$NES, decreasing = T),]$Description[1:5],
             top.kegg[order(top.kegg$NES, decreasing = F),]$Description[5:1])
p = gsea.plot(res.kegg, top.kegg, hub.genes$Gene[2])
p
ggsave("06.GSEA/02.DUSP1.png", p, width = 12, height = 6, units = "in", limitsize = 300)
ggsave("06.GSEA/02.DUSP1.pdf", p, width = 12, height = 6, units = "in", limitsize = 300)

hub = hub.genes$GeneID[3] %>% as.character()
hub.exp = df.exp[[hub]]
hub.cor = cor(df.exp, hub.exp) %>% as.data.frame %>% na.omit
hub.coreff = hub.cor[[1]]
names(hub.coreff) = rownames(hub.cor)
hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.05, seed = 1, pAdjustMethod = "BH", eps = 0)
top.kegg = res.kegg@result
top.kegg = c(top.kegg[order(top.kegg$NES, decreasing = T),]$Description[1:5],
             top.kegg[order(top.kegg$NES, decreasing = F),]$Description[5:1])
p = gsea.plot(res.kegg, top.kegg, hub.genes$Gene[3])
p
ggsave("06.GSEA/03.NFIL3.png", p, width = 12, height = 6, units = "in", limitsize = 300)
ggsave("06.GSEA/03.NFIL3.pdf", p, width = 12, height = 6, units = "in", limitsize = 300)

hub = hub.genes$GeneID[4] %>% as.character()
hub.exp = df.exp[[hub]]
hub.cor = cor(df.exp, hub.exp) %>% as.data.frame %>% na.omit
hub.coreff = hub.cor[[1]]
names(hub.coreff) = rownames(hub.cor)
hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.05, seed = 1, pAdjustMethod = "BH", eps = 0)
top.kegg = res.kegg@result
top.kegg = c(top.kegg[order(top.kegg$NES, decreasing = T),]$Description[1:5],
             top.kegg[order(top.kegg$NES, decreasing = F),]$Description[5:1])
p = gsea.plot(res.kegg, top.kegg, hub.genes$Gene[4])
p
ggsave("06.GSEA/04.FOS.png", p, width = 12, height = 6, units = "in", limitsize = 300)
ggsave("06.GSEA/04.FOS.pdf", p, width = 12, height = 6, units = "in", limitsize = 300)
