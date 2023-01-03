rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(WGCNA)
library(org.Hs.eg.db)

df.exp = read.delim2("01.Data/01.Expression.xls", row.names = 1) %>% lc.tableToNum()
df.pheno = read.delim2("01.Data/02.MetaInfo.xls")
df.exp = df.exp[df.pheno$geo_accession]
df.pheno$group = ifelse(str_detect(df.pheno$title, "Control"), "Control", "Stroke")
df.pheno$group = factor(df.pheno$group)
df.pheno = df.pheno[order(df.pheno$group, decreasing = F),]
df.exp = df.exp[df.pheno$geo_accession] %>% t %>% as.data.frame %>% lc.tableToNum 
df.exp = df.exp[which(colSums(df.exp)>0)]
df.annot = AnnotationDbi::select(org.Hs.eg.db, colnames(df.exp), columns = "GENETYPE", keytype = "ENTREZID")
df.annot = df.annot[!duplicated(df.annot$ENTREZID),]
identical(df.annot$ENTREZID, colnames(df.exp))
df.exp = df.exp[which(df.annot$GENETYPE == "protein-coding")]

goodSamplesGenes(df.exp, verbose = 3)
collectGarbage()

sampleTree = hclust(dist(df.exp), method = "average")
par(mfrow = c(1,1))
png("04.WGCNA/01.OutlierDetection.png", width = 12, height = 7, res = 1200, units = "in")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex = 0.8, labels = F)
abline(h = 70, col = "red", lty = 2)
dev.off()

pdf("04.WGCNA/01.OutlierDetection.pdf", width = 12, height = 7)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", labels = F)
abline(h = 70, col = "red", lty = 2)
dev.off()

out.cut = cutree(sampleTree, k = 2)
table(out.cut)
df.exp = df.exp[names(out.cut[out.cut==1]),]
df.pheno = df.pheno[match(names(out.cut[out.cut==1]),df.pheno$geo_accession),]

sampleTree2 = hclust(dist(df.exp), method = "average")
identical(df.pheno$geo_accession, rownames(df.exp))
traitColors = numbers2colors(df.pheno$group %>% as.numeric, signed = FALSE, colors = c("darkgreen","darkorange"))

png("04.WGCNA/02.SampleCluster.png", width = 12, height = 5, res = 300, units = "in")
plotDendroAndColors(sampleTree2, traitColors, dendroLabels = F,
                    groupLabels = "Group",
                    main = "Sample dendrogram and trait heatmap")
dev.off()

pdf("04.WGCNA/02.SampleCluster.pdf", width = 12, height = 5)
plotDendroAndColors(sampleTree2, traitColors, dendroLabels = F,
                    groupLabels = "Group",
                    main = "Sample dendrogram and trait heatmap")
dev.off()

powers = 1:20
sft = pickSoftThreshold(df.exp, powerVector = powers, verbose = 5, blockSize = ncol(df.exp))

png("04.WGCNA/03.SoftThred.png", width = 10, height = 6, units = "in", res = 1200)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

pdf("04.WGCNA/03.SoftThred.pdf", width = 10, height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

cor = WGCNA::cor
net = blockwiseModules(
  df.exp, power = 4, maxBlockSize = ncol(df.exp), networkType = "signed",
  minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 4, corType = "bicor",
  pamRespectsDendro = FALSE, reassignThreshold = 0, 
  robustY = FALSE, numericLabels = TRUE, verbose = 3)
table(net$colors)

mergedColors = labels2colors(net$colors)
png("04.WGCNA/04.Modules.png", width = 12, height = 6, res = 300, units = "in")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = F)
dev.off()
pdf("04.WGCNA/04.Modules.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = F)
dev.off()

df.im = read.delim2("03.Immune/01.ImmuneProfile.xls", row.names = 1) %>% lc.tableToNum()
df.diff = read.delim2("03.Immune/02.ImmuneDiff.xls") %>% lc.tableToNum()
df.diff = subset(df.diff, p.adj < 0.05)
df.im = df.im[df.diff$Immune,]
df.im = df.im[df.pheno$geo_accession]
design = t(df.im) %>% as.data.frame()
identical(rownames(design), rownames(df.exp))

moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(df.exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs = MEs[,-ncol(MEs)]
design = design[rownames(MEs),]

moduleTraitCor = cor(MEs, design)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(df.exp)) %>% signif(3)
moduleTraitPvalue = p.adjust(moduleTraitPvalue)

sizeGrWindow(15,10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

png("04.WGCNA/05.ModuleTraitCorrelation.png", width = 8, height = 10, units = "in", res = 300)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs) %>% str_remove("ME"),
               colorLabels = FALSE,
               keepLegendSpace = F,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               xLabelsAngle = 30,
               setStdMargins = T,
               x.adj.lab.y = 1,xLabelsAdj = 1,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

pdf("04.WGCNA/05.ModuleTraitCorrelation.pdf", width = 8, height = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs) %>% str_remove("ME"),
               colorLabels = FALSE,
               keepLegendSpace = F,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               xLabelsAngle = 30,
               setStdMargins = T,
               x.adj.lab.y = 1,xLabelsAdj = 1,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

genes.color = data.frame(names(net$colors),labels2colors(net$colors))
colnames(genes.color) = c("EntrezID", "Module")
table(genes.color$Module)

library(org.Hs.eg.db)
genes.color$Gene = AnnotationDbi::select(org.Hs.eg.db, genes.color$EntrezID, c("SYMBOL"), "ENTREZID")$SYMBOL
genes.color = genes.color[c(1,3,2)]
write.table(genes.color, "04.WGCNA/06.GeneModule.xls", sep = "\t", quote = F, col.names = T, row.names = F)

MEs = moduleEigengenes(df.exp, moduleColors)$eigengenes
MET = orderMEs(MEs)
png("04.WGCNA/07.EigengeneDendrogram.png", width = 10, height = 8, res = 300, units = "in", bg = "white")
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()
pdf("04.WGCNA/07.EigengeneDendrogram.pdf", width = 10, height = 8)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

png("04.WGCNA/08.EigengeneAdjacencyHeatmap.png", width = 8, height = 8, res = 300, units = "in", bg = "white")
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = F, xLabelsAngle = 90)
dev.off()
pdf("04.WGCNA/08.EigengeneAdjacencyHeatmap.pdf", width = 8, height = 8)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = F, xLabelsAngle = 90)
dev.off()

gene.sel = genes.color$EntrezID[genes.color$Module %in% c("brown","yellow","green","pink","salmon","red")]
gene.diff = read.delim2("02.DEG/05.Hypoxia.DEG.xls") %>% lc.tableToNum()
gene.diff = gene.diff$GeneID %>% as.character()
p.venn = lc.vennFromList(list(DEG = gene.sel, `Hypoxia\nDEG` = gene.diff), label.dist = -20)
png("04.WGCNA/09.Venn.png", width = 7, height = 7, res = 300, units = "in")
grid.newpage();grid.draw(p.venn)
dev.off()
pdf("04.WGCNA/09.Venn.pdf", width = 7, height = 7)
grid.newpage();grid.draw(p.venn)
dev.off()

gene.inter = intersect(gene.sel, gene.diff)
gene.inter = data.frame(GeneID = gene.inter, 
                        Gene = AnnotationDbi::select(org.Hs.eg.db, gene.inter %>% as.character(), "SYMBOL", "ENTREZID")$SYMBOL)
df.diff = read.delim2("02.DEG/01.DEG.xls") %>% lc.tableToNum()
gene.inter = subset(df.diff, GeneID %in% gene.inter$GeneID)
write.table(gene.inter, "04.WGCNA/10.WGCNA_DEG.xls", col.names = T, row.names = F, quote = F, sep = "\t")

rm(list = ls())
library(clusterProfiler)
library(ggplot2)
library(GOplot)
gene.inter = read.delim2("04.WGCNA/10.WGCNA_DEG.xls")
res.kegg = enrichKEGG(as.character(gene.inter$GeneID) %>% as.character, organism = "hsa", keyType = "ncbi-geneid")
res.go = enrichGO(as.character(gene.inter$GeneID), org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL")
res.kegg = setReadable(res.kegg, org.Hs.eg.db, "ENTREZID")
res.go = setReadable(res.go, org.Hs.eg.db, "ENTREZID")
df.kegg = res.kegg@result
df.go = res.go@result

df.deg = gene.inter[2:3] %>% lc.tableToNum()
colnames(df.deg) = c("ID","logFC")
df.go = df.go[c(1,2,3,9,7)]
colnames(df.go) = c("Category","ID","Term","genes","adj_pval")
df.go$genes = str_replace_all(df.go$genes, fixed("/"), ",")
circ = circle_dat(df.go, df.deg)
p = GOCircle(circ, label.size = 4.5)
ggsave("04.WGCNA/11.GOCircle.png", p, width = 11, height = 7, units = "in", dpi = 300)
ggsave("04.WGCNA/11.GOCircle.pdf", p, width = 11, height = 7, units = "in", dpi = 300)

process = df.go$Term[1:10]
chord = chord_dat(data=circ, process = process)
chord = chord[order(abs(chord[,ncol(chord)]), decreasing = T),]
p = GOChord(chord, space = 0.015, gene.order = 'logFC', lfc.min = -5, lfc.max = 5,
            gene.space = 0.25, gene.size = 3.2, process.label = 12,)
p$guides$size$ncol = 2
p = p + theme(aspect.ratio = 1)
p
ggsave("04.WGCNA/12.GOChord.png", p, width = 10, height = 10, units = "in", dpi = 300)
ggsave("04.WGCNA/12.GOChord.pdf", p, width = 10, height = 10, units = "in", dpi = 300)
write.table(df.go, "04.WGCNA/13.GO.xls", quote = F, sep = "\t", row.names = F)

df.kegg = df.kegg[c(1,2,8,6)]
df.kegg = cbind(Category = "KEGG", df.kegg)
colnames(df.kegg) = c("Category","ID","Term","genes","adj_pval")
df.kegg = subset(df.kegg, adj_pval < 0.05)
df.kegg$genes = str_replace_all(df.kegg$genes, fixed("/"), ",")
circ = circle_dat(df.kegg, df.deg)
cnetplot(res.kegg) + theme(legend.position = "none")
ggsave("04.WGCNA/14.KEGG.png",width = 10, height = 10, units = "in", dpi = 300)
ggsave("04.WGCNA/14.KEGG.pdf",width = 10, height = 10, units = "in", dpi = 300)

df.kegg = res.kegg@result
write.table(df.kegg, "04.WGCNA/15.KEGG.xls", quote = F, sep = "\t", row.names = F)






