

####1.1.mRNA####

library(rjson)
library(tidyverse)

setwd("")

json <- jsonlite::fromJSON("metadata.cart.2024-10-14.json")
#View(json)

sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)

count_file <- list.files('gdc_download_20241014_141123.497189/',
                         pattern = '*.tsv',recursive = TRUE)
#count_file[1:10]

count_file_name <- strsplit(count_file,split='/')

count_file_name <- sapply(count_file_name,function(x){x[2]})

matrix = data.frame(matrix(nrow=60660,ncol=0))

for (i in 1:length(count_file)){
  path = paste0('gdc_download_20241014_141123.497189//',count_file[i])   
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[6]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

path = paste0('gdc_download_20241014_141123.497189//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
#gene_name[1:10]
matrix0 <- cbind(gene_name,matrix)
#获取基因类型
gene_type <- data[-c(1:6),2]
#gene_type[1:10]
matrix0 <- cbind(gene_type,matrix0)
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)
table(gene_name)

matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
#table(gene_type)

rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
matrix1 = data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1) = gsub('[.]', '-', colnames(matrix1))

write.table(matrix1,'mRNA.txt', sep="\t", quote=F, row.names = F)

####1.2.lncRNA####

library(rjson)
library(tidyverse)

setwd("")

json <- jsonlite::fromJSON("metadata.cart.2024-10-14.json")

sample_id <- sapply(json$associated_entities,function(x){x[,1]})

file_sample <- data.frame(sample_id,file_name=json$file_name)  

count_file <- list.files('gdc_download_20241014_141123.497189/',
                         pattern = '*.tsv',recursive = TRUE)

count_file_name <- strsplit(count_file,split='/')

count_file_name <- sapply(count_file_name,function(x){x[2]})

matrix = data.frame(matrix(nrow=60660,ncol=0))

for (i in 1:length(count_file)){
  path = paste0('gdc_download_20241014_141123.497189//',count_file[i])  
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[6]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

path = paste0('gdc_download_20241014_141123.497189//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <- cbind(gene_name,matrix)
gene_type <- data[-c(1:6),2]
matrix0 <- cbind(gene_type,matrix0)

matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)
table(gene_name)

matrix0 <- subset(x = matrix0, gene_type == "lncRNA")

rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
matrix1 = data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1) = gsub('[.]', '-', colnames(matrix1))

write.table(matrix1,'lncRNA.txt', sep="\t", quote=F, row.names = F)

####1.3Telomerase Gene Expression####

library(limma)         
expFile="mRNA.txt"      
geneFile="Telomerase gene.txt"      

setwd("")

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

outTab=rbind(ID=colnames(geneExp),geneExp)
write.table(outTab, file="TelomeraseExp.txt", sep="\t", quote=F, col.names=F)

####1.4.WGCNA####

library("WGCNA")            
library("limma")            
expFile="TelomeraseExp.txt"     
normalCount=50              
tumorCount=374              
setwd("")      


rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=log2(data+1)
data=data[apply(data,1,sd)>0.5,]
datExpr0=t(data)


gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

abline(h = 20000, col = "red")
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]


traitData=data.frame(Normal=c(rep(1,normalCount),rep(0,tumorCount)),
                     Tumor=c(rep(0,normalCount),rep(1,tumorCount)))
row.names(traitData)=colnames(data)
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(traitData)
sameSample=intersect(fpkmSamples,traitSamples)
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]

sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2_sample_heatmap.pdf",width=15,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

enableWGCNAThreads()   
powers = c(1:20)       
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9   

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") 

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



sft 
softPower =sft$powerEstimate 
adjacency = adjacency(datExpr0, power = softPower)
softPower


TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



minModuleSize = 50     
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE, 
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="5_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="6_Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25                              
abline(h=MEDissThres, col = "red")
dev.off()



merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors(TCGA)")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("yellow", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs



nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="8_Module_trait.pdf",width=6,height=6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships(TCGA)"))
dev.off()



modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste("9_", trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}



probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)


for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("TCGA_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

####1.5DEGs####

library(limma)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library("BiocParallel")
register(MulticoreParam(5))

expFile="TelomeraseExp.txt"      
logFCfilter=1          
padjFilter=0.001          
setwd("")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data=round(data,0)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
normalData=data[,group==1]
tumorData=data[,group==0]
data=cbind(normalData, tumorData)
normalNum=length(group[group==1])      
tumorNum=length(group[group==0])       

group=c(rep("normal",normalNum), rep("tumor",tumorNum))
coldata=data.frame(condition=group, row.names=colnames(data))
dds=DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)
dds=DESeq(dds, parallel=TRUE)
diff=results(dds, parallel=TRUE)
diff=as.data.frame(diff)


diff=diff[is.na(diff$padj)==FALSE,]
diff=diff[order(diff$pvalue),]
diffOut=rbind(id=colnames(diff), diff)

write.table(diffOut, file="all.txt", sep="\t", quote=F, col.names=F)
diffSig=diff[(diff$padj<padjFilter & (diff$log2FoldChange>logFCfilter | diff$log2FoldChange<(-logFCfilter))),]

diffSigOut=rbind(id=colnames(diffSig), diffSig)
write.table(diffSigOut, file="diff.txt",sep="\t",quote=F,col.names=F)


newData=counts(dds, normalized=TRUE)
newData=log2(newData+1)
normalizeExp=rbind(id=colnames(newData), newData)
write.table(normalizeExp,file="normalize.txt",sep="\t",quote=F,col.names=F)


geneNum=30      
diffUp=diffSig[diffSig$log2FoldChange>0,]
diffDown=diffSig[diffSig$log2FoldChange<0,]
geneUp=row.names(diffUp)
geneDown=row.names(diffDown)
if(nrow(diffUp)>geneNum){geneUp=row.names(diffUp)[1:geneNum]}
if(nrow(diffDown)>geneNum){geneDown=row.names(diffDown)[1:geneNum]}
hmExp=newData[c(geneUp,geneDown),]
Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))
names(Type)=colnames(newData)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()


Significant=ifelse((diff$padj<padjFilter & abs(diff$log2FoldChange)>logFCfilter), ifelse(diff$log2FoldChange>logFCfilter,"Up","Down"), "Not")

p = ggplot(diff, aes(log2FoldChange, -log10(padj)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("blue", "grey", "red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf(file="vol.pdf",width=6,height=5)
print(p)
dev.off()

####1.6Venn1####


library(ggvenn)

diffFile="diff.txt"                  
moduleFile="TCGA_blue.txt"     
setwd("")     
geneList=list()


rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              
geneNames=gsub("^ | $","",geneNames)     
uniqGene=unique(geneNames)               
geneList[["DEG"]]=uniqGene


rt=read.table(moduleFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)   
uniqGene=unique(geneNames)               
geneList[["WGCNA"]]=uniqGene


pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C","#1E90FF"),
       set_name_color =c("#E41A1C","#1E90FF"),
       set_name_size=6, text_size=4.5)
dev.off()


interGenes=Reduce(intersect, geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)

####2.1Venn2####

library(ggvenn)     
scoreFile="score.csv"    
setwd("")     

rt=read.csv(scoreFile, header=T, sep=",", check.names=F, row.names=1)
colnames(rt)=c(colnames(rt)[2:ncol(rt)], "N")
scoreType=c("Degree", "EPC", "MCC", "MNC")

geneList=c()
for(i in scoreType){
  data=rt[order(rt[,i],decreasing=T),]
  hubGene=row.names(data)[1:20]
  geneList[[i]]=hubGene
}

pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C","#1E90FF","#4DAF4A","#984EA3"),
       set_name_color = c("#E41A1C","#1E90FF","#4DAF4A","#984EA3"),
       set_name_size=6, text_size=4.5)
dev.off()

interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="hubGenes.txt", sep="\t", quote=F, row.names=F, col.names=F)

####2.2GO####

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

wdPath <- ""
if (!dir.exists(wdPath)) dir.create(wdPath, recursive = TRUE)
setwd(wdPath)

updateProgress <- function(pb, value, message) {
  setTxtProgressBar(pb, value)
  cat(message, "\n")
}
totalSteps <- 9
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)

pvalThreshold <- 0.05
padjThreshold <- 0.05
colorParameter <- if (padjThreshold > 0.05) "pvalue" else "p.adjust"
updateProgress(pb, 1, "第1步：设置参数和工作环境完成")

geneFile <- "hubGenes.csv"
if (!file.exists(geneFile)) stop("错误：基因文件 不存在")
geneData <- read.table(geneFile, header = FALSE, sep = ",", check.names = FALSE)
updateProgress(pb, 2, "第2步：读取基因数据完成")

geneSymbols <- unique(as.vector(geneData[, 1]))
if (length(geneSymbols) == 0) stop("错误：未找到有效的基因符号")

entrezMapping <- mget(geneSymbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezMapping)
validGenes <- entrezIDs[!is.na(entrezIDs) & entrezIDs!="NA"]
if (length(validGenes) == 0) stop("错误：无有效的 Entrez 基因 ID")
updateProgress(pb, 3, "第3步：基因符号转换为 EntrezID完成")

goAnalysis <- enrichGO(
  gene = validGenes,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 1, qvalueCutoff = 1,
  ont = "all",
  readable = TRUE
)
goResult <- as.data.frame(goAnalysis)
if (nrow(goResult)==0) stop("警告：未检测到任何富集结果")
filteredGO <- goResult[goResult$pvalue < pvalThreshold & goResult$p.adjust < padjThreshold, ]
updateProgress(pb, 4, "第4步：GO 富集分析及结果过滤完成")

outputFile <- "GO_results.txt"
write.table(filteredGO, file = outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
updateProgress(pb, 5, "第5步：富集结果写入文件完成")

pdf("GO_barplot.pdf", width=7, height=6)
barPlot <- barplot(
  goAnalysis,
  drop=TRUE, showCategory=10, label_format=50,
  split="ONTOLOGY", color=colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_fill_gradientn(colors = c("#FF6666", "#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(barPlot)
dev.off()

pdf("GO_bubble.pdf", width=8, height=10)
bubblePlot <- dotplot(
  goAnalysis,
  showCategory = 10,
  orderBy = "GeneRatio",
  label_format = 50,
  split = "ONTOLOGY",
  color = colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_color_gradientn(colors = c("#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(bubblePlot)
dev.off()
updateProgress(pb, 6, "第6步：柱状图和气泡图生成完成")

topGO <- filteredGO %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
pdf("GO_grouped_barplot.pdf", width=11, height=8)
groupBarPlot <- ggbarplot(
  topGO,
  x="Description", y="Count", fill="ONTOLOGY", color="white",
  xlab="", palette="aaas",
  legend="right", sort.val="desc", sort.by.groups=TRUE,
  position=position_dodge(0.9)
) +
  rotate_x_text(75) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=10, color="black")) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  geom_text(
    aes(label=Count),
    position=position_dodge(0.9), vjust=-0.3, size=3
  )
print(groupBarPlot)
dev.off()
updateProgress(pb, 7, "第7步：分组条形图生成完成")

pdf("GO_chord_diagram.pdf", width=12, height=12)
go <- read.delim("GO_results.txt", header=TRUE, stringsAsFactors=FALSE)
top_terms <- go %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()
insert_linebreak <- function(text, line_length=35) {
  if(nchar(text) <= line_length) return(text)
  paste(strwrap(text, width=line_length), collapse = "\n")
}
top_terms$Description_new <- sapply(top_terms$Description, insert_linebreak)
mat <- table(
  factor(top_terms$ONTOLOGY, levels=c("BP","CC","MF")),
  factor(top_terms$Description_new, levels=unique(top_terms$Description_new))
)
n_ont <- 3
n_term <- ncol(mat)
gap.deg <- c(rep(1, n_ont-1), 10, rep(1, n_term-1), 10)
grid.col <- c(BP="#E69F00", CC="#56B4E9", MF="#009E73",
              setNames(rep("#BBBBBB", n_term), colnames(mat)))
circos.clear()
circos.par(gap.degree = gap.deg, start.degree = 90)
chordDiagram(
  mat,
  grid.col = grid.col,
  transparency = 0.4,
  annotationTrack = c("", "grid"),
  preAllocateTracks = list(track.height = 0.2)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
              facing="clockwise", niceFacing=TRUE,
              adj=c(0,0.5), cex=0.60)
}, bg.border=NA)
title("GO Ontology")
circos.clear()
dev.off()
updateProgress(pb, 8, "第8步：chord弦图生成完成")


library(dplyr)
library(ggplot2)
library(RColorBrewer)

filteredGO <- read.table("GO_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

bar_colors <- brewer.pal(7, "YlOrRd")
p_bar <- ggplot(topGO, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = bar_colors) +
  labs(x = "GO Term", y = "-log10(Adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11)
  )
ggsave("GO_barplot_custom.pdf", p_bar, width = 12, height = 12)

dot_colors <- brewer.pal(7, "Spectral")
p_dot <- ggplot(topGO, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_color_gradientn(colors = dot_colors) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11)
  )
ggsave("GO_dotplot_custom.pdf", p_dot, width = 12, height = 12)

filteredGO <- read.table("GO_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

bar_colors <- brewer.pal(7, "YlOrRd")
p_bar <- ggplot(topGO, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = bar_colors) +
  labs(x = "GO Term", y = "-log10(adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold") # 通路名字加粗
  )
ggsave("GO_barplot_custom.pdf", p_bar, width = 12, height = 12)

dot_colors <- brewer.pal(7, "Spectral")
p_dot <- ggplot(topGO, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_color_gradientn(colors = dot_colors) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold")
  )
ggsave("GO_dotplot_custom.pdf", p_dot, width = 10, height = 10)

cat("所有LASSO系数已保存为 LASSO_Coefficients_All.csv\n")

close(pb)
cat('全部分析及绘图完成！\n')

####2.3KEGG####

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(RColorBrewer)
library(circlize)

pvalueFilter <- 0.05
adjPvalFilter <- 1
colorSel <- "p.adjust"
if (adjPvalFilter > 1) {
  colorSel <- "pvalue"
}

setwd("")

rt <- read.table("hubGenes.csv", header = TRUE, sep = ",", check.names = FALSE)
colnames(rt)[1] <- "genes"  # Rename the first column to 'genes'

genes <- unique(as.vector(rt$genes))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt$entrezIDs <- entrezIDs
rt <- rt[rt$entrezIDs != "NA", ]
gene <- rt$entrezIDs  # Filtered Entrez IDs


kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)

kkx <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
KEGG <- as.data.frame(kkx)

KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < adjPvalFilter), ]

write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Define the number of pathways to display
showNum <- 30
if (nrow(KEGG) < showNum) {
  showNum <- nrow(KEGG)
}

topKEGG <- KEGG[order(KEGG$p.adjust), ][1:showNum, ]
topKEGG$Description <- factor(topKEGG$Description, levels = rev(topKEGG$Description))

p_bar <- ggplot(topKEGG, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradientn(colors = brewer.pal(7, "YlOrRd")) +
  theme_minimal(base_size = 14) +
  labs(x = "KEGG Pathway", y = "-log10(Adjusted p-value)", title = "KEGG Pathway Enrichment Analysis") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"))
ggsave("barplot_custom.pdf", p_bar, width = 8, height = 7)


if(all(grepl("/", topKEGG$GeneRatio))){
  topKEGG$GeneRatio_num <- sapply(topKEGG$GeneRatio, function(x) {
    parts <- unlist(strsplit(x, "/"))
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
} else {
  topKEGG$GeneRatio_num <- as.numeric(topKEGG$GeneRatio)
}

p_dot <- ggplot(topKEGG, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(7, "Spectral")) +
  theme_classic(base_size = 14) +
  labs(x = "Gene Ratio", y = "KEGG Pathway", title = "KEGG Dotplot") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"))
ggsave("dotplot_custom.pdf", p_dot, width = 8, height = 7)

pdf(file = "cnetplot_custom.pdf", width = 9, height = 5.25)
kkx <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
p_cnet <- cnetplot(kkx, circular = TRUE, showCategory = 5, colorEdge = TRUE)
print(p_cnet + 
        ggtitle("Gene-Pathway Network") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "purple")))
dev.off()

topKEGG_chord <- KEGG[1:min(10, nrow(KEGG)), ]
gene_path_df <- data.frame(Pathway = character(), Gene = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(topKEGG_chord)) {
  pathway <- topKEGG_chord$Description[i]
  genes <- unlist(strsplit(topKEGG_chord$geneID[i], "/"))
  temp_df <- data.frame(Pathway = rep(pathway, length(genes)), Gene = genes, stringsAsFactors = FALSE)
  gene_path_df <- rbind(gene_path_df, temp_df)
}
print(head(gene_path_df))

chord_matrix <- as.matrix(table(gene_path_df$Pathway, gene_path_df$Gene))
rownames(chord_matrix) <- as.character(rownames(chord_matrix))
colnames(chord_matrix) <- as.character(colnames(chord_matrix))
print(dim(chord_matrix))
print(chord_matrix)

allSectors <- union(rownames(chord_matrix), colnames(chord_matrix))
print(allSectors)

nPathways <- length(rownames(chord_matrix))
nGenes <- length(colnames(chord_matrix))
pathway_colors <- brewer.pal(n = min(nPathways, 8), name = "Set2")
if(nPathways > length(pathway_colors)){
  pathway_colors <- rep(pathway_colors, length.out = nPathways)
}
gene_colors <- brewer.pal(n = min(nGenes, 8), name = "Pastel1")
if(nGenes > length(gene_colors)){
  gene_colors <- rep(gene_colors, length.out = nGenes)
}

grid.col <- c(setNames(pathway_colors, rownames(chord_matrix)),
              setNames(gene_colors, colnames(chord_matrix)))
grid.col <- grid.col[allSectors]
print(grid.col)

pdf(file = "chordDiagram_custom.pdf", width = 8, height = 8)
chordDiagram(chord_matrix, grid.col = grid.col, transparency = 0.25, 
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05))
title("Gene-Pathway Chord Diagram")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] - mm_y(5),
              sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)
dev.off()

####2.4Hubgenes Expression####

library(limma)         
expFile="mRNA.txt"      
geneFile="hubGenes.txt"      

setwd("")

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

outTab=rbind(ID=colnames(geneExp),geneExp)
write.table(outTab, file="hubGenesExp.txt", sep="\t", quote=F, col.names=F)

####3.1Cor-expression####

library(limma)
corFilter=0.4            
pvalueFilter=0.001      

setwd("")

rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       
treatNum=length(group[group==0])    
sampleType=c(rep(1,conNum), rep(2,treatNum))

rt1=read.table("hubgenesExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
Telomerase=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
Telomerase=avereps(Telomerase)

group=sapply(strsplit(colnames(Telomerase),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
Telomerase=Telomerase[,group==0]

outTab=data.frame()
for(i in row.names(lncRNA)){
  if(sd(lncRNA[i,])>0.5){
    test=wilcox.test(data[i,] ~ sampleType)
    if(test$p.value<0.05){
      for(j in row.names(Telomerase)){
        x=as.numeric(lncRNA[i,])
        y=as.numeric(Telomerase[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((abs(cor)>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Telomerase=j,lncRNA=i,cor,pvalue))
        }
      }
    }
  }
}

write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

TelomeraseLncRNA=unique(as.vector(outTab[,"lncRNA"]))
TelomeraseLncRNAexp=data[TelomeraseLncRNA,]
TelomeraseLncRNAexp=rbind(ID=colnames(TelomeraseLncRNAexp), TelomeraseLncRNAexp)
write.table(TelomeraseLncRNAexp,file="TelomeraseLncExp.txt",sep="\t",quote=F,col.names=F)

library(dplyr)
library(ggplot2)
library(ggalluvial)

inputFile="corResult.txt"     
setwd("")     

rt=read.table(inputFile, header=T, sep="\t", check.names=F)

mycol=rainbow(length(unique(rt[,"Telomerase"])), s=0.8, v=0.8)
#mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

p1<-ggplot(data=rt, aes(axis1 =lncRNA , axis2 =Telomerase, y = 1))+
  geom_alluvium(aes(fill = Telomerase), width = 0.1, knot.pos = 0.1, reverse = F)+ 
  geom_stratum(fill=NA, color=NA, alpha= 0.5, width = 0.1)+
  scale_fill_manual(values = mycol) +
  geom_text(stat = 'stratum', size =1.5, color='black', label.strata = T)+
  scale_x_discrete(limits = c('lncRNA','Telomerase'), expand=c(0, 0))+
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(),axis.text.x = element_blank()) + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  coord_flip()+ggtitle("")

pdf(file="ggalluvial.pdf", width=10, height=5)
print(p1)
dev.off()

####3.2expTime####

library(limma)        
lncFile="TelomeraseLncExp.txt"      
cliFile="time.txt"    
setwd("")    

rt=read.table(lncFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
data=avereps(data)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)

out=cbind(id=row.names(out),out)
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)

####3.3Model####

library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

coxPfilter=0.05       
setwd("")     

set.seed(222)

rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)

bioForest=function(coxFile=null, forestFile=null, forestCol=null){

  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file=forestFile, width=9, height=16)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex = 10 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol[1],forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,10^a1)
  dev.off()
}

n=1000      
for(i in 1:n){
  inTrain<-createDataPartition(y=rt[,2], p=0.5, list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  outUniTab=data.frame()
  sigGenes=c("futime","fustat")
  for(i in colnames(train[,3:ncol(train)])){
    if(sd(train[,i])>0.1){
      cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
      coxSummary = summary(cox)
      coxP=coxSummary$coefficients[,"Pr(>|z|)"]
      if(coxP<coxPfilter){
        sigGenes=c(sigGenes,i)
        outUniTab=rbind(outUniTab,
                        cbind(id=i,
                              HR=coxSummary$conf.int[,"exp(coef)"],
                              HR.95L=coxSummary$conf.int[,"lower .95"],
                              HR.95H=coxSummary$conf.int[,"upper .95"],
                              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
        )
      }
    }	
  }
  uniSigExp=train[,sigGenes]
  uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
  if(ncol(uniSigExp)<6){next}
  
  x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
  y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
  lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
  if(nrow(geneCoef)<2){next}
  
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
  multiCox=step(multiCox, direction = "both")
  multiCoxSum=summary(multiCox)

  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  outMultiTab=outMultiTab[,1:2]
  
  riskScore=predict(multiCox,type="risk",newdata=train)      
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(riskScore)
  risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
  trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))

  riskScoreTest=predict(multiCox,type="risk",newdata=test)    
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))

  diff=survdiff(Surv(futime, fustat) ~risk,data = train)
  pValue=1-pchisq(diff$chisq, df=1)
  diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
  pValueTest=1-pchisq(diffTest$chisq, df=1)
  

  predictTime=1  
  roc=timeROC(T=train$futime, delta=train$fustat,
              marker=riskScore, cause=1,
              times=c(predictTime), ROC=TRUE)
  rocTest=timeROC(T=test$futime, delta=test$fustat,
                  marker=riskScoreTest, cause=1,
                  times=c(predictTime), ROC=TRUE)	
  
  if((n==1) | ((pValue<0.01) & (roc$AUC[2]>0.65) & (pValueTest<0.05) & (rocTest$AUC[2]>0.6))){

    write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)

    write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
    write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
    bioForest(coxFile="uni.trainCox.txt",forestFile="uni.foreast.pdf",forestCol=c("red","green"))

    write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("lasso.lambda.pdf")
    plot(fit, xvar = "lambda", label = F)
    dev.off()
    pdf("lasso.cvfit.pdf")
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
    dev.off()

    write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
    write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)

    allRiskOut=rbind(trainRiskOut, testRiskOut)
    write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
    break
  }
}

####3.4MultiCox####

library(survival)
library(survminer)

setwd("")               
rt=read.table("lasso.SigExp.txt",header=T,sep="\t",check.names=F,row.names=1)  
rt[,"futime"]=rt[,"futime"]/365

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

pdf(file="forest.pdf",
    width = 8,          
    height = 5,            
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

####3.5SHAP####

library("glmnet")
library("survival")
library(ggplot2)
library(kernelshap)
library(shapviz)

inputFile="lasso.SigExp.txt"    
setwd("")    

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

multiCox <- coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox, direction = "both")

fit=additive_shap(multiCox, rt[,-c(1, 2)])
shp <- shapviz(fit, X_pred = rt[,-c(1, 2)], X=rt[,-c(1, 2)], interactions=T)

important=sort(colMeans(abs(shp$S)), decreasing=T)
showVars=names(important)

pdf(file="barplot.pdf", width=6, height=6)
sv_importance(shp, kind="bar", show_numbers=TRUE)+theme_bw()
dev.off()

pdf(file="bee.pdf", width=7, height=6)
sv_importance(shp, kind = "bee", show_numbers=TRUE)+theme_bw()
dev.off()

pdf(file="waterfall.pdf", width=7, height=5)
sv_waterfall(shp, row_id = 1)
dev.off()

pdf(file="force.pdf", width=9, height=5)
sv_force(shp, row_id = 1)
dev.off()

multiCoxSum=summary(multiCox)
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
outMultiTab=outMultiTab[,1:2]
write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

####4.1Table 1####

trainFile="data.train.txt"    
testFile="data.test.txt"       
cliFile="clinical.1.txt"         
setwd("")     

train=read.table(trainFile, header=T, sep="\t", check.names=F, row.names=1)

test=read.table(testFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

trainCli=cli[row.names(train),]
trainCli=cbind(trainCli, Type="Train")
testCli=cli[row.names(test),]
testCli=cbind(testCli, Type="Test")
rt=rbind(trainCli, testCli)

cliStatOut=data.frame()
for(i in 1:(ncol(rt)-1)){
  nameStat=colnames(rt)[i]
  tableStat=table(rt[,c(nameStat,"Type")])
  tableStatSum=cbind(Total=rowSums(tableStat), tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatPaste=paste(tableStatSum,"(",tableStatRatio,"%)",sep="")
  tableStatOut=matrix(tableStatPaste,ncol=3,dimnames=dimnames(tableStatSum))
  pStat=chisq.test(tableStat[row.names(tableStat)!="unknow",])
  pValueStat=round(pStat$p.value,4)
  pValueCol=c(pValueStat,rep(" ",(nrow(tableStatOut)-1)) )
  tableStatOut=cbind(Variates=nameStat,Type=row.names(tableStatOut),tableStatOut,Pvalue=pValueCol)
  cliStatOut=rbind(cliStatOut,tableStatOut)
}

write.table(cliStatOut,file="cliStat.result.xls",sep="\t",quote=F,row.names=F)

####4.2clinicalCircos####

library(reshape2)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(grid)

riskFile="risk.all.txt"    
cliFile="clinical.txt"    
setwd("") 

data=read.table(cliFile, header=T, sep="\t", check.names=F)
data[,"Age"]=ifelse(data[,"Age"]=="unknow", "unknow", ifelse(data[,"Age"]>65,">65","<=65"))
riskdata = read.table(riskFile, header=T, sep="\t", check.names=F)
colnames(riskdata)[ncol(riskdata)]="Risk"

commonsamples=intersect(data[,1],riskdata[,1])
data = data[match(commonsamples,data[,1]),]
riskdata = riskdata[match(commonsamples,riskdata[,1]),]
data$Risk = riskdata$Risk
highnum = sum(data$Risk=="high")
lownum = sum(data$Risk=="low")

data$Risk[data$Risk=="high"]= sprintf("High\n(%s)",highnum)
data$Risk[data$Risk=="low"]= sprintf("Low\n(%s)",lownum)
piedata = melt(data[,-1],id="Risk")
piedata = piedata[piedata$value!="unknow",]
colnames(piedata) = c("Risk","type","index")
piedata$value = 1
melpiedata = melt(piedata,id=c('Risk','type','index'))
df = dcast(melpiedata,Risk+type+index~variable,length)
df <- df%>%group_by(Risk,type)%>%mutate(ymax = cumsum(value))%>%
  mutate(ymin = ymax - value)%>%
  mutate(ymin=ymin/max(ymax))%>%
  mutate(ymax=ymax/max(ymax))
dfout = df
dfout$Risk = sub('\\n','',dfout$Risk)
#write.table(file="stat.xls", dfout, col.names=T, row.names=F, sep="\t", quote=F)

unitype=levels(df$type)
column.color = rainbow(length(unitype), s=0.7, v=0.7)

p1 =ggplot(df,aes(fill=index,ymax=ymax,ymin=ymin,label=type,
                  xmax=5, xmin=3))+   
  facet_grid(Risk~type,switch="y")+
  theme(panel.spacing = unit(0, "lines"))+
  geom_rect(color="white")+
  coord_polar(theta = 'y')+
  xlim(c(0,5))+ 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="right",legend.box = "horizontal")+
  guides(fill = guide_legend(ncol = 2))+
  theme(legend.title=element_blank())+
  theme(panel.background = element_rect(
    fill = "white", color = "white", size=2)
  )+
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0.5, size = 20), # strip.text
        strip.background = element_rect(fill = "black", linetype = "dotted")) #strip.background
breaks = c()
values = c()
for(i in 1:length(unitype)){
  unitypei = unitype[i]
  dfi=df[df$type==unitypei,,drop=F];
  indexi= unique(dfi$index)
  indexi.col = colorRampPalette(c(column.color[i],"white"))(length(indexi)*2)[1:length(indexi)]
  breaks = c(breaks,indexi)
  values = c(values,rev(indexi.col))
}
p=p1+scale_fill_manual(breaks = breaks,values=values) 

pdf(file="clinicalCircos.pdf", width=10, height=4.5)
print(p)
a = 0.08
c = 0.16
b = (1-a-c)/length(unitype)

for(i in 1:length(unitype)){
  compare.data = spread(data=df[df$type==unitype[i],c('Risk','index','value')], key=Risk, value=value, fill = 0, convert = FALSE, drop = T, sep = NULL)
  chisq.data = chisq.test(compare.data[,-1])
  pi=chisq.data$p.value
  pi=ifelse(pi<0.001,'p<0.001',paste0("p=",sprintf("%.03f",pi)) )
  grid.text(pi, x=a+b*(2*i-1)/2, y=0.1) # middle
}
dev.off()

####4.3clinicalpheatmap####

library(pheatmap)

setwd("")

exp=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
exp = exp[,-c(1,2,c(length(colnames(exp))-1))]

expCluster= exp

cli=read.table("clinical.2.txt", header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

data=data[order(data$risk),]
Type=data[,((ncol(exp)):ncol(data))]
data=t(data[,1:c(ncol(exp)-1)])

bioCol=c("#0072B5","#BC3C29","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$risk)))]
names(prgCluCol)=levels(factor(Type$risk))
ann_colors[["risk"]]=prgCluCol

pdf("clinical.heatmap.pdf", width=5, height=3.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#BC3C29",5), "white", rep("#0072B5",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         show_colnames=F,
         scale="row",
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

####4.4clinicalsubgroup####

library(limma)
library(ggpubr)
setwd("")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

cli2=read.table("clinical.3.txt", header=T, sep="\t", check.names=F, row.names=1)
cli2[,"Age"]=ifelse(cli2[,"Age"]=="unknow", "unknow", ifelse(cli2[,"Age"]>65,">65","<=65"))
samesamples = intersect(rownames(risk),rownames(cli2))

risk = risk[samesamples,]
cli2 = cli2[samesamples,]

data2 = cbind(rownames(cli2),cli2,risk[,"riskScore",drop=F])


for (i in 2:c(length(colnames(data2)) -1)) {
  
  rt = data2[,c(1,i,length(colnames(data2)))]
  rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
  
  
  x=colnames(rt)[2]
  y=colnames(rt)[3]
  colnames(rt)=c("id","Type","Expression")
  
  group=levels(factor(rt$Type))
  rt$Type=factor(rt$Type, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

  boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                    xlab=x,
                    ylab=y,
                    legend.title=x,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  
  pdf(file=paste0("",x,".boxplot.pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
  
}

####4.5clinical survical####

library(survival)
library(survminer)

riskFile="risk.all.txt"    
cliFile="clinical(Age).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)

  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}
riskFile="risk.all.txt"    
cliFile="clinical(Gender).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}
riskFile="risk.all.txt"    
cliFile="clinical(Grade).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}
riskFile="risk.all.txt"    
cliFile="clinical(M0).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}
riskFile="risk.all.txt"    
cliFile="clinical(N0).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}
riskFile="risk.all.txt"    
cliFile="clinical(Stage).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}
riskFile="risk.all.txt"    
cliFile="clinical(T).txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette=c("#BC3C29", "#0072B5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 5.5,     
      height =4.8)    
  print(surPlot)
  dev.off()
}


####5.1Model####

library(survival)
library(survminer)
library(pheatmap) 

setwd("")      

bioSurvival=function(inputFile=null, outFile=null,title0=null){

  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
 
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("#BC3C29","#0072B5"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25,
                     title=paste0(title0))
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}

bioSurvival(inputFile="risk.train.txt", title0="Training set",outFile="surv.train.pdf")
bioSurvival(inputFile="risk.test.txt", title0="Testingg set",outFile="surv.test.pdf")
bioSurvival(inputFile="risk.all.txt", title0="All set",outFile="surv.all.pdf")

bioRiskPlot=function(inputFile=null, project=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1) 
  rt=rt[order(rt$riskScore),]   
  
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=3)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)",
       ylab="Risk score",
       col=c(rep("#0072B5",lowLength),rep("#BC3C29",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("#BC3C29","#0072B5"),cex=1.2)
  dev.off()
  
  color=as.vector(rt$fustat)
  color[color==1]="#BC3C29"
  color[color==0]="#0072B5"
  pdf(file=paste0(project, ".survStat.pdf"), width=7, height=3)
  plot(rt$futime, pch=19,ylim=c(0,18),
       xlab="Patients (increasing risk socre)",
       ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("#BC3C29","#0072B5"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  ann_colors=list()
  bioCol=c("#0072B5", "#BC3C29")
  names(bioCol)=c("low", "high")
  ann_colors[["Risk"]]=bioCol
  
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(Risk=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=paste0(project, ".heatmap.pdf"), width=4.5, height=2)
  pheatmap(rt1, 
           annotation=annotation,
           annotation_colors = ann_colors, 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = F,
           scale="row",
           color = colorRampPalette(c(rep("#0072B5",3.5), "white", rep("#BC3C29",3.5)))(50),
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}

bioRiskPlot(inputFile="risk.train.txt", project="train")
bioRiskPlot(inputFile="risk.test.txt", project="test")
bioRiskPlot(inputFile="risk.all.txt", project="all")

####5.2ROC####

library(survival)
library(survminer)
library(timeROC)

riskFile="risk.all.txt"    
cliFile="clinical.4.txt"     
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

bioCol=c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")

ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="all ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=3)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=3)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=3)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=3, bty = 'n')
dev.off()

predictTime=1   
aucText=c()
pdf(file="all cliROC.pdf", width=5, height=5)
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}

legend("bottomright", aucText,lwd=3,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()



riskFile="risk.train.txt"    
cliFile="clinical.4.txt"     
setwd("")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

bioCol=c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")

ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="train ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=3)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=3)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=3)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=3, bty = 'n')
dev.off()


predictTime=1     
aucText=c()
pdf(file="train cliROC.pdf", width=5, height=5)

i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}

legend("bottomright", aucText,lwd=3,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()



riskFile="risk.test.txt"    
cliFile="clinical.4.txt"     
setwd("")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

bioCol=c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")

ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="test ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=3)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=3)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=3)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=3, bty = 'n')
dev.off()

predictTime=1    
aucText=c()
pdf(file="test cliROC.pdf", width=5, height=5)

i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}

legend("bottomright", aucText,lwd=3,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

####6.1Prognosis factors####

library(survival)
library(survminer)

riskFile="risk.all.txt"         
cliFile="clinical.5.txt"    
setwd("")   

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk=risk[,"risk"])
rt$risk=factor(rt$risk, levels=c("low", "high"))

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab, file="cox.result.txt", sep="\t", row.names=F, quote=F)

pdf(file="forest.pdf", width=8, height=6, onefile = FALSE)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.18, 0.38), 
         fontsize = 0.8, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()

####6.2nomo####

library(survival)
library(regplot)
library(rms)
library(survcomp)

riskFile="risk.all.txt"      
cliFile="clinical.5.txt"       
setwd("")    


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)

res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),    
             dencol="#6EE2FFFF", boxcol="#F7C530FF", 
             clickable=F,
             title="",               
             points=TRUE,            
             droplines=TRUE,        
             observation=rt[1,],    
             rank="sd",
             failtime = c(1,3,5),   
             prfail = F)
dev.copy2pdf(file="Nomo.pdf", width=8, height=6, out.type="pdf")

nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

pdf(file="calibration.pdf", width=5, height=5)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

####6.3C-index####

library(dplyr)
library(survival)
library(rms)
library(pec)

riskFile="risk.all.txt"    
cliFile="clinical.4.txt"      
setwd("")    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
Gender=cph(Surv(futime,fustat)~Gender, data=rt, surv=TRUE)
Grade=cph(Surv(futime,fustat)~Grade, data=rt, surv=TRUE)
Stage=cph(Surv(futime,fustat)~Stage, data=rt, surv=TRUE)
c_index  <- cindex(list("Risk score"=riskScore, 
                        "Age"=Age,
                        "Gender"=Gender,
                        "Grade"=Grade,
                        "Stage"=Stage),
                   formula=Surv(futime,fustat)~ .,
                   data=rt,
                   eval.times=seq(0,10,1),
                   splitMethod="bootcv",
                   B=1000
)

pdf(file="C-index.pdf", width=5.5, height=5)
plot(c_index, 
     xlim=c(0,10), ylim=c(0.4,0.8), lwd=3,
     col=bioCol, xlab="Time (years)",
     legend.x=6, legend.y=0.82, legend.cex=1)
dev.off()

####6.4DCA####

library(survival)
library(ggDCA)

setwd("")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

cli=read.table("clinical.6.txt", header=T, sep="\t", check.names=F, row.names=1)
Nome=read.table("nomoRisk.txt", header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
Nome = Nome[samSample,"Nomogram",drop=F]
rt=cbind(risk1,Nome,cli)
rt[,"Age"]=ifelse(rt[,"Age"]>65, 1, 0)

predictTime=1    
Risk<-coxph(Surv(futime,fustat)~risk,rt)
Nomogram<-coxph(Surv(futime,fustat)~Nomogram,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
Gender<-coxph(Surv(futime,fustat)~Gender,rt)
Grade<-coxph(Surv(futime,fustat)~Grade,rt)
Stage<-coxph(Surv(futime,fustat)~Stage,rt)
Tstage<-coxph(Surv(futime,fustat)~Tstage,rt)

pdf(file=paste0("",predictTime,"year.DCA.pdf"), width=6.5, height=5.2)
d_train=dca(Risk,Nomogram,Age,Gender,Grade,Stage,Tstage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()

predictTime=3    
Risk<-coxph(Surv(futime,fustat)~risk,rt)
Nomogram<-coxph(Surv(futime,fustat)~Nomogram,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
Gender<-coxph(Surv(futime,fustat)~Gender,rt)
Grade<-coxph(Surv(futime,fustat)~Grade,rt)
Stage<-coxph(Surv(futime,fustat)~Stage,rt)
Tstage<-coxph(Surv(futime,fustat)~Tstage,rt)

pdf(file=paste0("",predictTime,"year.DCA.pdf"), width=6.5, height=5.2)
d_train=dca(Risk,Nomogram,Age,Gender,Grade,Stage,Tstage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()

predictTime=5   
Risk<-coxph(Surv(futime,fustat)~risk,rt)
Nomogram<-coxph(Surv(futime,fustat)~Nomogram,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
Gender<-coxph(Surv(futime,fustat)~Gender,rt)
Grade<-coxph(Surv(futime,fustat)~Grade,rt)
Stage<-coxph(Surv(futime,fustat)~Stage,rt)
Tstage<-coxph(Surv(futime,fustat)~Tstage,rt)

pdf(file=paste0("",predictTime,"year.DCA.pdf"), width=6.5, height=5.2)
d_train=dca(Risk,Nomogram,Age,Gender,Grade,Stage,Tstage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()

####6.5PCA####

library(limma)
library(scatterplot3d)
setwd("")   

myPCA=function(input=null, output=null){

  rt=read.table(input, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0.5,]
  
  type=sapply(strsplit(colnames(data),"\\-"),"[",4)
  type=sapply(strsplit(type,""),"[",1)
  type=gsub("2","1",type)
  data=t(data[,type==0])
  rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
  
  risk=read.table("risk.all.txt", header=T, sep="\t", row.names=1, check.names=F)
  sameSample=intersect(rownames(data),rownames(risk))
  data=data[sameSample,]
  risk=risk[sameSample,]
  group=as.vector(risk[,"risk"])
  
  data.class <- rownames(data)
  data.pca <- prcomp(data, scale. = TRUE)
  pcaPredict=predict(data.pca)
  
  color=ifelse(group=="low",4,2)
  pdf(file=output, width=7, height=7)
  par(oma=c(1,1,2.5,1))
  s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=50)
  legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
  dev.off()
}

myPCA(input="symbol.txt", output="PCA.allGene.pdf")
myPCA(input="TelomeraseExp.txt", output="PCA.TelomeraseGene.pdf")
myPCA(input="TelomeraseLncExp.txt", output="PCA.TelomeraseLncRNA.pdf")
myPCA(input="hubgenesExp.txt", output="PCA.hubgene.pdf")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])

data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

color=ifelse(group=="low",4,2)
pdf(file="PCA.riskLnc.pdf", width=6.5, height=6)
par(oma=c(1,1,2.5,1))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=45)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
dev.off()

####7.1mutation####

library(maftools)
library(dplyr)

setwd("")

files <- list.files(pattern = '*.gz',recursive = TRUE)
all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}

all_mut$Tumor_Sample_Barcode = substr(all_mut$Tumor_Sample_Barcode,1,12)

all_mut <- read.maf(all_mut)

a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))

sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}

mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))

colnames(gene_count)[1] = "Gene"
colnames(gene_count)[2] = "Num"
write.table(gene_count,'geneMut.txt', sep="\t", quote=F, row.names = F)
write.mafSummary(maf = all_mut,basename = "input")

tmb_table = tmb(maf = all_mut)   
tmb_table = tmb(maf = all_mut,logScale = F)  

tmb_table = tmb_table[,c(1,3)]
tmb_table = as.data.frame(tmb_table)
tmb_table[,1]=substr(tmb_table[,1],1,12)
colnames(tmb_table)
tmb_table <- aggregate( . ~ Tumor_Sample_Barcode,data=tmb_table, max)
colnames(tmb_table)[1] = "id"
colnames(tmb_table)[2] = "TMB"
write.table(tmb_table,'TMB.txt', sep="\t", quote=F, row.names = F)

####7.2waterfall plot####

library(maftools)     

setwd("")

all.maf = read.maf("input.maf")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

low.tab = subset(x = outTab,Risk == "low")
low.maf = subsetMaf(maf = all.maf, tsb = low.tab[,1])
write.mafSummary(maf = low.maf,basename = "low.maf")

high.tab = subset(x = outTab,Risk == "high")
high.maf = subsetMaf(maf = all.maf, tsb = high.tab[,1])
write.mafSummary(maf = high.maf,basename = "high.maf")

geneNum=15
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

ann_colors=list()
col=c("#0072B5", "#BC3C29")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.mafmaf", clinicalData="5treat/ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


####7.3TMB####


library(limma)
library(ggpubr)
library(survival)
library(survminer)
setwd("E:/10FS/48HCC")

tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB=log2(data$TMB+1)

data$risk=ifelse(data$risk=="high", "high", "low")
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggviolin(data, x="risk", y="TMB", fill="risk",
                 xlab="",
                 ylab="Tumor tmbation burden (log2)",
                 legend.title="",
                 palette = c("#BC3C29", "#0072B5"),
                 add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file="riskTMB.pdf", width=5, height=5)
print(boxplot)
dev.off()

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)  
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)    

sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  bioCol=c("#EE1D23","#19468B","#9932CC","#EE9A00","#008B8A","#f25f5c")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  pdf(file=outFile, onefile = FALSE, width=6.1, height=5.3)
  print(surPlot)
  dev.off()
}
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")

####7.4pre.TIDE####

setwd("")

rt=read.table("mRNA.txt", header=T, sep="\t", check.names=F,row.names = 1)
exp=rt
Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize.txt', sep="\t", quote=F, row.names = TRUE)

rt=read.table("mRNA.txt", header=T, sep="\t", check.names=F,row.names = 1)
exp=log2(exp+1)
Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize_log.txt', sep="\t", quote=F, row.names = TRUE)


####7.5TIDE####

setwd("")

library(limma)
library(ggpubr)
library(data.table)
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)

tide=fread("TIDE.csv")
tide = as.data.frame(tide)
rownames(tide) = tide[,1]

colnames(tide)
tide = tide[,c("TIDE","Dysfunction","Exclusion"),drop = F]

tide = t(tide)
group=sapply(strsplit(colnames(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide = tide[,group == 0,drop=F]
colnames(tide)=substr(colnames(tide),1,12)
tide = t(tide)
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)

for (l in c("TIDE","Dysfunction","Exclusion")) {
  #l = "TIDE"
  data=cbind(tide, risk)
  data = data[,c(l,"risk")]
  data$risk=ifelse(data$risk=="high", "high", "low")
  group=levels(factor(data$risk))
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  gg1=ggviolin(data, x="risk", y=l, fill = "risk", 
               xlab="", ylab=l,
               palette=c("#BC3C29", "#0072B5"),
               legend.title="risk",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0("",l,".pdf"), width=5, height=5)
  print(gg1)
  dev.off()
  
}

####7.6RNAss####

library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
setwd("")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

RNAss=read.table("StemnessScores_RNAexp_20170127.2.tsv", header=T, sep="\t",check.names=F, row.names=1)
RNAss=t(RNAss[1,,drop=F])
rownames(RNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(RNAss))
RNAss=avereps(RNAss)

sameSample=intersect(row.names(risk), row.names(RNAss))
risk=risk[sameSample,"riskScore",drop=F]
RNAss=RNAss[sameSample,,drop=F]
data=cbind(RNAss, risk)

xlab="riskScore"
ylab="RNAss"
outFile="RNAss.cor.pdf"
x=as.numeric(data[,xlab])
x[x>quantile(x,0.99)]=quantile(x,0.99)
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
  xlab("Risk score") + ylab(ylab)+ ylim(0,0.7)+
  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()

####8.1CIBERSORT####

library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
library(bseqsc)
library(ggplot2)
library(CIBERSORT)
library(corrplot)
library(vioplot)

setwd("")
data(LM22) 
data=read.table("lncRNA.txt", header=T, sep="\t", check.names=F,row.names = 1)
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

results <- cibersort(sig_matrix = LM22, mixture_file = data)

results=as.matrix(results[,1:(ncol(results)-3)])
results=rbind(id=colnames(results),results)
write.table(results, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

####8.2immuneCor####

library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
library(bseqsc)
library(ggplot2)
library(CIBERSORT)
library(corrplot)
library(vioplot)

setwd("")

data(LM22) 

immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=as.matrix(immune)
data=t(immune)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0,drop=F]
colnames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data = t(data)

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
lowSample=rownames(risk)[risk[,"risk"]=="low"]
highSample=rownames(risk)[risk[,"risk"]=="high"]

lowSameSample=intersect(row.names(data), lowSample)
highSameSample=intersect(row.names(data), highSample)
data=t(data[c(lowSameSample,highSameSample),])
conNum=length(lowSameSample)
treatNum=length(highSameSample)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
               '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', 
               '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', 
               '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', 
               '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', 
               '#3A6963','#968175')

pdf("immune.barplot.pdf", width=25, height=12)
col=rainbow(nrow(data),s=0.7,v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=my36colors,ylab="Relative Percent",xaxt="n",yaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="#0072B5")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="#BC3C29")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=my36colors,pch=15,bty="n",cex=1.3)
dev.off()

pdf(file="immune.corHeatmap.pdf", width=11, height=11)
par(oma=c(0.5,1,1,1.2))
corData=t(data)
corData=corData[,colMeans(corData)>0]
M=cor(corData)
corrplot(M,
         order="hclust",        
         method = "color",       
         diag = TRUE,            
         tl.col="black",         
         addCoef.col = "black",  
         number.cex=0.75,        
         col=colorRampPalette(c("#0072B5", "white", "#BC3C29"))(50))   
dev.off()


rt=t(data)
pdf("immune.vioplot.pdf", height=8, width=12)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.02),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  lowData=rt[1:conNum,i]
  highData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = "#0072B5")
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = "#BC3C29")
  wilcoxTest=wilcox.test(lowData, highData)
  p=wilcoxTest$p.value
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topleft", 
       c("Low risk", "High risk"),
       lwd=3.5, bty="n", cex=1.2,
       col=c("#0000FD","#FF0000"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()

####8.3immunecell####

library(limma)
library(scales)
library(ggplot2)
library(ggtext)
setwd("")

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

immune=read.csv("infiltration_estimation_for_tcga.csv", header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

x=as.numeric(risk)
outTab=data.frame()
for(i in colnames(immune)){
  y=as.numeric(immune[,i])
  corT=cor.test(x, y, method="spearman")
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<0.05){
    outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
  }
}
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

corResult=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))   
pdf(file="cor.pdf", width=10, height=12)      
ggplot(data=b, aes(x=cor, y=immune, color=Software))+
  labs(x="Correlation coefficient",y="Immune cell")+
  geom_point(size=4.1)+
  theme(panel.background=element_rect(fill="white",size=1,color="black"),
        panel.grid=element_line(color="grey75",size=0.5),
        axis.ticks = element_line(size=0.5),
        axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
dev.off()

####8.4ICGs####

library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"   
riskFile="risk.all.txt"   
geneFile="ICGs.txt"    
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
data=avereps(data)

risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]

sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
  if(sd(rt1[,i])<0.001){next}
  wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]

rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")

group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
                  xlab="",
                  ylab="Gene expression",
                  legend.title="Risk",
                  width=0.8,
                  palette = c("#0072B5", "#BC3C29") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.s=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

pdf(file="checkpoint.diff.pdf", width=10, height=4.5)
print(boxplot)
dev.off()

jco <- c("#0072B5","#BC3C29")

boxplot=ggplot(data = rt1,aes(x = Gene, y = Expression, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + 
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("Gene expression")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), 
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        
  )

pdf(file="checkpoint.diff1.pdf", width=10, height=4.5)
print(boxplot)
dev.off()

####8.5immFunction####

library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
expFile="symbol.txt"         
gmtFile="immune1.gmt"         
riskFile="risk.all.txt"          
socreFile="immFunScore.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)

risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)

sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(data, risk)

data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
            ylab="Score",add = "none",xlab="",palette = c("#0072B5", "#BC3C29") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.s=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

pdf(file="immFunction.pdf", width=10, height=5)
print(p)
dev.off()


jco <- c("#BC3C29", "#0072B5")

boxplot=ggplot(data = data,aes(x = Type, y = Score, fill = Risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + 
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("Score")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risk),method = "wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), 
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

pdf(file="immFunction.boxplot1.pdf", width=11, height=5.5)
print(boxplot)
dev.off()

####8.6immSubtype####

library(RColorBrewer)         
riskFile="risk.all.txt"                             
subtypeFile="Subtype_Immune_Model_Based.txt"    
setwd("")   

ChischisqTest <- function(tabledata){
  ct = chisq.test(tabledata)
  ct.pvalue = ct$p.value
  ct.pvalue = ifelse(ct.pvalue<0.001,0.001,round(ct.pvalue,3))
}

rect2text <- function(x1,y1,x2,y2,text,rect.col,text.cex,text.col="white",...){
  rect(x1,y1,x2,y2,col=rect.col,border=NA,...)
  text((x1+x2)/2,(y1+y2)/2,text,cex=text.cex,col=text.col,...)
}

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(risk)[ncol(risk)]="Risk"

subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))

sameSample=intersect(row.names(subtype), row.names(risk))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
risk=risk[sameSample,]
data=cbind(subtype, as.data.frame(risk))
typeTab=table(data$subtype)
typeName=names(typeTab[typeTab>6])
merge.data=data[which(data[,"subtype"] %in% typeName),]
cliName=colnames(merge.data)[1]
colnames(merge.data)[1]="Clinical"

groups = sort(unique(merge.data$Clinical))
groups.table = table(merge.data$Clinical)
groups.num = length(groups)
samples.num = nrow(merge.data)
groups.lownum = sum(merge.data$Risk=="low")
groups.highnum = sum(merge.data$Risk=="high")
pergrouptable = table(merge.data$Clinical, merge.data$Risk)
ct.pvalue = ChischisqTest(pergrouptable)

xlim1=0; xlim2=samples.num; ylim1=1; ylim2=10
space = dist.up = 0.1
space.x = 1/100*xlim2
groups.col = col = colorRampPalette(brewer.pal(9, "Paired"))(groups.num)
high_low.col = c(rgb(235/255,116/255,106/255),rgb(30/255,181/255,184/255))
bottom.unitwidth = xlim2 / (groups.num+1+1.3)
bottom.left.width = bottom.unitwidth*1.3

pdf(file="immSubtype.pdf", width=10, height=7)
plot(1,type="n",xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2),axes=F,xlab="",ylab="")

header.text.cex = 1.5
header.text = gettextf("%s TCGA patients",samples.num)
header.text.width = strwidth(header.text,cex=header.text.cex)
arrows.len = (xlim2 - header.text.width - space.x*2)/2
text(xlim2/2,9.5,header.text,cex=header.text.cex,font=2)
arrows(0,9.5,arrows.len,9.5,code=1,lwd=3,length=0.2)
arrows(xlim2-arrows.len,9.5,xlim2,9.5,code=2,lwd=3,length=0.2)
rect2text(x1=0,y1=4.5+space,x2=bottom.left.width,y2=7-space,text.cex=1.2,text = 
            "TRLs\n groups",rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)
#draw bottom header
bottom.header.text = paste0(cliName, " group (n=" , samples.num, ")")
rect2text(x1=bottom.left.width+space.x,y1=6+space/2,x2=xlim2,y2=7-space,text.cex=1.2,
          text= bottom.header.text,rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)

#draw bottom left 2
bottom.left2.text = gettextf("TRLs-low\n(n=%s)",groups.lownum)
rect2text(x1=0,y1=3+space/2,x2=bottom.left.width,y2=4.5-space/2,text.cex=1.2,text= 
            bottom.left2.text,rect.col=high_low.col[2],text.col="white",font=2)
#draw bottom left 3
bottom.left3.text = gettextf("TRLs-high\n(n=%s)",groups.highnum)
rect2text(x1=0,y1=1.5+space/2,x2=bottom.left.width,y2=3-space/2,text.cex=1.2,text= 
            bottom.left3.text,rect.col=high_low.col[1],text.col="white",font=2)

start = 0
for(i in 1:length(groups.table)){
  end = groups.table[i]+start
  namesi = names(groups.table)[i]
  # up rect
  rect2text(x1=start,y1=8+space,x2=end,y2=9-space,text.cex=1.1,text= namesi, 
            rect.col=groups.col[i],text.col="white")
  merge.datai = merge.data[merge.data$Clinical==namesi,,drop=F]
  high.num = sum(merge.datai$Risk=="high")
  low.num = sum(merge.datai$Risk=="low")
  # middle low rect
  rect(start,7+dist.up,start+low.num,8-dist.up,col=high_low.col[2],border=NA)
  # middle high rect
  rect(start+low.num,7+dist.up,end,8-dist.up,col=high_low.col[1],border=NA)
  bottom.starti = bottom.left.width+(i-1)* bottom.unitwidth
  bottom.endi = bottom.starti+bottom.unitwidth
  subheader.text = gettextf("%s\n(n=%s, %s)",namesi,nrow(merge.datai),paste0(round(nrow(merge.datai)/samples.num*100),"%"))
  rect2text(x1=bottom.starti+space.x,y1=4.5+space,x2=bottom.endi,y2=6-space,text.cex=1.1,text= subheader.text,rect.col=groups.col[i],text.col="white")
  low.texti = gettextf("%s(%s)",low.num,paste0(round(low.num/groups.lownum*100),"%"))
  rect2text(x1=bottom.starti+space.x,y1=3+space/2,x2=bottom.endi,y2=4.5-space/2,text.cex=1.1,text= low.texti,rect.col="grey90",text.col="black")
  high.texti = gettextf("%s(%s)",high.num,paste0(round(high.num/groups.highnum*100),"%"))
  rect2text(x1=bottom.starti+space.x,y1=1.5+space/2,x2=bottom.endi,y2=3-space/2,text.cex=1.1,text= high.texti,rect.col="grey70",text.col="black")
  start = end
}

rect2text(x1=bottom.endi+space.x,y1=4.5+space,x2=xlim2,y2=6-space,text.cex=1.1,text= "P-value",rect.col="grey70",text.col="black",font=3)
rect2text(x1=bottom.endi+space.x,y1=1.5+space/2,x2=xlim2,y2=4.5-space/2,text.cex=1.1,text= ct.pvalue,rect.col="grey70",text.col="black")

rect(0,8+space,xlim2,9-space,lwd=2,border="grey30")
rect(0,7+space,xlim2,8-space,lwd=2,border="grey30")
rect(0,1.5+space/2,xlim2,7-space,lwd=2,border="grey30")

legend("bottom",legend=c("TRLs-low","TRLs-high"),col=rev(high_low.col),bty="n",ncol=2,pch=15,pt.cex=1.3,cex=1.3)
dev.off()

####8.7oncopredict####

library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="symbol.txt"  
setwd("")   

rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

calcPhenotype(trainingExprData = GDSC2_Expr,   
              trainingPtype = GDSC2_Res,       
              testExprData = data,           
              batchCorrect = 'eb',            
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,   
              minNumSamples = 10,               
              printOutput = TRUE,               
              removeLowVaringGenesFrom = 'rawData')

####8.8plot####

library(limma)
library(oncoPredict)
library(parallel)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(ggrepel)
library(corrplot)
library(Hmisc)
library(ggExtra)

setwd("")

set.seed(999)

senstivity <- read.csv("DrugPredictions.csv", header = T, sep = ",", check.names = F, row.names = 1)
colnames(senstivity) <- gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

group <- sapply(strsplit(rownames(senstivity), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
senstivity <- senstivity[group == 0, ]

senstivity <- t(senstivity)
colnames(senstivity) <- substr(colnames(senstivity), 1, 12)
senstivity <- t(senstivity)

risk <- read.table("risk.all.txt", header = T, sep = "\t", check.names = F, row.names = 1)

sameSample <- intersect(row.names(risk), row.names(senstivity))
risk <- risk[sameSample, "risk", drop = F]
senstivity <- senstivity[sameSample, , drop = F]
senstivity[is.na(senstivity)] <- 0
senstivity <- log2(senstivity + 1)
rt <- cbind(risk, senstivity)

rt$risk <- factor(rt$risk, levels = c("low", "high"))
type <- levels(factor(rt[, "risk"]))
comp <- combn(type, 2)
my_comparisons <- list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]
}

sigGene <- c()
for (i in colnames(rt)[2:(ncol(rt))]) {
  if (sd(rt[, i]) < 0.05) next
  wilcoxTest <- wilcox.test(rt[, i] ~ rt[, "risk"])
  if (wilcoxTest$p.value < 0.05) {
    sigGene <- c(sigGene, i)
  }
}
sigGene <- c(sigGene, "risk")
rt <- rt[, sigGene]

jco <- c("#BC3C29", "#0072B5")

output_dir <- ""
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (drug in colnames(rt)[-1]) {  
  rt1 <- rt[, c("risk", drug)]
  colnames(rt1) <- c("risk", "sensitivity")
  
  if (nrow(rt1) == 0 || all(is.na(rt1$sensitivity))) {
    warning(paste("No valid data for drug:", drug))
    next
  }
  
  group_counts <- table(rt1$risk)
  cat("Drug:", drug, "Group counts:\n")
  print(group_counts)
  
  use_notch <- all(group_counts >= 5)  
  
  gg1 <- ggplot(data = rt1, aes(x = risk, y = sensitivity, fill = risk)) +
    scale_fill_manual(values = jco[2:1]) +
    geom_violin(alpha = 0.4, position = position_dodge(width = 0.75), size = 0.8, color = "black") +
    geom_boxplot(
      notch = use_notch,
      outlier.size = -1, 
      color = "black", 
      lwd = 0.8, 
      alpha = 0.7, 
      width = 0.3
    ) +
    geom_point(shape = 21, size = 2, position = position_jitterdodge(), color = "black", alpha = 1) +
    theme_classic() +
    ylab(paste0(drug, " sensitivity (IC50)")) +
    xlab("Risk group") +
    stat_compare_means(comparisons = my_comparisons) +
    theme(
      axis.ticks = element_line(size = 0.2, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.line = element_blank(),
      
      axis.title = element_text(face = "bold.italic", colour = "#441718", size = 16),
      axis.text = element_text(face = "bold.italic", colour = "#441718", size = 16),
      plot.title = element_text(face = "bold.italic", colour = "#441718", size = 16),
      
      legend.position = "none",
      legend.text = element_text(face = "bold.italic"),
      legend.title = element_text(face = "bold.italic", size = 13),
      
      panel.border = element_rect(fill = NA, color = "#35A79D", size = 1.5, linetype = "solid"),
      panel.background = element_rect(fill = "#F1F6FC"),
      panel.grid.major = element_line(color = "#CFD3D6", size = 0.5, linetype = "dotdash")
    )

  pdf_file <- file.path(output_dir, paste0(drug, ".notch.pdf"))
  pdf(file = pdf_file, width = 5, height = 4.5)
  print(gg1)
  dev.off()
  
  cat("yes:", pdf_file, "Notch:", use_notch, "\n")
}

####9.1TMP####

library(limma)            
inputFile="symbol.txt"  
setwd("")    

outTab=data.frame()
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)
data=t(data)

fpkmToTpm=function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(data, 2, fpkmToTpm)

tpmOut=rbind(ID=colnames(tpm), tpm)
write.table(tpmOut, file="TCGA.TPM.txt", sep="\t", col.names=F, quote=F)

####9.2CIBERSORT####

library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
library(bseqsc)
library(ggplot2)
library(CIBERSORT)
library(corrplot)
library(vioplot)

setwd("")

data(LM22) 

data=read.table("TCGA.TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

results <- cibersort(sig_matrix = LM22, mixture_file = data)

results=as.matrix(results[,1:(ncol(results)-3)])
results=rbind(id=colnames(results),results)
write.table(results, file="CIBERSORT-Results1.txt", sep="\t", quote=F, col.names=F)

####9.3score####

library(estimate)         
inputFile="TCGA.TPM.txt"   
setwd("")   

filterCommonGenes(input.f=inputFile, 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

scores=read.table("estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores=scores[,1:3]
out=rbind(ID=colnames(scores), scores)
write.table(out, file="scores.txt", sep="\t", quote=F, col.names=F)

####9.4cluster####

library(ConsensusClusterPlus)        
cellFile="CIBERSORT-Results1.txt"    =
workDir=""     
setwd(workDir)     

data=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


clusterNum=4      
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="ICIcluster.txt",sep="\t",quote=F,col.names=F)

####9.5immCor####

library(corrplot)       
cellFile="CIBERSORT-Results1.txt"   
scoreFile="scores.txt"          
setwd("")  

cell=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(cell), row.names(score))
cell=cell[sameSample,]
score=score[sameSample,]
data=t(cbind(cell, score[,1:2]))

data=t(data)   
M=cor(data)    
res1=cor.mtest(data, conf.level = 0.95)

pdf(file="immCor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "pie",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "blank",
         sig.level=0.05,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()

####9.6immCell####

library(reshape2)
library(ggpubr)
cellFile="CIBERSORT-Results1.txt"   
scoreFile="scores.txt"          
clusterFile="ICIcluster.txt"      
setwd("")  

cell=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(cell), row.names(score))
cell=cell[sameSample,]
score=score[sameSample,]
immData=cbind(cell, score[,1:2], ICIcluster=cluster[sameSample,])
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(immData$ICIcluster))
immData$ICIcluster=letter[match(immData$ICIcluster, uniqClu)]
immData[,-ncol(immData)]=scale(immData[,-ncol(immData)])

data=melt(immData, id.vars=c("ICIcluster"))
colnames(data)=c("ICIcluster", "Immune", "Fraction")
data$Fraction[data$Fraction>12]=12

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(immData[,"ICIcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="ICIcluster", 
            ylab="Scale of fraction",
            xlab="",
            legend.title="ICI cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf", width=8, height=6.5)                   
p+stat_compare_means(aes(group=ICIcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

####9.7survival####

library(survival)
library(survminer)
clusterFile="ICIcluster.txt"    
cliFile="time.txt"             
setwd("")   

cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,], ICIcluster=cluster[sameSample,])
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(rt$ICIcluster))
rt$ICIcluster=letter[match(rt$ICIcluster, uniqClu)]

length=length(levels(factor(rt$ICIcluster)))
diff=survdiff(Surv(futime, fustat) ~ ICIcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ ICIcluster, data = rt)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="ICI cluster",
                   legend.labs=levels(factor(rt[,"ICIcluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()

####9.8ICG####

#引用包
library(limma)
library(ggplot2)
library(ggpubr)
expFile="TCGA.TPM.txt"            
clusterFile="ICIcluster.txt"  
gene="IDO1"       #revisability           
showName="IDO1"   #revisability          
setwd("")  

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=rbind(data,gene=data[gene,])
data=t(data[c(gene,"gene"),])

cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,], ICIcluster=cluster[sameSample,])
data=as.data.frame(data)
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(data$ICIcluster))
data$ICIcluster=letter[match(data$ICIcluster, uniqClu)]

group=levels(factor(data$ICIcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"ICIcluster"])))]
pdf(file="IDO1.pdf", width=6, height=5)
ggviolin(data, x="ICIcluster", y="gene", fill = "ICIcluster",
         ylab=paste0(showName, " expression"),
         xlab="",
         legend.title="ICI cluster",
         palette=bioCol, 
         add="boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

####10.1pre.Single-cell####

packages <- c("Seurat", "dplyr", "magrittr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is not installed. Please install it before running this script."))
  }
}
lapply(packages, library, character.only = TRUE)

workDir <- "\\GSE242889_RAW"
setwd(workDir)

dirs <- list.dirs(workDir, recursive = FALSE, full.names = TRUE)
names(dirs) <- basename(dirs)

counts_list <- list()
gene_lists <- list()

message("Step 1/4: Reading 10X data for each sample...")
for (sample_name in names(dirs)) {
  sample_dir <- dirs[[sample_name]]
  message(paste("Reading data for sample:", sample_name))
  tryCatch({
    counts <- Read10X(data.dir = sample_dir)
    counts_list[[sample_name]] <- counts
    gene_lists[[sample_name]] <- rownames(counts)
  }, error = function(e) {
    warning(paste("Failed to read data for sample:", sample_name, ":", e$message))
  })
}

message("Step 2/4: Identifying common genes across all samples...")
if (length(gene_lists) < 2) stop("Not enough samples to find common genes.")
common_genes <- Reduce(intersect, gene_lists)
message(paste("Number of common genes across all samples:", length(common_genes)))

message("Step 3/4: Subsetting counts matrices and adjusting cell barcodes...")
for (sample_name in names(counts_list)) {
  counts <- counts_list[[sample_name]]
  counts <- counts[common_genes, , drop = FALSE]
  colnames(counts) <- paste(sample_name, colnames(counts), sep = "_")
  counts_list[[sample_name]] <- counts
}

message("Step 4/4: Combining counts matrices and creating Seurat object...")
counts_combined <- do.call(cbind, counts_list)
pbmc <- CreateSeuratObject(counts_combined, min.cells = 3, min.features = 100)

expression_matrix <- as.matrix(GetAssayData(object = pbmc, slot = "counts"))
expression_df <- data.frame(gene = rownames(expression_matrix), expression_matrix, check.names = FALSE)


gene_count <- nrow(expression_matrix)
cell_count <- ncol(expression_matrix)
info <- paste0("表达矩阵行数（基因数）: ", gene_count, 
               "\n表达矩阵列数（细胞数）: ", cell_count, "\n")
cat(info) 
writeLines(info, "matrix_info.txt")

expression_df <- data.frame(gene = rownames(expression_matrix), expression_matrix, check.names = FALSE)
colnames(expression_df)[1] <- "gene"
write.csv(expression_df, file = "expression_matrix.csv", quote = FALSE, row.names = FALSE)
message("Process completed successfully!")

####10.2Single-cell####

library(Seurat)       
library(dplyr)        
library(ggplot2)    
library(magrittr)    
library(RColorBrewer)
library(limma)      
library(tidyr)       
library(NMF)         
library(CellChat)     
library(ggalluvial)  
library(svglite)      
library(celldex)     
library(SingleR)     
library(monocle)     
library(patchwork)
library(Matrix)
library(monocle3)


workDir <- ""  
setwd(workDir)  

logFC_filter       <- 1        
p_adj_filter       <- 0.05     
min_cells_gene     <- 5       
min_genes_per_cell <- 200      
post_filter_cells  <- 300      
post_filter_mito   <- 20      
n_top_var_features <- 2500    
n_pcs              <- 22      
n_topmarker_heat   <- 10       
cluster_resolution <- 0.6      
neighbor_dims      <- 15      
qc_vln_width       <- 15       
qc_vln_height      <- 7       
input_expr_file    <- "expression_matrix.csv" 
output_dir         <- "analysis_results"     

target_genes <- c("AC026356.1", "AL512598.1", "CYTOR", "MIR210HG")   

if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)  
cat("输出目录设置为：", output_dir, "\n")  

cat("Step 1: 读取表达数据...\n")  
if(!file.exists(input_expr_file)) stop("表达数据文件不存在！")  
expr_table <- read.csv(input_expr_file, header=TRUE, sep=",", stringsAsFactors=FALSE, check.names=FALSE) 
if(nrow(expr_table)==0) stop("表达矩阵无内容！") 
rownames(expr_table) <- expr_table[,1]  
gene_names <- expr_table[,1]             
expr_values <- as.matrix(expr_table[,-1]) 
expr_matrix <- expr_values                

cat("Step 2: 创建Seurat对象...\n") # 流程提示
scObject <- CreateSeuratObject(
  counts = as.matrix(expr_matrix),         
  project = "scRNAseqProject",             
  min.cells = min_cells_gene,              
  min.features = min_genes_per_cell,       
  names.delim = "_"                       
)
scObject[["percent.mito"]] <- PercentageFeatureSet(scObject, pattern = "^MT-") 
pdf(file.path(output_dir, "QC_violin_basicMetrics_pubstyle.pdf"), width=12, height=6)
p <- VlnPlot(
  scObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  theme_classic(base_size = 16) +  
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 12, face = "bold"),
    strip.text   = element_text(size = 15, face = "bold", color = "black")
  )
print(p)  
dev.off()


cat("Step 3: 二次细胞质控...\n")
cell_num_before <- ncol(scObject) 
scObject <- subset(scObject, subset = nFeature_RNA > post_filter_cells & percent.mito < post_filter_mito) 
cell_num_after <- ncol(scObject) # 过滤后细胞数
cat(sprintf("细胞过滤: 原%d，剩%d。\n", cell_num_before, cell_num_after)) 
if(cell_num_after < 10) stop("过滤后过少细胞，检查阈值！") 

pdf(file.path(output_dir, "QC_scatter_metrics.pdf"), width=13, height=7) 
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 1.5)  
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)   
CombinePlots(plots = list(plot1, plot2)) 
dev.off()

cat("Step 4: 归一化+高变基因...\n") 
scObject <- NormalizeData(scObject, normalization.method="LogNormalize", scale.factor=10000)
scObject <- FindVariableFeatures(scObject, selection.method="vst", nfeatures=n_top_var_features) 
pdf(file.path(output_dir, "VarGenes_overview.pdf"), width=10, height=6) 
VariableFeaturePlot(scObject)
dev.off()

cat("Step 5: PCA降维...\n") 
scObject <- ScaleData(scObject) 
scObject <- RunPCA(scObject, npcs=n_pcs, features=VariableFeatures(scObject)) 
pdf(file.path(output_dir, "PCA_geneLoadings.pdf")) 
VizDimLoadings(scObject, dims=1:3, reduction="pca", nfeatures=25)
dev.off()
pdf(file.path(output_dir, "PCA_heatmap_topGenes.pdf"), width=8, height=7) 
DimHeatmap(scObject, dims=1:3, cells=400, balanced=TRUE)
dev.off()

cat("Step 6: 聚类与UMAP降维...\n")
scObject <- FindNeighbors(scObject, dims=1:neighbor_dims)      
scObject <- FindClusters(scObject, resolution=cluster_resolution) 
scObject <- RunUMAP(scObject, dims=1:neighbor_dims)            
nColors <- length(unique(scObject$seurat_clusters))             
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors) 
pdf(file.path(output_dir, "UMAP_clustered_samples.pdf"), width=7, height=5)
DimPlot(scObject, reduction="umap", label=TRUE, cols=myColors) 
dev.off()
write.csv(data.frame(Cell=colnames(scObject), Cluster=scObject$seurat_clusters), 
          file=file.path(output_dir, "CellCluster_UMAP_assignments.csv"), row.names=FALSE) 


cat("Step 7: 差异表达分析...\n") 
marker_df <- tryCatch({
  FindAllMarkers(scObject, only.pos=TRUE, min.pct=0.2, logfc.threshold=logFC_filter) 
}, error=function(e){cat("Marker分析异常\n"); NULL})
if(!is.null(marker_df) && nrow(marker_df)>0){
  sigMarkers <- marker_df[abs(marker_df$avg_log2FC)>logFC_filter & marker_df$p_val_adj<p_adj_filter,]
  write.csv(sigMarkers, file=file.path(output_dir, "Cluster_Markers_DEGs.csv"), row.names=FALSE)
  topmarker <- marker_df %>% group_by(cluster) %>% top_n(n_topmarker_heat, avg_log2FC)
  
  valid_genes <- topmarker$gene[topmarker$gene %in% rownames(scObject)]
  if(length(valid_genes) > 0){
    pdf(file.path(output_dir, "Markers_DoHeatmap.pdf"))
    print(DoHeatmap(scObject, features=valid_genes, size=4) + NoLegend())
    dev.off()
  } else {
    cat("未找到可用于热图的marker基因。\n")
  }
} else {
  cat('未检出显著marker。\n')
}

cat("Step 8: SingleR细胞类型自动注释...\n") 
expr4singler <- GetAssayData(scObject, slot = "data")             
clusters4singler <- scObject$seurat_clusters                     
ref_hpa <- celldex::HumanPrimaryCellAtlasData()                 
singler_result <- SingleR(
  test = expr4singler,
  ref = ref_hpa,
  labels = ref_hpa$label.main,
  clusters = clusters4singler
)  
scObject$SingleR_celltype <- singler_result$labels[as.numeric(scObject$seurat_clusters)+1]  
write.csv(
  data.frame(Cluster=rownames(singler_result), CellType=singler_result$labels),
  file = file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv"), row.names=FALSE
)  
if(!is.null(marker_df) && !is.null(singler_result)) {
  cluster_anno <- data.frame(
    cluster = as.character(rownames(singler_result)),
    CellType = singler_result$labels,
    stringsAsFactors = FALSE
  )  
  marker_df$cluster <- as.character(marker_df$cluster)
  marker_df_anno <- marker_df %>% left_join(cluster_anno, by = "cluster") 
  write.csv(marker_df_anno, file=file.path(output_dir,"Cluster_Markers_DEGs_with_CellType.csv"), row.names=FALSE)
}

nCellTypes <- length(unique(scObject$SingleR_celltype))  
celltypeColors <- colorRampPalette(brewer.pal(12, "Set3"))(nCellTypes) 
pdf(file.path(output_dir, "UMAP_celltype_autoAnnot_Set3.pdf"), width = 8, height = 6)
DimPlot(scObject, group.by = "SingleR_celltype", reduction = "umap", label = TRUE, cols = celltypeColors) 
dev.off()

cat("Step 9: 各种统计可视化输出...\n") 
meta_df <- scObject@meta.data
celltype_stats <- as.data.frame(table(
  Sample = meta_df$orig.ident, 
  CellType = meta_df$SingleR_celltype
))        
celltype_stats <- celltype_stats %>%
  group_by(Sample) %>%
  mutate(Ratio = Freq / sum(Freq)) %>%
  ungroup()  
celltype_stats$Sample <- factor(celltype_stats$Sample, levels = unique(celltype_stats$Sample))
celltype_stats$CellType <- factor(celltype_stats$CellType, levels = unique(celltype_stats$CellType)) 
cell_types <- levels(celltype_stats$CellType)
nColors <- length(cell_types) 
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors) 
names(myColors) <- cell_types 

pdf(file.path(output_dir, "Sample_CellType_Composition_barplot.pdf"), width=7, height=6) 
ggplot(celltype_stats, aes(x = Sample, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, color = "grey70") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent, name = "Proportion") +
  scale_fill_manual(values = myColors, name="Cell Type") +
  labs(x = "Sample") +
  theme_classic(base_size=16) +
  coord_flip() +  
  theme(
    legend.title=element_text(size=14),
    legend.text=element_text(size=13),
    axis.text = element_text(size=13),
    axis.title = element_text(size=14)
  ) 
dev.off()

cell_counts <- as.data.frame(table(Sample = meta_df$orig.ident, CellType = meta_df$SingleR_celltype)) 
write.csv(cell_counts, file.path(output_dir, "CellNumber_perSample_CellType.csv"), row.names = FALSE) 
celltype_matrix <- celltype_stats %>% select(Sample, CellType, Ratio) %>%
  tidyr::spread(CellType, Ratio, fill = 0)
write.csv(celltype_matrix, file.path(output_dir, "CellType_Proportion_Matrix.csv"), row.names=FALSE) 
cell_annot <- data.frame(CellBarcode = rownames(meta_df), Sample = meta_df$orig.ident,
                         Cluster = meta_df$seurat_clusters, CellType = meta_df$SingleR_celltype)
write.csv(cell_annot, file.path(output_dir, "Cell_Full_Annotation.csv"), row.names=FALSE) 

if (exists("marker_df_anno") && nrow(marker_df_anno) > 0) {
  gene_count_by_celltype <- marker_df_anno %>%
    group_by(CellType) %>%
    summarise(gene_number = n_distinct(gene)) %>%
    arrange(desc(gene_number))
  write.csv(gene_count_by_celltype, file = file.path(output_dir, "GeneNumber_per_CellType.csv"), row.names = FALSE)
  pdf(file.path(output_dir, "MarkerGeneCount_perCellType_bar.pdf"), width=8, height=5)
  p_bar <- ggplot(gene_count_by_celltype, aes(x = reorder(CellType, gene_number), y = gene_number, fill=CellType)) +
    geom_bar(stat="identity", width=0.7, color="grey40") +
    scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(nrow(gene_count_by_celltype))) +
    coord_flip() +
    labs(x="Cell Type", y="Number of Marker Genes", title="Number of Marker Genes per Cell Type") +
    theme_classic(base_size=16) +
    theme(legend.position = "none")
  print(p_bar)
  dev.off()
}


top_genes_per_CT <- marker_df_anno %>%
  group_by(CellType) %>% top_n(5, abs(avg_log2FC)) %>% ungroup() # 每类top5 marker
pdf(file.path(output_dir, "MarkerGene_CellType_Alluvial.pdf"), width=9, height=6)
ggplot(top_genes_per_CT, aes(axis1 = CellType, axis2 = gene, y = abs(avg_log2FC))) +
  geom_alluvium(aes(fill=CellType), width=1/12) +
  geom_stratum(width=1/6, fill="grey90", color="black") +
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
  scale_x_discrete(limits = c("CellType","target_genes"), expand = c(.05, .05)) +
  theme_minimal(base_size=14) +
  ggtitle("Top Marker Genes Alluvial Plot")
dev.off()


pca_stdev <- scObject[["pca"]]@stdev
pca_var_explained <- (pca_stdev^2) / sum(pca_stdev^2)
cumulative_var <- cumsum(pca_var_explained)
pca_table <- data.frame(
  PC = seq_along(pca_var_explained),
  Variance_Explained = pca_var_explained,
  Cumulative_Variance = cumulative_var
)
write.csv(pca_table, file = file.path(output_dir, "PCA_Variance_Explained.csv"), row.names=FALSE)

opt_cutoff <- 0.8 
best_n_pcs <- which(cumulative_var >= opt_cutoff)[1]
write.csv(
  data.frame(Recommended_n_PCs = best_n_pcs, Cumulative_Variance = cumulative_var[best_n_pcs]),
  file = file.path(output_dir, "PCA_Recommended_nPCs.csv"),
  row.names=FALSE
)

# 


for (target_gene in target_genes) {
  if(! target_gene %in% rownames(scObject)) {
    cat("警告：指定基因", target_gene, "不在表达矩阵中！\n")
    next
  }
  
  pdf(file.path(output_dir, paste0(target_gene, "_UMAP_FeaturePlot.pdf")), width=7, height=6)
  print(
    FeaturePlot(
      scObject, 
      features = target_gene, 
      reduction = "umap", 
      cols = c("lightgrey", "red"), 
      pt.size = 1,
      label = TRUE
    ) + 
      ggtitle(paste(target_gene, "Expression (UMAP)"))
  )
  dev.off()
  
  pdf(file.path(output_dir, paste0(target_gene, "_VlnPlot_byCluster.pdf")), width=8, height=6)
  p2 <- VlnPlot(
    scObject, 
    features = target_gene, 
    group.by = "seurat_clusters", 
    pt.size = 0,
    slot = "data"
  ) +
    geom_boxplot(width=0.15, outlier.size=0, fill=NA, color="black", lwd=0.7) +
    theme_classic(base_size=16) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
      axis.text.y  = element_text(size = 12, face = "bold"),
      strip.text   = element_text(size = 15, face = "bold", color = "black")
    ) +
    labs(x = NULL, y = "Normalized Expression") +
    ggtitle(paste(target_gene, "by Cluster"))
  print(p2)
  dev.off()
  
  avg_by_celltype <- AverageExpression(scObject, features = target_gene, group.by = "SingleR_celltype")$RNA
  avg_by_cluster  <- AverageExpression(scObject, features = target_gene, group.by = "seurat_clusters")$RNA
  
  write.csv(avg_by_celltype, file.path(output_dir, paste0(target_gene, "_Average_By_CellType.csv")))
  write.csv(avg_by_cluster,  file.path(output_dir, paste0(target_gene, "_Average_By_Cluster.csv")))
  
  gene_expr_vec <- FetchData(scObject, vars = target_gene)
  cell_anno <- data.frame(
    CellBarcode = colnames(scObject),
    CellType = scObject$SingleR_celltype,
    Cluster = scObject$seurat_clusters,
    Gene_Expression = as.vector(gene_expr_vec[,1])
  )
  write.csv(cell_anno, file.path(output_dir, paste0(target_gene, "_Expression_perCell.csv")), row.names=FALSE)
  
  anno_path <- file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv")
  anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
  cluster_col <- grep("cluster", colnames(anno_df), ignore.case = TRUE, value = TRUE)[1]
  celltype_col <- grep("cell.*type", colnames(anno_df), ignore.case = TRUE, value = TRUE)[1]
  if(is.na(cluster_col) | is.na(celltype_col)) stop("注释表中未找到cluster或celltype相关的列，请检查！")
  anno_df[[cluster_col]] <- as.character(anno_df[[cluster_col]])
  
  umap_df <- as.data.frame(Embeddings(scObject, "umap"))
  umap_cols <- grep("UMAP|Dim", colnames(umap_df), value = TRUE, ignore.case = TRUE)
  if(length(umap_cols) >= 2) {
    colnames(umap_df)[match(umap_cols[1:2], colnames(umap_df))] <- c("UMAP_1", "UMAP_2")
  } else if(ncol(umap_df) == 2) {
    colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  } else {
    stop("找不到UMAP坐标列，请检查umap_df的列名！")
  }
  umap_df$Cluster <- as.character(scObject$seurat_clusters)
  umap_df$CellBarcode <- rownames(umap_df)
  gene_found <- grep(paste0("^", target_gene, "$"), rownames(scObject), value = TRUE, ignore.case = TRUE)
  if(length(gene_found) == 0) stop(paste("指定基因", target_gene, "不在表达矩阵中！"))
  gene_to_use <- gene_found[1]
  umap_df$Expression <- FetchData(scObject, vars = gene_to_use)[,1]
  
  umap_df <- left_join(umap_df, anno_df, by = setNames(cluster_col, "Cluster"))
  if(!(celltype_col %in% colnames(umap_df))) stop(paste("合并后未找到细胞类型列", celltype_col, "！"))
  
  cluster_stats <- umap_df %>%
    group_by(Cluster, .data[[celltype_col]]) %>%
    summarise(AvgExpr = mean(Expression), .groups = "drop") %>%
    arrange(Cluster)
  colnames(cluster_stats)[2] <- "CellType"
  
  label_df <- umap_df %>%
    group_by(Cluster, .data[[celltype_col]]) %>%
    summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
  colnames(label_df)[2] <- "CellType"
  
  cluster_levels <- sort(unique(umap_df$Cluster))
  cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(cluster_levels)), cluster_levels)
  
  p_umap <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes(color=factor(Cluster)), size=1, alpha=0.7) +
    scale_color_manual(values=cluster_colors, name="Cluster") +
    geom_text(data=label_df,
              aes(label=paste0("c-", Cluster, "\n", CellType)),
              color="black", size=4, fontface="bold") +
    theme_classic(base_size=16) +
    ggtitle(target_gene) +
    theme(
      plot.title = element_text(size=22, face="bold", hjust=0.5),
      legend.position = "none"
    )
  
  cluster_stats$ClusterLabel <- paste0("c-", cluster_stats$Cluster)
  cluster_stats$ClusterLabel <- factor(cluster_stats$ClusterLabel, levels = paste0("c-", cluster_levels))
  cluster_stats$CellTypeLabel <- cluster_stats$CellType
  max_expr <- max(cluster_stats$AvgExpr)
  offset <- max_expr * 0.03 
  
  cluster_stats$Cluster <- as.character(cluster_stats$Cluster)
  cluster_stats$ClusterLabel <- paste0("c-", cluster_stats$Cluster)
  cluster_stats$ClusterLabel <- factor(cluster_stats$ClusterLabel, levels = paste0("c-", cluster_levels))
  cluster_stats$CellTypeLabel <- cluster_stats$CellType
  
  p_legend <- ggplot(cluster_stats, aes(x=ClusterLabel, y=AvgExpr, fill=Cluster)) +
    geom_bar(stat="identity", width=0.7) +
    scale_fill_manual(values=cluster_colors, guide="none") +
    geom_text(aes(x=ClusterLabel, y=AvgExpr + offset, label=CellTypeLabel),
              hjust=0, vjust=0.5, size=5, fontface="plain", color="black") +
    coord_flip(clip = "off") +
    labs(x=NULL, y="Avg Expression", title="") +
    theme_minimal(base_size=14) +
    theme(
      axis.text.y = element_text(size=13, face="bold"),
      axis.text.x = element_text(size=12),
      axis.title.x = element_text(size=14),
      plot.margin = margin(5, 80, 5, 5)
    )
  
  p_final <- p_umap + p_legend + patchwork::plot_layout(widths = c(2.5, 1.2))
  
  pdf(file.path(output_dir, paste0(target_gene, "_UMAP_withLegend.pdf")), width=12, height=6)
  print(p_final)
  dev.off()
  cat("已输出自定义UMAP表达图：", file.path(output_dir, paste0(target_gene, "_UMAP_withLegend.pdf")), "\n")
}

cat("所有结果/表格/图形均已输出到文件夹：", output_dir, "\n") # 流程结束提示

expr_matrix <- GetAssayData(scObject, assay = "RNA", slot = "counts")
cell_metadata <- scObject@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expr_matrix))
rownames(gene_metadata) <- rownames(expr_matrix)

cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds, reduction_method = "UMAP")


pdf("Trajectory_Pseudotime.pdf", width=7, height=6)
print(plot_cells(cds, color_cells_by = "pseudotime"))
dev.off()

if("SingleR_celltype" %in% colnames(colData(cds))) {
  pdf("Trajectory_CellType.pdf", width=7, height=6)
  print(plot_cells(cds, color_cells_by = "SingleR_celltype"))
  dev.off()
}


for (gene in target_genes) {
  if (!(gene %in% rownames(cds))) {
    cat("警告：基因", gene, "不在cds对象中，跳过。\n")
    next
  }
  
  expr_vec <- Matrix::t(assay(cds)[gene, ]) 
  expr_vec <- as.numeric(expr_vec)
  names(expr_vec) <- colnames(cds)
  
  med <- median(expr_vec)
  expr_group <- ifelse(expr_vec > med, "High", "Low")
  
  colData(cds)[[paste0(gene, "_Group")]] <- factor(expr_group, levels = c("Low", "High"))
  
  pdf(paste0("Trajectory_", gene, "_HighLow.pdf"), width=7, height=6)
  print(
    plot_cells(
      cds,
      color_cells_by = paste0(gene, "_Group"),
      show_trajectory_graph = TRUE,
      label_groups_by_cluster = FALSE,
      label_leaves = TRUE,
      label_branch_points = TRUE,
      graph_label_size = 4,
      cell_size = 1.5,
      reduction_method = "UMAP"
    ) +
      scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
      ggtitle(paste("Trajectory -", gene, "(High/Low by median)"))
  )
  dev.off()
  cat("已输出：Trajectory_", gene, "_HighLow.pdf\n")
}

pseudotime <- pseudotime(cds)
pseudotime <- as.numeric(pseudotime)
names(pseudotime) <- colnames(cds)

for (gene in target_genes) {
  if (!(gene %in% rownames(cds))) {
    cat("警告：基因", gene, "不在cds对象中，跳过。\n")
    next
  }
  
  expr_vec <- Matrix::t(assay(cds)[gene, ])
  expr_vec <- as.numeric(expr_vec)
  names(expr_vec) <- colnames(cds)
  
  med <- median(expr_vec)
  expr_group <- ifelse(expr_vec > med, "High", "Low")
  
  plot_df <- data.frame(
    Pseudotime = pseudotime,
    Group = factor(expr_group, levels = c("Low", "High"))
  )
  plot_df <- plot_df[!is.na(plot_df$Pseudotime), ]  # 去除NA
  
  p <- ggplot(plot_df, aes(x = Pseudotime, fill = Group, color = Group)) +
    geom_density(alpha = 0.35, size = 1.5, adjust = 1.1) +
    scale_fill_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
    scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
    labs(
      title = bquote(.(gene)~" density along pseudotime"),
      x = "Pseudotime",
      y = "Density",
      fill = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(face = "bold", size = 16),
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 16),
      legend.key = element_blank(),
      panel.grid = element_blank()
    )
  
  pdf(paste0("Pseudotime_Density_", gene, "_HighLow.pdf"), width=6, height=5)
  print(p)
  dev.off()
  cat("已输出：Pseudotime_Density_", gene, "_HighLow.pdf\n")
}
umap_df <- as.data.frame(Embeddings(scObject, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$Cell <- rownames(umap_df)
umap_df$CellType <- scObject$SingleR_celltype[umap_df$Cell]

for (gene in target_genes) {
  cat("当前基因：", gene, "\n")
  if (!(gene %in% rownames(scObject))) {
    cat("警告：基因", gene, "不在scObject中，跳过。\n")
    next
  }
  expr_vec <- FetchData(scObject, vars = gene)[,1]
  med <- median(expr_vec)
  expr_group <- ifelse(expr_vec > med, "High", "Low")
  
  plot_df <- umap_df
  plot_df$Group <- factor(expr_group, levels = c("Low", "High"))
  plot_df$Expression <- expr_vec
  
  gene_dir <- file.path(output_dir, paste0(gene, "_Celltype_UMAP"))
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  for (ct in sort(unique(na.omit(plot_df$CellType)))) {
    sub_df <- subset(plot_df, CellType == ct)
    cat("  注释类型：", ct, "，细胞数：", nrow(sub_df), "\n")
    if (nrow(sub_df) == 0) next
    
    sub_df$GroupLabel <- ifelse(sub_df$Group == "High",
                                paste0("High ", gene, " ", ct),
                                paste0("Low ", gene, " ", ct))
    color_map <- setNames(
      c("#4682B4", "#CD2626"),
      c(paste0("Low ", gene, " ", ct), paste0("High ", gene, " ", ct))
    )
    
    p <- ggplot(sub_df, aes(x = UMAP_1, y = UMAP_2, color = GroupLabel)) +
      geom_point(size = 1.2, alpha = 0.8) +
      scale_color_manual(values = color_map) +
      labs(
        title = paste0(gene, " expression in ", ct, " (by median)"),
        color = NULL
      ) +
      theme_classic(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 14),
        axis.line = element_line(size = 1.1),
        axis.ticks = element_line(size = 1.1),
        legend.position = "top",
        legend.text = element_text(face = "bold", size = 14),
        legend.key = element_blank(),
        panel.grid = element_blank()
      )
    
    pdf_path <- file.path(gene_dir, paste0(gene, "_", gsub("[/\\:*?\"<>| ]", "_", ct), "_UMAP.pdf"))
    cat("  输出：", pdf_path, "\n")
    pdf(pdf_path, width=6, height=6)
    print(p)
    dev.off()
  }
  cat("已输出：", gene_dir, "下所有注释细胞类型UMAP图\n")
}

anno_path <- file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv")
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") 

csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  cat("正在处理基因：", gene_name, "\n")
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) {
    cat("文件", csv_file, "缺少 Cluster/CellType/Gene_Expression 列，跳过。\n")
    next
  }
  anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
  cluster_levels <- as.character(anno_df$Cluster)
  cluster_labels <- paste0("c-", anno_df$Cluster, " (", anno_df$CellType, ")")
  names(cluster_labels) <- cluster_levels
  nClusters <- length(cluster_levels)
  cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
  names(cluster_colors) <- cluster_levels
  
  bubble_df <- expr_df %>%
    group_by(Cluster) %>%
    summarise(
      avg_expr = mean(Gene_Expression, na.rm = TRUE),
      pct_pos = mean(Gene_Expression > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  bubble_df$Cluster <- factor(bubble_df$Cluster, levels = cluster_levels)
  
  bubble_df$Y <- " "  
  p <- ggplot(bubble_df, aes(x = Cluster, y = Y, size = pct_pos, color = avg_expr)) +
    geom_point() +
    scale_x_discrete(labels = cluster_labels, expand = expansion(add = 0.5)) +  
    scale_y_discrete(expand = expansion(add = 0.5)) +  
    scale_size(range = c(5, 20), name = "Percent Positive") +
    scale_color_gradientn(colors = c("blue", "white", "red"), name = "Avg Expression") +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold", size = 16),  
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 40)
    ) +
    labs(title = paste0(gene_name, " DotPlot by Cluster (Main Cell Type)"))
  
  pdf_name <- file.path(output_dir, paste0(gene_name, "_DotPlot_byCluster_withMainCellType.pdf"))
  pdf(pdf_name, width=14, height=7)
  print(p)
  dev.off()
  cat("已输出：", pdf_name, "\n")
}

output_dir <- "analysis_results"  
anno_path <- file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv")

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") # 确保列名一致
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)
cluster_labels <- paste0("c-", anno_df$Cluster, " (", anno_df$CellType, ")")
names(cluster_labels) <- cluster_levels
nClusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
names(cluster_colors) <- cluster_levels

csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

all_bubble_df <- data.frame()
for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) next
  expr_df$Gene <- gene_name
  all_bubble_df <- rbind(all_bubble_df, expr_df[, c("Cluster", "Gene_Expression", "Gene")])
}

if (nrow(all_bubble_df) == 0) stop("没有有效的表达数据可用于绘图！")

bubble_df <- all_bubble_df %>%
  group_by(Cluster, Gene) %>%
  summarise(
    avg_expr = mean(Gene_Expression, na.rm = TRUE),
    pct_pos = mean(Gene_Expression > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

bubble_df$Cluster <- factor(bubble_df$Cluster, levels = cluster_levels)
gene_order <- unique(bubble_df$Gene)
bubble_df$Gene <- factor(bubble_df$Gene, levels = gene_order)

bubble_df$Y <- bubble_df$Gene 
p <- ggplot(bubble_df, aes(x = Cluster, y = Y, size = pct_pos, color = avg_expr)) +
  geom_point() +
  scale_x_discrete(labels = cluster_labels, expand = expansion(add = 0.5)) + 
  scale_y_discrete(expand = expansion(add = 0.5)) + 
  scale_size(range = c(5, 20), name = "Percent Positive") +
  scale_color_gradientn(colors = c("blue", "white", "red"), name = "Avg Expression") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 40)  
  ) +
  labs(title = "Gene Expression DotPlot by Cluster (Main Cell Type)")

pdf_name <- file.path(output_dir, "AllGenes_DotPlot_byCluster_withMainCellType.pdf")
pdf(pdf_name, width = 15, height = 4 + 0.7 * length(gene_order))
print(p)
dev.off()
cat("已输出：", pdf_name, "\n")

output_dir <- "analysis_results"  
anno_path <- file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv")

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") 
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)
cluster_labels <- paste0("c-", anno_df$Cluster, "\n", anno_df$CellType)
names(cluster_labels) <- cluster_levels
nClusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
names(cluster_colors) <- cluster_levels

csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  cat("正在处理基因：", gene_name, "\n")
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) {
    cat("文件", csv_file, "缺少 Cluster/CellType/Gene_Expression 列，跳过。\n")
    next
  }
  expr_df$Cluster <- factor(expr_df$Cluster, levels = cluster_levels)
  expr_df$ClusterLabel <- factor(expr_df$Cluster, levels = cluster_levels, labels = cluster_labels)
  p_vln <- ggplot(expr_df, aes(x = ClusterLabel, y = Gene_Expression, fill = Cluster)) +
    geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.9) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
    scale_fill_manual(values = cluster_colors) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(
        angle = 45,      
        hjust = 1,      
        vjust = 1,       
        face = "bold",
        size = 12       
      ),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 16),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 60, l = 40)  
    ) +
    labs(y = "Expression", title = paste0(gene_name, " Expression by Cluster (Main Cell Type)"))
  
  
  pdf_name_vln <- file.path(output_dir, paste0(gene_name, "_VlnPlot_byCluster_withMainCellType.pdf"))
  pdf(pdf_name_vln, width = 14, height = 7)
  print(p_vln)
  dev.off()
  cat("已输出：", pdf_name_vln, "\n")
}

output_dir <- "analysis_results" 
anno_path <- file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv")

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") 
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)

cluster_labels <- paste0("c-", anno_df$Cluster, "\n", anno_df$CellType)
names(cluster_labels) <- cluster_levels
nClusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
names(cluster_colors) <- cluster_levels

csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

all_vln_df <- data.frame()
for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) next
  expr_df$Gene <- gene_name
  all_vln_df <- rbind(all_vln_df, expr_df[, c("Cluster", "Gene_Expression", "Gene")])
}

if (nrow(all_vln_df) == 0) stop("没有有效的表达数据可用于绘图！")

all_vln_df$Cluster <- factor(all_vln_df$Cluster, levels = cluster_levels)
all_vln_df$ClusterLabel <- factor(all_vln_df$Cluster, levels = cluster_levels, labels = cluster_labels)
gene_order <- unique(all_vln_df$Gene)
all_vln_df$Gene <- factor(all_vln_df$Gene, levels = gene_order)

p_vln_merge <- ggplot(all_vln_df, aes(x = ClusterLabel, y = Gene_Expression, fill = Cluster)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.9) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(
      angle = 45,     
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 60, l = 40)
  ) +
  labs(y = "Expression", title = "All Genes Expression by Cluster (Main Cell Type)")

pdf_name_vln_merge <- file.path(output_dir, "AllGenes_VlnPlot_byCluster_withMainCellType.pdf")
pdf(pdf_name_vln_merge, width = 14, height = 4 + 2 * length(gene_order))
print(p_vln_merge)
dev.off()
cat("已输出：", pdf_name_vln_merge, "\n")

output_dir <- "analysis_results"
anno_path <- file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv")
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType")
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)

pseudotime_vec <- as.numeric(monocle3::pseudotime(cds))
names(pseudotime_vec) <- colnames(cds)

for (gene in target_genes) {
  if (!(gene %in% rownames(cds))) {
    cat("警告：基因", gene, "不在cds对象中，跳过。\n")
    next
  }
  gene_dir <- file.path(output_dir, paste0(gene, "_Trajectory_ClusterGroup_HighLow"))
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  cell_meta <- as.data.frame(colData(cds))
  if (!"seurat_clusters" %in% colnames(cell_meta)) {
    stop("cds对象中没有seurat_clusters列，请检查聚类信息是否写入。")
  }
  cell_meta$Cluster <- as.character(cell_meta$seurat_clusters)
  
  for (i in seq_len(nrow(anno_df))) {
    cluster_id <- as.character(anno_df$Cluster[i])
    celltype <- as.character(anno_df$CellType[i])
    cells_in_cluster <- rownames(cell_meta)[cell_meta$Cluster == cluster_id]
    n_cells <- length(cells_in_cluster)
    if (n_cells == 0) next
    
    expr_vals <- assay(cds)[gene, cells_in_cluster]
    med <- median(expr_vals, na.rm = TRUE)
    expr_group <- ifelse(expr_vals > med, "High", "Low")
    names(expr_group) <- cells_in_cluster
    
    highlight_vec <- rep("Other", nrow(cell_meta))
    names(highlight_vec) <- rownames(cell_meta)
    highlight_vec[cells_in_cluster] <- expr_group
    colData(cds)$ExprGroup_HighLow <- factor(highlight_vec, levels = c("Other", "Low", "High"))
    
    p_traj <- plot_cells(
      cds,
      color_cells_by = "ExprGroup_HighLow",
      show_trajectory_graph = TRUE,
      label_groups_by_cluster = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      cell_size = 1.5,
      reduction_method = "UMAP"
    ) +
      scale_color_manual(
        values = c("Other" = "grey80", "Low" = "#079EDF", "High" = "#D377A9"),
        name = "Expression Group",
        labels = c("Other", "Low (Blue)", "High (Red)")
      ) +
      labs(
        title = paste0("Gene: ", gene, " | Cluster: ", cluster_id),
        subtitle = paste0("CellType: ", celltype, " | n=", n_cells)
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "top"
      )
    
    pdf_path_traj <- file.path(
      gene_dir,
      paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_Trajectory_HighLow.pdf")
    )
    pdf(pdf_path_traj, width = 7, height = 6)
    print(p_traj)
    dev.off()
    
    plot_df_expr <- data.frame(
      Expression = expr_vals,
      ExprGroup = factor(expr_group, levels = c("Low", "High"))
    )
    p_expr_density <- ggplot(plot_df_expr, aes(x = Expression, fill = ExprGroup, color = ExprGroup)) +
      geom_density(alpha = 0.35, size = 1.5, adjust = 1.1) +
      scale_fill_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
      scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
      labs(
        title = paste0("Gene: ", gene, " | Cluster: ", cluster_id),
        subtitle = paste0("CellType: ", celltype, " | n=", n_cells, " | AvgExpr=", signif(mean(expr_vals, na.rm=TRUE), 3)),
        x = "Expression",
        y = "Density",
        fill = NULL,
        color = NULL
      ) +
      theme_classic(base_size = 18) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 14)
      )
    pdf_path_expr_density <- file.path(
      gene_dir,
      paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_ExpressionDensity_HighLow.pdf")
    )
    pdf(pdf_path_expr_density, width = 7, height = 5)
    print(p_expr_density)
    dev.off()
    
    pseudotime_vec_cluster <- pseudotime_vec[cells_in_cluster]
    valid_idx <- which(!is.na(pseudotime_vec_cluster) & !is.na(expr_vals))
    if (length(valid_idx) > 0) {
      pt_valid <- pseudotime_vec_cluster[valid_idx]
      expr_group_valid <- expr_group[valid_idx]
      expr_vals_valid <- expr_vals[valid_idx]
      plot_df_pt <- data.frame(
        Pseudotime = pt_valid,
        ExprGroup = factor(expr_group_valid, levels = c("Low", "High")),
        Expression = expr_vals_valid
      )
      p_pt_density <- ggplot(plot_df_pt, aes(x = Pseudotime, fill = ExprGroup, color = ExprGroup)) +
        geom_density(alpha = 0.35, size = 1.5, adjust = 1.1) +
        scale_fill_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
        scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
        labs(
          title = paste0("Gene: ", gene, " | Cluster: ", cluster_id),
          subtitle = paste0("CellType: ", celltype, " | n=", n_cells, " | Pseudotime cells: ", length(valid_idx)),
          x = "Pseudotime",
          y = "Density",
          fill = NULL,
          color = NULL
        ) +
        theme_classic(base_size = 18) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(face = "bold", size = 14)
        )
      pdf_path_pt_density <- file.path(
        gene_dir,
        paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_PseudotimeDensity_HighLow.pdf")
      )
      pdf(pdf_path_pt_density, width = 7, height = 5)
      print(p_pt_density)
      dev.off()
      
      expr_table_path <- file.path(
        gene_dir,
        paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_Pseudotime_Expression.csv")
      )
      write.csv(
        data.frame(Cell = cells_in_cluster[valid_idx],
                   Pseudotime = pt_valid,
                   Expression = expr_vals_valid,
                   ExprGroup = expr_group_valid),
        file = expr_table_path,
        row.names = FALSE
      )
    }
  }
  cat("已输出：", gene_dir, "yes\n")
}

####11.1AC026356.1 cli.heatmat####

library(limma)
library(ComplexHeatmap)
setwd("")

data=read.table("AC026356.1.txt", header=T, sep="\t", check.names=F,row.names = 1)

Type = data$risk
Type=factor(Type, levels=c("low","high"))
data=cbind(as.data.frame(data), Type)

data=data[,c("AC026356.1","Type")]
data=data[order(data[,"AC026356.1"]),] 

cli=read.table("clinical.3.txt", header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,"Type",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

sigVec=c("AC026356.1")
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("Type", clinical)]
  colnames(data)=c("Type", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))
}
colnames(rt)=sigVec

bioCol=c("#BC3C29", "#0072B5","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
colorList[["AC026356.1"]]=c("low"="#0072B5", "high"="#BC3C29")
j=0
for(cli in colnames(rt[,2:ncol(rt)])){
  cliLength=length(levels(factor(rt[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(rt[,cli]))
  cliCol["unknow"]="grey75"
  colorList[[cli]]=cliCol
}

ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

pdf(file="cli.heatmap.pdf", width=8, height=5)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

####11.2AC026356.1 survival####

library(survival)       
library(survminer)       
library(timeROC)         
library(ggplot2)         
library(viridis)       
library(ggsci)          
library(MetBrewer)      
library(ComplexHeatmap)  
library(circlize)        


setwd("")

rt = read.table("risk.csv", sep=",", header=TRUE, row.names=1, check.names=FALSE)
rt = rt[order(rt$riskScore),]
n = nrow(rt)
riskClass = rt[,"risk"]
lowLength = length(riskClass[riskClass == "low"])
highLength = length(riskClass[riskClass == "high"])

rt$riskScore[rt$riskScore > 10] = 10
rt$Patient = 1:n
rt$RiskGroup = factor(rt$risk, levels = c("low", "high"))

median_score <- median(rt$riskScore)
pdf(file="riskScore.pdf", width = 10, height = 3.5)
print(
  ggplot(rt, aes(x=Patient, y=riskScore, color=RiskGroup)) +
    geom_point(size=2.2, alpha=0.85) +
    scale_color_manual(values = c("low" = "#0072B5", "high" = "#BC3C29")) +
    geom_vline(xintercept=lowLength, linetype="dashed", color="grey50", size=1) +
    geom_hline(yintercept=median_score, linetype="dashed", color="grey50", size=1) +
    annotate("text", 
             x = n/2, 
             y = median_score, 
             label = paste0("Median = ", round(median_score, 2)), 
             vjust = -1, 
             hjust = 0.5, 
             size = 5, 
             color = "grey30") +
    theme_minimal(base_size=16) +
    labs(x="Patients (increasing risk score)", y="Risk score") +
    theme(legend.position="top") +
    guides(color=guide_legend(title="Risk Group"))
)
dev.off()

rt$Status = factor(rt$fustat, levels = c(0,1), labels = c("Alive", "Dead"))
rt$RiskGroup = factor(rt$risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))
risk_colors = c("Low risk" = "#0072B5", "High risk" = "#BC3C29")
status_shapes = c("Alive" = 16, "Dead" = 17)

pdf(file="survStat.pdf", width = 10, height = 3.5)
print(
  ggplot(rt, aes(x=Patient, y=futime, color=RiskGroup, shape=Status)) +
    geom_point(size=3, alpha=0.85) +
    scale_color_manual(values = risk_colors) +
    scale_shape_manual(values = status_shapes) +
    geom_vline(xintercept=length(rt$RiskGroup[rt$RiskGroup=="Low risk"]), 
               linetype="dashed", color="grey50", size=1) +
    theme_minimal(base_size=16) +
    labs(
      x="Patients (increasing risk score)", 
      y="Survival time (years)",
      color="Risk group",
      shape="Status"
    ) +
    theme(
      legend.position="top",
      legend.title=element_text(size=14, face="bold"),
      legend.text=element_text(size=12)
    ) +
    guides(
      color=guide_legend(override.aes = list(size=4)),
      shape=guide_legend(override.aes = list(size=4))
    )
)
dev.off()

expr_cols = sapply(rt, is.numeric)
exclude_cols = c("riskScore", "futime", "fustat", "Patient")
expr_cols[names(expr_cols) %in% exclude_cols] = FALSE
rt1 = log2(rt[, expr_cols] + 0.01)
rt1 = t(rt1) 

annotation = data.frame(type = rt$risk)
rownames(annotation) = rownames(rt)
ha = HeatmapAnnotation(df = annotation, col = list(type = c("low"="#0072B5", "high"="#BC3C29")))

pdf(file="heatmap.pdf", width = 10, height = 3)
Heatmap(rt1,
        name = "log2(Expression)",
        top_annotation = ha,
        col = colorRamp2(
          c(min(rt1, na.rm=TRUE), median(rt1, na.rm=TRUE), max(rt1, na.rm=TRUE)),
          MetBrewer::met.brewer("Hokusai1", 3)
        ),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 3))
dev.off()
data = read.table("risk.csv", header=TRUE, sep=",", check.names=FALSE)
data$risk <- factor(data$risk, levels = c("low", "high"))
fit = survfit(Surv(futime, fustat) ~ risk, data=data)
coxFit = coxph(Surv(futime, fustat) ~ risk, data=data)
hr = summary(coxFit)$coefficients[1,2]
hr_confint = summary(coxFit)$conf.int[1, c("lower .95", "upper .95")]
hr_p = summary(coxFit)$coefficients[1,5]
hr_text = paste0("HR = ", sprintf("%.2f", hr), " (95% CI: ", 
                 sprintf("%.2f", hr_confint[1]), "-", sprintf("%.2f", hr_confint[2]), ")")
if (hr_p < 0.001) {
  p_text = "p < 0.001"
} else {
  p_text = paste0("p = ", sprintf("%.3f", hr_p))
}
show_text = paste(hr_text, p_text, sep = "\n")
surPlot = ggsurvplot(
  fit, 
  data = data, 
  conf.int = TRUE, 
  pval = show_text, 
  pval.size = 5, 
  risk.table = TRUE, 
  legend.labs = c("Low risk", "High risk"), 
  legend.title = "Risk", 
  xlab = "Time (years)", 
  break.time.by = 1, 
  risk.table.title = "", 
  palette = c("#0072B5", "#BC3C29"), 
  risk.table.height = .25,
  ggtheme = theme_bw(base_size = 18)
)
pdf(file = "survival.pdf", onefile = FALSE, width = 9, height = 7)
print(surPlot)
dev.off()

inputFile = "risk.csv"
rocFile = "ROC.pdf"
rt = read.table(inputFile, header=TRUE, sep=",")
ROC_rt = timeROC(T=rt$futime, delta=rt$fustat, 
                 marker=rt$riskScore, cause=1, 
                 weighting='aalen', 
                 times=c(1, 3, 5), ROC=TRUE)

roc_df = data.frame(
  FPR = c(ROC_rt$FP[,1], ROC_rt$FP[,2], ROC_rt$FP[,3]),
  TPR = c(ROC_rt$TP[,1], ROC_rt$TP[,2], ROC_rt$TP[,3]),
  Time = factor(rep(c("1 year", "3 years", "5 years"), each=nrow(ROC_rt$FP)))
)
auc_labels = paste0("AUC at ", c(1,3,5), " years: ", round(ROC_rt$AUC, 3))

pdf(file=rocFile, width=5, height=5)
print(
  ggplot(roc_df, aes(x=FPR, y=TPR, color=Time)) +
    geom_line(size=2) +
    scale_color_manual(values=MetBrewer::met.brewer("Hokusai1", 3)) +
    theme_minimal(base_size=16) +
    labs(x="False Positive Rate", y="True Positive Rate", color="Time") +
    annotate("text", x=0.3, y=seq(0.2,0.05,length.out=3), label=auc_labels, hjust=0, size=5, color=MetBrewer::met.brewer("Hokusai1", 3)) +
    theme(legend.position="top")
)
dev.off()

write.csv(rt, file="risk_group_info.csv")

cox_summary <- summary(coxFit)
cox_table <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[,2],
  Lower95 = cox_summary$conf.int[,3],
  Upper95 = cox_summary$conf.int[,4],
  p.value = cox_summary$coefficients[,5]
)
write.csv(cox_table, file="cox_result.csv", row.names=FALSE)

roc_auc_df <- data.frame(
  Time = c(1, 3, 5),
  AUC = ROC_rt$AUC
)
write.csv(roc_auc_df, file="roc_auc.csv", row.names=FALSE)
write.csv(rt1, file="log2_expr_matrix.csv")

####11.3AC026356.1 Independent prognostic factors####

library(survival)  
setwd("")    

bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file=forestFile, width=6.6, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.8)
  axis(1)
  dev.off()
}

indep=function(riskFile=null, cliFile=null, project=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)  
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)    
  
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  uniCoxFile=paste0(project,".uniCox.txt")
  uniCoxPdf=paste0(project,".uniCox.pdf")
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="#0072B5")
  
  
  multiCoxFile=paste0(project,".multiCox.txt")
  multiCoxPdf=paste0(project,".multiCox.pdf")
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
  bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="#BC3C29")
}

indep(riskFile="AC026356.1.txt", cliFile="clinical.4.txt", project="all")

####12.1hsa-mir-126 expression####

library("limma")
library("ggpubr")

setwd("")  

data=read.table("miR.txt", header=T, sep="\t", check.names=F,row.names = 1)
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
normal=data[,group!=0]
tumor=data[,group==0]


gene = "hsa-mir-126"

normal=rbind(normal,gene=normal[gene,])
normal=as.matrix(t(normal[c("gene",gene),]))
rownames(normal)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(normal))

tumor=rbind(tumor,gene=tumor[gene,])
tumor=as.matrix(t(tumor[c("gene",gene),]))
rownames(tumor)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(tumor))

samSample=intersect(row.names(normal),row.names(tumor))
normal=normal[samSample,]
tumor=tumor[samSample,]
data=cbind(normal,tumor)
data=as.data.frame(data[,c(1,3)])
data = log2(data+1)
colnames(data)=c("Normal","Tumor")

pdf(file="pairDiff.hsa-mir-126.pdf",width=5.5,height=5)
ggpaired(data, cond1 = "Normal", cond2 = "Tumor",fill = c("#0072B5", "#BC3C29"), palette = "jco",
         xlab="",ylab = paste0(gene," expression"))+
  stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.35)
dev.off()

####13.1GEO normalize####

library(limma)

inputFile="geneMatrix.txt"  
conFile="S1.txt"            
treatFile="S2.txt"             
geoID="GSE14520"            
setwd("\\GSE14520")      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

inputFile="geneMatrix.txt"  
conFile="S1.txt"            
treatFile="S2.txt"             
geoID="GSE17856"            
setwd("\\GSE17856")      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

inputFile="geneMatrix.txt"  
conFile="S1.txt"            
treatFile="S2.txt"             
geoID="GSE54236"            
setwd("\\GSE54236")      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

inputFile="geneMatrix.txt"  
conFile="S1.txt"            
treatFile="S2.txt"             
geoID="GSE76427"            
setwd("\\GSE76427")      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

inputFile="geneMatrix.txt"  
conFile="S1.txt"            
treatFile="S2.txt"             
geoID="GSE121248"            
setwd("\\GSE121248")      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

inputFile="geneMatrix.txt"  
conFile="S1.txt"            
treatFile="S2.txt"             
geoID="GSE174570"            
setwd("\\GSE174570")      

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

####13.2GEO expression####

library(limma)
library(reshape2)
library(ggpubr)
library(pROC)

expFile="GSE174570.normalize.txt"   
geneFile="TargetGenes.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
            xlab="", ylab="Gene expression", legend.title="Type",
            palette = c("#0072B5", "#BC3C29"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="GSE174570.boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()


expFile="GSE121248.normalize.txt"   
geneFile="TargetGenes.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
            xlab="", ylab="Gene expression", legend.title="Type",
            palette = c("#0072B5", "#BC3C29"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="GSE121248.boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()


expFile="GSE14520.normalize.txt"   
geneFile="TargetGenes.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
            xlab="", ylab="Gene expression", legend.title="Type",
            palette = c("#0072B5", "#BC3C29"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="GSE14520.boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()


expFile="GSE17856.normalize.txt"   
geneFile="TargetGenes.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
            xlab="", ylab="Gene expression", legend.title="Type",
            palette = c("#0072B5", "#BC3C29"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="GSE17856.boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()


expFile="GSE54236.normalize.txt"   
geneFile="TargetGenes.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
            xlab="", ylab="Gene expression", legend.title="Type",
            palette = c("#0072B5", "#BC3C29"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="GSE54236.boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()


expFile="GSE76427.normalize.txt"   
geneFile="TargetGenes.txt"   
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
            xlab="", ylab="Gene expression", legend.title="Type",
            palette = c("#0072B5", "#BC3C29"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="GSE76427.boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()

