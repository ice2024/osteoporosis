###OP###
geo=read.delim("GSE56815_series_matrix.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
anno=read.delim("GPL96-57554.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
geo1=merge(anno,geo,by="ID")
length(unique(geo1$Gene))
geo1$ID=NULL
expr.matrix <- aggregate(geo1[, -1], by = list( symbols = (geo1[, 1])), mean) 
write.csv(expr.matrix,"GSE56815_exp.csv",row.names = FALSE)
#
ExprMatrix=read.csv("GSE59867_exp.csv",header=T,row.names=1,check.names=FALSE)
fenzu=read.table("group.GSE59867.txt",header=T,row.names=1,stringsAsFactors = FALSE)
all(colnames(ExprMatrix)==row.names(fenzu))
library(limma)
group=factor(fenzu$group)
design <- model.matrix(~-1+group)
colnames(design)<-levels(group)
contrasts <- makeContrasts(OP-Control,levels=design)
fit<-lmFit(ExprMatrix, design)
fit1<-contrasts.fit(fit, contrasts)
fit2<-eBayes(fit1)
diff.matrix <-topTable(fit2,coef=1,number= nrow(ExprMatrix), adjust.method="BH",sort.by="B",resort.by="M")
write.csv(diff.matrix, file = "GSE56815_diff.matrix.csv")
diff.matrix=read.csv("GSE56815_diff.matrix.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
s=which(diff.matrix$adj.P.Val<0.05&abs(diff.matrix$logFC)>0)
s1=diff.matrix[s,]
write.csv(s1, file = "GSE56815_diff.adjp0.05.csv")

library(ggplot2)
library(ggrepel)
data=read.csv("GSE56815_diff.matrix.csv",header=T,row.names=1,stringsAsFactors = FALSE)
threshold <- as.factor(ifelse(data$adj.P.Val < 0.05 &
                                
                                abs(data$logFC) >= 0 ,
                              
                              ifelse(data$logFC >= 0 ,'Up','Down'),'Not'))
pdf("volcano.pdf",height=8,width=8,onefile = FALSE)
ggplot(data,aes(x=logFC,y=-log10(adj.P.Val),colour=threshold)) +
  
  xlab("log2(Fold Change)")+ylab("-log10(adj.P.Value)") +
  
  geom_point(size = 2,alpha=1) +
  
  scale_color_manual(values=c("#2b8cbe","grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(1.3), lty = 2,colour="#000000")+
  theme(
    axis.text=element_text(size=20),
    axis.title=element_text(size=20)
  )
dev.off()

uni_matrix=read.delim("GSE56815_exp.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
diff=read.csv("GSE56815_diff.adjp0.05.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
s=which(row.names(uni_matrix)%in%row.names(diff))
diff_exp=uni_matrix[s,]
col_anno=read.table("group.GSE56815.txt",header=T,row.names=1,stringsAsFactors = FALSE)
diff_exp=diff_exp[,row.names(col_anno)]
all(colnames(diff_exp)==row.names(col_anno))
ann_colors = list(group=c(OP="orange",Control="blue"))
library(pheatmap)
pdf("pheatmap.pdf",height=6,width=6,onefile=FALSE)
pheatmap(diff_exp, scale = "row", color = colorRampPalette(c("darkgreen","white","red"))(50),fontsize=9, fontsize_row=9,show_colnames = F,show_rownames = FALSE,annotation_col=col_anno,annotation_colors =ann_colors, annotation_legend=T,annotation_names_col=T,cluster_rows = T,cluster_cols = F)
dev.off()
###WGCNA#
library(WGCNA)  
dataExpr=read.csv("GSE56815_exp.csv",header=T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
group=read.delim("group.wgcna.txt",header=T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
dataExprVar=dataExpr
diff_exp=dataExprVar
mydata <- as.data.frame(t(dataExprVar))
nGenes = ncol(mydata)
nSamples = nrow(mydata)
mytree=hclust(dist(mydata), method="complete") 
pdf(file="1.sampleClustering.pdf",width = 18, height = 9)
par(mar = c(5,5,5,5))
plot(mytree, main = "Sample clustering", sub="", xlab="")
dev.off()
powers = seq(1,20)

sft = pickSoftThreshold(mydata, powerVector=powers, verbose=5)  
pdf(file = "2.softPower.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylim = c(-1,1), ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")
abline(h=0.85,col="blue")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.6,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.6,col="red")
dev.off()
softPower=sft$powerEstimate
adjacency=adjacency(mydata,power=softPower)
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize =70
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize) 

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(mydata, colors = dynamicColors)
MEs = MEList$eigengenes 

MEDiss = 1-cor(MEs)

METree = hclust(as.dist(MEDiss), method = "average")

pdf(file = "4.eigengeneClustering.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "") 
dev.off()

pdf("5.Eigengene_adjacency_heatmap.pdf")
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)
dev.off()

MEDissThres=0.2

merge=mergeCloseModules(mydata,dynamicColors,cutHeight=MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs=merge$newMEs

pdf("6.ModuleTree.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf("7.mergeModuleTree.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree,mergedColors, "Mergeddynamic",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
write.table(paste(colnames(mydata), moduleColors, sep = "\t"), file = "4.netcolor2gene.xls", row.names=FALSE, col.names=FALSE, quote=FALSE)

gene_color=read.table("4.netcolor2gene.xls",row.names=1)
all(row.names(gene_color)==colnames(mydata))
color=unique(gene_color[,1])
color=color[color!="grey"]

traitColors = numbers2colors(group, signed = TRUE,centered=TRUE);
pdf("9.Sample_dendrogram_and_trait_heatmap.pdf",width = 18, height = 9)
plotDendroAndColors(mytree, 
                    traitColors, 
                    groupLabels = names(group), 
                    rowTextAlignment = "right-justified",
                    addTextGuide = TRUE ,
                    hang = 0.03,
                    dendroLabels = NULL, 
                    addGuide = FALSE, 
                    guideHang = 0.05,
                    main = "Sample dendrogram and trait heatmap")
dev.off()

all(colnames(diff_exp)==row.names(group))
moduleTraitCor_noFP <- cor(mergedMEs, group, use = "p");
moduleTraitPvalue_noFP = corPvalueStudent(moduleTraitCor_noFP, nSamples); 
textMatrix_noFP <-paste(signif(moduleTraitCor_noFP, 2), "\n(", signif(moduleTraitPvalue_noFP, 1), ")", sep = ""); 
dim(textMatrix_noFP) = dim(moduleTraitCor_noFP)
pdf("10.Module-trait_relationships.pdf",width = 5, height = 7)
par(mar = c(4, 6, 5, 3)); 
labeledHeatmap(Matrix = moduleTraitCor_noFP, 
               xLabels = colnames(group), 
               yLabels = names(mergedMEs), 
               ySymbols = names(mergedMEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix_noFP,
               setStdMargins = FALSE, 
               cex.text = 0.65, 
               zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off()
###
library(clusterProfiler)
dif_enrich=read.delim("jVenn.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
gene=dif_enrich$node
gene.id <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
egoBP<-enrichGO(gene = gene.id$ENTREZID,  pAdjustMethod ="none",    
                OrgDb = 'org.Hs.eg.db',          
                ont = "ALL",                              
                pvalueCutoff = 0.05, 
                qvalueCutoff=1,
                minGSSize = 2,
                readable= TRUE)  
go.all=egoBP
write.csv(go.all@result,'go.all.result.csv',row.names=F)
dotplot(go.all,showCategory = 10,color = "pvalue", size = NULL,font.size = 10, title = "TOP10 GO enrichment", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
####
library(tidyverse)
library(glmnet)
source('msvmRFE.R')
library(e1071)
library(caret)
library(randomForest)
train <- read.csv("lasso_svm_input.csv",row.names = 1, 
                  as.is = F)
dim(train)
x <- as.matrix(train[,-1])  
y <- ifelse(train$group == "Control", 0,1) 
library(glmnet)
set.seed(2023)
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
pdf("A_lasso_model.pdf", width = 5, height = 5)
plot(fit, xvar = "dev", label = TRUE)
dev.off()
cvfit = cv.glmnet(x, y, 
                  nfold=10,
                  family = "binomial", type.measure = "class")


pdf("A_lasso_cvfit.pdf", width = 6, height = 4)
plot(cvfit)
dev.off()
cvfit$lambda.min
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
write.csv(lasso_fea,"feature_lasso_0.004827212.csv")
predict <- predict(cvfit, newx = x[1:nrow(x),], s = "lambda.min", type = "class")
table(predict,y)
input=train
input$group=ifelse(input$group == "Normal", 0,1)
set.seed(2023)
svmRFE(input, k = 5, halve.above = 100) 
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) 

top.features = WriteFeatures(results, input, save=F) 
head(top.features)
write.csv(top.features,"feature_svm.csv")
featsweep = lapply(1:25, FeatSweep.wrap, results, input) 
save(featsweep,file = "featsweep.RData")
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
#dev.new(width=4, height=4, bg='white')
pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()
#dev.new(width=4, height=4, bg='white')
pdf("B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) 
dev.off()
library(randomForest)
df_train <- read.csv("lasso_input.csv",sep=",",header=T,check.names=F,row.names = 1)
set.seed(2023)
df_train$group <- factor(df_train$group)
train.forest <- randomForest(group ~ ., 
                             data = df_train, 
                             ntree = 500, 
                             importance = TRUE,
                             proximity = TRUE) 
train.forest
pdf("RF_trees.pdf",width = 5,height = 5)
plot(train.forest,main="Random Forest")
dev.off()
summary(train.forest)
importance_otu <- train.forest$importance
head(importance_otu)
importance_otu <- data.frame(importance(train.forest))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy,decreasing = T),]
head(importance_otu)
write.csv(importance_otu,"RF_importance_otu.csv")
pdf("RF_top10.pdf",width = 6.5,height = 6.5)
varImpPlot(train.forest,n.var = min(10, nrow(train.forest$importance)), main = 'Top 10 - variable importance')
dev.off()
####
library(pROC)
df <- read.csv("input.csv",header = T,check.names = F,row.names = 1,stringsAsFactors = FALSE)
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2","red")
pdf("ROC.pdf",height=6,width=6)
auc.out <- c()
x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=FALSE, 
              ci=TRUE,
              col=mycol[2],
              lwd=2, 
              legacy.axes=T)
ci.lower <- round(as.numeric(x$ci[1]),3) 
ci.upper <- round(as.numeric(x$ci[3]),3)

auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)
for (i in 3:ncol(df)){
  x <- plot.roc(df[,1],df[,i],
                add=T,
                smooth=FALSE,
                ci=TRUE,
                col=mycol[i],
                lwd=2,
                legacy.axes=T)
  
  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)
  
  auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
}
legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out[,2],sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:length(df)],
       lwd = 2,
       bty="n")
dev.off()
#######
library(rms)
bc=read.csv("GSE56815_group.csv",header = T,stringsAsFactors = FALSE,row.names = 1)
bc$group=ifelse(bc$group == "Control", 0,1)
bc <- na.omit(bc)
dd <- datadist(bc)
options(datadist="dd")
colnames(bc)
formula1<-as.formula(group~.)
fit1<-lrm(formula1,data = bc,x=T,y=T)
summary(fit1)
nom1<-nomogram(fit1,
               fun=function(x)1/(1+exp(-x)),
               lp=F,
               fun.at = c(0.1,0.5,0.9),
               funlabel = "group")
pdf("nomogram.pdf")
plot(nom1)
dev.off() 

cal1<-calibrate(fit1,method = "boot",B=1000)
pdf("nomogram_predict.pdf")
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Nomogram Predicted Probability", ylab = "Actual Probability")
dev.off()
###
library("GSVA")
##1-ssGSEA##
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
uni_matrix=read.delim("GSE56815_exp.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
gene_immune=read.csv("cell.csv",header=T,stringsAsFactors = FALSE)
list<- split(as.matrix(gene_immune)[,1], gene_immune[,2])
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
immune_score=as.data.frame(t(gsva_matrix))
phenotype_file=read.delim("group.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
immune_score$sample=row.names(immune_score)
imm=merge(phenotype_file,immune_score,by="sample")
row.names(imm)=imm$sample
imm$sample=NULL
write.csv(imm,"immune.cells.score.group.csv")
data=imm
test=as.list(data[,-1])
height<-stack(test)
Group=rep(data$group,ncol(data)-1)     
df=as.data.frame(cbind(Group,height))
colnames(df)=c("group","Infiltration_level","cells")
library(ggplot2)
library(ggpubr)
library(ggsignif)
pdf("immune.cells_boxplot.pdf",height=6,width=14)
ggplot(data=df,aes(x=cells,y=Infiltration_level,fill=group))+
  geom_boxplot(width=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#4dbbd5", "#e64b35"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method = "wilcox.test",
                     hide.ns = TRUE)+
  theme(axis.text.x = element_text(size=12,colour="black",angle=45,hjust = 1),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))

dev.off()

