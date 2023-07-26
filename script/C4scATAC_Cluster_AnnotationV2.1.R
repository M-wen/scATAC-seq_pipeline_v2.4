### Get the parameters
parser = argparse::ArgumentParser(description="Script to clustering and celltype annotation scATAC data")
parser$add_argument('-I','--countpeak', help='input peak count matrix')
parser$add_argument('-P','--countpromoter', help='input promoter count matrix')
parser$add_argument('-G','--promoter', help='input promoter bed file')
parser$add_argument('-Q','--qc', help='input qc file')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-C','--chrmt', help='chrmt')
args = parser$parse_args()

###library packages
library(data.table)
library(Matrix)
library(Signac)
library(Seurat)
library(scMCA)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(scMCA)
library(stringr)
library(IRanges)
library(GenomicRanges)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a=getPalette(8)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
b=getPalette(8)

setwd(args$out)

raw_count=as.data.frame(fread(args$countpeak))
raw_meta=read.table(args$qc,header=T,comment.char="")

### 
#peak=as.data.frame(raw_count$peak)
#peak$tag1=unlist(lapply(strsplit(as.character(peak[,1]),"_"),"[",1)) #Feb 21,2022 modified 
#peak$tag2=unlist(lapply(strsplit(as.character(peak[,1]),"_"),"[",2))
#peak$tag3=unlist(lapply(strsplit(as.character(peak[,1]),"_"),"[",3))
#peak$Tag1=paste(peak$tag1,peak$tag2,sep = ":")
#peak$Tag=paste(peak$Tag1,peak$tag3,sep = "-")
rownames(raw_count)=raw_count[,1]
raw_count=raw_count[,-1]
raw_count <- raw_count[grep(paste0("^",args$chrmt),rownames(raw_count),invert=T),]
my_counts=as(as.matrix(raw_count), "dgCMatrix")

###
rownames(raw_meta)=raw_meta[,1]
raw_meta=raw_meta[,-1]
my_meta=raw_meta
colnames(my_meta)[2]="uniqueNuclearFrags" #0604 modified
my_meta$FRIP=colSums(raw_count)/2/my_meta$uniqueNuclearFrags
paste0(100*round(median(my_meta$FRIP),5),"%")
###------------------------------------------------------------------------###

scATAC <- CreateSeuratObject(
  counts = my_counts,
  assay = 'peaks',
  project = 'C4scATAC',
  min.cells = 10,
  meta.data = my_meta
)

###QC
scATAC@meta.data$log10_uniqueFrags=log10(scATAC@meta.data$uniqueNuclearFrags)
scATAC@meta.data$batch="1"
QC_plot=list()
QC_plot[[1]]=VlnPlot(
  object = scATAC,
  features = c('log10_uniqueFrags'),
  pt.size = 0,
  ncol = 3,group.by = "batch") + NoLegend()+
  scale_fill_manual(values = b[1]) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())  
QC_plot[[2]]=VlnPlot(
  object = scATAC,
  features = c('tssProportion'),
  pt.size = 0,
  ncol = 3,group.by = "batch") + NoLegend()+
  scale_fill_manual(values = b[2]) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
QC_plot[[3]]=VlnPlot(
  object = scATAC,
  features = c('FRIP'),
  pt.size = 0,
  ncol = 3,group.by = "batch") + NoLegend()+
  scale_fill_manual(values = b[6]) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

png("plot4_QC.png",width = 523,height = 428)
do.call(grid.arrange,c(QC_plot,ncol=3))
dev.off()

svg("plot4_QC.svg",width = 6,height = 4)
do.call(grid.arrange,c(QC_plot,ncol=3))
dev.off()

###output report file

qc5=data.frame(qc="Fraction of fragments overlapping TSS",num=paste0(100*round(sum(scATAC$tssProportion*scATAC$uniqueNuclearFrags)/sum(scATAC$uniqueNuclearFrags),5),"%"),stringsAsFactors = FALSE)
qc5[2,1]="Called peak number"
qc5[2,2]=prettyNum(nrow(raw_count),big.mark = ",")
qc5[3,1]="Fraction of fragments overlapping called peaks"
qc5[3,2]=paste0(100*round(sum(scATAC$FRIP*scATAC$uniqueNuclearFrags)/sum(scATAC$uniqueNuclearFrags),5),"%")
qc5[4,1]="Percent duplicates"
qc5[4,2]=paste0(100*round(1-(sum(scATAC$uniqueNuclearFrags)/sum(scATAC$totalFrags)),5),"%")
write.table(qc5,"5.library.QC.csv",sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)

qc1=data.frame(qc="Estimated number of cells",num=prettyNum(ncol(raw_count),big.mark = ","),stringsAsFactors = FALSE)
qc1[2,1]="Median fragments per cell"
qc1[2,2]=prettyNum(median(scATAC$uniqueNuclearFrags),big.mark = ",")
qc1[3,1]="Median fraction of fragments overlapping peaks"
qc1[3,2]=paste0(100*round(median(scATAC$FRIP),5),"%")
qc1[4,1]="Median fraction of fragments overlapping TSSs"
qc1[4,2]=paste0(100*round(median(scATAC$tssProportion),5),"%")
write.table(qc1,"1.cell_report.csv",sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)

###Number of beads per droplet

name=as.data.frame(rownames(my_meta))
colnames(name)="DropBarcode"
name$Num=unlist(lapply(strsplit(as.character(name$DropBarcode),"N0"),"[",2))
table=as.data.frame(table(name$Num))
colnames(table)=c("Num","Count")

png("plot3_DropBeadsnum.png",width = 458,height = 377)
p=ggplot(data = table, mapping = aes(x = factor(Num), y = Count, fill = Num)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = 'Set2',labels = paste(levels(table$Num)," ",table$Count))+
  theme_bw()+
  xlab("Number of beads per droplet")+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  ggtitle(paste("Total cell number",nrow(my_meta)))
print(p)
dev.off()

svg("plot3_DropBeadsnum.svg",width = 5,height = 4)
print(p)
dev.off()

###Cluster

scATAC <- subset(scATAC, subset = uniqueNuclearFrags > 3 & tssProportion > 0.1)
scATAC <- RunTFIDF(scATAC)
scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q0')
scATAC <- RunSVD(
  object = scATAC,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 1:30)
scATAC <- FindNeighbors(object = scATAC, reduction = 'lsi', dims = 1:30)
scATAC <- FindClusters(object = scATAC, verbose = FALSE)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a=getPalette(length(unique(scATAC@meta.data$seurat_clusters)))

c1=DimPlot(object = scATAC, label = TRUE) + NoLegend()+
  scale_color_manual(values = a)

png("plot7_Cluster_peak.png",width = 524,height = 488)
print(c1)
dev.off()

svg("plot7_Cluster_peak.svg",width = 5,height = 5)
print(c1)
dev.off()

c2=FeaturePlot(scATAC,features = "log10_uniqueFrags")+
  scale_color_viridis(direction = -1,option = "D",name="log10 Frags")+
  ggtitle("")

png("plot8_Cluster_depth.png",width = 524,height = 488)
print(c2)
dev.off()

svg("plot8_Cluster_depth.svg",width = 6,height = 5)
print(c2)
dev.off()

###Cell annotation
input=as.data.frame(fread(args$countpromoter))
tss=read.table(args$promoter,header = F)
tss$Name=unlist(lapply(strsplit(as.character(tss$V4),":"),"[",1))
tss$Name=make.unique(tss$Name)
tss$Name=str_to_title(tss$Name)
rownames(input)=make.unique(tss$Name)
input=input[,-1]
scATAC_result <- scMCA(scdata = input, numbers_plot = 3)
if (length(scATAC_result$scMCA) != 0){
  out=as.data.frame(unlist(scATAC_result$scMCA))
  out$`unlist(scATAC_result$scMCA)`=as.character(out$`unlist(scATAC_result$scMCA)`)
  
  scATAC@meta.data$cell_type=out[match(rownames(scATAC@meta.data),rownames(out)),1]
  
  out_meta=scATAC@meta.data
  
  table_list=list()
  
  for(i in 0:(length(unique(scATAC$seurat_clusters))-1)){
    a=i+1
    sub=subset(out_meta,out_meta$seurat_clusters==i)
    tab=as.data.frame(table(sub$cell_type))
    tab_order=tab[order(tab[,2],decreasing = T),]
    tab_order$Cluster=i
    table_list[[a]]=tab_order[1,]
  }
  
  sum=do.call(rbind,table_list)
  use=out_meta[,c(15,16)]
  use$seurat_clusters=as.character(use$seurat_clusters)
  use$ID=rownames(out_meta)
  colnames(sum)=c("predicated.cell.type","Freq","seurat_clusters")
  res=merge(use,sum,by="seurat_clusters")
  scATAC@meta.data$predicated.cell.type=res[match(rownames(out_meta),res$ID),4]
  
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  a=getPalette(length(unique(scATAC@meta.data$predicated.cell.type)))
  
  c3=DimPlot(object = scATAC, label = FALSE, group.by = "predicated.cell.type") +
    scale_color_manual(values = a)
  
  png("plot9_Cluster_annotation.png",width = 990,height = 479)
  print(c3)
  dev.off()
  
  svg("plot9_Cluster_annotation.svg",width = 11,height = 5)
  print(c3)
  dev.off()
  
  table_cell.type=as.data.frame(table(as.character(scATAC@meta.data$predicated.cell.type)))
  colnames(table_cell.type)=c("predicated.cell.type","number")
  table_cell.type$ratio=table_cell.type$number/colSums(table_cell.type[2])
  order_table_cell.type=table_cell.type[order(table_cell.type$number,decreasing = T),]
  write.table(order_table_cell.type,"Table2.celltypecount.csv",sep = ",",quote = FALSE,row.names = FALSE)

} else{
  table_cell.type=data.frame(predicated.cell.type="None",number=dim(scATAC)[2],ratio=1)
  write.table(table_cell.type,"Table2.celltypecount.csv",sep = ",",quote = FALSE,row.names = FALSE)
}


###Find differentially accessible peaks between clusters
da_peaks=FindAllMarkers(scATAC,
                        min.pct = 0.2,
                        test.use = 'LR')
sub_da_peak=subset(da_peaks,da_peaks$p_val_adj<0.01 & da_peaks$avg_log2FC > log(exp(1)))
sub_da_peak$chr=unlist(lapply(strsplit(sub_da_peak$gene,":"),"[",1))
sub_da_peak$tag=unlist(lapply(strsplit(sub_da_peak$gene,":"),"[",2))
sub_da_peak$start=unlist(lapply(strsplit(sub_da_peak$tag,"-"),"[",1))
sub_da_peak$end=unlist(lapply(strsplit(sub_da_peak$tag,"-"),"[",2))
sub_da_peak=sub_da_peak[,-9]
use=sub_da_peak[,c(8,9,10,7)]
colnames(use)[4]="peak ID"
tss=tss[,c(1:3,5)]
colnames(tss)=c("chr","start","end","gene symbol")
GRange_TSS=makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)
GRange_Peak=makeGRangesFromDataFrame(use, keep.extra.columns = TRUE)
distance=distanceToNearest(GRange_Peak,GRange_TSS)
out_Peak=as.data.frame(GRange_Peak[queryHits(distance)])[,c(6,4)]
out_TSS=as.data.frame(GRange_TSS[subjectHits(distance)])[6]
summary=cbind(out_Peak,out_TSS)
summary$distance=distance@elementMetadata$distance
sub_da_peak=sub_da_peak[,c(1:7)]
colnames(sub_da_peak)[7]="peak.ID"
result_da_peak=unique(merge(summary,sub_da_peak,by="peak.ID"))
result_da_peak=result_da_peak[order(result_da_peak$cluster),]
write.table(result_da_peak,"Table1.dapeak.csv",sep = ",",quote = FALSE,row.names = FALSE)
