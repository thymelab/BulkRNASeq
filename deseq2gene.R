library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("dplyr")
library("regionReport")
library("pathview")
library("gage")
library("GenomicAlignments")
library("biomaRt")
library("pheatmap")
library("dplyr")
library("EnhancedVolcano")
library("biomaRt")
library("clusterProfiler")
library("tidyverse")
library("data.table")
library("ReactomePA")
library("org.Dr.eg.db")
library("formattable")
library("stringr")
library("recount")
library("gageData")
library("tidyr")
library(optparse)

option_list = list(
	make_option(c("-d", "--directory"), type="character", default=NULL, 
		help="working directory", metavar="character"),
#	make_option(c("-o", "--out"), type="character", default="out.txt",
#		help="output file name [default= %default]", metavar="character"),
	make_option(c("-p", "--prefix"), type="character", default="defaultname",
		help="prefix name for all graphs [default= %default]", metavar="character"),
	make_option(c("-s", "--sampletable"), type="character", default="defaultname",
		help="sample table csv file [default= %default]", metavar="character")
#	make_option(c("-c", "--mincell"), type="integer", default=3,
#		help="minimum number of cells [default= %default]", metavar="integer"),
#	make_option(c("-d", "--minfeat"), type="integer", default=200, 
#		help="minimum number of molecules [default= %default]", metavar="integer")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#if (is.null(opt$file)){
#	print_help(opt_parser)
#	stop("At least one argument must be supplied (input file)", call.=FALSE)
#}

outputPrefix <- opt$prefix
directory <- opt$directory
sampleTablefile <- opt$sampletable

#sampleTable <- read.csv("input.csv")
sampleTable <- read.csv(sampleTablefile)
# CHANGE BELOW
treatments = c("wt","hom") ##genotype again

ddsHTSeq <-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=directory,
                                      design=~condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels=treatments)
dds <-DESeq(ddsHTSeq)

# filtering might need to change depending on sample number
keep <- rowSums(counts(dds) == 0) < 4
dds <- dds[keep,]

# CHANGE BELOW
res <-results(dds,contrast=c("condition","hom","wt")) ##contrast= sets the order of analysis to hom vs wt. If obmitted, the order will be alphabetical
res<-res[order(res$pvalue),]
rlog<-rlog(dds)
vst <-vst(dds)


gene_name <-read.csv("llgeneid_genename.csv") ##this file has LLgeneID as column 1 and LLgeneAbbrev/gene name as column 2
dataframe_res <-as.data.frame(res)
dataframe_res$LLgeneID<-row.names(dataframe_res)
dataframe_res<-dataframe_res[c(7,1:6)]
res_gene <-inner_join(dataframe_res,gene_name,by="LLgeneID")
res_gene<-res_gene[!is.na(res_gene$padj),]
res_gene_p<-subset(res_gene,padj<0.05)

res_gene_p_up<-subset(res_gene_p,log2FoldChange>0)
res_gene_p_down<-subset(res_gene_p,log2FoldChange<0)


write.csv(res_gene,file=paste0(outputPrefix,"_allresults_wt_hom-with-normalized.csv")) ##this output file is normalized change in expression of all genes
write.csv(res_gene_p,file=paste0(outputPrefix,"_topp_wt_hom-with-normalized.csv")) ##this output file is normalized changed in expression of genes with an adjusted p values of <0.05
write.csv(res_gene_p_up,file=paste0(outputPrefix,"_upreg_topp_wt_hom-with-normalized.csv")) ##this output file is normalized changed in expression of genes with an adjusted p values of <0.05
write.csv(res_gene_p_down,file=paste0(outputPrefix,"_downreg_topp_wt_hom-with-normalized.csv")) ##this output file is normalized changed in expression of genes with an adjusted p values of <0.05



{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("drerio_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoDRerio <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version","entrezgene_id", 
                                     "zfin_id_symbol", "description", "external_gene_name"), mart=ensembl)
}

anno_gene_list_up <- inner_join(res_gene_p_up,annoDRerio,by=c("LLgeneAbbrev" ="external_gene_name")) 

keggPA_up = as.data.frame(enrichKEGG(gene=anno_gene_list_up$entrezgene_id, organism="dre",pvalueCutoff=0.05))
write.csv(keggPA_up,file=paste0(outputPrefix,"_up_pathways.csv")) ##this output file give kegg pathways IDs, names, and statistical values 




{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("drerio_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoDRerio <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version","entrezgene_id", 
                                     "zfin_id_symbol", "description", "external_gene_name"), mart=ensembl)
}

anno_gene_list_down <- inner_join(res_gene_p_down,annoDRerio,by=c("LLgeneAbbrev" ="external_gene_name")) 


keggPA_down = as.data.frame(enrichKEGG(gene=anno_gene_list_down$entrezgene_id, organism="dre",pvalueCutoff=0.05))
write.csv(keggPA_down,file=paste0(outputPrefix,"_down_pathways.csv")) ##this output file give kegg pathways IDs, names, and statistical values 


# pathview <-pathview(gene.data=res_gene_p[,3],pathway.id="dre04110",species="dre")



## the next few plots will be saved in your working directory

##MA plot
jpeg(file=paste0(outputPrefix,"_MA_basemeans.jpeg"))
plotMA(dds) ##MA plot of base means
dev.off()

jpeg(file=paste0(outputPrefix,"_MA_logfoldchange.jpeg"))
plotMA(res) ##MA plot of logfold change
dev.off()


##PCA
jpeg(file=paste0(outputPrefix,"_PCA.jpeg"))
plotPCA(rlog,"condition")
dev.off()

#print(res_gene)
#print(res_gene$LLgeneAbbrev)
#head(res)
#head(res_gene)
#rownames(res_gene) <- res_gene$LLgeneAbbrev
# there are non-unique values here, need to solve this
#head(res_gene)
##volcano
pdf(file=paste0(outputPrefix,"_volcano.pdf"))
EnhancedVolcano(res, x = 'log2FoldChange', lab = rownames(res), labSize = 6.0,
                y = 'padj',ylab=bquote(~-Lot[10] ~ italic(Padj)),col=c("grey","grey","magenta","magenta"),pCutoff=0.05,
                legendPosition = 'none',cutoffLineType="blank",xlim=c(-3,3),ylim=0,5)
dev.off()

##each gene (with llgene ID, gene name, and entrez ID) and then normalized counts for each sample
countsdata <-counts(dds,normalized=TRUE)
namelist <-read.csv("LLgeneID_entrezID.csv")
##want to merge countsdata with namelist and gene_name
rownames(namelist) <- namelist[,1]
rownames(gene_name) <- gene_name[,1]
genecounts<-merge(gene_name,namelist,by=0)
rownames(genecounts)<-genecounts[,1]
genecounts<-merge(genecounts,countsdata,by=0)
genecounts<-genecounts[,-(1:3)]

write.csv(genecounts,file=paste0(outputPrefix,"_normalized_reads_gene_list.csv")) ##this output file is normalized change in expression of all genes


##dot plot of only significant genes (in res_gene_p)

#make counts data with 1st column is LLgeneID, then merge with res_gene_p
countsdata<-as.data.frame(countsdata)
countsdata1 <-setDT(countsdata, keep.rownames = "LLgeneID")
countsdata1<-countsdata1[match(res_gene_p$LLgeneID,countsdata1$LLgeneID,)]
countsdata1_good<-as.data.frame(t(countsdata1))
colnames(countsdata1_good) <-countsdata1_good[1,]
countsdata1_good<-countsdata1_good[-c(1),]
 

#for(i in 1:ncol(countsdata1_good)) {
#  jpeg(file=paste(outputPrefix,"_expression_",i,".jpeg",sep=""))  
#  print(ggplot(countsdata1_good,aes(x=rownames(countsdata1_good),y=countsdata1_good[,i]))+
#          geom_point()+labs(x="sample",y="normalized counts",title=colnames(countsdata1_good[i])))
#  dev.off()
#}



##dot plot of every gene---this is really a lot so I would recommend only using it if you really want it or altering this code to get a specific gene graph
#genecountst<-as.data.frame(t(genecounts))
#genecounts_good <-genecountst[-c(1,3),]
#colnames(genecounts_good) <-genecounts_good[1,]
#genecounts_good<-genecounts_good[-c(1),]

#for(i in 1:ncol(genecounts_good)) { 
 # jpeg(file=paste0(outputPrefix,"_expression_[i].jpeg"))
  #print(ggplot(genecounts_good,aes(x=rownames(genecounts_good),y=genecounts_good[,i]))+geom_point())
 # dev.off()
#}


##pathway graph

keggPA_down$NegLogPAdj <-log10(keggPA_down$p.adjust)

keggPA_up$NegLogPAdj <--log10(keggPA_up$p.adjust)

keggPA_forplot <-rbind(keggPA_down,keggPA_up) ##only use this for making a plot because I flipped the log adjusted p values in order to force the down regulated pathways to be on the left....



color<-ifelse(keggPA_forplot$NegLogPAdj<0,"blue","yellow")
jpeg(file=paste0(outputPrefix,"_pathways.jpeg"))
ggplot(keggPA_forplot, aes(y=reorder(Description,NegLogPAdj), x=NegLogPAdj))+geom_bar(stat="identity",fill=color)+
  labs(x="log 10 padj",y="pathway")
dev.off()



jpeg(file=paste0(outputPrefix,"_pathways_up.jpeg"))
ggplot(keggPA_up, aes(y=Description, x=NegLogPAdj))+geom_bar(stat="identity")
dev.off()

jpeg(file=paste0(outputPrefix,"_pathways_down.jpeg"))
ggplot(keggPA_down, aes(y=Description, x=-NegLogPAdj))+geom_bar(stat="identity")
dev.off()



###grouping of samples in tree and in heatmap
sampleDists <-dist(t(assay(vst)))
plot(hclust(sampleDists))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$condition, vst$type, sep="-")
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
