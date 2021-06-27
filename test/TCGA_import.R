##############################################################
## 使用TCGAbiolinks向 https://github.com/BioinfoFungi/bioinfo_server 导入TCGA数据
## reference:http://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html
##
## author: wangyang
## data: 2021.6.27
##############################################################


library(TCGAbiolinks)
library(BioinfoR)
initParam(baseUrl = "http://8.140.164.151:8080")
#initParam(baseUrl = "http://localhost:8080")
showParam()

project_df <- (function(){
  projects <- TCGAbiolinks:::getGDCprojects()$project_id
  projects_full <- projects[grepl('^TCGA',projects,perl=T)]
  projects_ID <- sapply(projects_full,function(x){stringi::stri_sub(x,6,10)})
  res <- data.frame(cancer=projects_ID,TCGA_ID=projects_full)
  return(res)
})()

## 添加癌症英文名称
res <- apply(project_df, 1, function(x){
  message(x[1])
  addCancer(name = x[1],enName = x[1])
})
## 添加研究类型
addStudy(name = "mRNA的Count数据",enName = "mRNA_counts")
addStudy(name = "mRNA的FPKM数据",enName = "mRNA_fpkm")

addStudy(name = "lncRNA的Count数据",enName = "lncRNA_counts")
addStudy(name = "lncRNA的FPKM数据",enName = "lncRNA_fpkm")

addStudy(name = "miRNA数据",enName = "miRNA_count")
addStudy(name = "临床数据",enName = "clinical")


## 添加数据来源
addDataOrigin(name = "The Cancer Genome Atlas",enName = "TCGA")


match.file.cases.all <- NULL
for(project in  project_df[,2]){
  #project <- project_df[5,2]
  query<- GDCquery(project = project,
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts")
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- project
  match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
}
readr::write_tsv(match.file.cases.all, path =  "filename_case.txt")


query<- GDCquery(project = "TCGA-ESCA",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "HTSeq - Counts")
GDCdownload(query,directory="TCGA")

#################################################################
data <- GDCprepare(query,directory = "TCGA",summarizedExperiment = F,
                   save = F,)
library(SummarizedExperiment)
expr <- assay(data)
dim(expr)
class(expr)
readr::write_tsv(as.data.frame(expr),file = "TCGA_ESCA_mRNA_count.tsv")

df <- readr::read_tsv("TCGA_ESCA_mRNA_count.tsv")

gene_v22 <- read.csv("~/Downloads/gencode.gene.info.v22 (1).tsv",sep = "\t")
dim(gene_v22)
gene_id <- stringr::str_sub(gene_v22$gene_id,1,15)
length(intersect(rownames(expr),gene_id))


uploadCancerStudy(cancer = "ESCA",study = "mRNA_counts",dataOrigin = "mRNA_counts",path = "TCGA_ESCA_mRNA_count.csv")
df <- readFile("ESCA","mRNA_counts","TCGA",isLocalPath = F)
showParam()
colData(data)









