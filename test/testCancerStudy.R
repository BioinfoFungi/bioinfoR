library(BioinfoR)
####################################################################
# 添加
####################################################################
## 添加term
addCancer(name = "食管癌",enName = "ESCA")
addCancer(name = "结直肠癌",enName = "COAD")
addCancer(name = "乳腺癌",enName = "BRAC")
addStudy(name = "count",enName = "count")
addDataOrigin(name = "TCGA",enName = "TCGA")
addAnalysisSoftware(name = "Deseq2",enName = "Deseq2")
addExperimentalStrategy(name = "WGS",enName = "WGS")
addExperimentalStrategy(name = "RNA-Seq",enName = "RNA-Seq")
addCancerStudy(cancer = "ESCA",
               study = "Counts",
               dataOrigin = "TCGA",
               analysisSoftware="Deseq2",
               experimentalStrategy="WGS",
               relativePath = "data/TCGA_ESCA_clinical.tsv",
               absolutePath =  "/home/wy/Downloads/TCGA_ESCA_clinical.tsv")
addCancerStudy(cancer = "ESCA",
               study = "Counts",
               dataOrigin = "TCGA",
               absolutePath =  "/home/wy/Downloads/TCGA_ESCA_clinical.tsv")
addCancerStudy(cancer = "BRAC",
               study = "count",
               dataOrigin = "TCGA",
               relativePath = "data/TCGA_ACC_Counts.tsv.gz",
               absolutePath =  "/home/wy/Downloads/TCGA_ACC_Counts.tsv.gz")

## 添加organizeFile
addOrganizeFile(enName = "132456",
                absolutePath = "/home/wy/Downloads/TCGA_ESCA_clinical.tsv",
                relativePath = "Downloads/TCGA_ACC_Counts.tsv")

addOrganizeFile(enName = "init_lncRNA",
                absolutePath = "/home/wy/Downloads/init_lncRNA.tsv.gz",
                relativePath = "Downloads/init_lncRNA.tsv.gz")
    ## 保存 Attachment
addAttachment(projectId = 55,absolutePath = "/home/wy/Downloads/TCGA_ACC_Counts.tsv")



####################################################################
# 读取测试
####################################################################
#ALIOSS,LOCAL;
(function(){
  cancer="ESCA"
  study="count"
  dataOrigin="TCGA"
  isLocalPath=F
  global_env <- new.env()
  global_env$remote <- "http://localhost:8080/api/cancer_study/download"
  global_env$host <-"http://localhost:8080"
  global_env$pathPrefix <- NULL
  global_env$isAbsolutePath <- F

  res <- getCancerStudyFile(cancer,study,dataOrigin)
  if(isLocalPath){
    path<- NULL
    # 如果服务器在本地
    if(grepl("localhost_",global_env$host) | global_env$isAbsolutePath){
      path <- res$absolutePath
      if(!file.exists(path)){
        message("服务器上不存在该文件!")
      }
    }else{
      path <- paste0(c(global_env$pathPrefix,res$relativePath),collapse =   "/")
      message("Load file from: ",path)
      if(!file.exists(path)){
        downloadById(id = res$id,toPath = path)

      }
    }
  }else{
    if(res$location=="LOCAL"){
      path <- paste0(c(global_env$host,res$relativePath),collapse = "/")
    }else{
      message("oss文件加载暂不支持")
    }
    message("Load file from: ",path)
  }

  if(res$fileType=="csv"){
    df <- readr::read_csv(path)
  }else if(res$fileType=="tsv"){
    df <- readr::read_tsv(path)
  }else{
    return(message("文件类型不支持！"))
  }
})()
initParam(host = "http://127.0.0.1:8080")
readCancerFile(cancer = "ESCA",
               study = "count",
               dataOrigin = "TCGA",experimentalStrategy="WGS",isLocalPath = T,location = "LOCAL")
dim(df)


gff_v22 <- readOrganizeFile("132456")
res <- getCancerStudyFile(cancer = "ESCA",
                   study = "count",
                   dataOrigin = "TCGA")
sapply(res$content, function(x){
return(c(analysisSoftware=x$fileName,uuid=x$uuid,relativePath=x$relativePath))
})

res <- getCancerStudyFile(fileName = "TCGA_ESCA_DESeq2_202162637891")
sapply(res$content, function(x){
  return(c(fileName=x$fileName,uuid=x$uuid,relativePath=x$relativePath))
})

readCancerFile(uuid = "c752940e-9351-4f96-8880-76232ce9a327",location = "ALIOSS")
res <- getCancerStudyFile(uuid = "c752940e-9351-4f96-8880-76232ce9a327")
sapply(res$content, function(x){
  return(c(fileName=x$fileName,uuid=x$uuid,relativePath=x$relativePath))
})
readCancerFile(keyword = "DESeq2")
res <- getCancerStudyFile(keyword = "DESeq2")
msg_df <- sapply(res$content, function(x){
  return(c(fileName=x$fileName,uuid=x$uuid,relativePath=x$relativePath))
})
as.data.frame(t(msg_df))

res <- getCancerStudyFile(cancer = "ESCA",
                          study = "count",
                          dataOrigin = "TCGA")
sapply(res$content, function(x){
  return(c(fileName=x$fileName,uuid=x$uuid,relativePath=x$relativePath))
})


getCancerStudyByUUID("285144cb-397f-43a2-94e9-eaa51eeec686")
getCancerStudyById(5)
readFileById(id=1)
readFileByEnName(enName="TCGA_ACC_Counts",type = "attachment")

(function(){
  res <- http_get(paste0("/base_file/findById/1"))
  return(res)
})()

####################################################################
# 下载测试
####################################################################
(function(){
  headers <- c("Authorization_SDK"= "")
  download.file("http://localhost:8080/api/base_file/downloadById/1",
                destfile="a.txt",headers=headers)
})
downloadFileById(id = 1,toPath = ".")
downloadFileById(id = 1,toPath = "abc")
downloadFileById(id = 1)
####################################################################
# 上传测试
####################################################################
(function(){
  library(httr)
  projectId=55
  path<- "TCGA_ACC_Counts.tsv"
  fileName=NULL
  fileType=NULL
  body <- list(projectId=projectId,
               file=upload_file(path),
               fileName=fileName,
               enName= "aaaaa",
               fileType=fileType)
  res <- http_post("/attachment/upload",body = body,encode = "multipart")
  return(res)
})
uploadAttachment(  projectId=55,
                   path<- "data/TCGA_ESCA_clinical.tsv")

  (function(){
  library(httr)
  path<- "TCGA_ACC_Counts.tsv"
  enName<- "aaaaa"
  fileName=NULL
  fileType=NULL
  body <- list(file=upload_file(path),
               enName=enName,
               fileName=fileName,
               fileType=fileType)
  res <- http_post("/organize_file/upload",body = body,encode = "multipart")
  return(res)
})
uploadOrganizeFile( path="data/TCGA_ESCA_clinical.tsv",
                    enName= "aaaaa")
uploadOrganizeFile( path="TCGA_ACC_Counts.tsv")

(function(){
  cancer = "ESCA"
  study = "count"
  dataOrigin = "TCGA"
  path="TCGA_ACC_Counts.tsv"
  fileType=NULL
  fileName=NULL
  width=NULL
  height=NULL

  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               file=upload_file(path),
               fileType=fileType,
               fileName=fileName,
               width=width,
               height=height)

  res <- http_post("/cancer_study/upload",body = body,encode = "multipart")
})
uploadCancerStudy( cancer = "COAD",
                   study = "count",
                   dataOrigin = "TCGA",
                   path="data/TCGA_ESCA_clinical.tsv")
####################################################################
# 真正远程测试
####################################################################
initParam(host  = "http://8.140.164.151:8080")
showParam()
#ALIOSS,LOCAL;
readCancerFile(cancer = "ESCA",
               study = "clinical",
               dataOrigin = "TCGA",isLocalPath = F,location="ALIOSS")

library(tools)
#md5sum(dir(R.home(), pattern = "^COPY", full.names = TRUE))
md5sum("data/TCGA_ESCA_clinical.tsv")

url <- "http://localhost:8080/api/base_file/downloadById/1?authorize=wangyang1749748955"
readr::read_tsv("http://127.0.0.1:8080/api/base_file/downloadById/13?authorize=wangyang1749748955")
read.csv("http://127.0.0.1:8080/api/base_file/downloadById/13?authorize=wangyang1749748955",sep = "\t")
x <- read.csv(url, header=FALSE, stringsAsFactors=FALSE, fileEncoding="latin1")
readr::read_tsv("http://localhost:8080/gzip")
con <- gzcon(url("http://localhost:8080/gzip"))
txt <- readLines(con)
read.csv(textConnection(txt),sep = "\t")

data.table::fread("http://localhost:8080/gzip")
initParam(host = "http://127.0.0.1:8080")
readOrganizeFile("aaaaa",isLocalPath = T,location = "ALIOSS")
df <- readGzip("http://localhost:8080/api/base_file/downloadById/13?authorize=wangyang1749748955")
dim(df)





###################################################################
df <- readCancerFile("OV","Counts","TCGA",isLocalPath = F)
uploadCancerStudy( cancer = "OV",
                   study = "miRNA_ISO",
                   dataOrigin = "GEO",
                   enName = "GSE83693",
                   path="/home/wangyang/workspace/OvarianCancer/data/GSE83693_sample.tsv")
