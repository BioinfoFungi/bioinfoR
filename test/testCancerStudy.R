library(BioinfoR)
####################################################################
# 添加
####################################################################
## 添加term
addCancer(name = "食管癌",enName = "ESCA")
addStudy(name = "count",enName = "count")
addDataOrigin(name = "TCGA",enName = "TCGA")
addCancerStudy(cancer = "ESCA",
               study = "count",
               dataOrigin = "TCGA",
               absolutePath =  "/home/wy/Downloads/TCGA_ACC_Counts.tsv")

## 添加organizeFile
addOrganizeFile(enName = "test123",
                absolutePath = "/home/wy/Downloads/TCGA_ACC_Counts.tsv",
                relativePath = "Downloads/TCGA_ACC_Counts.tsv")
  ## 保存 Attachment
addAttachment(projectId = 55,absolutePath = "/home/wy/Downloads/TCGA_ACC_Counts.tsv")



####################################################################
# 读取测试
####################################################################
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
initParam(host = "http://127.0.0.1")
readCancerFile(cancer = "ESCA",
               study = "count",
               dataOrigin = "TCGA")

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
                   path<- "TCGA_ACC_Counts.tsv")

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
uploadOrganizeFile( path="TCGA_ACC_Counts.tsv",
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
uploadCancerStudy( cancer = "ESCA",
                   study = "count",
                   dataOrigin = "TCGA",
                   path="TCGA_ACC_Counts.tsv")






