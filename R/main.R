global_env <- new.env()
global_env$authorize <- ""
global_env$headers <- c(
  "Authorization_SDK"= ""
  #"Content-Type"="application/json"
)
global_env$host <-"http://localhost:8080"
global_env$baseUrl <-  paste0(global_env$host,"/api")
global_env$isLocalPath <- T
global_env$isAbsolutePath <- F
global_env$pathPrefix <- NULL
{
  if(file.exists("~/.bioinfo/authorize")){
    authorize <- read.table("~/.bioinfo/authorize")
    global_env$headers["Authorization_SDK"] <-  authorize$V1
    global_env$authorize<-authorize$V1
  }
}




#' @export
showParam <- function(){
  return(list(headers=global_env$headers,
              host=global_env$host,
              baseUrl=global_env$baseUrl,
              pathPrefix=global_env$pathPrefix,
              isLocalPath=global_env$isLocalPath))
}


#' usage custom url and authorization token
#'
#' @export
initParam <- function(host=NULL,authorization=NULL,pathPrefix=NULL,isLocalPath=NULL){
  if(!is.null(host)){
    global_env$host <- host
    global_env$baseUrl <-  paste0(host,"/api")
  }
  if(!is.null(authorization)){
    global_env$headers["Authorization_SDK"] <- authorization
    global_env$authorize<-authorization
  }
  if(!is.null(isLocalPath)){
    global_env$isLocalPath=isLocalPath
  }
  if(!is.null(pathPrefix)){
    global_env$pathPrefix <-pathPrefix
  }
}
#' @export
readGzip <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con)
  return(read.csv(textConnection(txt),sep = "\t"))
}

###########################GET##################################

#' http get method
#'
#' @param url url
  #' @importFrom httr GET add_headers content
#' @return list
#'
#' @export
http_get <- function(url,query=NULL,showStatus=F){
  resp <- GET(paste0(global_env$baseUrl,url), add_headers(global_env$headers),query=query)
  resp  <- content(resp)

  if(resp$status!=200){
    message(resp$message)
  }

  data <- resp[["data"]]
  return(data)
}

#' globalConfig
#'
#' API: /global
#' @return attachment and cancer study path
#' @export
globalConfig <- function(){
  res <- http_get(url="/global",showStatus=T)
  return(res)
}

getDownloadPath <- function(id,type,location=NULL){
  path <- paste0(global_env$baseUrl,"/",type,"/downloadById/",id,
                 "?authorize=",global_env$authorize)
  if(!is.null(location)){
    path<- paste0(path,"&location=",location)
  }
  return(path)
}


#' @export
getFileById <-function(Id){
  res <- http_get(paste0("/base_file/findById/",Id))
  return(res)
}





#' @export
downloadById <- function(id,type,location=NULL,toPath){
  options(timeout = max(600, getOption("timeout")))
  if(!is.null(toPath) && !dir.exists(dirname(toPath))){
    dir.create(dirname(toPath),recursive = T)
  }
  path <-  getDownloadPath(id,type,location)
  message(toPath," not found, start downloading from :",path)
  download.file(path,destfile=toPath)
  message("Save the file to: ",toPath)
}

#' @export
downloadFileById <- function(id,toPath=NULL){
  res <- getFileById(id)
  basename <- basename(res$absolutePath)
  downloadById(id,paste0(c(toPath,basename),collapse = "/"))
}







readFileByData <-function(data,type,location=NULL,isLocalPath=global_env$isLocalPath){
  path<- NULL
  if(isLocalPath){
    # 如果服务器在本地
    if(grepl("localhost",global_env$host) | global_env$isAbsolutePath){
      path <- data$absolutePath
      if(!file.exists(path)){
        message("服务器上不存在该文件!")
      }
    }else{
      path <- paste0(c(global_env$pathPrefix,data$relativePath),collapse ="/")
      message("local: ",tools::md5sum(path))
      message("serve: ",data$md5)
      if(!file.exists(path)| tools::md5sum(path)!=data$md5){

        downloadById(id =  data$id,type = type,location=location,toPath = path)
      }
    }
  }else{
    path <-  getDownloadPath(data$id,type,location)
    if(data$fileType=="tsv.gz"){
      df <- readGzip(path)
      message("Load network file from: ",path)
      return(df)
    }
  }
  message("Load file from: ",path)
  if(data$fileType=="csv" | data$fileType=="csv.gz"){
    df <- readr::read_csv(path)
  }else if(data$fileType=="tsv" | data$fileType=="tsv.gz"){
    df <- readr::read_tsv(path)
  }else{
    return(message("文件类型不支持！"))
  }
  return(df)
}



#' @export
readOrganizeFile <- function(enName,location=NULL,isLocalPath=global_env$isLocalPath){
  res <- http_get(paste0("/organize_file/findName/",enName))
  readFileByData(data = res,type = "organize_file",location = location,isLocalPath = isLocalPath)
}



#' @export
readFileById <-function(id,location=NULL,isLocalPath=global_env$isLocalPath){
  res <- getFileById(id)
  readFileByData(data = res,type = "",location =location ,isLocalPath = isLocalPath)
}



#' find one cancer study
#'
#' API: /cancer_study/findOne
#' @param cancer A character(1) canser type
#' @param study A character(1) study type
#' @param dataOrigin A character(1) database source
#'
#' @examples
#'
#' getCancerStudyFile("BRAC","transcript","TCGA")
#'
#' @export
getCancerStudyFile <- function(cancer=NULL,study=NULL,dataOrigin=NULL,uuid=NULL,keyword=NULL,fileName=NULL,analysisSoftware=NULL,experimentalStrategy=NULL){
  query <- list( cancer = cancer,
                 study= study,
                 dataOrigin=dataOrigin,
                 fileName=fileName,
                 keyword=keyword,
                 uuid=uuid,
                 analysisSoftware=analysisSoftware,
                 experimentalStrategy=experimentalStrategy)
  res <- http_get("/cancer_study/findByCategory",query = query)
  return(res)
}

#' @export
getCancerStudyByUUID <- function(uuid){
  res <- http_get(paste0("/cancer_study/findOne/",uuid))
  return(res)
}
#' @export
getCancerStudyById <- function(id){
  res <- http_get(paste0("/cancer_study/findById/",id))
  return(res)
}


#' @export
readCancerFile <-function(... , location=NULL,isLocalPath=global_env$isLocalPath){
  res <- getCancerStudyFile(...)
  if(length(res$content)==1){
    readFileByData(data = res$content[[1]],type = "cancer_study",location = location,isLocalPath = isLocalPath)
  }else if(length(res$content)==0){
    message("没有找到Cancer数据!")
  }else{
    message("文件不唯一!")
    msg_df <- sapply(res$content, function(x){
      return(c(fileName=x$fileName,uuid=x$uuid,relativePath=x$relativePath))
    })
    as.data.frame(t(msg_df))
  }
}




#' @export
listAllCancer <- function(){
  res <- http_get("/cancer/listAll")
  return(res)
}

#' @export
listAllStudy <- function(){
  res <- http_get("/study/listAll")
  return(res)
}

#' @export
listAllDataOrigin <- function(){
  res <- http_get("/data_origin/listAll")
  return(res)
}

###########################POST##################################
#' http post method
#'
#' @param url url
#' @importFrom httr POST add_headers content
#' @return list
#'
#' @export
http_post <- function(url,body,encode = "json",showStatus=F){
  resp <- POST(paste0(global_env$baseUrl,url),
               add_headers(global_env$headers),
               encode=encode,
               body = body)
  resp  <- content(resp)

  if(resp$status!=200){
    message(resp$message)
  }
  data <- resp[["data"]]
  return(data)
}

#(r <- POST(paste0(baseUrl,"/attachment"), add_headers(headers),encode = "json",body = body))

#' @export
addAttachment <- function(projectId,absolutePath,enName=NULL,relativePath=NULL,fileName=NULL,fileType=NULL){
  body <- list(projectId=projectId,
               absolutePath=absolutePath,
               relativePath=relativePath,
               enName=enName,
               fileName=fileName,
               fileType=fileType)
  res <- http_post("/attachment",body = body)
  return(res)
}
#' @export
addCancerStudy <- function(cancer,study,dataOrigin,experimentalStrategy=NULL,analysisSoftware=NULL,absolutePath,relativePath=NULL){

  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               analysisSoftware=analysisSoftware,
               experimentalStrategy=experimentalStrategy,
               absolutePath=absolutePath,
               relativePath=relativePath)
  res <- http_post("/cancer_study",body = body)
  return(res)
}

#' @export
addOrganizeFile <- function(enName,absolutePath,relativePath=NULL,fileType=NULL,fileName=NULL){
  body <- list(enName = enName,
               absolutePath= absolutePath,
               relativePath=relativePath,
               fileType=fileType,
               fileName=fileName)
  res <- http_post("/organize_file",body = body)
  return(res)
}

#' @export
addCancer <- function(name,enName){
  body <- list(name=name,
               enName=enName)
  res <- http_post("/cancer",body = body)
  return(res)
}

#' @export
addStudy <- function(name,enName){
  body <- list(name=name,
               enName=enName)
  res <- http_post("/study",body = body)
  return(res)
}


#' @export
addDataOrigin <- function(name,enName){
  body <- list(name=name,
               enName=enName)
  res <- http_post("/data_origin",body = body)
  return(res)
}

#' @export
addAnalysisSoftware <- function(name,enName){
  body <- list(name=name,
               enName=enName)
  res <- http_post("/AnalysisSoftware",body = body)
  return(res)
}
#' @export
addExperimentalStrategy<- function(name,enName){
  body <- list(name=name,
               enName=enName)
  res <- http_post("/ExperimentalStrategy",body = body)
  return(res)
}
#################################################################


#' @importFrom httr upload_file
#' @export
uploadAttachment <- function(projectId,path,enName=NULL,fileName=NULL,fileType=NULL){
  body <- list(projectId=projectId,
               file=upload_file(path),
               enName=enName,
               fileName=fileName,
               fileType=fileType)
  res <- http_post("/attachment/upload",body = body,encode = "multipart")
  return(res)
}

#' @importFrom httr upload_file
#' @export
uploadCancerStudy<- function(cancer,study,dataOrigin,path,analysisSoftware=NULL,experimentalStrategy=NULL){

  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               file=upload_file(path),
               analysisSoftware=analysisSoftware,
               experimentalStrategy=experimentalStrategy)

  res <- http_post("/cancer_study/upload",body = body,encode = "multipart")
  return(res)
}

#' @importFrom httr upload_file
#' @export
uploadOrganizeFile<- function(path,enName=NULL,fileType=NULL,fileName=NULL){
  body <- list(file=upload_file(path),
               enName=enName,
               fileName=fileName,
               fileType=fileType)
  res <- http_post("/organize_file/upload",body = body,encode = "multipart")
  return(res)
}


#' @export
updateProject <- function(projectId,jupyterUrl,projectStatus=0){
  body <- list(jupyterUrl=jupyterUrl,
               projectStatus=projectStatus)
  res <- http_post(paste0("/project/updateSDK/",projectId),body = body)
  return(res)
}


#(r <- POST(paste0(baseUrl,"/cancer_study"), add_headers(headers),encode = "json",body = body))




#(r <- POST(paste0(baseUrl,"/cancer_study/upload"), add_headers(headers),encode = "multipart",body = body))


















