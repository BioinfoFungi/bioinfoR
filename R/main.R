global_env <- new.env()
global_env$headers <- c(
  "Authorization_SDK"= ""
  #"Content-Type"="application/json"
)
global_env$host <-"http://localhost:8080"
global_env$baseUrl <-  paste0(global_env$host,"/api")
global_env$remote <-  NULL
global_env$isLocalPath <- T
global_env$isAbsolutePath <- F
global_env$pathPrefix <- NULL
{
  if(file.exists("~/.bioinfo/authorize")){
    authorize <- read.table("~/.bioinfo/authorize")
    global_env$headers["Authorization_SDK"] <-  authorize$V1
  }
}




#' @export
showParam <- function(){
  return(list(headers=global_env$headers,
              baseUrl=global_env$baseUrl,
              remote=global_env$remote,
              isLocalPath=global_env$isLocalPath))
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

#' @export
getFileById <-function(Id){
  res <- http_get(paste0("/base_file/findById/",Id))
  return(res)
}

#' @export
download <- function(fromurl,toPath){
  full_url <- paste0(global_env$baseUrl,fromurl)
  message(toPath," not found, start downloading from :",full_url)
  download.file(full_url,
                destfile=toPath,headers=global_env$headers)
  message("Save the file to: ",toPath)
}


#' @export
downloadById <- function(id,toPath){
  if(!is.null(toPath) && !dir.exists(dirname(toPath))){
    dir.create(dirname(toPath))
  }
  download(paste0("/base_file/downloadById/",id),toPath)
}

#' @export
downloadFileById <- function(id,toPath=NULL){
  res <- getFileById(id)
  basename <- basename(res$absolutePath)
  downloadById(id,paste0(c(toPath,basename),collapse = "/"))
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


#' usage custom url and authorization token
#'
#' @export
initParam <- function(host=NULL,authorization=NULL,remote=NULL,isLocalPath=NULL){
  if(!is.null(host)){
    global_env$host <- host
  }
  if(!is.null(authorization)){
    global_env$headers["Authorization_SDK"] <- authorization
  }
  if(!is.null(remote)){
    global_env$remote <-remote
  }
  if(!is.null(isLocalPath)){
    global_env$isLocalPath=isLocalPath
  }
}
#initParam()


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
getCancerStudyFile <- function(cancer,study,dataOrigin){
  query <- list( cancer = cancer, study=study,dataOrigin=dataOrigin)
  res <- http_get("/cancer_study/findOne",query = query)
  return(res)
}



readFileByData <-function(data,isLocalPath=global_env$isLocalPath){
  path<- NULL
  if(isLocalPath){
    # 如果服务器在本地
    if(grepl("localhost",global_env$host) | global_env$isAbsolutePath){
      path <- data$absolutePath
      if(!file.exists(path)){
        message("服务器上不存在该文件!")
      }
    }else{
      path <- paste0(c(global_env$pathPrefix,data$relativePath),collapse =   "/")

      if(!file.exists(path)){
        downloadById(id =  data$id,toPath = path)
      }
    }
  }else{
    if(data$location=="LOCAL"){
      path <- paste0(c(global_env$host,data$relativePath),collapse = "/")
    }else{
      message("oss文件加载暂不支持")
    }

  }
  message("Load file from: ",path)
  if(data$fileType=="csv"){
    df <- readr::read_csv(path)
  }else if(data$fileType=="tsv"){
    df <- readr::read_tsv(path)
  }else{
    return(message("文件类型不支持！"))
  }
  return(df)
}

#' organize_file cancer_study attachment base_file
#'
#' @export
getFileByEnName <-function(enName,type){
  if(is.null(type)){
    type<- "base_file"
  }
  res <- http_get(paste0("/",type,"/findOne/",enName))
  return(res)
}

#' @export
readFileByEnName <- function(enName,type=NULL,isLocalPath=global_env$isLocalPath){
  res <- getFileByEnName(enName,type)
  readFileByData(res,isLocalPath)
}



#' @export
readFileById <-function(id,isLocalPath=global_env$isLocalPath){
  res <- getFileById(id)
  readFileByData(res,isLocalPath)
}

#' @export
readCancerFile <-function(cancer,study,dataOrigin,isLocalPath=global_env$isLocalPath){
  res <- getCancerStudyFile(cancer,study,dataOrigin)
  readFileByData(res,isLocalPath)
}


#' @export
readFileByName <-function(name,isLocalPath=global_env$isLocalPath){
  res <- http_get(paste0("/organize_file/findByEnName/",name))
  if(isLocalPath){
    path <- res$localPath
    if(!file.exists(path)){
      return(message(path,"不存在！"))
    }
  }else{
    path <- paste0(global_env$remote,"/", res$networkPath)
    if(substr(path, 1, 4)!="http" && !file.exists(path)){
      return(message(path," is not found in your computer!"))
    }
  }

  message("File loading path: ",path)
  if(res$fileType=="csv"){
    df <- read.csv(path,row.names = 1)
  }else if(res$fileType=="tsv"){
    df <- readr::read_tsv(path)
  }else{
    return(message("文件类型不支持！"))
  }
  return(df)
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
addAttachment <- function(projectId,absolutePath,relativePath=NULL,enNamefileName=NULL,fileType=NULL){
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
addCancerStudy <- function(cancer,study,dataOrigin,absolutePath,relativePath=NULL,fileType=NULL,fileName=NULL){
  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               absolutePath=absolutePath,
               relativePath=relativePath,
               fileType=fileType,
               fileName=fileName)
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
uploadCancerStudy<- function(cancer,study,dataOrigin,path,fileType=NULL,fileName=NULL){
  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               file=upload_file(path),
               fileType=fileType,
               fileName=fileName)

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


















