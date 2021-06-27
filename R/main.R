global_env <- new.env()
global_env$headers <- c(
  "Authorization_SDK"= ""
  #"Content-Type"="application/json"
)
global_env$baseUrl <-  "http://localhost:8080/api"
global_env$remote <-  "http://localhost:8080"
global_env$isLocalPath <- T
{
  if(file.exists("~/.bioinfo/authorize")){
    authorize <- read.table("~/.bioinfo/authorize")
    global_env$headers["Authorization_SDK"] <-  authorize$V1
  }
}


#' usage custom url and authorization token
#'
#' @export
initParam <- function(baseUrl=NULL,authorization=NULL,remote=NULL,isLocalPath=NULL){
  if(!is.null(baseUrl)){
    global_env$baseUrl <- paste0(baseUrl,"/api")
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

#' globalConfig
#'
#' API: /global
#' @return attachment and cancer study path
#' @export
globalConfig <- function(){
  res <- http_get(url="/global",showStatus=T)
  return(res)
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
#' getFile("BRAC","transcript","TCGA")
#'
#' @export
getFile <- function(cancer,study,dataOrigin){
  query <- list( cancer = cancer, study=study,dataOrigin=dataOrigin)
  res <- http_get("/cancer_study/findOne",query = query)
  return(res)
}


#' @export
readFile <-function(cancer,study,dataOrigin,isLocalPath=global_env$isLocalPath){
  res <- getFile(cancer,study,dataOrigin)
  if(isLocalPath){
    path <- res$localPath
    if(!file.exists(path)){
      return(message("文件不存在！"))

    }
  }else{
    path <- paste0(global_env$remote,"/", res$networkPath)
  }

  message("File loading path: ",path)
  if(res$fileType=="csv"){
    df <- read.csv(path)
  }else{
    return(message("文件类型不支持！"))
  }
  return(df)
}


listAllCancer <- function(){
  res <- http_get("/cancer/listAll")
  return(res)
}

listAllStudy <- function(){
  res <- http_get("/study/listAll")
  return(res)
}

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
addAttachment <- function(projectId,path,fileName=NULL,fileType=NULL){
  body <- list(projectId=projectId,
               path=path,
               fileName=fileName,
               fileType=fileType)
  res <- http_post("/attachment",body = body)
  return(res)
}

#' @importFrom httr upload_file
#' @export
uploadAttachment <- function(projectId,path,fileName=NULL,fileType=NULL){
  body <- list(projectId=projectId,
               file=upload_file(path),
               fileName=fileName,
               fileType=fileType)
  res <- http_post("/attachment/upload",body = body,encode = "multipart")
  return(res)
}

# (r <- POST(paste0(baseUrl,"/project/updateSDK/56"), add_headers(headers),encode = "json",body = list(jupyterUrl = "testets", projectStatus= 2)))

#' @export
updateProject <- function(projectId,jupyterUrl,projectStatus=0){
  body <- list(jupyterUrl=jupyterUrl,
               projectStatus=projectStatus)
  res <- http_post(paste0("/project/updateSDK/",projectId),body = body)
  return(res)
}


#(r <- POST(paste0(baseUrl,"/cancer_study"), add_headers(headers),encode = "json",body = body))

#' @export
addCancerStudy <- function(cancer,study,dataOrigin,path,fileType=NULL,fileName=NULL,width=NULL,height=NULL){
  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               path=path,
               fileType=fileType,
               fileName=fileName,
               width=width,
               height=height)
  res <- http_post("/cancer_study",body = body)
  return(res)
}



#(r <- POST(paste0(baseUrl,"/cancer_study/upload"), add_headers(headers),encode = "multipart",body = body))


#' @importFrom httr upload_file
#' @export
uploadCancerStudy<- function(cancer,study,dataOrigin,path,fileType=NULL,fileName=NULL,width=NULL,height=NULL){
  body <- list(cancer = cancer,
               study= study,
               dataOrigin=dataOrigin,
               file=upload_file(path),
               fileType=fileType,
               fileName=fileName,
               width=width,
               height=height)

  res <- http_post("/cancer_study/upload",body = body,encode = "multipart")
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



