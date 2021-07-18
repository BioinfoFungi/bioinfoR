
#initParam(baseUrl = "http://8.140.164.151:8080",remote = "/home/wangyang/workspace/www/data/TCGADOWNLOAD",isLocalPath = F)
initParam(host  = "http://8.140.164.151:8080",isLocalPath = T,pathPrefix  = "myData")
#initParam(baseUrl = "http://localhost:8080",isLocalPath = F,remote = "/home/wy/Documents/bioinfoR")
showParam()

#############################################################
## 生存分析
#############################################################
readCancerFile("CHOL","FPKM","TCGA",isLocalPath = T)

fit <- tcgaSurvival("CHOL","TP53")
tcgaGGsurvplot(fit)

expr <- tcgaMiRNA("UCEC")
expr@expr
expr@metadata

expr <- tcgaExpr("UCEC")
expr@expr$gene_type
expr@metadata


expr <- tcgaMRNA("UCEC")
expr@expr$gene_type
expr@metadata



### test download
gff_v22 <- readFileByName("init_mRNA")
#http://8.140.164.151/data/gff_v22.tsv
download.file("http://8.140.164.151/data/gff_v22.tsv",destfile = "gff_v22.tsv")
basename("http://8.140.164.151/data/gff_v22.tsv")

checkFile <- function(path,isCheck=T){

  #download.file(path,destfile =localPath)
  message("download file:",path," to ",basename(path))
  res <- suppressWarnings(tryCatch(download.file(url=path,
                                                 destfile="a.txt"),
                                   error=function(e) e[[1]]))
  return(res)

}
checkFile("http://localhost:8080/data/aaa")

headers <- c("Authorization_SDK"= "wangyang1749748955")
download.file("http://localhost:8080/api/cancer_study/download/testFile",
              destfile="a.txt",headers=headers)

file.exists("http://8.140.164.151s:8080/data/gff_v22.tsv")
suppressWarnings(tryCatch(download.file("http://localhost:8080/data/gff_v522.tsv",
                                               destfile="a.txt"),
                                 error=function(e) e[[1]]))
suppressWarnings(
  tryCatch(
    download.file("http://8.140.164.151:8080/data/gff_v522.tsv",
                  destfile="a.txt",
                  method="auto"),
    error=function(e) e[[1]]
  )
)


