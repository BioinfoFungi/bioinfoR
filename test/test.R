
library(BioinfoR)
#authorize <- read.table("~/.bioinfo/authorize")
#initParam(authorization = authorize$V1,isLocalPath = T)
initParam(baseUrl = "http://8.140.164.151:8080",authorization = NULL)
showParam()
headers <- c(
  "Authorization_SDK"= ""
)
res <- GET("http://localhost:8080/api/lncRNA", add_headers(headers))
content(res)
resp <- POST("http://localhost:8080/api/miRNA",
             add_headers(headers),
             encode="json",
             body = list(name="miRNAGAPDH"))
res <- apply(project_df, 1, function(x){
  message(x[1])
  resp <- POST("http://localhost:8080/api/lncRNA",
               add_headers(headers),
               encode="json",
               body = list(name=x[1]))
})


  content(res)

global <- globalConfig()

global
global$attachment
global$cancerStudy

#query <- list(cancer = "BRAC", study="transcript",dataOrigin="TCGA")
#http_get("/cancer_study/findOne",query = query)
df <- readFile("ESCA","transcript","TCGA")
head(df)

#body <- list(projectId=55,
#             name="testets",
#             path="http://baidu.com/",
#             mediaType='svg')
addAttachment(projectId = 4,path="http://baidu.com/test.csv")
addAttachment(projectId = 56,path="http://baidu.com/test.csv",fileName = "testFileName",fileType = "png")

uploadAttachment(projectId = 4,path = "~/Downloads/clinical_arrange.csv")
uploadAttachment(projectId = 4,path = "README.md",fileName = "testFileName",fileType = "png")

updateProject(projectId = 56,jupyterUrl = "123456")
updateProject(projectId = 56,jupyterUrl = "123456",projectStatus = 1)

#body <- list(cancer = "BRAC",
 #            study= "transcript",
  #           dataOrigin="TCGA",
   #          filename="123456.txt")
addCancerStudy(cancer = "COAD",study = "transcript",dataOrigin = "TCGA",path = "1234567.txt")
addCancerStudy(cancer = "COAD",study = "transcript",dataOrigin = "TCGA",path = "123456.txt",fileType = "svg",fileName = "testFilename")
addCancerStudy(cancer = "COAD",study = "transcript",dataOrigin = "TCGA",path = "123456.csv",width = 300,height = 300)


df <- read.csv("~/Downloads/clinical_arrange.csv")
dim(df)

uploadCancerStudy(cancer = "COAD",study = "Methylation",dataOrigin = "GEO",path = "~/Downloads/clinical_arrange.csv")
uploadCancerStudy(cancer = "BRAC",study = "transcript",dataOrigin = "TCGA",path = "~/Downloads/clinical_arrange.csv",fileType = "csv",fileName = "testFile")
uploadCancerStudy(cancer = "BRAC",study = "transcript",dataOrigin = "TCGA",path = "~/Downloads/clinical_arrange.csv",width = 300,height = 300)
df <- readFile("COAD","transcript","TCGA")
head(df)

addCancer(name = "食管癌",enName = "ESCA")
listAllCancer()

addStudy(name = "基因表达count数据",enName = "HTSeq-Counts")
addDataOrigin(name = "The Cancer Genome Atlas",enName = "TCGA")


listAllStudy()
listAllDataOrigin()
#
initParam(baseUrl = "http://8.140.164.151:8080",remote = "/home/wangyang/workspace/www/data/TCGADOWNLOAD",isLocalPath = F)
showParam()
getFilePath("init_miRNA")

gff_v22 <- readFileByName("gff_v22",isLocalPath = F)

readr::write_tsv(gff_v22,file = "/home/wangyang/workspace/www/data/TCGADOWNLOAD/data/gff_v22.tsv")


readFile("CHOL","FPKM","TCGA",isLocalPath = F)

fit <- tcgaSurvival("CHOL","TP53")
tcgaGGsurvplot(fit)

library(survminer)
ggsurvplot(fit$fit, pval=T, risk.table=F,
           risk.table.height = 0.3,
           data = fit$expr,
           xlab="Time(Year)",
           pval.size=7,
           font.x = c(18),
           font.y = c(18),
           font.legend=c(18),
           legend.title = "",
           legend.labs = c("high expression", "low expression"),
           title="TP53")








