#' @importFrom tibble remove_rownames column_to_rownames rownames_to_column
#' @importFrom dplyr mutate '%>%' left_join
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom survival  Surv
#'
#' @export
tcgaSurvival <- function(cancer,gene,dataType = "FPKM",isLocalPath=global_env$isLocalPath){
  #     expr <- readFile(cancer = "CHOL",study = "FPKM",dataOrigin = "TCGA")
  #     gene <- "ARID1A"
  gff_v22 <- readFileByName("gff_v22",isLocalPath = isLocalPath)

  expr <- readFile(cancer,dataType,"TCGA",isLocalPath = isLocalPath)%>%
    tibble::column_to_rownames("X1")%>%
    {.[gff_v22$gene_id,]}%>%
    dplyr::mutate(symbol=gff_v22$gene_name)%>%
    {.[!duplicated(.$symbol),]}%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames("symbol")

  clinical <- readFile(cancer,"clinical","TCGA",isLocalPath = isLocalPath)%>%
    dplyr::mutate(sample_id=submitter_id,
           days_to_last_followup = ifelse(vital_status=='Alive',days_to_last_follow_up,days_to_death),
           Overall_Survival_Status = factor(vital_status,levels = c("Alive","Dead"),labels = c(0,1)),
           days_to_last_followup = days_to_last_followup/365)%>%
    dplyr::select(sample_id,Overall_Survival_Status,days_to_last_followup)

  t(expr[gene,])%>%
    as.data.frame()%>%
    rownames_to_column("sample_id")%>%
    mutate(sample_id = stringr::str_sub(sample_id, 1, 12))%>%
    dplyr::left_join(clinical,by="sample_id")-> expr_gene
  low_high <- ifelse(expr_gene[,gene]<=median(expr_gene[,gene]),"Low","Hight")

  expr_gene$Overall_Survival_Status <-  as.numeric(expr_gene$Overall_Survival_Status)

  #     diff <- survdiff(Surv(days_to_last_followup,Overall_Survival_Status)~low_high,data = expr_gene)
  #     pVal <- 1 -pchisq(diff$chisq,df =1)
  #     pValue <- signif(pVal,4)
  #     cat(pValue)

  fit <- surv_fit(Surv(days_to_last_followup,Overall_Survival_Status) ~ low_high, data = expr_gene)
  return(list(fit=fit,expr=expr_gene))
}

#' @importFrom survminer ggsurvplot
#'
#' @export
tcgaGGsurvplot <- function(fit,title=NULL,legend.labs=c("high expression", "low expression"),xlab="Time(Year)"){
  ggsurvplot(fit$fit, pval=T, risk.table=F,
             risk.table.height = 0.3,
             data = fit$expr,
             xlab=xlab,
             pval.size=7,
             font.x = c(18),
             font.y = c(18),
             font.legend=c(18),
             legend.title = "",
             legend.labs = legend.labs,
             title=title)
}


