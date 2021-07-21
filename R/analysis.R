#' @importFrom tibble tibble remove_rownames column_to_rownames rownames_to_column
#'
#' @export
tcgaExpr <-function(cancer,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath){
  gff_v22 <- readOrganizeFile("gff_v22",isLocalPath = isLocalPath,location = location)
  expr <- readCancerFile(cancer,dataType,"TCGA",location = location)%>%
    tibble::column_to_rownames("X1")%>%
    {.[gff_v22$gene_id,]}%>%
    dplyr::mutate(symbol=gff_v22$gene_name,gene_type=gff_v22$gene_type)%>%
    {.[!duplicated(.$symbol),]}%>%
    remove_rownames()%>%
    tibble::column_to_rownames("symbol")

  tibble::tibble(
    TCGA_id_full=colnames(expr),
    TCGA_id = stringr::str_sub(TCGA_id_full, 1, 16),
    patient_id = stringr::str_sub(TCGA_id, 1, 12),
    tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
    tissue_type = sapply(tissue_type_id, function(x) {
      switch(x,
             "01" = "Primary Solid Tumor",
             "02" = "Recurrent Solid Tumor",
             "03" = "Primary Blood Derived Cancer - Peripheral Blood",
             "05" = "Additional - New Primary",
             "06" = "Metastatic",
             "07" = "Additional Metastatic",
             "11" = "Solid Tissue Normal")}),
    group = ifelse(tissue_type_id == "11", "Normal", "Tumor")
  )->metadata

  obj <- new("Expr", expr=expr,metadata=metadata)
  return(obj)
}

#' @export
tcgaMRNA <- function(cancer,gene,dataType = "FPKM",isLocalPath=global_env$isLocalPath){
  expr <- tcgaExpr(cancer,dataType,isLocalPath)%>%
    filter(gene_type=="protein_coding")%>%
    dplyr::select(-gene_type)
  return(expr)
}



#' @importFrom dplyr rename_at vars contains select
#' @importFrom tibble tibble
#'
#' @export
tcgaMiRNA <- function(cancer,isLocalPath=global_env$isLocalPath){
  expr <- readCancerFile(cancer,"miRNA","TCGA",isLocalPath = isLocalPath)%>%
    select("miRNA_ID",starts_with("read_count"))%>%
    rename_at(vars(contains("read_count")), ~ substr(.,12,length(.)))%>%
    column_to_rownames("miRNA_ID")
  tibble::tibble(
    TCGA_id_full=colnames(expr),
    TCGA_id = stringr::str_sub(TCGA_id_full, 1, 16),
    patient_id = stringr::str_sub(TCGA_id, 1, 12),
    tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
    tissue_type = sapply(tissue_type_id, function(x) {
      switch(x,
             "01" = "Primary Solid Tumor",
             "02" = "Recurrent Solid Tumor",
             "03" = "Primary Blood Derived Cancer - Peripheral Blood",
             "05" = "Additional - New Primary",
             "06" = "Metastatic",
             "07" = "Additional Metastatic",
             "11" = "Solid Tissue Normal")}),
    group = ifelse(tissue_type_id == "11", "Normal", "Tumor")
  )->metadata
  obj <- new("Expr", expr=expr,metadata=metadata)
  return(obj)
}


###############################################################################
## Conjoint analysis
###############################################################################

#' @importFrom tibble remove_rownames column_to_rownames rownames_to_column
#' @importFrom dplyr mutate '%>%' left_join
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom survival  Surv
#'
#' @export
tcgaSurvival <- function(cancer,gene,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath){
  #     expr <- readFile(cancer = "CHOL",study = "FPKM",dataOrigin = "TCGA")
  #     gene <- "ARID1A"
  gff_v22 <- readOrganizeFile("gff_v22",isLocalPath = isLocalPath)

  expr <- readCancerFile(cancer,dataType,"TCGA",isLocalPath = isLocalPath,location = location)%>%
    tibble::column_to_rownames("X1")%>%
    {.[gff_v22$gene_id,]}%>%
    dplyr::mutate(symbol=gff_v22$gene_name)%>%
    {.[!duplicated(.$symbol),]}%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames("symbol")

  clinical <- readCancerFile(cancer,"clinical","TCGA",isLocalPath = isLocalPath,location = location)%>%
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
  return(list(fit=fit,expr=expr_gene,legend.labs=c("high expression", "low expression"),title=gene))
}


#' @importFrom tidyr  pivot_wider
#' @importFrom dplyr right_join
#'
#' @export
tcgaMutationWider <- function(cancer,location=NULL,isLocalPath=global_env$isLocalPath){
  tcga_mutation <- readCancerFile(cancer = cancer ,study = "mutation_varscan2",dataOrigin = "TCGA",isLocalPath = isLocalPath,location = location)%>%
    mutate(Tumor_Sample_Barcode = stringr::str_sub(Tumor_Sample_Barcode, 1, 12))
  gene_sample <- tcga_mutation%>%
    dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode)
  gene_sample <- gene_sample[!duplicated(gene_sample),]
  gene_sample%>%
    mutate(type="Mutation")%>%
    pivot_wider(names_from=Tumor_Sample_Barcode,values_from = c(type))->gene_sample_wider
  gene_sample_wider[is.na(gene_sample_wider)] <- "Wild"
  gene_sample_wider <- column_to_rownames(gene_sample_wider,"Hugo_Symbol")
  return(gene_sample_wider)
}


#' @importFrom tibble remove_rownames column_to_rownames rownames_to_column
#' @importFrom dplyr mutate '%>%' left_join
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom survival  Surv
#'
#' @export
tcgaMutationSurvival <- function(cancer,gene,location=NULL,isLocalPath=global_env$isLocalPath){
  tcga_clinical <- readCancerFile(cancer = cancer,study = "clinical",dataOrigin = "TCGA",location = location,isLocalPath = isLocalPath)%>%
    plyr::rename(c(bcr_patient_barcode="Tumor_Sample_Barcode"))%>%
    mutate(days_to_last_followup = ifelse(vital_status=='Alive',days_to_last_follow_up,days_to_death),
           Overall_Survival_Status = factor(vital_status,levels = c("Alive","Dead"),labels = c(0,1)),
           days_to_last_followup = days_to_last_followup/365)%>%
    dplyr::select(Tumor_Sample_Barcode,days_to_last_followup,Overall_Survival_Status)
  gene_sample_wider <- tcgaMutationWider(cancer = cancer,location = location,isLocalPath=isLocalPath)

  gene_sample_wider[gene,]%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    right_join(tcga_clinical,by="Tumor_Sample_Barcode")%>%
    na.omit()%>%
    dplyr::select(Tumor_Sample_Barcode,mutation=2,days_to_last_followup,Overall_Survival_Status)-> mutation_survival
  mutation_survival$Overall_Survival_Status <-  as.numeric(mutation_survival$Overall_Survival_Status)
  fit <- surv_fit(Surv(days_to_last_followup,Overall_Survival_Status) ~ mutation, data = mutation_survival)
  return(list(fit=fit,expr=mutation_survival,legend.labs=c("Mutation","Wild"),title=gene))
}



#' Correlation between mutation and expression
#'
#' Correlation between mutation and expression
#'
#' @param cancer The type of TCGA cancer, such as CHOL
#' @param gene gene symbol, such as ARID1A
#' @return a data frame
#' @author YangWang <1749748955@qq.com>
#' @examples
#' plot_data <- tcgaMutationExpr(cancer="CHOL",gene="ARID1A")
#'
#' boxplot(expr~mutation,data =plot_data,col=c("red","blue"))
#'
#' ggplot(plot_data,aes(x=mutation,y=expr,color=mutation))+
#'   stat_boxplot()
#'
#' @export
tcgaMutationExpr <- function(cancer,gene,location=NULL,isLocalPath=global_env$isLocalPath){
  gene_sample_wider <- tcgaMutationWider(cancer = cancer,location = location,isLocalPath=isLocalPath)
  gene_sample_wider[gene,]%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    dplyr::rename("mutation"=2) ->gene_sample_select
  tcga_expr <- tcgaExpr(cancer = cancer,location = location,isLocalPath = isLocalPath)
  tcga_expr@expr[gene,]%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    filter(grepl("TCGA",Tumor_Sample_Barcode))%>%
    mutate(Tumor_Sample_Barcode = stringr::str_sub(Tumor_Sample_Barcode, 1, 12))%>%
    dplyr::rename("expr"=2) -> tcga_expr_select
  plot_data <- inner_join(tcga_expr_select,gene_sample_select,by="Tumor_Sample_Barcode")
  plot_data$expr <- as.numeric(plot_data$expr)
  return(plot_data)
}








