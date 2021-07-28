#' @importFrom tibble tibble remove_rownames column_to_rownames rownames_to_column
#'
#' @export
tcgaExpr <-function(cancer,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath){
  gff_v22 <- readOrganizeFile("gff_v22",isLocalPath = isLocalPath,location = location)
  expr <- readCancerFile(cancer,dataType,"TCGA",location = location,isLocalPath = isLocalPath)%>%
    tibble::column_to_rownames("X1")%>%
    {.[gff_v22$gene_id,]}%>%
    dplyr::mutate(symbol=gff_v22$gene_name,gene_type=gff_v22$gene_type)%>% #
    {.[!duplicated(.$symbol),]}%>%
    remove_rownames()%>%
    tibble::column_to_rownames("symbol")

  TCGA_id_full <-  colnames(expr)
  TCGA_id_full <- TCGA_id_full[which(TCGA_id_full!="gene_type")]
  tibble::tibble(
    TCGA_id_full=TCGA_id_full,
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
tcgaExprNoType <-function(cancer,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath){
  expr_obj <- tcgaExpr(cancer = cancer,dataType = dataType,location = location ,isLocalPath = isLocalPath)
  expr_obj@expr%>%
    dplyr::select(-gene_type) -> expr
  return(expr)
}

#' @export
tcgaMRNA <- function(cancer,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath){
  expr_obj <- tcgaExpr(cancer = cancer,dataType = dataType,location = location,isLocalPath = isLocalPath)
  expr_obj@expr%>%
    filter(gene_type=="protein_coding")%>%
    dplyr::select(-gene_type) -> expr
  obj <- new("Expr", expr=expr,metadata=expr_obj@metadata)
  return(obj)
}


#' @export
tcgaLncRNA <- function(cancer,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath){
  expr_obj <- tcgaExpr(cancer = cancer,dataType = dataType,location = location,isLocalPath = isLocalPath)
  lncRNA_types <- "3prime_overlapping_ncrna, antisense, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic, sense_overlapping"
  lncRNA_types <- unlist(str_split(lncRNA_types, pattern = ", "))
  expr_obj@expr%>%
    filter(gene_type %in% lncRNA_types)%>%
    dplyr::select(-gene_type) -> expr
  obj <- new("Expr", expr=expr,metadata=expr_obj@metadata)
  return(obj)
}


#' @importFrom dplyr rename_at vars contains select
#' @importFrom tibble tibble
#'
#' @export
tcgaMiRNA <- function(cancer,location=location,isLocalPath=global_env$isLocalPath){
  expr <- readCancerFile(cancer,"miRNA","TCGA",location = location,isLocalPath = isLocalPath)%>%
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


#' @export
tcgaClinical<- function(cancer,location=NULL,isLocalPath=global_env$isLocalPath,time=365){
  tcga_clinical <- readCancerFile(cancer = cancer,study = "clinical",dataOrigin = "TCGA",location = location,isLocalPath = isLocalPath)%>%
    plyr::rename(c(bcr_patient_barcode="Tumor_Sample_Barcode"))%>%
    mutate(days_to_last_followup = ifelse(vital_status=='Alive',days_to_last_follow_up,days_to_death),
           Overall_Survival_Status = ifelse(vital_status=="Alive",0,1),
           days_to_last_followup = days_to_last_followup/time)%>%
    dplyr::select(Tumor_Sample_Barcode,days_to_last_followup,Overall_Survival_Status)
  #tcga_clinical$Overall_Survival_Status <-  as.numeric(tcga_clinical$Overall_Survival_Status)
  return(tcga_clinical)
}

checkGeneExist <- function(existGene,inputGene){
  gene_intersect<- intersect(existGene,inputGene)
  gene_diff <- setdiff(inputGene,gene_intersect)
  if(length(gene_diff)!=0){
    message("remove this gene: ",gene_diff)
  }
  return(gene_intersect)
}
###############################################################################
## Conjoint analysis
###############################################################################

#' Survival analysis of gene expression
#'
#' @examples
#' library(survival)
#' library(survminer)
#' genes <- c("TP53","ELF3","ARID1A")
#' survival_df <- tcgaSurvival("CHOL",genes)
#'
#' df_result <- NULL
#' for(gene in  intersect(colnames(survival_df),genes)){
#'  formula <-as.formula(paste0("Surv(days_to_last_followup,Overall_Survival_Status)~",paste0(gene,collapse = "+")))
#'   diff <- survdiff(formula,data = survival_df)
#'   pVal <- 1 -pchisq(diff$chisq,df =1)
#'   pValue <- signif(pVal,4)
#'   df_result <- rbind(df_result,data.frame(gene=gene,pValue=pValue))
#' }
#'
#' lapply(genes, function(gene) {
#'   formula <-as.formula(paste0("Surv(days_to_last_followup,Overall_Survival_Status)~",paste0(gene,collapse = "+")))
#'   fit <- surv_fit(formula, data = survival_df)
#'   return(list(fit=fit,expr=survival_df,legend=c("Mutation","Wild"),title=gene))
#' })%>%
#'   lapply(function(fit){
#'     dev.new()
#'    res <- ggsurvplot(fit$fit, pval=T, risk.table=F,
#'                      risk.table.height = 0.3,
#'                      data = fit$expr,
#'                      xlab="",
#'                       pval.size=7,
#'                       font.x = c(18),
#'                       font.y = c(18),
#'                       font.legend=c(18),
#'                       legend.title = fit$legend,
#'                       legend.labs = fit$legend,
#'                       title=fit$title)
#'     print(res)
#'   })
#'
#' @importFrom tibble remove_rownames column_to_rownames rownames_to_column
#' @importFrom dplyr mutate '%>%' left_join mutate_at inner_join
#'
#' @export
tcgaSurvival <- function(cancer,genes,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath,time=365){
  #     expr <- readFile(cancer = "CHOL",study = "FPKM",dataOrigin = "TCGA")
  #     gene <- "ARID1A"
  # gff_v22 <- readOrganizeFile("gff_v22",isLocalPath = isLocalPath)

  tcga_expr <- tcgaExprNoType(cancer = cancer,dataType = dataType,location = location,isLocalPath = isLocalPath)
  gene_intersect <- checkGeneExist(rownames(tcga_expr),genes)
  tcga_clinical <- tcgaClinical(cancer,location = location,isLocalPath = isLocalPath,time=time)

  apply(tcga_expr[gene_intersect,], 1,function(row){
    low_hight <- ifelse(row<=median(row),"Low","Hight")
    return(low_hight)})%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    mutate(Tumor_Sample_Barcode = stringr::str_sub(Tumor_Sample_Barcode, 1, 12)) -> tcga_expr_select

  plot_data <- inner_join(tcga_expr_select,tcga_clinical,by="Tumor_Sample_Barcode")

  return(plot_data)
}

#' @export
tcgaMiRNASurvival <- function(cancer,genes,location=NULL,isLocalPath=global_env$isLocalPath,time=365){
  #     expr <- readFile(cancer = "CHOL",study = "FPKM",dataOrigin = "TCGA")
  #     gene <- "ARID1A"
  # gff_v22 <- readOrganizeFile("gff_v22",isLocalPath = isLocalPath)

  tcga_obj <- tcgaMiRNA(cancer = cancer,location = location,isLocalPath = isLocalPath)
  tcga_expr <- tcga_obj@expr
  gene_intersect <- checkGeneExist(rownames(tcga_expr),genes)
  tcga_clinical <- tcgaClinical(cancer,location = location,isLocalPath = isLocalPath,time=time)
  apply(tcga_expr[gene_intersect,], 1,function(row){
    low_hight <- ifelse(row<=median(row),"Low","Hight")
    return(low_hight)})%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    mutate(Tumor_Sample_Barcode = stringr::str_sub(Tumor_Sample_Barcode, 1, 12))-> tcga_expr_select

  plot_data <- inner_join(tcga_expr_select,tcga_clinical,by="Tumor_Sample_Barcode")

  return(plot_data)
}

tcgaSurvival <- function(cancer,genes,dataType = "FPKM",location=NULL,isLocalPath=global_env$isLocalPath,time=365){
  #     expr <- readFile(cancer = "CHOL",study = "FPKM",dataOrigin = "TCGA")
  #     gene <- "ARID1A"
  # gff_v22 <- readOrganizeFile("gff_v22",isLocalPath = isLocalPath)

  tcga_expr <- tcgaExprNoType(cancer = cancer,location = location,isLocalPath = isLocalPath,dataType = dataType)
  gene_intersect <- checkGeneExist(rownames(tcga_expr),genes)
  tcga_clinical <- tcgaClinical(cancer,location = location,isLocalPath = isLocalPath,time=time)
  apply(tcga_expr[gene_intersect,], 1,function(row){
    low_hight <- ifelse(row<=median(row),"Low","Hight")
    return(low_hight)})%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    mutate(Tumor_Sample_Barcode = stringr::str_sub(Tumor_Sample_Barcode, 1, 12))-> tcga_expr_select

  plot_data <- inner_join(tcga_expr_select,tcga_clinical,by="Tumor_Sample_Barcode")

  return(plot_data)
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

#' The correlation between mutation and survival
#'
#'
#' @examples
#' library(survival)
#' library(survminer)
#' genes <- c("TP53","ELF3","ARID1A")
#'
#' df_result <- NULL
#' for(gene in  intersect(colnames(survival_df),genes)){
#'  formula <-as.formula(paste0("Surv(days_to_last_followup,Overall_Survival_Status)~",paste0(gene,collapse = "+")))
#'   diff <- survdiff(formula,data = survival_df)
#'   pVal <- 1 -pchisq(diff$chisq,df =1)
#'   pValue <- signif(pVal,4)
#'   df_result <- rbind(df_result,data.frame(gene=gene,pValue=pValue))
#' }
#'
#' survival_df <- tcgaMutationSurvival("CHOL",genes)
#' lapply(genes, function(gene) {
#'   formula <-as.formula(paste0("Surv(days_to_last_followup,Overall_Survival_Status)~",paste0(gene,collapse = "+")))
#'   fit <- surv_fit(formula, data = survival_df)
#'   return(list(fit=fit,expr=survival_df,legend=c("Mutation","Wild"),title=gene))
#' })%>%
#'   lapply(function(fit){
#'     dev.new()
#'    res <- ggsurvplot(fit$fit, pval=T, risk.table=F,
#'                      risk.table.height = 0.3,
#'                      data = fit$expr,
#'                      xlab="",
#'                       pval.size=7,
#'                       font.x = c(18),
#'                       font.y = c(18),
#'                       font.legend=c(18),
#'                       legend.title = fit$legend,
#'                       legend.labs = fit$legend,
#'                       title=fit$title)
#'     print(res)
#'   })
#'
#'
#' @importFrom tibble remove_rownames column_to_rownames rownames_to_column
#' @importFrom dplyr mutate '%>%' left_join
#'
#' @export
tcgaMutationSurvival <- function(cancer,genes,location=NULL,isLocalPath=global_env$isLocalPath,time=365){
  tcga_clinical <- tcgaClinical(cancer,location = location,isLocalPath = isLocalPath,time=time)
  gene_sample_wider <- tcgaMutationWider(cancer = cancer,location = location,isLocalPath=isLocalPath)
  gene_intersect <- checkGeneExist(rownames(gene_sample_wider),genes)
  gene_sample_wider[gene_intersect,]%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    right_join(tcga_clinical,by="Tumor_Sample_Barcode")%>%
    na.omit()%>%
    dplyr::select(Tumor_Sample_Barcode,gene_intersect,days_to_last_followup,Overall_Survival_Status)-> mutation_survival
  return(mutation_survival)
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
#' boxplot(expr~mutation,data =plot_data,col=c("red","blue"),main="ARID1A")
#'
#' plot_data <- tcgaMutationExpr(cancer="CHOL",gene=c("ARID1A","PTEN"))
#' lapply(plot_data,function(x){
#'   dev.new()
#'   boxplot(expr~mutation,data =x$expr,col=c("red","blue"),main=x$gene)
#' })
#'
#' @export
tcgaMutationExpr <- function(cancer,genes,location=NULL,isLocalPath=global_env$isLocalPath){
  gene_sample_wider <- tcgaMutationWider(cancer = cancer,location = location,isLocalPath=isLocalPath)
  gene_intersect <- checkGeneExist(rownames(gene_sample_wider),genes)
  gene_sample_wider[gene_intersect,]%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")->tcga_mutation_select


  tcga_expr <- tcgaExpr(cancer = cancer,location = location,isLocalPath = isLocalPath)
  gene_intersect <- checkGeneExist(rownames(tcga_expr@expr),genes)
  tcga_expr@expr[gene_intersect,]%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tumor_Sample_Barcode")%>%
    mutate(Tumor_Sample_Barcode = stringr::str_sub(Tumor_Sample_Barcode, 1, 12))%>%
    mutate_at(gene_intersect, as.numeric) -> tcga_expr_select


  if(length(gene_intersect)>1){
    rename_at(tcga_expr_select,gene_intersect, ~ paste0(., "_expr"))-> tcga_expr_select
    rename_at(tcga_mutation_select,gene_intersect, ~ paste0(., "_mutation")) -> tcga_mutation_select
  }else{
    rename_at(tcga_expr_select,gene_intersect, ~ paste0("expr"))-> tcga_expr_select
    rename_at(tcga_mutation_select,gene_intersect, ~ paste0("mutation")) -> tcga_mutation_select
  }

  plot_data <- inner_join(tcga_expr_select,tcga_mutation_select,by="Tumor_Sample_Barcode")
  if(length(gene_intersect)>1){
    plot_data <- lapply(gene_intersect,function(x){
      expr <- paste0(x,"_expr")
      mutation <- paste0(x,"_mutation")
      mutation_expr <- plot_data[,c(expr,mutation)]
      colnames(mutation_expr) <- c("expr","mutation")
      return(list(expr=mutation_expr,gene=x))
    })
  }else{
    plot_data <- subset(plot_data,select=c("expr","mutation"))
  }
  return(plot_data)
}








