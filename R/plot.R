#' @importFrom survminer ggsurvplot
#'
#' @export
tcgaGGsurvplot <- function(fit,title=NULL,legend.labs=NULL,xlab="Time(Year)"){
  if(!is.null(fit$legend.labs)){
    legend.labs<-fit$legend.labs
  }
  if(!is.null(fit$title)){
    title <-fit$title
  }
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
