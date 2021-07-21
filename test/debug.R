library(BioinfoR)
initParam(host="http://8.140.164.151:8080",authorization = "wangyang1749748955")
cancer="CHOL"
fit <- tcgaMutationSurvival("CHOL","ARID1A",location="ALIOSS")
tcgaGGsurvplot(fit)
