initParam(host="http://8.140.164.151:8080",authorization = "wangyang1749748955")

genes <- c("TP53","ARID1A")
survival_df <- tcgaSurvival("CHOL",genes,time=30)
