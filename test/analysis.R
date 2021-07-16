
initParam(baseUrl = "http://8.140.164.151:8080",remote = "/home/wangyang/workspace/www/data/TCGADOWNLOAD",isLocalPath = F)
expr <- tcgaMiRNA("UCEC")
expr@expr
expr@metadata

expr <- tcgaExpr("UCEC")
expr@expr$gene_type
expr@metadata


expr <- tcgaMRNA("UCEC")
expr@expr$gene_type
expr@metadata
