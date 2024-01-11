#Annotation in R 

library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("Homo sapiens","EnsDb"))
require("ensembldb")
ensembl.104 <- ah[["AH95744"]]
