matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}


X<-structure(list(gene = c("AT1G01040", "AT1G01270", "AT1G01471", "AT1G01680"), log2.fold_change._Mer7_2.1_Mer7_2.2 = c(0, 0, 0, 0), log2.fold_change._Mer7_1.2_W29_S226A_1 = c(0, 0, -1.14, 0 ), log2.fold_change._Mer7_1.2_W29_1 = c(0, 0, 0, 0)), .Names = c("gene", "log2.fold_change._Mer7_2.1_Mer7_2.2", "log2.fold_change._Mer7_1.2_W29_S226A_1", "log2.fold_change._Mer7_1.2_W29_1"), row.names = c(NA, 4L), class = "data.frame")

M <- matrix.please(X)
str(M)


a <- c(62.3, 55.3, 65.3, 59.3, 67.3)
b <- c(2.2, 5.4, 1.3, 2.8, 5.4)
c <- c(0.1, 1.5, 1.6, 2.1, 0.3)
data <- cbind(a, b, c)
data

dimnames(data)
rownames(data) <- c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5")


a <- matrix(rexp(200, rate=.1), ncol=20)
rownames(a) <- c(rep("a", 5), rep("g", 5))
a2 <- limma::avereps(a, ID = rownames(a))






a <- as.data.frame(matrix(rexp(200, rate=.1), ncol=20))
a <- data.frame(
  geneID =c(rep("a", 5), rep("g", 5)),
  a)
limma::avereps(a[,2:ncol(a)], a$geneID)




limma::avereps(eset_df2[,2:ncol(eset_df2)], eset_df2$GB_ID)







