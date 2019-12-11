cosmic=read.table("./data-raw/cosmic_signatures.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
cosmic_mat=cosmic[,4:ncol(cosmic)]
rownames(cosmic_mat)=paste(cosmic$Substitution.Type, cosmic$Trinucleotide, sep="_")
colnames(cosmic_mat)=paste("Signature", seq_len(30), sep="")
cosmic_result <- new("Result", signatures = as.matrix(cosmic_mat), samples = matrix(), type = "Cosmic")
usethis::use_data(cosmic_result, internal = TRUE)
