source("https://bioconductor.org/biocLite.R")
biocLite(c("MotifDb",  "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10",
           "org.Sc.sgd.db", "TxDb.Mmusculus.UCSC.mm10.knownGene",
           "motifStack", "seqLogo", "Biostrings"))
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10 = BSgenome.Mmusculus.UCSC.mm10
AP1motif = DNAString("TGACTCA")
gene_trans <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")
promoter.seqs <- getPromoterSeq(gene_trans, mm10, upstream=2000, downstream=0)
# multiple entries per gene correspond to different tss
# get the first one
promoter.seqs.unique = DNAStringSet()
for(i in names(promoter.seqs)){
	seq = promoter.seqs[[i]]
	promoter.seqs.unique = append(promoter.seqs.unique, seq[1])
}
AP1matches = sapply(promoter.seqs.unique, function(x) countPattern(AP1motif, x, max.mismatch = 0, min.mismatch = 0, with.indels = FALSE) +  countPattern(AP1motif, reverseComplement(x), max.mismatch = 0, min.mismatch = 0, with.indels = FALSE))

ensembl=useMart("ensembl")
ensembl=useDataset("mmusculus_gene_ensembl",mart=ensembl)
entrez2gene_name = getBM(attributes= c("entrezgene", "external_gene_name"), mart = ensembl)
length(which(names(AP1matches) %in% entrez2gene_name$entrezgene))
AP1matches = AP1matches[which(names(AP1matches) %in% entrez2gene_name$entrezgene)]
names(AP1matches) = entrez2gene_name$external_gene_name[match(names(AP1matches), entrez2gene_name$entrezgene)]
AP1matches = AP1matches[which(names(AP1matches) %in% colnames(datExpr))]
AP1matches.mean_abundance = data.frame(abundance.mean = c(mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[1:4]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[5:8]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[9:12]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[13:16]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[1:4]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[5:8]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[9:12]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[13:16]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[1:4]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[5:8]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[9:12]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[13:16]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[1:4]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[5:8]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[9:12]), mean(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[13:16])), abundance.sd = c(sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[1:4]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[5:8]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[9:12]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 0)], ], 2, mean)[13:16]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[1:4]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[5:8]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[9:12]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 1)], ], 2, mean)[13:16]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[1:4]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[5:8]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[9:12]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 2)], ], 2, mean)[13:16]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[1:4]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[5:8]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[9:12]), sd(apply(neural_diff_txi$abundance[names(AP1matches)[which(AP1matches == 3)], ], 2, mean)[13:16])), days = rep(c(0, 2, 5, 12), times = 4), matches = rep(c(0, 1, 2, 3), each = 4))
AP1matches.mean_abundance.barplot = ggplot(AP1matches.mean_abundance, aes(x = factor(matches), y = abundance.mean, fill = factor(days)), color = factor(days)) + stat_summary(fun.y = mean, position = position_dodge(), geom = "bar")
AP1matches.mean_abundance.barplot + geom_linerange(aes(ymax = abundance.mean + 2.16*abundance.sd, ymin = abundance.mean - 2.16*abundance.sd), position = position_dodge(width = 0.9)) + theme_bw()

x = neural_diff_kallisto_deseq2.day2$log2FoldChange[match(names(AP1matches)[which(AP1matches == 0)], rownames(neural_diff_kallisto_deseq2.day2))]
y = neural_diff_kallisto_deseq2.day2$log2FoldChange[match(names(AP1matches)[which(AP1matches > 0)], rownames(neural_diff_kallisto_deseq2.day2))]
AP1nonmatches_day0vs2_log2foldchange_gaussian_kernel = density(x, kernel = "gaussian")
AP1matches_day0vs2_log2foldchange_gaussian_kernel = density(y, kernel = "gaussian")

x = neural_diff_kallisto_deseq.day2vsday5.deseq.results$log2FoldChange[match(names(AP1matches)[which(AP1matches == 0)], rownames(neural_diff_kallisto_deseq.day2vsday5.deseq.results))]
x = x[!is.na(x)]
y = neural_diff_kallisto_deseq.day2vsday5.deseq.results$log2FoldChange[match(names(AP1matches)[which(AP1matches > 0)], rownames(neural_diff_kallisto_deseq.day2vsday5.deseq.results))]
y = y[!is.na(y)]
AP1nonmatches_day2vs5_log2foldchange_gaussian_kernel = density(x, kernel = "gaussian")
AP1matches_day2vs5_log2foldchange_gaussian_kernel = density(y, kernel = "gaussian")

# CRE motif matches
CREmotif = DNAString("TGATGTCA")
CREmatches = sapply(promoter.seqs.unique, function(x) countPattern(CREmotif, x, max.mismatch = 0, min.mismatch = 0, with.indels = FALSE) +  countPattern(CREmotif, reverseComplement(x), max.mismatch = 0, min.mismatch = 0, with.indels = FALSE))
CREmatches = CREmatches[which(names(CREmatches) %in% entrez2gene_name$entrezgene)]
names(CREmatches) = entrez2gene_name$external_gene_name[match(names(CREmatches), entrez2gene_name$entrezgene)]
CREmatches = CREmatches[which(names(CREmatches) %in% colnames(datExpr))]


x = neural_diff_kallisto_deseq2.day2$log2FoldChange[match(names(CREmatches)[which(CREmatches == 0)], rownames(neural_diff_kallisto_deseq2.day2))]
y = neural_diff_kallisto_deseq2.day2$log2FoldChange[match(names(CREmatches)[which(CREmatches > 0)], rownames(neural_diff_kallisto_deseq2.day2))]
CREnonmatches_day0vs2_log2foldchange_gaussian_kernel = density(x, kernel = "gaussian")
CREmatches_day0vs2_log2foldchange_gaussian_kernel = density(y, kernel = "gaussian")

x = neural_diff_kallisto_deseq.day2vsday5.deseq.results$log2FoldChange[match(names(CREmatches)[which(CREmatches == 0)], rownames(neural_diff_kallisto_deseq.day2vsday5.deseq.results))]
x = x[!is.na(x)]
y = neural_diff_kallisto_deseq.day2vsday5.deseq.results$log2FoldChange[match(names(CREmatches)[which(CREmatches > 0)], rownames(neural_diff_kallisto_deseq.day2vsday5.deseq.results))]
y = y[!is.na(y)]
CREnonmatches_day2vs5_log2foldchange_gaussian_kernel = density(x, kernel = "gaussian")
CREmatches_day2vs5_log2foldchange_gaussian_kernel = density(y, kernel = "gaussian")
