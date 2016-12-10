library(c("pathview", "gage", "gageData", "GenomicAlignments", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
data(kegg.gs)
kegg.mouse = kegg.gsets("mouse")
kegg.gs = kegg.mouse$kg.sets[kegg.mouse$sigmet.idx]
neural_diff_kallisto_deseq_var_stab = varianceStabilizingTransformation(neural_diff_kallisto_deseq2, blind = FALSE)
neural_diff_kallisto_deseq_var_stab = assay(neural_diff_kallisto_deseq_var_stab)
neural_diff_kallisto_deseq_var_stab = neural_diff_kallisto_deseq_var_stab[which(!is.na(entrez.id$entrezgene[match(rownames(neural_diff_kallisto_deseq_var_stab), entrez.id$external_gene_name)])),]
rownames(neural_diff_kallisto_deseq_var_stab) = entrez.id$entrezgene[match(rownames(neural_diff_kallisto_deseq_var_stab), entrez.id$external_gene_name)]
neural_diff_kallisto_deseq2.day2.kegg = gage(neural_diff_kallisto_deseq_var_stab, gsets = kegg.gs, ref = 1:4, samp = 5:8, compare = "paired")
neural_diff_kallisto_deseq2.day5.kegg = gage(neural_diff_kallisto_deseq_var_stab, gsets = kegg.gs, ref = 5:8, samp = 9:12, compare = "paired")

neural_diff_kallisto_deseq2.day5.kegg$greater= neural_diff_kallisto_deseq2.day5.kegg$greater[which(!is.na(neural_diff_kallisto_deseq2.day5.kegg$greater[ ,"p.geomean"])), ]
neural_diff_kallisto_deseq2.day5.kegg$less = neural_diff_kallisto_deseq2.day5.kegg$less[which(!is.na(neural_diff_kallisto_deseq2.day5.kegg$less[ ,"p.geomean"])), ]

neural_diff_kallisto_deseq2.day2.kegg$greater= neural_diff_kallisto_deseq2.day2.kegg$greater[which(!is.na(neural_diff_kallisto_deseq2.day2.kegg$greater[ ,"p.geomean"])), ]
neural_diff_kallisto_deseq2.day2.kegg$less = neural_diff_kallisto_deseq2.day2.kegg$less[which(!is.na(neural_diff_kallisto_deseq2.day2.kegg$less[ ,"p.geomean"])), ]

neural_diff_kallisto_deseq2.day2.kegg.padj = apply(cbind(2*neural_diff_kallisto_deseq2.day2.kegg$greater[order(neural_diff_kallisto_deseq2.day2.kegg$greater[,"stat.mean"], decreasing =TRUE),"p.val"], 2*neural_diff_kallisto_deseq2.day2.kegg$less[order(neural_diff_kallisto_deseq2.day2.kegg$greater[,"stat.mean"], decreasing = FALSE),"p.val"]), 1, min)
neural_diff_kallisto_deseq2.day2.kegg.padj = p.adjust(neural_diff_kallisto_deseq2.day2.kegg.padj)
neural_diff_kallisto_deseq2.day2.kegg.padj = data.frame(p.adj = neural_diff_kallisto_deseq2.day2.kegg.padj, direction = apply(cbind(neural_diff_kallisto_deseq2.day2.kegg$greater[order(neural_diff_kallisto_deseq2.day2.kegg$greater[,"stat.mean"], decreasing =TRUE),"p.val"], neural_diff_kallisto_deseq2.day2.kegg$less[order(neural_diff_kallisto_deseq2.day2.kegg$greater[,"stat.mean"], decreasing = FALSE),"p.val"]), 1, function(x) which.min(x)))
neural_diff_kallisto_deseq2.day2.kegg.padj = neural_diff_kallisto_deseq2.day2.kegg.padj[order(neural_diff_kallisto_deseq2.day2.kegg.padj[,1], decreasing = FALSE), ]
neural_diff_kallisto_deseq2.day2.kegg.padj$direction = as.factor(neural_diff_kallisto_deseq2.day2.kegg.padj$direction)
levels(neural_diff_kallisto_deseq2.day2.kegg.padj$direction) = c("greater", "less")


neural_diff_kallisto_deseq2.day5.kegg.padj = apply(cbind(2*neural_diff_kallisto_deseq2.day5.kegg$greater[order(neural_diff_kallisto_deseq2.day5.kegg$greater[,"stat.mean"], decreasing =TRUE),"p.val"], 2*neural_diff_kallisto_deseq2.day5.kegg$less[order(neural_diff_kallisto_deseq2.day5.kegg$greater[,"stat.mean"], decreasing = FALSE),"p.val"]), 1, min)
neural_diff_kallisto_deseq2.day5.kegg.padj = p.adjust(neural_diff_kallisto_deseq2.day5.kegg.padj)
neural_diff_kallisto_deseq2.day5.kegg.padj = data.frame(p.adj = neural_diff_kallisto_deseq2.day5.kegg.padj, direction = apply(cbind(neural_diff_kallisto_deseq2.day5.kegg$greater[order(neural_diff_kallisto_deseq2.day5.kegg$greater[,"stat.mean"], decreasing =TRUE),"p.val"], neural_diff_kallisto_deseq2.day5.kegg$less[order(neural_diff_kallisto_deseq2.day5.kegg$greater[,"stat.mean"], decreasing = FALSE),"p.val"]), 1, function(x) which.min(x)))
neural_diff_kallisto_deseq2.day5.kegg.padj = neural_diff_kallisto_deseq2.day5.kegg.padj[order(neural_diff_kallisto_deseq2.day5.kegg.padj[,1], decreasing = FALSE), ]
neural_diff_kallisto_deseq2.day5.kegg.padj$direction = as.factor(neural_diff_kallisto_deseq2.day5.kegg.padj$direction)
levels(neural_diff_kallisto_deseq2.day5.kegg.padj$direction) = c("greater", "less")


neural_diff_kallisto_deseq2.day2andday5.keggintersect = cbind(neural_diff_kallisto_deseq2.day2.kegg.padj[intersect(rownames(neural_diff_kallisto_deseq2.day5.kegg.padj)[which(neural_diff_kallisto_deseq2.day5.kegg.padj$p.adj < 0.05)], rownames(neural_diff_kallisto_deseq2.day2.kegg.padj)[which(neural_diff_kallisto_deseq2.day2.kegg.padj$p.adj < 0.05)]), ], neural_diff_kallisto_deseq2.day5.kegg.padj[intersect(rownames(neural_diff_kallisto_deseq2.day5.kegg.padj)[which(neural_diff_kallisto_deseq2.day5.kegg.padj$p.adj < 0.05)], rownames(neural_diff_kallisto_deseq2.day2.kegg.padj)[which(neural_diff_kallisto_deseq2.day2.kegg.padj$p.adj < 0.05)]), ])

cJunKEGGpathways = data.frame(day2padj = c(neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04010 MAPK signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04210 Apoptosis", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04310 Wnt signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04510 Focal adhesion", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04668 TNF signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04722 Neurotrophin signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04912 GnRH signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04921 Oxytocin signaling pathway", "p.adj"]), day2direction = c(neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04010 MAPK signaling pathway", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04210 Apoptosis", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04310 Wnt signaling pathway", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04510 Focal adhesion", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04668 TNF signaling pathway", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04722 Neurotrophin signaling pathway", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04912 GnRH signaling pathway", "direction"], neural_diff_kallisto_deseq2.day2.kegg.padj["mmu04921 Oxytocin signaling pathway", "direction"]), day5padj = c(neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04010 MAPK signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04210 Apoptosis", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04310 Wnt signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04510 Focal adhesion", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04668 TNF signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04722 Neurotrophin signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04912 GnRH signaling pathway", "p.adj"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04921 Oxytocin signaling pathway", "p.adj"]), day5direction = c(neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04010 MAPK signaling pathway", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04210 Apoptosis", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04310 Wnt signaling pathway", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04510 Focal adhesion", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04668 TNF signaling pathway", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04722 Neurotrophin signaling pathway", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04912 GnRH signaling pathway", "direction"], neural_diff_kallisto_deseq2.day5.kegg.padj["mmu04921 Oxytocin signaling pathway", "direction"]))
rownames(cJunKEGGpathways) = c("mmu04010 MAPK signaling pathway", "mmu04210 Apoptosis", "mmu04310 Wnt signaling pathway", "mmu04510 Focal adhesion", "mmu04668 TNF signaling pathway", "mmu04722 Neurotrophin signaling pathway", "mmu04912 GnRH signaling pathway", "mmu04921 Oxytocin signaling pathway")
cJunKEGGpathways$day2direction = as.factor(cJunKEGGpathways$day2direction)
cJunKEGGpathways$day5direction = as.factor(cJunKEGGpathways$day5direction)
levels(cJunKEGGpathways$day2direction) = c("greater", "less")
levels(cJunKEGGpathways$day5direction) = c("greater", "less")


dat = data.frame(day0mean = apply(neural_diff_kallisto_deseq_var_stab[ ,1:4], 1, mean), day2mean = apply(neural_diff_kallisto_deseq_var_stab[ ,5:8], 1, mean), day5mean = apply(neural_diff_kallisto_deseq_var_stab[ ,9:12], 1, mean))
mapk_genes = kegg.mouse$kg.sets['mmu04010 MAPK signaling pathway']
mapk_genes = entrez.id$external_gene_name[match(noquote(mapk_genes[[1]]), entrez.id$entrezgene)]
mapk_genes = mapk_genes[which(!is.na(mapk_genes))]
mapk_means = data.frame(day0mean = apply(neural_diff_kallisto_deseq_var_stab[ mapk_genes,1:4], 1, mean), day2mean = apply(neural_diff_kallisto_deseq_var_stab[ mapk_genes,5:8], 1, mean), day5mean = apply(neural_diff_kallisto_deseq_var_stab[ mapk_genes,9:12], 1, mean))
all_day0vs2.lm = lm(dat$day2mean ~ dat$day0mean)
mapk_day0vs2.lm = lm(mapk_means$day2mean ~ mapk_means$day0mean)

ggplot(dat, aes(x = day0mean, y = day2mean)) + geom_point(alpha = 0.1) + geom_abline(slope = all_day0vs2.lm$coefficients[2], intercept = all_day0vs2.lm$coefficients[1]) + geom_point(data = mapk_means, aes(x = day0mean, y = day2mean), colour = "red") + geom_abline(slope = mapk_day0vs2.lm$coefficients[2], intercept = mapk_day0vs2.lm$coefficients[1], size = 1, linetype = 2, colour = "red")

all_day2vs5.lm = lm(dat$day5mean ~ dat$day2mean)
mapk_day2vs5.lm = lm(mapk_means$day5mean ~ mapk_means$day2mean)
ggplot(dat, aes(x = day2mean, y = day5mean)) + geom_point(alpha = 0.1) + geom_abline(slope = all_day2vs5.lm$coefficients[2], intercept = all_day2vs5.lm$coefficients[1]) + geom_point(data = mapk_means, aes(x = day0mean, y = day2mean), colour = "red") + geom_abline(slope = mapk_day2vs5.lm$coefficients[2], intercept = mapk_day2vs5.lm$coefficients[1], size = 1, linetype = 2, colour = "red")

