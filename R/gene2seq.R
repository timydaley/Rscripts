# load required libraries
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("org.Hs.eg.db")
biocLite("GenomicRanges")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
library(biomaRt)
library(org.Hs.eg.db)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

get_top_tss <- function(gene_names){
	mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
	tss_map = getBM(attributes = c("hgnc_symbol", "chromosome_name", "strand", "transcript_start"), mart = mart)
	# remove genes without hugo symbols
	hugo_gene_names = gene_names[gene_names %in% tss_map$hgnc_symbol]
	if(length(hugo_gene_names) != length(gene_names)){
		print("The following genes are not in HGNC and TSS cannot be determined:")
		print(setdiff(gene_names, hugo_gene_names))
	}
	return(tss_map[match(hugo_gene_names, tss_map$hgnc_symbol), ])
}

get_all_tss <- function(gene_names){
	mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
	tss_map = getBM(attributes = c("hgnc_symbol", "chromosome_name", "strand", "transcript_start"), mart = mart)
	# remove genes without hugo symbols
	hugo_gene_names = gene_names[gene_names %in% tss_map$hgnc_symbol]
	if(length(hugo_gene_names) != length(gene_names)){
		print("The following genes are not in HGNC and TSS cannot be determined:")
		print(setdiff(gene_names, hugo_gene_names))
	}
	return(getBM(attributes = c("hgnc_symbol", "chromosome_name", "strand", "transcript_start"), filters = "hgnc_symbol", values = hugo_gene_names, mart = mart))
}

get_seq <- function(start_relative2tss, end_relative2tss, tss_pos, strand, chrom){
	start_pos = tss + strand*start_relative2tss + 1;
	end_pos = tss + strand*end_relative2tss;
	chrom = paste0("chr", unique(na.omit(as.numeric(unlist(strsplit(unlist(chrom), "[^0-9]+")))))[1])
	wanted_range = GRanges(chrom, IRanges(min(start_pos, end_pos), max(start_pos, end_pos)))
    return(getSeq(Hsapiens, wanted_range))
}

