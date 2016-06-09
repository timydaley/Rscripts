# usage:  source("gene2seq.R); gene_names2target_regions("input.txt", "output", -50, -400);
# input is a single column file of gene names
# ouput is a fasta file containing the sequence of the target region defined by (start_relative2tss, end_relative2tss) in hg38 (latest human genome build)


# load required libraries
install.packages("seqinr", repos = "http://cran.r-project.org")
library(seqinr)
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

chrom_name <- function(chrom){
	if(grepl("CHR", chrom)){
		chrom = unique(na.omit(as.numeric(unlist(strsplit(unlist(chrom), "[^0-9]+")))))[1]
	}
	return(paste0("chr", chrom))
}

get_seqs <- function(start_relative2tss, end_relative2tss, tss_pos, strand, chrom, gene_names){
	start_pos = tss_pos + as.numeric(strand)*start_relative2tss + 1;
	end_pos = tss_pos + as.numeric(strand)*end_relative2tss;
	chroms = sapply(chrom, function(x) chrom_name(x));
	wanted_ranges = GRanges(chroms, IRanges(apply(cbind(start_pos, end_pos), 1, min) , apply(cbind(start_pos, end_pos), 1, max)))
	seqs = c()
	for(i in 1:length(start_pos)){
		seqs = c(seqs, getSeq(Hsapiens, wanted_ranges[i], as.character=TRUE))
	}
    return(list(genes = gene_names, seqs = seqs, tss = tss_pos))
}

write_seqs <- function(seqs, gene_names, tss, filename){
	stopifnot(dim(seqs)[1] == length(gene_names))
	write.fasta(file.out = filename, sequences = seqs[1], names = gene_names[1], open = "w", nbchar = 80, as.string = TRUE)
	if(length(gene_names) > 1){
	    for(i in 2:length(gene_names)){
	        write.fasta(file.out = filename, sequences = seqs[i], names = paste(gene_names[i], "\t", tss[i]), open = "a", nbchar = 80, as.string = TRUE)
	   	}
	}
}

gene_names2target_regions <- function(input_filename, output_filename, start_relative2tss, end_relative2tss){
	gene_names = as.vector(t(read.table(file = input_filename)))
	tss = get_top_tss(gene_names)
	gene_plus_seq = get_seqs(start_relative2tss, end_relative2tss, tss$transcript_start, tss$strand, tss$chromosome_name, tss$hgnc_symbol)
	n_splits = ceiling(length(gene_plus_seq$genes)/10)
	for(i in 1:n_splits){
	    write_seqs(gene_plus_seq$seqs[((i - 1)*25 + 1):min(i*25, length(gene_plus_seq$seqs))], gene_plus_seq$genes[((i - 1)*25 + 1):min(i*25, length(gene_plus_seq$seqs))], gene_plus_seq$tss[((i - 1)*25 + 1):min(i*25, length(gene_plus_seq$seqs))], paste0(output_filename, i, ".txt"))
	 }
}

