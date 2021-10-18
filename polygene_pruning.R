#### STAGE 1: CLUMPING ####

#### run clumping with PLINK
# here we run it chromosome-wise

# plink --bfile /your/bfile/prefix_chr --clump bp.assoc --clump-range genes4plink.dat --clump-range-border 5000 --clump-p1 0.0001 --clump-p2 0.05 --clump-r2 0.2 --clump-kb 5000 --clump-allow-overlap --clump-best --out plink_chr

#### manage clumping results

### select unique snps and sort by P values

library(data.table)

ch0 <- paste0('chr', 1:22)

allsnps.all <- c()

for (ch in ch0) {

print(ch)
cl <- fread(paste0('plink_', ch, '.clumped'), h = T, data.table = F) # PLINK output file
cl$SNP <- as.character(cl$SNP)
cl$SP2 <- as.character(cl$SP2)

allsnps <- cl$SNP
for (i in 1:dim(cl)[1]) {
	if (!i%%1000) cat(i, '')
	rs <- cl$SP2[i]
	rs <- gsub('(1)', '', rs, fixed = T)
	rs <- strsplit(rs, ',')
	rs <- unlist(rs)
	allsnps <- c(allsnps, rs)
}
allsnps <- unique(allsnps)

allsnps.all <- c(allsnps.all, allsnps)

}
save(allsnps.all, file = paste0('allsnps.RDa'))

### subset clumped SNPs from GWAS results data

ch0 <- paste0('chr', 1:22)

gwas <- fread('summary.statistics.txt', h = T, data.table = F) # a text file in the format of the prep.score.files() input file 
gwas <- gwas[, c('CHROM', 'POS', 'ID', 'P')]
table(gwas$P <= 0.05)

allsnps <- get(load('allsnps.RDa'))
length(allsnps)
gwas0 <- gwas
gwas <- gwas[gwas$ID %in% allsnps, ]

### for each top SNP, create a list of clumped SNPs with P > PtopSNP

clumping.list0 <- list()

for (ch in ch0) {

	print(ch)

	cl <- fread(paste0('plink_', ch, '.clumped'), h = T, data.table = F)
	cl$SNP <- as.character(cl$SNP)
	cl$SP2 <- as.character(cl$SP2)

	clumping.list <- list()

	for (i in 1:dim(cl)[1]) {
		if (!i%%1000) cat(i, '')
		topsnp <- cl$SNP[i]
		rs <- cl$SP2[i]
		rs <- gsub('(1)', '', rs, fixed = T)
		rs <- strsplit(rs, ',')
		rs <- unlist(rs)
		ptopsnp <- gwas[gwas[, 'ID'] == topsnp, 'P']
		gwas.tmp <- gwas[gwas[, 'ID'] %in% rs, ]
		rs <- gwas.tmp[gwas.tmp[, 'P'] > ptopsnp, 'ID']
		clumping.list[[i]] <- rs
	}
	
	names(clumping.list) <- cl$SNP
	clumping.list0 <- c(clumping.list0, clumping.list)

}

save(clumping.list0, file = paste0('clumping.list0.RDa'))

### merge clumped.ranges for all chromosomes

system('cat plink_chr1.clumped.ranges > plink.clumped.ranges')
for (ch in ch0[2:22]) {
	system(paste0('tail -n +2 plink_', ch, '.clumped.ranges >> plink.clumped.ranges'))
}


#### STAGE 2: PRUNING ####

library(data.table)
library(sumFREGAT)
corpath <- '/path/to/cor.files' # LD matrices
geneFile <-  system.file("testfiles/NCBI37.3.geneFile.txt.gz", package = "sumFREGAT") # for build 37
#geneFile <-  system.file("testfiles/ensembl.hg38.txt", package = "sumFREGAT") # for build 38

# read genefile
df0 <- read.table(geneFile)
df0$V2 <- df0$V4 <- NULL
colnames(df0)[1:4] <- c('gene', 'chr', 'start', 'stop')
df0$chr <- gsub('chr', '', df0$chr)
df0$gene <- as.character(df0$gene)
# 5e6 bp threshold:
df0$start2 <- df0$start - 5e6
df0$stop2 <- df0$stop + 5e6

clumping.list <- get(load('clumping.list0.RDa')) # file from clumping.R
ran <- read.table('plink.clumped.ranges', h = T) # file from plink
ran$SNP <- as.character(ran$SNP)

#type0 <- c('.nsyn', '.exon', '.intron', '')
meth0 <- c('SKATO11', 'PCA11', 'ACAT11')

# read a reference file for UKBB LD matrices, needed for prep.score.files()
ref <- get(load('path/to/reference/file/ref.ldstore265.all.RData'))

#for (type in type0) { # cycle for annotation scenarios
type <- ''

cat(type, '')

tmp <- gsub('.', '', type, fixed = T)
if (type == '') tmp <- 'all'
fn <- paste0('rerun.genes', type, '.RDa')

# read GWAS stats
gwas <- fread(paste0('/path/to/summary.statistics', type, '.txt'), h = T, data.table = F)# variable file name, specific for a scenario 
# keep only SNPs present in 1KG (only if 1KG reference is used)
#map1kg <- fread('/path/to/g1000p3_EUR/origin/g1000p3_EUR.bim', data.table = F) # 1KG snp list
#v <- gwas$ID %in% map1kg$V2
#table(v)
#gwas <- gwas[v, ]
# remove SNPs monomorphic in 1KG
#freq <- fread('/home/common/Neurotism_Nagel_2018/sumFREGAT/res/r2/g1000p3_EUR.frq', data.table = F, h = T)
#freq <- freq[freq$MAF > 0, ]
#v <- gwas$ID %in% freq$SNP
#table(v)
#gwas <- gwas[v, ]

######## read input files for gene-based analysis

stats1 <- paste0('/path/to/summary.statistics', type, '.txt') # input name, same as 'gwas' for type = ''
stats2 <- paste0('../MV_Back_Disc_gwas_stats', type, '.clump.txt') # output name
stats3 <- paste0('../scores.stats', type, '.clump') # name for a tempopary file

prep0 <- fread(stats1, h = T, data.table = F)

# keep only genes passed 2.5e-5 in prior analysis

rerun.genes <- read.table(paste0('path/to/results/', tmp, '.combined.2.5e-5.csv'), h = T, sep = ';', as.is = T)$gene
df <- df0[df0$gene %in% rerun.genes, ]

# pruning

write.table('gene chrom start end markers filtered.markers pvalue', file = 'ACAT11.pruned.out', col = F, row = F, qu = F)
write.table('gene chrom start end markers filtered.markers pvalue ncomp varcomp', file = 'PCA11.pruned.out', col = F, row = F, qu = F)
write.table('gene chrom start end markers filtered.markers pvalue', file = 'SKATO11.pruned.out', col = F, row = F, qu = F)

rerun.genes <- c()
snps2keep <- c()

for (i in 1:dim(df)[1]) {

	# this part is equal for all scenarios, except for genes listed in 'df'
	# it can be done once for all scenarios if the same genes are analysed
	g <- df[i, 'gene']
	cat(g, '')
	topsnps <- ran[grep(paste0('\\b', g, '\\b'), ran[, 'RANGES']), 'SNP']
	vvv <- gwas$CHR == df$chr[i] & gwas$POS >= df$start[i] & gwas$POS <= df$stop[i]
	genesnps <- gwas[vvv, 'ID']
	topsnps <- topsnps[!topsnps %in% genesnps] # outer topsnps only
	snps2waste <- c()
	for (topsnp in topsnps) {
		snps2waste <- c(snps2waste, clumping.list[[topsnp]]) # correlated snps from clumps
	}
	if (sum(genesnps %in% snps2waste) == 0) next
	# if there are snps to exclude:
	cat('!!! ')
	rerun.genes <- c(rerun.genes, g)
	v <- !genesnps %in% snps2waste
	genesnps <- genesnps[v]
	if (sum(v) == 0) { # all snps excluded
		write.table(t(c(g, NA, NA, NA, 0, 0, NA)), file = 'ACAT11.pruned.out', col = F, row = F, qu = F, app = T)
		write.table(t(c(g, NA, NA, NA, 0, 0, NA)), file = 'SKATO11.pruned.out', col = F, row = F, qu = F, app = T)
		write.table(t(c(g, NA, NA, NA, 0, 0, NA, NA, NA)), file = 'PCA11.pruned.out', col = F, row = F, qu = F, app = T)
		next
	}
	snps2keep <- rbind(snps2keep, cbind(g, genesnps))
	# end of the part equal for all scenarios, except for genes list in 'df'

	######## write new files with pruned snps only
	prep <- prep0[prep0$rs %in% genesnps, ]
	write.table(prep, file = stats2, row = F, qu = F)

	######## run prep.score.files()
	prep.score.files(stats2, reference = ref, output.file.prefix = stats3)

	reg <- ACAT(score.file = stats3, gene.file = geneFile, genes = g, write.file = 'ACAT11.tmp.out', quiet = T)
	system('tail -1 ACAT11.tmp.out >> ACAT11.pruned.out')
	reg <- SKAT(score.file = stats3, cor.path = corpath, gene.file = geneFile, beta.par = c(1, 1), method = "kuonen", rho = TRUE, genes = g, write.file = 'SKATO11.tmp.out', quiet = T)
	system('tail -1 SKATO11.tmp.out >> SKATO11.pruned.out')
	reg <- PCA(score.file = stats3, cor.path = corpath, gene.file = geneFile, n = 380506, genes = g, write.file = 'PCA11.tmp.out', quiet = T)
	system('tail -1 PCA11.tmp.out >> PCA11.pruned.out')

}

save(rerun.genes, snps2keep, file = fn)

system(paste0('mv ACAT11.pruned.out ACAT11', type, '.out'))
system(paste0('mv SKATO11.pruned.out SKATO11', type, '.out'))
system(paste0('mv PCA11.pruned.out PCA11', type, '.out'))

#}