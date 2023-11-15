# Get DNA sequence of genes of interest for CpG density Calculation!
# Its needed to recalculate CpG density like observed to expected CpG ratio:
# To do this read following paper: "https://www.sciencedirect.com/science/article/pii/S0002929720302445"
# and consider following statement: CpG density was calculated as the observed to expected ratio:
"(OE = [number of CpGs/{number of Cs × number of Gs}] × length of the region in nucleotides).
The red line indicates the threshold in OE used in the standard definition of CpG islands 
(Gardiner-Garden and Frommer 1987)."

# Install libraries

library("BSgenome.Mmusculus.UCSC.mm10")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.ensGene")
library("GenomicFeatures")
library("data.table")
library("stringr")
library("org.Mm.eg.db")


# Get gene types from NCBI
# NCBI > Gene > Download/FTP > DATA > gene_info.gz
# Download from: "https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz"
# Note that this file contains info for all genes/species
# Taxanomy ID of Mus musculus is 10090
# Taxanomy ID of Homo sapiens is 9606

gene_info <- read.delim("/path_to_directory/gene_info")

gene_info_Mus_musculus <- subset(gene_info, gene_info$X.tax_id == "10090")
gene_info_Mus_musculus_protein_coding_genes <- subset(gene_info_Mus_musculus,
                                                      gene_info_Mus_musculus$type_of_gene == "protein-coding")

# In NCBI, there is 26340 Mus musculus protein-coding genes!

# Get UCSC mm10 gene coordinates and gene IDS as GRanges object

genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene) # to get NCBI-based gene IDs
genes2 <- genes(TxDb.Mmusculus.UCSC.mm10.ensGene) # to get ENSEMBL-based gene IDs

# Get coordinates of Mus Musculus protein-coding genes!

Mm_protein_coding_genes_coordinates <- subset(genes, genes$gene_id %in% gene_info_Mus_musculus_protein_coding_genes$GeneID)

# I can't get coordinates of all 26340 protein coding genes because of multi isoforms which makes hard to get
# true coordinates of all genes, so I obtained coordinates of only 20595 out of 26340 protein coding genes!

# Filter out unconventional chromosomes names like "chrUn_GL456396"...

Mm_protein_coding_genes_coordinates <- subset(Mm_protein_coding_genes_coordinates,
                                              Mm_protein_coding_genes_coordinates@seqnames %in% c("chr1", "chr2",
                                                                                                  "chr3", "chr4",
                                                                                                  "chr5","chr6",
                                                                                                  "chr7", "chr8",
                                                                                                  "chr9", "chr10",
                                                                                                  "chr11", "chr12",
                                                                                                  "chr13", "chr14",
                                                                                                  "chr15","chr16",
                                                                                                  "chr17", "chr18",
                                                                                                  "chr19", "chrX", "chrY"))

Mm_protein_coding_genes_coordinates_as_data_frame_without_chromosome_names <- as.data.frame(Mm_protein_coding_genes_coordinates@ranges)
Mm_protein_coding_genes_chromosome_coordinates <- as.data.frame(Mm_protein_coding_genes_coordinates@seqnames)
colnames(Mm_protein_coding_genes_chromosome_coordinates) <- "chr"
Mm_protein_coding_genes_strand_info <- as.data.frame(Mm_protein_coding_genes_coordinates@strand)
colnames(Mm_protein_coding_genes_strand_info) <- c("strand_info_of_genes")

# Get DNA sequence of Mus musculus protein-coding genes

Mus_musculus_DNA_sequence_of_protein_coding_genes <- as.character(Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, Mm_protein_coding_genes_chromosome_coordinates$chr,
                                                                                     Mm_protein_coding_genes_coordinates_as_data_frame_without_chromosome_names$start,
                                                                                     Mm_protein_coding_genes_coordinates_as_data_frame_without_chromosome_names$end, strand = Mm_protein_coding_genes_strand_info$strand_info_of_genes))

Mus_musculus_DNA_sequence_of_protein_coding_genes[1] # See DNA sequence of first coordinate
Mus_musculus_DNA_sequence_of_protein_coding_genes[2] # See DNA sequence of second coordinate
Mus_musculus_DNA_sequence_of_protein_coding_genes[3] # See DNA sequence of third coordinate

# Export 20585 Mus musculus protein coding genes' DNA sequence (UCSC mm10 genome)

Mus_musculus_DNA_sequence_of_protein_coding_genes2 <- as.data.frame(Mus_musculus_DNA_sequence_of_protein_coding_genes)

# Add chromosome, start, and end columns to each gene. Each row is gene and its DNA sequence

Mus_musculus_DNA_sequence_of_protein_coding_genes2$chr <- Mm_protein_coding_genes_chromosome_coordinates$chr
Mus_musculus_DNA_sequence_of_protein_coding_genes2$start <- Mm_protein_coding_genes_coordinates_as_data_frame_without_chromosome_names$start
Mus_musculus_DNA_sequence_of_protein_coding_genes2$end <- Mm_protein_coding_genes_coordinates_as_data_frame_without_chromosome_names$end
Mus_musculus_DNA_sequence_of_protein_coding_genes2$geneID <- Mm_protein_coding_genes_coordinates_as_data_frame_without_chromosome_names$names
Mus_musculus_DNA_sequence_of_protein_coding_genes2$strand <- Mm_protein_coding_genes_strand_info$strand_info_of_genes

# Re-assign column names

colnames(Mus_musculus_DNA_sequence_of_protein_coding_genes2) <- c("DNA_sequences_gene_of_interest_from_mm10_UCSC", "chr", "start", "end", "geneID", "strand")

# Convert Entrez geneIDs to gene names

write.table(Mus_musculus_DNA_sequence_of_protein_coding_genes2$geneID, "/path_to_directory/Mus_musculus_DNA_sequence_of_protein_coding_genes2$geneID.txt", row.names = F, quote = F)
file_path <- '/path_to_directory/Mus_musculus_DNA_sequence_of_protein_coding_genes2$geneID.txt'
gene_ids <- readLines(file_path) # read rows in text file containing gene IDs. Every row has one geneID.

gene_names <- mapIds(org.Mm.eg.db, gene_ids, "SYMBOL", 'ENTREZID')
gene_names <- as.data.frame(gene_names)

# Create a column which has gene names

Mus_musculus_DNA_sequence_of_protein_coding_genes2$gene_names <- gene_names$gene_names

write.csv(Mus_musculus_DNA_sequence_of_protein_coding_genes2, "/path_to_directory/Mus_musculus_UCSC_mm10_coordinates_DNA_sequence_of_protein_coding_genes.csv", quote = F, row.names = F)

# Note/ To get two different coordinates' DNA sequences at once
# my.dnastring <- as.character(Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, c("chr15", "chr16"), c(81585397, 4081328), c(81652077, 4213997)))

# Split genes according to strand info

Forward_strand_genes <- subset(Mus_musculus_DNA_sequence_of_protein_coding_genes2,
                               Mus_musculus_DNA_sequence_of_protein_coding_genes2$strand == "+")

# 10,279 genes are in forward strand of DNA

Reverse_strand_genes <- subset(Mus_musculus_DNA_sequence_of_protein_coding_genes2,
                               Mus_musculus_DNA_sequence_of_protein_coding_genes2$strand == "-")

# 10,306 genes are in reverse strand of DNA

# Calculate CpG density of forward genes:

Forward_strand_genes$number_of_Cytosines <- lengths(regmatches(Forward_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC, gregexpr("C", Forward_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC)))
Forward_strand_genes$number_of_Guanine <- lengths(regmatches(Forward_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC, gregexpr("G", Forward_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC)))
Forward_strand_genes$number_of_CpGs <- lengths(regmatches(Forward_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC, gregexpr("CG", Forward_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC)))
Forward_strand_genes$ratio_of_CpGs <- (Forward_strand_genes$number_of_CpGs/(Forward_strand_genes$number_of_Cytosines + Forward_strand_genes$number_of_Guanine))*100

# Calculate CpG density of reverse genes:

# function that reverses a string by characters
# This function for reversing DNA sequence of genes that are in minus strand to count CpGs properly.

reverse_chars <- function(string)
{
  # split string by characters
  string_split = strsplit(string, split = "")
  # reverse order
  rev_order = nchar(string):1
  # reversed characters
  reversed_chars = string_split[[1]][rev_order]
  # collapse reversed characters
  paste(reversed_chars, collapse = "")
}

Reverse_strand_genes$reversed_DNA_sequence <- ""
Reverse_strand_genes$reversed_DNA_sequence <- lapply(Reverse_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC, reverse_chars)

# Calculate CpG density of reverse genes:

Reverse_strand_genes$number_of_Cytosines <- lengths(regmatches(Reverse_strand_genes$reversed_DNA_sequence, gregexpr("C", Reverse_strand_genes$reversed_DNA_sequence)))
Reverse_strand_genes$number_of_Guanine <- lengths(regmatches(Reverse_strand_genes$reversed_DNA_sequence, gregexpr("G", Reverse_strand_genes$reversed_DNA_sequence)))
Reverse_strand_genes$number_of_CpGs <- lengths(regmatches(Reverse_strand_genes$reversed_DNA_sequence, gregexpr("CG", Reverse_strand_genes$reversed_DNA_sequence)))
Reverse_strand_genes$ratio_of_CpGs <- (Reverse_strand_genes$number_of_CpGs/(Reverse_strand_genes$number_of_Cytosines + Reverse_strand_genes$number_of_Guanine))*100

# Combine reverse and forward genes

# Replace un-reversed DNA sequences to reversed DNA sequences for reverse genes
Reverse_strand_genes$DNA_sequences_gene_of_interest_from_mm10_UCSC <- Reverse_strand_genes$reversed_DNA_sequence
Reverse_strand_genes <- Reverse_strand_genes[, c(1:7,9:12)]

CpG_ratio_for_all_Mus_musculus_protein_coding_genes <- rbind(Forward_strand_genes, Reverse_strand_genes)

# Check type of columns if they are list, convert them to characters and then dataframe:

apply(CpG_ratio_for_all_Mus_musculus_protein_coding_genes, 2, typeof)
CpG_ratio_for_all_Mus_musculus_protein_coding_genes <- apply(CpG_ratio_for_all_Mus_musculus_protein_coding_genes, 2, as.character)
CpG_ratio_for_all_Mus_musculus_protein_coding_genes <- as.data.frame(CpG_ratio_for_all_Mus_musculus_protein_coding_genes)
write.csv(CpG_ratio_for_all_Mus_musculus_protein_coding_genes, "/path_to_directory/CpG_ratio_for_all_Mus_musculus_protein_coding_genes_regarding_strand_info.csv", quote = F, row.names = F)

# Calculate ten deciles of CpG ratio to find which genes are in 1th, 2th, 3rd ..... and 10th deciles?

CpG_ratio_for_all_Mus_musculus_protein_coding_genes <- fread("/path_to_directory/CpG_ratio_for_all_Mus_musculus_protein_coding_genes_regarding_strand_info.csv")

# Keep only gene name and CpG ratio on columns for decile calculation!

CpG_ratio_for_all_Mus_musculus_protein_coding_genes2 <- CpG_ratio_for_all_Mus_musculus_protein_coding_genes[,c(7, 11)]

# Calculate deciles
CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$Decile <- cut(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$ratio_of_CpGs, quantile(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$ratio_of_CpGs, probs = 0:10/10), include.lowest = TRUE, labels = FALSE)

# Order genes within each decile
CpG_ratio_for_all_Mus_musculus_protein_coding_genes2 <- CpG_ratio_for_all_Mus_musculus_protein_coding_genes2[order(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$Decile, CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$ratio_of_CpGs), ]

# Print the resulting data frame
print(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2)

write.csv(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2, "/path_to_directory/Ten_decile_of_CpG_ratio_for_all_Mus_musculus_protein_coding_genes.csv", quote = F, row.names = F)

# To see how differentially expressed genes of Kmt2a distributed into deciles!
# They are mostly in higher deciles or lower? What is mean of it?

Kmt2a_DEGs_with_adjusted_pvalue_less_than_005 <- fread("/path_to_directory/Kmt2a_DEGs_with_adjusted_pvalue_less_than_005.csv")


Kmt2a_DEGs_and_Deciled_CpG_ratio <- subset(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2,
                                           CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$gene_names %in% Kmt2a_DEGs_with_adjusted_pvalue_less_than_005$V1)

write.table(Kmt2a_DEGs_and_Deciled_CpG_ratio, "/path_to_directory/Kmt2a_DEGs_and_Deciled_CpG_ratio.csv", quote = F, row.names = F)

# 319 of 335 Kmt2a DEGs are protein-coding!

# Calculate number of DEGs in each decile divided by total number of genes in each decile

Number_of_DEGs_in_1st_decile <- 40
Number_of_DEGs_in_2nd_decile <- 36
Number_of_DEGs_in_3rd_decile <- 21
Number_of_DEGs_in_fourth_decile <- 29
Number_of_DEGs_in_fifth_decile <- 21
Number_of_DEGs_in_sixth_decile <- 41
Number_of_DEGs_in_seventh_decile <- 53
Number_of_DEGs_in_eight_decile <- 22
Number_of_DEGs_in_nineth_decile <- 21
Number_of_DEGs_in_tenth_decile <- 35

# Total number of genes in each decile

table(CpG_ratio_for_all_Mus_musculus_protein_coding_genes2$Decile)

Total_Number_of_genes_in_1st_decile <- 2059
Total_Number_of_genes_in_2nd_decile <- 2058
Total_Number_of_genes_in_3rd_decile <- 2059
Total_Number_of_genes_in_fourth_decile <- 2058
Total_Number_of_genes_in_fifth_decile <- 2059
Total_Number_of_genes_in_sixth_decile <- 2058
Total_Number_of_genes_in_seventh_decile <- 2058
Total_Number_of_genes_in_eight_decile <- 2060
Total_Number_of_genes_in_nineth_decile <- 2057
Total_Number_of_genes_in_tenth_decile <- 2059

# Ratio of Number of DEGs in each decile/Total number of genes in each decile

first_decile_ratio <- (40/2059)*100
second_decile_ratio <- (36/2058)*100
third_decile_ratio <- (21/2059)*100
fourth_decile_ratio <- (29/2058)*100
fifth_decile_ratio <- (21/2059)*100
sixth_decile_ratio <- (41/2058)*100
seventh_decile_ratio <- (53/2058)*100
eigth_decile_ratio <- (22/2060)*100
nineth_decile_ratio <- (21/2057)*100
tenth_decile_ratio <- (35/2059)*100

data2 <- data.frame(Vectors = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                    Values = c(first_decile_ratio, second_decile_ratio, third_decile_ratio, fourth_decile_ratio, fifth_decile_ratio, sixth_decile_ratio, seventh_decile_ratio, eigth_decile_ratio, nineth_decile_ratio, tenth_decile_ratio))

# Create a PDF device
pdf("/path_to_directory/Percent_of_Kmt2a_DEGs_in_each_CpG_ratio_decile2.pdf", width = 8, height = 6)

# Create the dot plot
par(mar = c(5, 7, 4, 2))  # Adjust margins
plot(data2$Values, axes = FALSE, xlab = "", ylab = "",
     pch = 16, col = "dodgerblue", cex = 2, ylim = c(0, 4.5))
axis(1, at = 1:10, labels = data2$Vectors)
axis(2, las = 1, at = c(0, 1.5, 3, 4.5), labels = c("0", "1.5", "3", "4.5"))

# Add a title
title(main = "Percentage of Kmt2a DEGs in each decile of CpG ratio")

# Close the PDF device
dev.off()

# pLI scores of DEGs vs CpG density of DEGs

# Download pLI scores from gnomAD:
# "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

gnomAD_pLI_scores <- fread("/path_to_directory/gnomad.v2.1.1.lof_metrics.by_gene.txt")

Kmt2a_DEGs_and_Deciled_CpG_ratio_and_pLI_scores <- subset(gnomAD_pLI_scores,
                                                          gnomAD_pLI_scores$gene %in% toupper(Kmt2a_DEGs_and_Deciled_CpG_ratio$gene_names))

write.table(Kmt2a_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, "/path_to_directory/Kmt2a_DEGs_and_Deciled_CpG_ratio_and_pLI_scores.csv", quote = F, row.names = F)

# Scatter plot of 319 Kmt2a DEGs pLI vs CpG density

Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores <- fread("/path_to_directory/Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores.csv")

ggplot(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, aes(x = Decile, y = pLI_score)) +
  geom_point() +
  labs(x = "Decile CpG Density", y = "pLI Score", title = "Scatter Plot of CpG Density vs. pLI Score") +
  theme_minimal()

# Save the plot as a PDF file
ggsave("/path_to_directory/Kmt2a DEGs CpG Density Decile vs pLI Score.pdf", width = 8, height = 6)

# Average pLI score in each CpG decile

first_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 1)$pLI_score, na.rm = TRUE)
second_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 2)$pLI_score, na.rm = TRUE)
third_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 3)$pLI_score, na.rm = TRUE)
fourth_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 4)$pLI_score, na.rm = TRUE)
fifth_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 5)$pLI_score, na.rm = TRUE)
sixth_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 6)$pLI_score, na.rm = TRUE)
seventh_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 7)$pLI_score, na.rm = TRUE)
eigth_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 8)$pLI_score, na.rm = TRUE)
nineth_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 9)$pLI_score, na.rm = TRUE)
tenth_decile_avarage_pLI <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 10)$pLI_score, na.rm = TRUE)

data3 <- data.frame(Vectors = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                    Values = c(first_decile_avarage_pLI, 
                               second_decile_avarage_pLI, 
                               third_decile_avarage_pLI, 
                               fourth_decile_avarage_pLI, 
                               fifth_decile_avarage_pLI, 
                               sixth_decile_avarage_pLI, 
                               seventh_decile_avarage_pLI, 
                               eigth_decile_avarage_pLI, 
                               nineth_decile_avarage_pLI, 
                               tenth_decile_avarage_pLI))

# Create a PDF device
pdf("/path_to_directory/average_pLI_score_in_each_CpG_decile_for_Kmt2a_DEGs.pdf", width = 8, height = 6)

# Create the dot plot
par(mar = c(5, 7, 4, 2))  # Adjust margins
plot(data3$Values, axes = FALSE, xlab = "", ylab = "",
     pch = 16, col = "dodgerblue", cex = 2)
axis(1, at = 1:10, labels = data2$Vectors)
axis(2, las = 1)

# Add a title
title(main = "average pLI score in each CpG decile for Kmt2a DEGs")

# Close the PDF device
dev.off()

# Average LOUEF score in each CpG decile

first_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 1)$LOUEF_score, na.rm = TRUE)
second_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 2)$LOUEF_score, na.rm = TRUE)
third_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 3)$LOUEF_score, na.rm = TRUE)
fourth_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 4)$LOUEF_score, na.rm = TRUE)
fifth_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 5)$LOUEF_score, na.rm = TRUE)
sixth_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 6)$LOUEF_score, na.rm = TRUE)
seventh_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 7)$LOUEF_score, na.rm = TRUE)
eigth_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 8)$LOUEF_score, na.rm = TRUE)
nineth_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 9)$LOUEF_score, na.rm = TRUE)
tenth_decile_avarage_LOUEF_score <- mean(subset(Kmt2a_protein_coding_DEGs_and_Deciled_CpG_ratio_and_pLI_scores, Decile == 10)$LOUEF_score, na.rm = TRUE)

data4 <- data.frame(Vectors = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                    Values = c(first_decile_avarage_LOUEF_score, 
                               second_decile_avarage_LOUEF_score, 
                               third_decile_avarage_LOUEF_score, 
                               fourth_decile_avarage_LOUEF_score, 
                               fifth_decile_avarage_LOUEF_score, 
                               sixth_decile_avarage_LOUEF_score, 
                               seventh_decile_avarage_LOUEF_score, 
                               eigth_decile_avarage_LOUEF_score, 
                               nineth_decile_avarage_LOUEF_score, 
                               tenth_decile_avarage_LOUEF_score))

# Create a PDF device
pdf("/path_to_directory/average_LOUEF_score_in_each_CpG_decile_for_Kmt2a_DEGs.pdf", width = 8, height = 6)

# Create the dot plot
par(mar = c(5, 7, 4, 2))  # Adjust margins
plot(data4$Values, axes = FALSE, xlab = "", ylab = "",
     pch = 16, col = "dodgerblue", cex = 2)
axis(1, at = 1:10, labels = data2$Vectors)
axis(2, las = 1)

# Add a title
title(main = "average LOUEF score in each CpG decile for Kmt2a DEGs")

# Close the PDF device
dev.off()
