# Preprocessing 
# Files pre-processed using ScPipe package

library(scPipe)
# Barcode annotation
barcode_annotation <- "barcode_annotation.csv"

# Parameters
bc.len <- 8
UMI.len <- 10
UMI.cor <- 1
threads <- 16

# Genome, annotation file, and Index file
genome <- "GCF_000001405.40_GRCh38.p14_genomic.fna"
genome.index <- "Index/GRCH_38"
gff3 <- "GCF_000001405.40_GRCh38.p14_genomic.gff"
gtf <- "GCF_000001405.40_GRCh38.p14_genomic.gtf"

# Output folder and files
out_dir <- ""
combined_fastq <- file.path(out_dir, "trimmed.fastq.gz")
aligned_bam <- file.path(out_dir, "aligned.bam")
mapped_bam <- file.path(out_dir, "aligned.mapped.bam")

# Raw reads: R2 has the transcript sequence, R1 has the barcode and UMI
fq.transcript <- "MultiplexRNASeq_S1_R2_001.fastq.gz"
fq.barcode <- "MultiplexRNASeq_S1_R1_001.fastq.gz"


# Analysis
# Read structure; UMI Barcode bs2 = 0 bl2 = 8 us = 8 ul = 10
sc_trim_barcode(
  outfq = combined_fastq,
  r1 = fq.transcript,
  r2 = fq.barcode,
 read_structure = list(
    bs1 = -1, bl1 = 0,
    bs2 = 0, bl2 = 8,
   us = 8, ul = 10
  )
  )

# Align reads
Rsubread::align(
 index = genome.index,
 readfile1 = combined_fastq,
 output_file = aligned_bam,
 nthreads = threads
)

# Assigning reads to annotated exons
sc_exon_mapping(
  inbam = aligned_bam,
  outbam = mapped_bam,
  annofn = gff3,
  bc_len = bc.len,
  UMI_len = UMI.len,
  nthreads = threads
)

# Demultiplex
sc_demultiplex(
  inbam = mapped_bam,
  outdir = out_dir,
  bc_anno = barcode_annotation,
  has_UMI = TRUE,
  nthreads = threads
)

# Gene counting
sc_gene_counting(
  outdir = out_dir,
  bc_anno = barcode_annotation,
  UMI_cor = 1
)

