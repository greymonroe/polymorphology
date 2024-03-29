% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tss_tts.calc.R
\name{tss_tts.calc}
\alias{tss_tts.calc}
\title{calculate minor allele frequencies of polymorphisms in relation to gene bodies}
\usage{
tss_tts.calc(
  gff,
  vcf,
  type = "both",
  chrs = "all",
  num = "all",
  feature = "gene",
  variable
)
}
\arguments{
\item{gff}{directory and file name of reference gff OR data.table object of gff data}

\item{vcf}{directory and file name vcf data OR data.table object of variant data (data.table with CHROM, POS, MAF)}

\item{type}{"snp", "indel", or "both", default is "both"}

\item{chrs}{selection specific chromosomes, default is "all"}

\item{num}{set number of genes, default is "all}

\item{feature}{name of the feature in the gff data that you want to look at, default is "gene"}

\item{variable}{name of variable in the vcf you want the calculated value of in relation to TSS TTS}
}
\value{
a vector of the relative position of polymorphsims in relation to feature
}
\description{
calculate minor allele frequencies of polymorphisms in relation to gene bodies
}
