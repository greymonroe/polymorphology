% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_features.R
\name{variant_features}
\alias{variant_features}
\title{annotate variants with features from gff}
\usage{
variant_features(vcf, gff, type, ann = "info", window = 0)
}
\arguments{
\item{vcf}{data.table object of variant data (data.table with CHROM, POS, REF, and ALT columns)}

\item{gff}{data.table object of gff data (data.table with "chr","source","type","start","stop","V1","direction","V2","info" columns)}

\item{type}{type of feature (e.g, "CDS","gene", etc)}

\item{ann}{which column to annotate with (values will be concatenated to track the names of the features a varint is found in), default is "info" column of gff}

\item{window}{number of base pairs to expand (e.g. to include regions +- 3000 base pairs from feature you set 3000)}
}
\value{
a vector of pasted (" AND " delimited) annotations (based on column provided by ) for each feature the variant is found in
}
\description{
annotate variants with features from gff
}
