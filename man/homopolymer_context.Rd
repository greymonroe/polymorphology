% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homopolymer_context.R
\name{homopolymer_context}
\alias{homopolymer_context}
\title{homopolymer_context}
\usage{
homopolymer_context(vars, genome, homopolymers, size)
}
\arguments{
\item{vars}{mutations with CHROM REF ALT POS columns}

\item{genome}{genome sequence (from seqinr)}

\item{homopolymers}{homopolymer object from make_homopolymer()}

\item{size}{size of homopolymer to consider}
}
\value{
a data table
}
\description{
homopolymer_context
}
