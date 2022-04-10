# Counts in gff object

#' add counts of variants to gff data object
#' @param regions data.table object with columns chr, start, stop
#' @param vcf data.table object of variant data (data.table with CHROM, POS, REF, and ALT columns)
#' @return a vector of variant counts for each row of input gff
#' @export

region_counts<-function(regions, vcf){

   if(!all(colnames(regions)[1:3] == c("chr","start","stop" ))){
    stop('First 3 column names of regions must be "chr","start","stop"')
  }
  if(!"CHROM" %in% colnames(vcf) | !"POS" %in% colnames(vcf)){
    stop('CHROM and POS must be colnames in vcf')
  }
  chromosomes_vcf<-unique(vcf$CHROM)
  chromosomes_regions<-unique(regions$chr)
  if(!all(chromosomes_vcf %in% chromosomes_regions) | !all(chromosomes_regions %in% chromosomes_vcf)){
    warning("chromosomes do not match between regions and vcf")
    warning(paste("chromosomes unique to vcf:", chromosomes_vcf[!chromosomes_vcf %in% chromosomes_regions]))
    warning(paste("chromosomes unique to regions:", chromosomes_regions[!chromosomes_regions %in% chromosomes_vcf]))
  }

  vcf_split<-split(vcf, by="CHROM")
  vcf_split[]<-lapply(vcf_split, function(x) as.numeric(unique(x$POS)))

  counts<-apply(regions, 1, function(x){
    chrom<-as.character(x[1])
    sum(x[2]:x[3] %in% snp_split[[counts]])
  })
  return(counts)
}
