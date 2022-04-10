# Annotate variants with features they are found in

#' annotate variants with features from gff
#' @param gff data.table object of gff data (data.table with "chr","source","type","start","stop","V1","direction","V2","info" columns)
#' @param vcf data.table object of variant data (data.table with CHROM, POS, REF, and ALT columns)
#' @param type type of feature (e.g, "CDS","gene", etc)
#' @param ann which column to annotate with (values will be concatenated to track the names of the features a varint is found in), default is "info" column of gff
#' @param window number of base pairs to expand (e.g. to include regions +- 3000 base pairs from feature you set 3000)
#' @return a vector of pasted (" AND " delimited) annotations (based on column provided by ) for each feature the variant is found in
#' @export
#'

variant_features<-function(vcf, gff, type, ann='info', window=0){

if(!all(colnames(gff)[1:9] == c("chr","source","type","start","stop","V1","direction","V2","info" ))){
  stop('First 9 column names of gff must be "chr","source","type","start","stop","V1","direction","V2","info" ')
}
if(!all(colnames(vcf)[1:2] == c("CHROM","POS"))){
  stop('CHROM and POS must be first 2 colnames in vcf')
}

annotations<-apply(vcf, 1, function(x){
  CHROM<-x[1]
  if(!CHROM %in% gff$chr){
    stop(paste(CHROM, "from vcf not found in gff"))
  }
  POS<-as.numeric(x[2])
  features<-gff[chr==CHROM & POS>=start-window & POS<=stop+window]
  if(nrow(features)==0){
    return(NA)
  } else {
  return(paste(features[,ann, with=F], collapse=" AND "))
  }
})

return(annotations)

}
