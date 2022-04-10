# get TSS and TTS variant distributions
# input:

#' extract relative position of polymorphisms in relation to gene bodies
#' @param gff directory and file name of reference gff OR data.table object of gff data
#' @param vcf directory and file name vcf data OR data.table object of variant data (data.table with CHROM, POS, REF, and ALT columns)
#' @param type "snp", "indel", or "both", default is "both"
#' @param chrs selection specific chromosomes, default is "all"
#' @param num set number of genes, default is "all
#' @param feature name of the feature in the gff data that you want to look at, default is "gene"
#' @return a vector of the relative position of polymorphsims in relation to feature
#' @import data.table
#' @import vcfR
#' @export

tss_tts.variants<-function(gff, vcf, type="both", chrs="all", num="all", feature="gene"){

  if (is.character(gff)){
    gff<-fread(gff)
    colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
  }


  if (is.character(vcf)){
    vcf<-read.vcfR(vcf)
    vcf<-data.table(vcf@fix)
  }

  genes<-gff[type==feature]
  if (chrs!="all"){genes<-genes[chr %in% chrs]}
  if (num!="all"){genes<-genes[1:num]}

  if(type=="snp"){var_split<-split(vcf[nchar(REF)==nchar(ALT)], by="CHROM")}
  if(type=="indel"){var_split<-split(vcf[nchar(REF)!=nchar(ALT)], by="CHROM")}
  if(type=="both"){var_split<-split(vcf[], by="CHROM")}
  var_split[]<-lapply(var_split, function(x) unique(x$POS))


  # extract positions of variants relative to transcription start sites (tss)
  tss<-rbindlist(apply(genes, 1, function(x){
    if(x[7]=="-"){
      out<- c(-3000:3000)[rev((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% var_split[[x[1]]])]
    } else {
      out<-c(-3000:3000)[((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% var_split[[x[1]]])]
    }
    if(length(out)==0){
      return(NULL)
    }
    dt<-data.table(pos=out, loc="TSS", gene=x[9])
    return(dt)
  }
  ))

  # extract positions of variants relative to transcription termination sites (tts)
  tts<-rbindlist(apply(genes, 1, function(x){
    if(x[7]=="-"){
      out<- c(-3000:3000)[rev((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% var_split[[x[1]]])]
    } else {
      out<-c(-3000:3000)[((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% var_split[[x[1]]])]
    }
    if(length(out)==0){
      return(NULL)
    }
    dt<-data.table(pos=out, loc="TTS", gene=x[9])
    return(dt)
  }
  ))
  # merge into single data.table of tss and tts
  results<-rbind(tss, tts)

  return(results[!is.na(pos)])

}
