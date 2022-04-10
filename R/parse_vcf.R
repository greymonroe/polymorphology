# Parse VCF

#' parse through vcf file so it can be used in R
#' @param vcf_file file name vcf data from gatk haplotypecaller
#' @param snpeff boolean: the vcf is annotated with snpeff, default = F
#' @return a data.table object with parsed vcf, new columns for various quality meaures, statistics, annotations, etc
#' @import data.table
#' @import vcfR
#' @export
parse_vcf<-function(vcf_file, snpeff=F){

  v<-read.vcfR(vcf_file)
  fix<-data.table(v@fix)
  gt<-data.table(v@gt)
  fix$file<-x
  fix$calls<-gt[,2]
  fix$BaseQRankSum<-as.numeric(gsub(".+BaseQRankSum=(.+);DP.+","\\1",fix$INFO ))
  fix$FS<-as.numeric(gsub(".+FS=(.+);MLEAC.+","\\1",fix$INFO ))
  fix$MQRankSum<-as.numeric(gsub(".+MQRankSum=(.+);QD.+","\\1",fix$INFO ))
  fix$ReadPosRankSum<-as.numeric(gsub(".+ReadPosRankSum=(.+);SOR.+","\\1",fix$INFO ))
  fix$MQ<-as.numeric(gsub(".+MQ=(.+);MQRankSum.+","\\1",fix$INFO ))
  fix$QD<-as.numeric(gsub(".+QD=(.+);ReadPosRankSum.+","\\1",fix$INFO ))
  fix$SOR<-as.numeric(gsub(".+SOR=(.+)","\\1",fix$INFO ))
  fix$DP<-as.numeric(sapply(fix$calls, function(x) unlist(strsplit(x, split=":"))[3]))
  fix$RD<-as.numeric(sapply(fix$calls, function(x) unlist(strsplit(unlist(strsplit(x, split=":"))[2], split=","))[1]))
  fix$AD<-as.numeric(sapply(fix$calls, function(x) unlist(strsplit(unlist(strsplit(x, split=":"))[2], split=","))[2]))
  fix$GT<-sapply(fix$calls, function(x) unlist(strsplit(x, split=":"))[1])
  fix$MULTIALLELIC<-grepl(",", fix$ALT)
  fix$TYPE<-ifelse(nchar(fix$REF)==1 & nchar(fix$ALT)==1, "SNP","InDel")
  fix$TYPE[MULTIALLELIC]<-NA
  fix$SIZE<-nchar(fix$ALT)-nchar(fix$REF)
  fix$SIZE[MULTIALLELIC]<-NA

  if(snpeff){
    fix$intron<-grepl("intron", fix$INFO)
    fix$nonsynonymous<-grepl("missense", fix$INFO)
    fix$synonymous<-grepl("synonymous", fix$INFO)
    fix$stop_gained<-grepl("stop_gained", fix$INFO)
    fix$inframe<-grepl("inframe", fix$INFO)
    fix$frameshift<-grepl("frameshift", fix$INFO)
  }
  return(fix)
}
