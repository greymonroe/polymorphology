# get TSS and TTS minor allele frequency distributions
# input:

#' calculate minor allele frequencies of polymorphisms in relation to gene bodies
#' @param gff directory and file name of reference gff OR data.table object of gff data
#' @param vcf directory and file name vcf data OR data.table object of variant data (data.table with CHROM, POS, MAF)
#' @param type "snp", "indel", or "both", default is "both"
#' @param chrs selection specific chromosomes, default is "all"
#' @param num set number of genes, default is "all
#' @param feature name of the feature in the gff data that you want to look at, default is "gene"
#' @param variable name of variable in the vcf you want the calculated value of in relation to TSS TTS
#' @return a vector of the relative position of polymorphsims in relation to feature
#' @import data.table
#' @import vcfR
#' @export

tss_tts.calc<-function(gff, vcf, type="both", chrs="all", num="all", feature="gene", variable){

  if (is.character(gff)){
    gff<-fread(gff)
    colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
  }


  genes<-gff[type==feature]
  if (chrs!="all"){genes<-genes[chr %in% chrs]}
  if (num!="all"){genes<-genes[1:num]}

  freqs_split<-split(vcf, by="CHROM")

  maf_tss<-rbindlist(apply(genes, 1, function(x){

    if(x[7]=="-"){
      which<-(((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% freqs_split[[x[1]]]$POS))
      nums<-c(as.numeric(x[5])-3000):(as.numeric(x[5])+3000)
      pos<- (-3000:3000)[rev(which)]
      vals<-rev(freqs_split[[x[1]]][POS %in% nums[which]][,variable, with=F])

    } else {
      which<-((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% freqs_split[[x[1]]]$POS)
      nums<-c(as.numeric(x[4])-3000):(as.numeric(x[4])+3000)
      pos<- (-3000:3000)[which]
      vals<-freqs_split[[x[1]]][POS %in% nums[which]][,variable, with=F]
    }
    return(data.table(pos, vals, loc="TSS", info=x["info"]))
  }
  ))


  maf_tts<-rbindlist(apply(genes, 1, function(x){
    if(x[7]=="-"){
      which<-((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% freqs_split[[x[1]]]$POS)
      nums<-c(as.numeric(x[4])-3000):(as.numeric(x[4])+3000)
      pos<- (-3000:3000)[rev(which)]
      vals<-rev(freqs_split[[x[1]]][POS %in% nums[which]][,variable, with=F])
    } else {
      which<-((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% freqs_split[[x[1]]]$POS)
      nums<-c(as.numeric(x[5])-3000):(as.numeric(x[5])+3000)
      pos<- (-3000:3000)[which]
      vals<-freqs_split[[x[1]]][POS %in% nums[which]][,variable, with=F]  }
    return(data.table(pos, vals, loc="TTS", info=x["info"]))
  }
  ))

  results<-rbind(maf_tts, maf_tss)

  return(results[!is.na(pos)])

}
