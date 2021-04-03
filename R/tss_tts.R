# get TSS and TTS variant distributions
# input:

#' extract relative position of polymorphisms in relation to gene bodies
#' @param gff directory and file name of reference gff OR data.table object of gff data
#' @param vcf directory and file name vcf data OR data.table object of vcf data (from package vcfR, and only the @fix data)
#' @param loc "tss" or "tts", default is "tss"
#' @param type "snp", "indel", or "both"
#' @param chrs selection specific chromosomes, default is "all"
#' @param num set number of genes, default is "all
#' @param feature name of the feature in the gff data that you want to look at, default is "gene"
#' @return a vector of the relative position of polymorphsims in relation to feature
#' @export

tss_tts<-function(gff=NA, vcf=NA, loc="tss", type="snp", chrs="all", num="all", feature="gene"){

  if (is.character(gff)){
    gff<-fread(gff)
    colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
  }


  if (is.character(vcf)){
    vcf<-read.vcfR(vcf)
    vcf<-vcf@fix
    }

  genes<-gff[type==feature]
  if (chrs!="all"){genes<-genes[chr %in% chrs]}
  if (num!="all"){genes<-genes[1:num]}

  if(type=="snp"){var_split<-split(vcf[nchar(REF)!=nchar(ALT)], by="CHROM")}
  if(type=="indel"){var_split<-split(vcf[nchar(REF)==nchar(ALT)], by="CHROM")}
  if(type=="both"){var_split<-split(vcf[], by="CHROM")}
  var_split[]<-lapply(var_split, function(x) unique(x$POS))


  if(loc=="tss"){
    positions<-unlist(apply(genes[1:maxgene], 1, function(x){
    if(x[7]=="-"){
      out<- c(-3000:3000)[rev((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% var_split[[x[1]]])]
    } else {
      out<-c(-3000:3000)[((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% var_split[[x[1]]])]
    }
    return(out)
  }
  ))
  }

  if(loc=="tts"){
    positions<-unlist(apply(genes, 1, function(x){
      if(x[7]=="-"){
        out<- c(-3000:3000)[rev((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% var_split[[x[1]]])]
      } else {
        out<-c(-3000:3000)[((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% var_split[[x[1]]])]
      }
      return(out)
    }
    ))
  }

}
