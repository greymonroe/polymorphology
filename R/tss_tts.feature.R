# get TSS and TTS feature distributions
# input:

#' extract relative position of polymorphisms in relation to gene bodies
#' @param gff directory and file name of reference gff OR data.table object of gff data
#' @param chrs selection specific chromosomes, default is "all"
#' @param num set number of genes, default is "all
#' @param feature name of the feature in the gff data that you want to look in relation to, default is "gene"
#' @param feature2 name of the feature in the gff data that you want test for obverlap default is "gene"
#' @param loc "tss" or "tts"
#' @return a vector of the probablity of overlapping in relation to feature
#' @import data.table
#' @import vcfR
#' @export

tss_tts.feature<-function(gff, chrs="all", num="all", feature="gene", feature2="gene", loc="tss"){


  if (is.character(gff)){
    gff<-fread(gff)
    colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
  }

  genes<-gff[type==feature]
  if (chrs!="all"){genes<-genes[chr %in% chrs]}
  if (num!="all"){genes<-genes[1:num]}

  features<-gff[type==feature2]
  if (chrs!="all"){features<-features[chr %in% chrs]}
  if (num!="all"){features<-features[1:num]}

  if(type=="snp"){var_split<-split(vcf[nchar(REF)==nchar(ALT)], by="CHROM")}
  if(type=="indel"){var_split<-split(vcf[nchar(REF)!=nchar(ALT)], by="CHROM")}
  features_split<-split(features[], by="chr")
  features_split[]<-lapply(features_split, function(x) unique(unlist(apply(x, 1, function(y) y[4]:y[5]))))


  if(loc=="tss"){
    positions<-unlist(apply(genes, 1, function(x){
      if(x[7]=="-"){
        out<- c(-3000:3000)[rev((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% features_split[[x[1]]])]
      } else {
        out<-c(-3000:3000)[((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% features_split[[x[1]]])]
      }
      return(out)
    }
    ))
  }

  if(loc=="tts"){
    positions<-unlist(apply(genes, 1, function(x){
      if(x[7]=="-"){
        out<- c(-3000:3000)[rev((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% features_split[[x[1]]])]
      } else {
        out<-c(-3000:3000)[((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% features_split[[x[1]]])]
      }
      return(out)
    }
    ))
  }
  return(positions)
}
