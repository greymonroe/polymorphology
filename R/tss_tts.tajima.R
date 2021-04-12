# get TSS and TTS Tajima's D distributions

#' extract relative position of polymorphisms in relation to gene bodies
#' @param gff directory and file name of reference gff OR data.table object of gff data
#' @param tajima directory and file name tajima data OR data.table object of tajima data (data.table with CHROM, BIN_START, N_SNPS, and TajimaD columns)
#' @param chrs selection specific chromosomes, default is "all"
#' @param num set number of genes, default is "all
#' @param feature name of the feature in the gff data that you want to look in relation to, default is "gene"
#' @param loc "tss" or "tts"
#' @return data.table with pos and D columns
#' @import data.table
#' @import vcfR
#' @export

tss_tts.tajima<-function(gff, tajima, type="both", chrs="all", num="all", loc="tss", feature="gene"){

  if (is.character(gff)){
    gff<-fread(gff)
    colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
  }

  genes<-gff[type==feature]
  if (chrs!="all"){genes<-genes[chr %in% chrs]}
  if (num!="all"){genes<-genes[1:num]}

  tajima_split<-split(tajima, by="CHROM")

  #get TSS snp coordinates
  if(loc=="tss"){
    D<-apply(genes[], 1, function(x){
      if(x[7]=="-"){
        stop<-as.numeric(x[4])
        start<-as.numeric(x[5])
        window<-(start-3000):(start+3000)

        bins<-window[window %in% tajima_split[[x[1]]]$BIN_START]
        #nums<-((as.numeric(x[5])-3000):(as.numeric(x[5])+3000))[which(which)]
        pos<- (start-bins)
        taj<-tajima_split[[x[1]]][(BIN_START) %in% bins]$TajimaD

      } else {
        which<-((as.numeric(x[4])-3000):(as.numeric(x[4])+3000) %in% (tajima_split[[x[1]]]$BIN_START))
        nums<-c(as.numeric(x[4])-3000):(as.numeric(x[4])+3000)
        pos<- (-3000:3000)[which]
        taj<-tajima_split[[x[1]]][(BIN_START) %in% nums[which(which)]]$TajimaD

      }
      return(data.table(pos, D=taj))
    }
    )
  }

  if(loc=="tts"){
    D<-apply(genes[], 1, function(x){
      if(x[7]=="-"){
        stop<-as.numeric(x[4])
        start<-as.numeric(x[5])
        window<-(stop-3000):(stop+3000)

        bins<-window[window %in% tajima_split[[x[1]]]$BIN_START]
        #nums<-((as.numeric(x[5])-3000):(as.numeric(x[5])+3000))[which(which)]
        pos<- (stop-bins)
        taj<-tajima_split[[x[1]]][(BIN_START) %in% bins]$TajimaD
      } else {
        which<-((as.numeric(x[5])-3000):(as.numeric(x[5])+3000) %in% (tajima_split[[x[1]]]$BIN_START))
        nums<-c(as.numeric(x[5])-3000):(as.numeric(x[5])+3000)
        pos<- (-3000:3000)[which]
        taj<-tajima_split[[x[1]]][(BIN_START) %in% nums[which(which)]]$TajimaD
      }
      return(data.table(pos=pos, D=taj))
    }
    )
  }
  return(rbindlist(D))
}
