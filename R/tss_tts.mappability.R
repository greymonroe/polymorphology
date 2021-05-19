# get TSS and TTS mappbility

#' extract relative position of polymorphisms in relation to gene bodies
#' @param gff directory and file name of reference gff OR data.table object of gff data
#' @param mapp directory and file name of gen map object
#' @param chrs selection specific chromosomes, default is 1
#' @param num set number of genes, default is "all
#' @return data.table with pos and D columns
#' @import data.table
#' @import vcfR
#' @import readr
#' @export

tss_tts.mappability<-function(gff, mapp, chrs=1, num="all"){

mappability1<-read_lines("~/Dropbox/Research/Collaborations/Pablo/vitis/data/mapp/mapp100_2.txt", skip=1,n_max=2)
mappability1<-mappability1[1]
mappability1<-strsplit(mappability1, split=" ")
mappability1<-unlist(mappability1)
mappability1<-as.numeric(mappability1)

if (is.character(gff)){
  gff<-fread(gff)
  colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
}


genes<-gff[type=="gene" & chr %in% "chr00" & direction=="+"]

mapp<-c()
for(i in -3000:3000){
  cat(i, "\n")
  mapp<-c(mapp, mean(mappability1[genes$start-i], na.rm=T))

}

plot(mapp, type="l");abline(v=3000)
plot(mapp, ylim=c(0,1), type="l", col="green");abline(v=3000, col="gray", lty=3)

tss<-mapp

mapp<-c()
for(i in -3000:3000){
  cat(i, "\n")
  mapp<-c(mapp, mean(mappability1[genes$stop-i], na.rm=T))

}

plot(mapp, type="l");abline(v=3000)
plot(mapp, ylim=c(0,1), type="l", col="green");abline(v=3000, col="gray", lty=3)

tts<-mapp

tss_tts_mapp<-data.table(mapp=c(tss, tts), pos=-3000:3000, loc=rep(c("TSS", "TTS"), each=6001))


}
