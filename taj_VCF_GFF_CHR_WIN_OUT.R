#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
vcf=args[1]
head=args[2]
gff=args[3]
chrom=args[4]
win=args[5]
out=args[6]


setwd(out)

library(data.table)

gff<-fread(gff)
colnames(gff)<-c("chr","source","type","start","stop","V1","direction","V2","info")
#gff$chr<-gsub("Chr","", gff$chr)
# looking only at genes
genes<-gff[type=="gene"]

i<-0
taj1000<-rbindlist(apply(genes[chr==chrom][1:10], 1, function(x){
  i<<-i+1
  cat(i,",")
  start<-as.numeric(x[4])
  stop<-as.numeric(x[5])
  if(x[7]=="+"){
    upstream_from<-start-win
    upstream_to<-start

    fiveprime_from<-start
    fiveprime_to<-start+win

    threeprime_from<-stop-win
    threeprime_to<-stop

    downstream_from<-stop
    downstream_to<-stop+win
  }

  if(x[7]=="-"){
    upstream_from<-stop
    upstream_to<-stop+win

    fiveprime_from<-stop-win
    fiveprime_to<-stop

    threeprime_from<-start
    threeprime_to<-start+win

    downstream_from<-start-win
    downstream_to<-start
  }


  system(paste0("tabix ",vcf," ",chrom,":",upstream_from,"-",upstream_to," > testwindow.vcf"))
  system(paste0("cat ",head," testwindow.vcf > test.vcf"))
  system("vcftools --vcf test.vcf --out test --TajimaD 10000", ignore.stdout = T, ignore.stderr  = T)
  upstream<-fread("test.Tajima.D")
  if(nrow(upstream)==0){cat("nothing")}
  upstream$region<-"upstream"

  system(paste0("tabix ",vcf," ",chrom,":",fiveprime_from,"-",fiveprime_to," > testwindow.vcf"))
  system(paste0("cat ",head," testwindow.vcf > test.vcf"))
  system("vcftools --vcf test.vcf --out test --TajimaD 10000", ignore.stdout = T, ignore.stderr  = T)
  five<-fread("test.Tajima.D")
  if(nrow(five)==0){cat("nothing")}
  five$region<-"five"

  system(paste0("tabix ",vcf," ",chrom,":",threeprime_from,"-",threeprime_to," > testwindow.vcf"))
  system(paste0("cat ",head," testwindow.vcf > test.vcf"))
  system("vcftools --vcf test.vcf --out test --TajimaD 10000", ignore.stdout = T, ignore.stderr  = T)
  three<-fread("test.Tajima.D")
  if(nrow(three)==0){cat("nothing")}
  three$region<-"three"

  system(paste0("tabix ",vcf," ",chrom,":",downstream_from,"-",downstream_to," > testwindow.vcf"))
  system(paste0("cat ",head," testwindow.vcf > test.vcf"))
  system("vcftools --vcf test.vcf --out test --TajimaD 10000", ignore.stdout = T, ignore.stderr  = T)
  downstream<-fread("test.Tajima.D")
  if(nrow(downstream)==0){cat("nothing")}
  downstream$region<-"downstream"

  out<-rbindlist(list(upstream, five, three, downstream))
  out$gene<-x[9]
  return(out)
}))
fwrite(taj1000, "taj1000.csv")


