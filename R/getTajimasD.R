# get TajimasD for each row in a GFF
# input:

#' extract Tajimas D for features
#' @param gff directory and file name of reference gff
#' @param vcf directory and file name vcf data
#' @return writes out a file with tajimnas D values for each feature of a GFF file
#' @import data.table
#' @import vcfR
#' @export

getTajimasD<-function(vcf, gff){

  #we then create a file to paste the output
  system(paste0("touch tajimas_",vcf,".txt"))

  gff<-read.table(gff)

#w we then loop through the gff
# for each row we create a temporary VCF file containing the coordinates of that gff item
# we then run vcftools --TajiamsD with a huge window size (so there is only 1 window per feature)
# we concatinate ths to the file where we paste the output

for( i in 1:nrow(gff)){
  chr<-gff$V1[i]
  start<-gff$V2[i]
  stop<-gff$V3[i]

  system(paste0("tabix -h ",vcf," ", chr,":",start,"-",stop," > ",vcf,"_tmp_region.vcf"))
  system(paste0("vcftools --vcf ",vcf,"_tmp_region.vcf --out ",vcf," --TajimaD 100000000"))
  system(paste0("tail -1 ",vcf,".Tajima.D >> tajimas_",vcf,".txt"))
}

}
