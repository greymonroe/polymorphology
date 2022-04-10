# Extract values for epigenomic features in a given region

#' Extract values for epigenomic features in a given region. Used in Monroe et al 2022 Not a very generalizable function as it depends on formatting of datasets.
#' @param regions data.table object with columns columns chr, start, stop
#' @param peaks_file file with ATACseq peaks
#' @param meth_file file with methylation data
#' @param genome_file reference genome fasta file
#' @param histone_dir directory of files with histone modifications
#' @param tmp_dir directory for temporary files
#' @param multiBigwigSummary system location for multiBigwigSummary function
#' @return a data.table annotated with genomic and epigenomic features
#' @export
#'

extract_epigenome<-function(regions, peaks_file, meth_file, genome_file, histone_dir, tmp_dir, multiBigwigSummary){

  if(!all(colnames(regions)[1:3] == c("chr","start","stop" ))){
    stop('First 3 column names of regions must be "chr","start","stop"')
  }
  if(typeof(regions$chr)!="character"){
    stop('regions$chr should be of type character')
  }
  if(!file.exists(tmp_dir)){
    stop('tmp_dir provided does not exist.')
  }
  if(!file.exists(peaks_file)){
    stop('peaks_file dprovided oes not exist.')
  }
  if(!file.exists(meth_file)){
    stop('meth_file provided does not exist.')
  }
  if(!file.exists(genome_file)){
    stop('genome_file provided does not exist.')
  }
  if(!file.exists(histone_dir)){
    stop('histone_dir provided does not exist.')
  }
  if(!file.exists(multiBigwigSummary)){
    stop('multiBigwigSummary provided does not exist.')
  }

  cat("It is always recommended to check your data files to see the formatting of the chromosome names. You may need to modify the code in this function (lines where you see ##CHROM NAME FORMATTING##) to make all chromosome names in the same format for it to work properly and avoid esoteric errors...")

  cat("Reading reference genome")

  genome<-read.fasta(genome_file)
  chr_lengths<-lapply(genome, length)
  ingenome<-apply(regions, 1, function(x){
    chr=as.character(x[1])
    stop=as.numeric(3)
    stop<chr_lengths[chr]
  })
  if(!all(ingenome)){
    stop("regions outside bounds of reference genome provided")
  }

  regions$length<-regions$stop-regions$start+1

  # add histone marks  ------------------------------------------------------
  cat("Extracting histone marks...")

  tmp_bed<-data.table(paste0("chr", regions$chr), regions$start, regions$stop) ##CHROM NAME FORMATTING##
  tmp_bed$V2<-as.integer(tmp_bed$V2)
  tmp_bed$V3<-as.integer(tmp_bed$V3)

  cat(paste0("writing ", tmp_dir, "/tmp_bed.bed"))
  write.table(tmp_bed, paste0(tmp_dir, "/tmp_bed.bed"), row.names = F, col.names = F, quote = F, sep="\t")

  cat("using files with ^H3K.+.bw format...")
  histone_files<-list.files(histone_dir, pattern = "^H3K.+.bw")

  # Monroe et al 2022 uses data from http://systemsbiology.cau.edu.cn/chromstates
  for (f in histone_files){
    cat("\n", f,"...")
  system(paste0(multiBigwigSummary, " BED-file -b ", histone_dir,f, " --outRawCounts ", tmp_dir,f,".tab -o ",tmp_dir ,f,".sum --BED ",tmp_dir, "tmp_bed.bed"))
  }

  for(f in histone_files){
    bwdata<-fread(paste0(tmp_dir,f,".tab"))
    colnames(bwdata)<-c("chr", "start","stop",f)
    bwdata$chr<-as.character(gsub("chr","",bwdata$chr))
    bwdata<-unique(bwdata)
    regions<-merge(regions, bwdata, all.x=T,all.y=F, by=c("chr","start","stop"))
  }

  # add ATACseq data ----------------------------------------------------

  cat("Extracting ATACseq peaks...")
  peaks<-fread(peaks_file)
  peaks$CHR<-gsub("Chr","",peaks$CHR) ##CHROM NAME FORMATTING##
  peaks_split<-split(peaks, by="CHR")
  peaks_split[]<-lapply(peaks_split, function(x) unlist(apply(x,1, function(y) y[3]:y[4])))

  atac<-apply(regions, 1, function(x) {
    chrom<-x[1]
    f_start<-as.numeric(x[2])
    f_stop<-as.numeric(x[3])
    sum(f_start:f_stop %in% peaks_split[[chrom]])
  })

  regions$atac<-atac
  regions$atac_pct<-regions$atac/regions$length

  # add methylation data ----------------------------------------------------
  cat("Extracting C_methylatiton data...")

  meth<-fread(meth_file)
  colnames(meth)<-c("chr","start","stop","direction","meth","unmeth","type","seq")
  meth$chr<-gsub("Chr","",meth$chr) ##CHROM NAME FORMATTING##
  meth_split<-split(meth, by=c("chr"))
  meth_split<-lapply(meth_split, function(x) split(x, by="type"))

  CG<-apply(regions, 1, function(x) {
    chrom<-as.character(x[1])
    f_start<-as.numeric(x[2])
    f_stop<-as.numeric(x[3])
    sum<-nrow(meth_split[[chrom]][["CG"]][start>=f_start & start<=f_stop])
  })

  CHG<-apply(regions, 1, function(x) {
    chrom<-as.character(x[1])
    f_start<-as.numeric(x[2])
    f_stop<-as.numeric(x[3])
    sum<-nrow(meth_split[[chrom]][["CHG"]][start>=f_start & start<=f_stop])
  })

  CHH<-apply(regions, 1, function(x) {
    chrom<-as.character(x[1])
    f_start<-as.numeric(x[2])
    f_stop<-as.numeric(x[3])
    sum<-nrow(meth_split[[chrom]][["CHH"]][start>=f_start & start<=f_stop])
  })

  regions$CG<-CG
  regions$CHG<-CHG
  regions$CHH<-CHH
  regions$CG_pct<-regions$CG/regions$length
  regions$CHG_pct<-regions$CHG/regions$length
  regions$CHH_pct<-regions$CHH/regions$length

  # add GC content data ----------------------------------------------------
  cat("Extracting CG content...")

  genome<-read.fasta(genome_file)

  GC_content<-apply(regions, 1, function(x) {
    f_start<-as.numeric(x[2])
    f_stop<-as.numeric(x[3])
    seq<-(genome[[as.character(x[1])]][f_start:f_stop])
    return(sum(grepl("c|g", seq)))
  })

  regions$GC_content<-GC_content
  regions$GC_content_pct<-regions$GC_content/regions$length
return(regions)
}

