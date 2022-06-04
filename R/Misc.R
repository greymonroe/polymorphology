

homopolymer_read<-function(seq, size){
  s=rle(seq)
  v=cumsum(rle(seq)$lengths)
  runs<-data.table('var'=s$values,'start'=v+1-s$lengths,'end'=v)
  runs$length<-runs$end-runs$start+1
  runs$var<-toupper(runs$var)
  homopolymers<-runs[length>=size]
}


make_homopolymer<-function(genome, size){
  seqs<-lapply(genome, function(x) unlist(x[1:length(x)]))
  homopolymers<-lapply(seqs, function(x) homopolymer_read(x, size))

}

neighbor<-function(vars,fasta){
  rbindlist(apply(vars, 1, function(x){
    chr<-x["CHROM"]
    pos<-as.numeric(x["POS"])
    ref_vcf<-x["REF"]
    alt<-x["ALT"]
    upstream<-toupper(paste0(fasta[[chr]][pos+1], collapse=""))
    downstream<- toupper(paste0(fasta[[chr]][pos-1], collapse=""))
    ref_genome<- toupper(fasta[[chr]][pos])
    neighbor_same<-(alt==upstream | alt==downstream)

    return(data.table(alt, downstream, upstream,ref_genome, ref_vcf, neighbor_same))
  }))
}
