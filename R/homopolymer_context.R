# get homopolymer contexts of mutations


#' homopolymer_context
#' @param vars  mutations with CHROM REF ALT POS columns
#' @param genome  genome sequence (from seqinr)
#' @param homopolymers  homopolymer object from make_homopolymer()
#' @param size  size of homopolymer to consider
#' @return a data table
#' @export


homopolymer_context<-function(vars, genome, homopolymers, size){

  homopolymers[]<-lapply(homopolymers, function(x) x[length>=size])

  rbindlist(apply(vars, 1, function(x){
    chr<-x["CHROM"]
    pos<-as.numeric(x["POS"])
    ref<-x["REF"]
    alt<-x["ALT"]
    homo<-homopolymers[[chr]][(start-1==pos | end+1==pos) & var==alt]
    nexttohomopolymer<-nrow(homo)
    context<-toupper(paste0(genome[[chr]][(pos-10):(pos+10)], collapse=""))
    return(data.table(alt, ref, context, nexttohomopolymer))
  }))
}

