# get contexts of mutations


#' trinucleotide context
#' @param vars  mutations with CHROM REF ALT POS columns
#' @param fasta  genome sequence (from seqinr)
#' @return a vector of contexts for each mutation
#' @export


contexts<-function(vars, fasta){
  SBS<-c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  apply(vars, 1, function(x){
    cat(x)
    chr<-x["CHROM"]
    pos<-as.numeric(x["POS"])
    ref<-x["REF"]
    alt<-x["ALT"]
    if(nchar(alt)!=1 | nchar(ref)!=1){
      return(NA)
    }
    up<-toupper(paste0(fasta[[chr]][(pos-1)], collapse=""))
    down<-toupper(paste0(fasta[[chr]][(pos+1)], collapse=""))
    context<-paste0(up,"[",ref,">",alt,"]", down)
    mutation<-paste0(up,ref,down,">", up,alt,down)
    if(!context %in% SBS){
      ref<-toupper(comp(ref))
      alt<-toupper(comp(alt))
      up<-toupper(comp(up))
      down<-toupper(comp(down))
      context<-paste0(down,"[",ref,">",alt,"]", up)
      mutation<-paste0(up,ref,down,">", up,alt,down)
    }
    return(mutation)
  })
}
