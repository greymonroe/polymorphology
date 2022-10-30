# plot contexts of mutations


#' plot trinucleotide context
#' @param contexts  results from contexts()
#' @return a list with context_table and a ggplot
#' @export

contexts_plot<-function(contexts, full=T){
  SBS<-c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  mutations<-paste0(substr(SBS,1,1),substr(SBS,3,3),substr(SBS,7,7),">", substr(SBS,1,1),substr(SBS,5,5),substr(SBS,7,7))


  if(sum(grepl("NA",contexts))>0){
    warning("NA values found in contexts")
  }
  contexts<-contexts[!grepl("NA",contexts)]
if (full ==T){
context_table<-data.table(table(context=c(contexts, mutations)))
context_table$context_only<-substr(context_table$context, 1, 3)
context_table$mut<-paste(substr(context_table$context, 2,2),substr(context_table$context, 6,6), sep=">")
context_table$N<-context_table$N-1
plot<-ggplot(context_table, aes(x=context_only, y=(N), fill=mut))+
  geom_bar(stat="identity", width=0.5)+
  facet_grid(~mut, scales = "free")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Context")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(x = 0.3, units = "line"))+
  scale_fill_manual(values=c("cyan3","black","red4","gray","green3","pink3"), guide="none")
return(list(context_table, plot))
} else {
  context_table<-data.table(table(context=paste(substr(contexts, 2,2),substr(contexts, 6,6), sep=">")))
  plot<-ggplot(context_table, aes(x=context, y=(N), fill=context))+
    geom_bar(stat="identity", width=0.5)+
    theme_classic(base_size = 6)+
    scale_x_discrete(name="Mutation")+
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(x = 0.3, units = "line"))+
    scale_fill_manual(values=c("cyan3","black","red4","gray","green3","pink3"), guide="none")
  return(list(context_table, plot))
}
}
