# plot contexts of mutations


#' plot trinucleotide context
#' @param contexts  results from contexts()
#' @return a list with context_table and a ggplot
#' @export

contexts_plot<-function(contexts){
context_table<-data.table(table(context=contexts))
context_table$context_only<-substr(context_table$context, 1, 3)
context_table$mut<-paste(substr(context_table$context, 2,2),substr(context_table$context, 6,6), sep=">")
plot<-ggplot(context_table, aes(x=context_only, y=(N), fill=mut))+
  geom_bar(stat="identity", width=0.5)+
  facet_grid(~mut, scales = "free")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Context")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(x = 0.3, units = "line"))+
  scale_fill_manual(values=c("cyan3","black","red4","gray","green3","pink3"), guide="none")
return(list(context_table, plot))
}
