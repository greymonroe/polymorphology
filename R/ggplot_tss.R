# plot values around tss tts

#' ggplot_tss
#' @param tss_means_out  results from tss_means
#' @param variable  which variable are means calulated for
#' @return a vector of contexts for each mutation
#' @export

ggplot_tss<-function(tss_means_out, variable){
  ggplot(tss_means_out, aes(x=bins, y=mean))+
    geom_line()+
    facet_grid(~loc)+
    scale_y_continuous(name=variable)
}
