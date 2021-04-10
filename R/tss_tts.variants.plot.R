
#' plot relative position of polymorphisms in relation to gene bodies
#' @param positions object creatd by tss_tts.variants() with columns pos and loc
#' @param window window size for counting variants, default is 10
#' @return a ggplot object plotting the tss and tts data
#' @import data.table
#' @import ggplot2
#' @export

tss_tts.variants.plot<-function(positions, window=10, color="green4"){
  
  positions$bins<-as.numeric(as.character(cut(positions$pos, breaks=seq(-3000, 3000, by=window), labels=seq(-3000, 3000, by=window)[-1])))
  positions_means<-positions[,.N, by=.(bins, loc)]
  plot<-ggplot(positions_means, aes(x=bins, y=N))+
    geom_line(group=1, col=color)+
    facet_grid(~loc)+
    theme_classic(base_size = 6)+
    scale_x_continuous(name="Relative genomic positions")
  return(plot)

}
