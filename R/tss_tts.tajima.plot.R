
#' plot relative position of tajimas D in relation to gene bodies
#' @param tajima object creatd by tss_tts.tajima() with columns pos, D, and loc
#' @param bins window size for counting variants, default is 10
#' @return a ggplot object plotting the tss and tts data
#' @import data.table
#' @import ggplot2
#' @export
tss_tts.tajima.plot<-function(tajima, window=300, color="green4"){

  tsstaj_means<-tajima[,.(D=mean(D, na.rm=T)), by=.(bins=as.numeric(cut(pos, 300)), loc)]
  plot<-ggplot(tsstaj_means, aes(x=bins, y=D))+
    geom_line(aes(group=1), col=color, size=0.25)+
    facet_grid(~loc, scales = "free")+
    theme_classic(base_size = 6)+
    scale_x_continuous(name="Relative genomic positions")+
    geom_vline(xintercept=0, linetype="dashed", size=0.25)+
    scale_y_continuous(name="Polymorphisms")
  return(plot)

}