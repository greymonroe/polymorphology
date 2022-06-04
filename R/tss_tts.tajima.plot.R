
#' plot relative position of tajimas D in relation to gene bodies
#' @param tajima object creatd by tss_tts.tajima() with columns pos, D, and loc
#' @param bins window size for counting variants, default is 10
#' @return a ggplot object plotting the tss and tts data
#' @import data.table
#' @import ggplot2
#' @export

tss_tts.tajima.plot<-function(tajima, window=10, color="green4", y="D"){

  if (!y %in% c("D","SNPS")){stop(" y= must be either 'D' or 'SNPS'")}

  tajima$bins <- as.numeric(as.character(cut(tajima$pos,
                                             breaks = seq(-3000, 3000, by = window), labels = seq(-3000,
                                                                                              3000, by = window)[-1])))
  tsstaj_means <- tajima[, .(D = mean(D, na.rm = T), SNPS=mean(N_SNPS)), by = .(bins = bins, loc)]

  if(y=="D"){
   plot <- ggplot(tsstaj_means[!is.na(pos)], aes(x = bins, y = D)) + geom_line(aes(group = 1),
                                                                 col = color, size = 0.25) + facet_grid(~loc) +
    theme_classic(base_size = 6) + scale_x_continuous(name = "Relative genomic positions") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.25) +
    scale_y_continuous(name = "Tajima's D")

  }
  if(y=="SNPS"){
    plot <- ggplot(tsstaj_means[!is.na(pos)], aes(x = bins, y = SNPS)) + geom_line(aes(group = 1),
                                                                   col = color, size = 0.25) + facet_grid(~loc) +
      theme_classic(base_size = 6) + scale_x_continuous(name = "Relative genomic positions") +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.25) +
      scale_y_continuous(name = "Polymorphisms")

  }

  return(plot)

}
