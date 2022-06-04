# get calculate means in windows relative to tss tss

#' tss_means
#' @param tss_out  results from tss_tts.calc
#' @param variable  which variable are means calulated for
#' @param window  window size to summarize in
#' @param na_action what to do with NA values in mean()
#' @return a vector of contexts for each mutation
#' @export

tss_means<-function(tss_out=out, variable="dist", window=100, na_action=F){
  tss_out$bins<-as.numeric(as.character(cut(tss_out$pos, breaks=seq(-3000, 3000, by=window), labels=seq(-3000, 3000, by=window)[-1])))
  tss_out_means<-tss_out[!is.na(bins),.(mean=mean(get(variable), na.rm=na_action)), by=.(bins, loc)]
  return(tss_out_means)
}
