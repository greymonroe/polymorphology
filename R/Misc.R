#' @export


make_feature_windows<-function(encode_data, deciles=10, gene=F){

  library(seqinr)
  library(openxlsx)
  library(data.table)

  windows<-rbindlist(apply(encode_data, 1, function(x) {

    chr=x["chr"]
    body_starts=round(seq(as.numeric(x["start"]), as.numeric(x["stop"]), length.out=deciles+1)[-(deciles+1)])
    body_stops<-round(seq(as.numeric(x["start"]), as.numeric(x["stop"]), length.out=deciles+1)[-1])
    upstream_starts<-seq(as.numeric(x["start"])-2000, as.numeric(x["start"]), length.out=deciles+1)[-(deciles+1)]
    upstream_stops<-seq(as.numeric(x["start"])-2000, as.numeric(x["start"]), length.out=deciles+1)[-1]
    downstream_starts<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+2000, length.out=deciles+1)[-(deciles+1)]
    downstream_stops<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+2000, length.out=deciles+1)[-1]

    out<-data.table(chr=x["chr"],
                    start=c(upstream_starts, body_starts, downstream_starts),
                    stop=c(upstream_stops, body_stops, downstream_stops),
                    region=c(rep("upstream", length(upstream_starts)),rep("gene body", length(body_starts)),rep("downstream", length(downstream_starts))),
                    ID=as.numeric(x["ID"]))

    out$pos<-1:nrow(out)
    out$length<-out$stop-out$start

    if(gene == T){
      direction=x["direction"]
      if(direction=="-"){
        out$pos<-rev(out$pos)
        out$region<-rev(out$region)
      }
    }

    return(out)

  }))

  setkey(windows, chr, start, stop)
  return(windows)
}

read_encode<-function(file){
  encode_data<-fread(file)
  encode_data$V1<-gsub("chr","Chr", encode_data$V1)
  encode_data$chr<-gsub("chr","Chr", encode_data$V1)
  encode_data$ID<-1:nrow(encode_data)
  encode_data$start<-encode_data$V2
  encode_data$stop<-encode_data$V3
  encode_data$length<-encode_data$stop-encode_data$start+1

  setkey(encode_data, chr, start, stop)

  return(encode_data)
}

dimerfrequency<-function(regions, genome){

  if(!"chr" %in% colnames(regions)){
    stop("'chr' column not found in region object")
  }
  if(!"start" %in% colnames(regions)){
    stop("'start' column not found in region object")
  }
  if(!"stop" %in% colnames(regions)){
    stop("'stop' column not found in region object")
  }

  out<-rbindlist(apply(regions, 1, function(x){

    chrom<-x["chr"]

    if(!chrom %in% names(genome)){
      stop(paste(chrom,"from region object not found in names of genome object"))
    }

    start<-as.numeric(x["start"])
    if(start<0){start<-1}
    stop<-as.numeric(x["stop"])
    if(stop>length(genome[[chrom]])){stop<-length(genome[[chrom]])}
    seq<-toupper(paste0(unlist(genome[[chrom]][start:stop]), collapse=""))

    motifs<-expand.grid(c("A","T","C","G"), c("A","T","C","G"))
    motifs<-paste0(motifs$Var1, motifs$Var2)
    motif<-"YY"
    motif<-motifs[1]
    allhits<-rbindlist(lapply(motifs, function(motif){
      hits<-str_locate_all(seq, motif)
      hits<-data.table(hits[[1]])
      return(data.table(motif, count=nrow(hits)))
    }))
    counts<-prop.table(allhits$count)
    motif_props<-data.table(t(matrix(counts)))
    colnames(motif_props)<-motifs
    return(motif_props)
  }))

  return(out)
}

make_gene_windows<-function(data, window=150){
  deciles<-3000/window
  data$ID<-1:nrow(data)
  windows<-rbindlist(apply(data, 1, function(x) {

    chr=x["CHROM"]
    body_starts=seq(as.numeric(x["start"]), as.numeric(x["stop"]), by=3000/deciles);body_starts<-body_starts[-length(body_starts)]
    body_stops<-seq(as.numeric(x["start"]), as.numeric(x["stop"]),  by=3000/deciles)[-1]
    upstream_starts<-seq(as.numeric(x["start"])-3000, as.numeric(x["start"]), length.out=deciles+1)[-(deciles+1)]
    upstream_stops<-seq(as.numeric(x["start"])-3000, as.numeric(x["start"]), length.out=deciles+1)[-1]
    downstream_starts<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+3000, length.out=deciles+1)[-(deciles+1)]
    downstream_stops<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+3000, length.out=deciles+1)[-1]

    out<-data.table(chr=x["chr"],
                    start=c(upstream_starts, body_starts, downstream_starts),
                    stop=c(upstream_stops, body_stops, downstream_stops),
                    region=c(rep("upstream", length(upstream_starts)),rep("gene body", length(body_starts)),rep("downstream", length(downstream_starts))),
                    ID=x["ID"],
                    gene=x["locus"])
    out$pos<-1:nrow(out)
    out$length<-out$stop-out$start
    return(out)

  }))
  windows$window_ID<-1:nrow(windows)
  setkey(windows, chr, start, stop)
  gene_windows<-windows
  gene_windows$H3K4me1<-encode_overlap(H3K4me1, gene_windows)
  gene_windows$H3K9me1<-encode_overlap(H3K9me1, gene_windows)
  gene_windows$H3K4me3<-encode_overlap(H3K4me3, gene_windows)
  gene_windows$H3K36me3<-encode_overlap(H3K36me3, gene_windows)
  gene_windows$H3K9me2<-encode_overlap(H3K9me2, gene_windows)
  gene_windows$H3K27me3<-encode_overlap(H3K27me3, gene_windows)
  gene_windows$H3K27ac<-encode_overlap(H3K27ac, gene_windows)
  gene_windows$PII<-encode_overlap(PII, gene_windows)
  gene_windows$H3K4ac<-encode_overlap(H3K4ac, gene_windows)
  gene_windows$H3K12ac<-encode_overlap(H3K12ac, gene_windows)
  gene_windows$H3K9ac<-encode_overlap(H3K9ac, gene_windows)

  return(gene_windows)
}

make_genome_windows<-function(genome, window, overlap=T){
  lengths<-lapply(genome, length)
  names(lengths)<-names(genome)
  chrs<-rbindlist(lapply(names(genome)[1:12], function(c){
    starts<-seq(1, lengths[[c]], by=window);starts<-starts[-length(starts)]
    stops<-c(starts[-1],lengths[[c]])
    return(data.table(chr=c, start=starts, stop=stops))
  }))
  gene_windows<-chrs
  gene_windows$window_ID<-1:nrow(gene_windows)
  setkey(gene_windows, chr, start, stop)
  if(overlap==T){
    gene_windows$H3K4me1<-encode_overlap(H3K4me1, gene_windows)
    gene_windows$H3K9me1<-encode_overlap(H3K9me1, gene_windows)
    gene_windows$H3K4me3<-encode_overlap(H3K4me3, gene_windows)
    gene_windows$H3K36me3<-encode_overlap(H3K36me3, gene_windows)
    gene_windows$H3K9me2<-encode_overlap(H3K9me2, gene_windows)
    gene_windows$H3K27me3<-encode_overlap(H3K27me3, gene_windows)
    gene_windows$H3K27ac<-encode_overlap(H3K27ac, gene_windows)
    gene_windows$PII<-encode_overlap(PII, gene_windows)
    gene_windows$H3K4ac<-encode_overlap(H3K4ac, gene_windows)
    gene_windows$H3K12ac<-encode_overlap(H3K12ac, gene_windows)
    gene_windows$H3K9ac<-encode_overlap(H3K9ac, gene_windows)
    gene_windows$snp<-add_vars_to_gene_windows(gene_windows, snps)
    gene_windows$snp_noCT<-add_vars_to_gene_windows(gene_windows, snps[!Single.base.substitution %in% c("C>T","G>A")])
    gene_windows$snp_homoz<-add_vars_to_gene_windows(gene_windows, snps[Genotype=="Homozygous"])
    gene_windows$indel<-add_vars_to_gene_windows(gene_windows, indels)
    gene_windows$chr<-factor(gene_windows$chr, levels=paste0("Chr",1:12))
  } else {
    gene_windows$H3K4me1<-encode_hits(H3K4me1, gene_windows)
    gene_windows$H3K9me1<-encode_hits(H3K9me1, gene_windows)
    gene_windows$H3K4me3<-encode_hits(H3K4me3, gene_windows)
    gene_windows$H3K36me3<-encode_hits(H3K36me3, gene_windows)
    gene_windows$H3K9me2<-encode_hits(H3K9me2, gene_windows)
    gene_windows$H3K27me3<-encode_hits(H3K27me3, gene_windows)
    gene_windows$H3K27ac<-encode_hits(H3K27ac, gene_windows)
    gene_windows$PII<-encode_hits(PII, gene_windows)
    gene_windows$H3K4ac<-encode_hits(H3K4ac, gene_windows)
    gene_windows$H3K12ac<-encode_hits(H3K12ac, gene_windows)
    gene_windows$H3K9ac<-encode_hits(H3K9ac, gene_windows)
    gene_windows$snp<-add_vars_hits_to_gene_windows(gene_windows, snps)
    gene_windows$snp_noCT<-add_vars_hits_to_gene_windows(gene_windows, snps[!Single.base.substitution %in% c("C>T","G>A")])
    gene_windows$snp_homoz<-add_vars_hits_to_gene_windows(gene_windows, snps[Genotype=="Homozygous"])
    gene_windows$indel<-add_vars_hits_to_gene_windows(gene_windows, indels)
    gene_windows$chr<-factor(gene_windows$chr, levels=paste0("Chr",1:12))

  }
  setkey(gene_windows, chr, start, stop)

  return(gene_windows)
}

encode_overlap<-function(encode_data, gene_windows){
  encode_overlap<-foverlaps(gene_windows, encode_data,type="any")
  marked<-encode_overlap[,.(marked=ifelse(sum(!is.na(ID))==0, "unmarked","marked")), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  return(as.factor(marked$marked))

}

encode_hits<-function(encode_data, gene_windows, out="marked"){
  encode_overlap<-foverlaps(gene_windows, encode_data,type="any")
  marked<-encode_overlap[,.(marked=sum(!is.na(ID)), length=sum(length,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  if(out=="length"){
    return(as.numeric(marked$length))
  } else  return(as.numeric(marked$marked))
}

plot_peaks<-function(encode_data, ggtitle, xtitle, deciles, var_data, gene=F, marky=F){

  windows<-make_feature_windows(encode_data =  encode_data, deciles, gene=gene)
  windows_vars<-foverlaps(windows, var_data,type="any")
  windows_means<-windows_vars[,.(mut=sum(!is.na(start)), N=.N, length=mean(length, na.rm=T)),by=.(pos, region, i.ID)]
  windows_means<-windows_means[,.(pct=sum(mut)/sum(length), mut=sum(mut), length=sum(length)), by=.(pos, region)]
  if(marky==T){windows_means$pct<-windows_means$mut}
  plot<-ggplot(windows_means, aes(x=pos, y=pct, col=region=="gene body", group=region))+
    geom_vline(xintercept = c(deciles, deciles*2), linetype="dashed", size=0.25)+
    geom_line()+
    scale_color_manual(values=c("gray75","green3"), guide="none")+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="Mutations/bp")+
    ggtitle(ggtitle)+
    scale_x_continuous(breaks=c(0,max(windows_means$pos)/3, max(windows_means$pos)/3*2, max(windows_means$pos)), labels=c("-2kb","0%","100%","+2kb"), name=xtitle)
  return(list(plot, windows_means))
}

add_vars_to_gene_windows<-function(gene_windows, var_object){
  vars_overlap<-foverlaps(gene_windows, var_object,type="any")
  vars<-vars_overlap[,.(mutations=sum(!is.na(Mutant.ID))), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  vars$mutated<-vars$mutations>0
  return(vars$mutated)
}

add_vars_hits_to_gene_windows<-function(gene_windows, var_object){
  vars_overlap<-foverlaps(gene_windows, var_object,type="any")
  vars<-vars_overlap[,.(mutations=sum(!is.na(Mutant.ID))), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  vars$mutated<-vars$mutations>0
  return(vars$mutations)
}

plot_model<-function(model_sum, ggtitle){
  ggplot(model_sum[predictor!="(Intercept)"], aes(x=predictor, y=y, fill=log10(`P`)))+
    geom_bar(stat="identity", col="black",size=0.25)+
    theme_classic(base_size = 6)+
    scale_fill_gradientn(colors=c("dodgerblue","white"), name="-log10(P)")+
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.key.size = unit(.3,"line"), legend.key=element_rect(color="black"))+
    scale_x_discrete(name="Predictor")+
    ggtitle(ggtitle)
}

log_model<-function(gene_windows, variable){
  form<-formula(paste0("as.numeric(",variable,")~H3K4me1+H3K9me1+H3K4me3+H3K27me3+H3K9me2+H3K27ac+H3K36me3+PII+H3K4ac+H3K12ac+H3K9ac"))
  model<-summary(glm(form, gene_windows, family="binomial"))
  model_sum<-data.table(model$coefficients)
  model_sum$predictor<-gsub("unmarked","",row.names(model$coefficients))
  model_sum$Estimate<--model_sum$Estimate
  model_sum$P<-model_sum$`Pr(>|z|)`
  model_sum$predictor<-factor(model_sum$predictor, levels=model_sum$predictor[order(-model_sum$`z value`)])
  model_sum$y<--model_sum$`z value`
  return(model_sum)
}

lm_model<-function(gene_windows, variable, aic=F){
  form<-formula(paste0("as.numeric(",variable,")~H3K4me1+H3K4me3+H3K9me1+H3K27me3+H3K9me2+H3K27ac+H3K36me3+PII+H3K4ac+H3K12ac+H3K9ac"))
  model<-lm(form, gene_windows)
  if(aic==T){model<-MASS::stepAIC(model, direction="both", trace = F)}
  model_sum<-data.table(summary(model)$coefficients)
  model_sum$predictor<-gsub("unmarked","",row.names(summary(model)$coefficients))
  model_sum$Estimate<-model_sum$Estimate
  model_sum$P<-model_sum$`Pr(>|t|)`
  model_sum$predictor<-factor(model_sum$predictor, levels=model_sum$predictor[order(model_sum$`t value`)])
  model_sum$y<-model_sum$`t value`
  return(model_sum)
}

chip_overlaps<-function(bedfile, featureobject){
  cat("\nreading ");cat(bedfile)
  in1<-fread(bedfile)
  colnames(in1)<-c("chr","start","stop","depth")
  in1$length<-as.numeric(in1$stop-in1$start)
  in1$depth<-as.numeric(in1$depth)
  setkey(in1, chr, start, stop)
  out<-unlist(lapply(1:5, function(c) {
    cat(" chr ");cat(c)
    CDS_input_overlap<-foverlaps(featureobject[chr==c], in1[chr==c],type="any")
    CDS_input<-CDS_input_overlap[,.(len=sum(length,na.rm=T), dep=sum(depth,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, gene)]
    CDS_input$input<-CDS_input$len*CDS_input$dep
    rm("CDS_input_overlap")
    input<-CDS_input$input
    return(input)
  }))
  return(out)
}

peaks_randomized<-function(featureobject){
  rand<-rbindlist(apply(rbindlist(list(featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject)), 1, function(x){
    chr<-x["CHROM"]
    start<-sample(1:max(featureobject[CHROM==chr]$stop), 1)
    stop<-start+as.numeric(x["length"])
    return(data.table(chr=chr, start=start, stop=stop))
  }))
  rand$length<-rand$stop-rand$start
  rand$ID<-1:nrow(rand)
  setkey(rand, chr, start, stop)
}

plot_bars_rice<-function(sumstable, yvar, xlab, ylab, ggtitle){
  ggplot(gene_annotations_all_sums, aes_string(x="grp", y=yvar))+
    geom_bar(stat="identity", fill="dodgerblue4", col="black", width=0.5)+
    theme_classic(base_size = 6)+
    scale_x_discrete(name=xlab)+
    scale_y_continuous(name=ylab)+
    ggtitle(ggtitle)

}

homopolymer_read<-function(seq, size){
  s=rle(seq)
  v=cumsum(rle(seq)$lengths)
  runs<-data.table('var'=s$values,'start'=v+1-s$lengths,'end'=v)
  runs$length<-runs$end-runs$start+1
  runs$var<-toupper(runs$var)
  homopolymers<-runs[length>=size]
}

homopolymer_var_annotate<-function(vars, homopolymers, size, dist){

  homopolymers<-rbindlist(lapply(1:length(homopolymers), function(x) {
    out<-homopolymers[[x]][length>=size]
    out$CHROM<-names(homopolymers)[x]
    return(out)
  }))

  homopolymers$startminus<-homopolymers$start-dist
  homopolymers$endplus<-homopolymers$end+dist
  setkey(homopolymers, CHROM, startminus, endplus)
  setkey(vars, chr, start, stop)
  overlap<-foverlaps(vars, homopolymers)
  overlap_sum<-overlap[ALT==var, .(homopolymer_neighbor=sum(!is.na(var))>0),by=unique]
  vars$homopolymer_neighbor=overlap_sum$homopolymer_neighbor[match(vars$unique, overlap_sum$unique)]
  return(vars$homopolymer_neighbor)
}

long_context<-function(vars,dist){
  apply(vars, 1, function(x){
  chr<-x["CHROM"]
  pos<-as.numeric(x["POS"])
  ref<-x["REF"]
  alt<-x["ALT"]
  long_context<-toupper(paste0(genome[[chr]][(pos-dist):(pos+dist)], collapse=""))
})
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

pctile<-function(raw_object, types, variable, char=F, x_name, title){
  if(char){
    tmp<-raw_object[type%in%types,c("mutations","sumlength",variable), with=F]
    tmp$vardata<-tmp[,variable,with=F]
    tmp<-tmp[!is.na(vardata) & vardata!=""]
    pcts<-(tmp)[, .(raw=sum(mutations), length=sum(sumlength), N=.N), by=.(cut=vardata)]
    pcts$pct<-pcts$raw/pcts$length
    pcts$variable<-variable
    chi<-chisq.test(pcts[,2:3])
    plot<-ggplot(pcts, aes(x=cut, y=pct))+
      geom_bar(stat="identity", fill="dodgerblue4", col="black", width=0.5)+
      scale_y_continuous(name="Mutations/bp")+
      theme_classic(base_size = 6)+
      theme(axis.text.x = element_text(angle=45, hjust=1))+
      scale_x_discrete(name=x_name)+
      ggtitle(title, subtitle=paste0("X-squared=",round(chi$statistic, digits = 1),"\n",ifelse(chi$p.value<0.05,"p<0.05*","n.s.")))

  } else{
    tmp<-raw_object[type%in%types,c("mutations","sumlength",variable), with=F]
    tmp$vardata<-tmp[,variable,with=F]
    tmp<-tmp[is.finite(vardata)]
    pcts<-tmp[, .(raw=sum(mutations), length=sum(sumlength), N=.N), by=.(cut=Hmisc::cut2(vardata, g = 2))]
    pcts$pct<-pcts$raw/pcts$length
    pcts$variable<-variable
  }
  bootpop<-rep(unlist(apply(pcts, 1, function(x){
    rep(x[1], times=as.numeric(x[2]))
  })), times=1000)

  pct_boot<-rbindlist(lapply(1:200, function(i){
    cut<-sample(bootpop, size=sum(pcts$raw))
    out<-data.table(table(cut))
    out<-merge(out,pcts, by="cut")
    out$pct<-out$N.x/out$length
    return(out)
  }))

  CI<-rbindlist(lapply(unique(pcts$cut), function(s){
    lvl<-pct_boot[cut==s][order(pct)][-c(1:5, 195:200)]
    out<-data.table(cut=s, min=min(lvl$pct), max=max(lvl$pct))
  }))
  pcts<-merge(pcts, CI)

  chi<-chisq.test(pcts[,2:3])
  plot<-ggplot(pcts, aes(x=cut, y=pct))+
    geom_point(col="dodgerblue4")+
    geom_errorbar(aes(ymin=min, ymax=max), width=0)+
    scale_y_continuous(name="Mutations/bp")+
    theme_classic(base_size = 6)+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_x_discrete(name=x_name)

  return(plot)
}

mutations_in_features<-function(features, mutations){
  # if(length(setdiff(c("CHROM", "START","STOP"), colnames(features)))>0){
  #   stop("features object needs to hace CHROM START and STOP columns")
  # }
  # if(length(setdiff(c("CHROM", "POSITION"), colnames(mutations)))>0){
  #   stop("mutations object needs to hace CHROM and POSITION columns")
  # }

  features$CHROM<-features$chr
  features$STOP<-features$stop
  features$START<-features$start
  if("POSITION" %in% colnames(mutations)){
    mutations$START<-mutations$POSITION
    mutations$STOP<-mutations$POSITION
  }
  if("POS" %in% colnames(mutations)){
    mutations$START<-mutations$POS
    mutations$STOP<-mutations$POS
    mutations$POSITION<-mutations$POS

  }


  if("start" %in% colnames(features)){
    features$START<-features$start
    features$STOP<-features$stop
  }

  if("chr" %in% colnames(features)){
    features$CHROM<-features$chr
  }

  features$CHROM<-as.character(features$CHROM)
  mutations$CHROM<-as.character(mutations$CHROM)

  features$feature_ID<-1:nrow(features)
  setkey(features, CHROM, START, STOP)
  setkey(mutations, CHROM, START, STOP)

  overlaps<-foverlaps(features, mutations)
  muts<-overlaps[,.(muts=sum(!is.na(POSITION))), by=.(feature_ID)]
  return(muts$muts)
}

features_overlap_mutation<-function(features, mutations){
  if(length(setdiff(c("CHROM", "START","STOP"), colnames(features)))>0){
    stop("features object needs to hace CHROM START and STOP columns")
  }
  if(length(setdiff(c("CHROM", "POSITION"), colnames(mutations)))>0){
    stop("mutations object needs to hace CHROM and POSITION columns")
  }

  mutations$START<-mutations$POSITION
  mutations$STOP<-mutations$POSITION

  if("start" %in% colnames(features)){
    features$START<-features$start
    features$STOP<-features$stop
  }

  features$CHROM<-as.character(features$CHROM)
  mutations$CHROM<-as.character(mutations$CHROM)

  features$feature_ID<-1:nrow(features)
  setkey(features, CHROM, START, STOP)
  setkey(mutations, CHROM, START, STOP)
  mutations$mutation_ID<-1:nrow(mutations)

  overlaps<-foverlaps(mutations, features)
  muts<-overlaps[,.(overlaps=sum(!is.na(feature_ID))>0), by=.(mutation_ID)][order(mutation_ID)]
  return(muts$overlaps)
}


