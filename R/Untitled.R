
genes_dimers<-dimerfrequency(genes, genome)
genes_dimers$gene<-genes$gene

setkey(genes, chr, start, stop)

gene_signal<-rbindlist(lapply(1:5, function(x){
  cat(x)
  overlap<-foverlaps(bw_dt[chr==x], genes[chr==x])
  overlap[!is.na(gene)]
  gene_signal<-overlap[,.(score=mean(score)), by="gene"]
}))


genes$score<-gene_signal$score[match(genes$gene, gene_signal$gene)]
genes_dimers$dir<-genes$direction[match(genes_dimers$gene, genes$gene)]
genes_dimers$score<-genes$score[match(genes_dimers$gene, genes$gene)]

all<-genes_dimers[, c(1:16, 19), with=F]
allMod<-lm(score~., all)

forward<-data.table(genes_dimers[direction=="+"])[, c(1:16, 19), with=F]
forwardMod<-lm(score~., forward)
summary(forwardMod)

reverse<-data.table(genes_dimers[direction=="-"])[, c(1:16, 19), with=F]
reverseMod<-lm(score~., reverse)
summary(reverseMod)


plot_mod_table<-function(mod){
  sum<-summary(mod)
  coef<-data.table(sum$coefficients)
  coef$predictor<-row.names(sum$coefficients)
  coef<-coef[-1]s

  ggplot(coef, aes(x=predictor, y=`t value`, fill=-log10(`Pr(>|t|)`)))+
    geom_bar(stat="identity", col="black")+
    theme_classic(base_size = 6)
}

plot_mod_table(reverseMod)+scale_x_discrete(name="Dimer")
plot_mod_table(forwardMod)+scale_x_discrete(name="Dimer")
plot_mod_table(allMod)+scale_x_discrete(name="Dimer")




dimerData<-forward
dimer_single<-function(dimerData){
  motifs<-expand.grid(c("A","T","C","G"), c("A","T","C","G"))
  motifs<-paste0(motifs$Var1, motifs$Var2)

  cors<-rbindlist(lapply(motifs, function(m){
    dt<-data.table(motif=unlist(dimerData[,m,with=F]), score=dimerData$score)
    test<-cor.test(dt$motif, dt$score)
    p=test$p.value
    estimate=test$estimate
    return(data.table(m, p, estimate))
  }))


  ggplot(cors, aes(x=m, y=estimate, fill=-log10(p)))+
    geom_bar(stat="identity", col="black")+
    theme_classic(base_size = 6)

}
dimer_single(forward)

