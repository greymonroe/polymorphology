# get estimate null nonsynonymous/synonymous ratio
# input:

#' get estimate null nonsynonymous/synonymous ratio
#' @param mutations  file with columns 'ref' 'alt' and 'N' counts for each mutation
#' @param composition  file where each row is coding regions and columns 2,4,6,8, have counts of A, T, C, G
#' @return a list non synonymous/synonymous and stop/synonymous ratios
#' @export

Null_ns_s<-function(mutations, composition){

code = c(
  'ATA'='I', 'ATC'='I', 'ATT'='I', 'ATG'='M',
  'ACA'='T', 'ACC'='T', 'ACG'='T', 'ACT'='T',
  'AAC'='N', 'AAT'='N', 'AAA'='K', 'AAG'='K',
  'AGC'='S', 'AGT'='S', 'AGA'='R', 'AGG'='R',
  'CTA'='L', 'CTC'='L', 'CTG'='L', 'CTT'='L',
  'CCA'='P', 'CCC'='P', 'CCG'='P', 'CCT'='P',
  'CAC'='H', 'CAT'='H', 'CAA'='Q', 'CAG'='Q',
  'CGA'='R', 'CGC'='R', 'CGG'='R', 'CGT'='R',
  'GTA'='V', 'GTC'='V', 'GTG'='V', 'GTT'='V',
  'GCA'='A', 'GCC'='A', 'GCG'='A', 'GCT'='A',
  'GAC'='D', 'GAT'='D', 'GAA'='E', 'GAG'='E',
  'GGA'='G', 'GGC'='G', 'GGG'='G', 'GGT'='G',
  'TCA'='S', 'TCC'='S', 'TCG'='S', 'TCT'='S',
  'TTC'='F', 'TTT'='F', 'TTA'='L', 'TTG'='L',
  'TAC'='Y', 'TAT'='Y', 'TAA'='_', 'TAG'='_',
  'TGC'='C', 'TGT'='C', 'TGA'='_', 'TGG'='W'
)

mutdif<-function(a,b){
  a.<-strsplit(a,split="")[[1]]
  b.<-strsplit(b,split="")[[1]]
  mut<- a.!= b.
  pos<-which(mut==T)
  # return(c(a.[pos],b.[pos]))
  if(length(pos)>1){
    return(NA)
  }else{
    return(paste0(a.[pos],b.[pos]))
  }
}
aadif<-function(a,b){
  a!=b
}
aadif<-function(a,b){
  if(a==b){"S"}
  else if(a=="_" & b!="_"){"STOPLOSS"}
  else if(a!="_" & b=="_"){"STOPGAIN"}
  else{"NS"}
}

# Create tables of NS and S transitions based on codon structure
Q<-matrix(nrow = length(code),ncol = length(code))
NS<-matrix(nrow = length(code),ncol = length(code))
for(i in 1:nrow(Q)){
  for(j in 1:ncol(Q)){
    if(i==j){
      next
    }else{
      Q[i,j] <- mutdif(names(code)[i],names(code)[j])
      NS[i,j] <- aadif(code[i],code[j])
    }
  }
}

#### Read base content genome
message("Read nucleotide frequencies")
bp<-read.table(composition,fill=T)
# bp[,]
bp2<-na.omit(bp) # there are some exons so small that do not contain one base. their effect in the total is very small so remove
bp<-bp2
bp<-head(bp,1000)
as<-sum(bp$V2)
cs<-sum(bp$V4)
gs<-sum(bp$V6)
ts<-sum(bp$V8)

relprop<-c(as,cs,gs,ts)/sum(as,cs,gs,ts)
relprop

#### Read mutation spectrum
message("Read SNP frequencies")
snpfreq<-read.table(mutations,header=T)
snpfreq

# Correcting those by the relative proportion of SNPs
# Create translator of frequency of each mutation
snpfreqtr<-snpfreq$N/sum(snpfreq$N)
names(snpfreqtr)<-paste0(snpfreq$ref,snpfreq$alt)
snpfreqtr<-c(snpfreqtr,"NA"="NA")
# Create new transition matrix based on the SNP frequencies in Ath
Q.<-Q
Q.[]<-as.numeric(snpfreqtr[as.character(Q.)])
Q.<-apply(Q.,2,as.numeric)
# dim(Q.)
# Q.

# create new transition matrix for
Q2<-apply(Q,2,function(x)substr(x,start = 1,stop = 1))
# Q2
bptr<-relprop
names(bptr)<-c('A','C','G','T')
Q2[]<-as.numeric(bptr[as.character(Q2)])
Q2<-apply(Q2,2,as.numeric)

Freqchanges<-(Q.*Q2)[!is.na(Q)]
Freqchanges<-Freqchanges/sum(Freqchanges)

NSS<- NS[!is.na(Q)]

message("Ratios of neutral mutations corrected for base frequencies and SNP transition frequency")
synonymous=sum(Freqchanges[NSS=="S"])
# corrected for base composition of genes and frequency of mutations
message("Nonsynonymous / Synonymous")
ns_s<-sum(Freqchanges[NSS=="NS"]) /synonymous
# corrected for base composition of genes and frequency of mutations
message("STOPGAIN / Synonymous")
stop_s<-sum(Freqchanges[NSS=="STOPGAIN"]) / synonymous
return(list(ns_s, stop_s))
}
