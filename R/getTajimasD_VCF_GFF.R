
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# when we run this we give the argument vcf which is the VCF file we want to analyze
# and gff which is the gff file for example:
# Rscript getTajimasD_for_gff.R TAIR10.gff 1001G.vcf.gz
vcf<-args[1]
gff<-read.table(args[2])

getTajimasD(vcf, gff)


