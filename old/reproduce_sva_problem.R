library(tidyverse)
library(recount3)
library(sva)
library(limma)

#from https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
#filter for humans, sample size >= 50. 
# studyMetdata = read.csv("recount3_selection_2022-08-09 16_01_21.csv") 
# #get rid of single cell experiments. 
# studyMetdata = studyMetdata[-grep("SINGLE", toupper(studyMetdata$study_abstract)) ,]

human_projects <- available_projects()
n_sv_agg = numeric(4)
n_sv = numeric(4)
#1. 
#Overall design: Matched high-grade (GS=7(4+3)) prostate tumor and adjacent
#normal specimens from 16 patients (8 AAM and 8 EAM) were subjected to two 
#replicate runs of RNA-sequencing.

proj_info <- subset(
  human_projects,
  project == "SRP118614"
)

rse_gene_SRP118614 <- create_rse(proj_info)
phenotype = data.frame(colData(rse_gene_SRP118614))
genotype = data.frame(rowData(rse_gene_SRP118614))

phenotype$phenotype_tumor = grepl("prostate tumor", phenotype$sra.sample_attributes)
mod_SRP118614 = model.matrix(~phenotype_tumor, data = phenotype)

geneCounts_SRP118614 = assays(rse_gene_SRP118614)$raw_counts
geneCounts_SRP118614 = geneCounts_SRP118614 / colSums(geneCounts_SRP118614) #normalize by bases sequenced in genes

n_sv[1] = num.sv(geneCounts_SRP118614, mod_SRP118614, method = "be")
n_sv_agg[1] = NA

#2.
#Overall design: Strand specific total RNA seq was performed using frozen 
#patient matched prostate cancer tissue in biological duplicates.
#Purpose: The goal of present study is to compare transcript level changes
#between normal and tumor of same individuals

proj_info <- subset(
  human_projects,
  project == "SRP212704"
)

rse_gene_SRP212704 <- create_rse(proj_info)
phenotype = data.frame(colData(rse_gene_SRP212704))
genotype = data.frame(rowData(rse_gene_SRP212704))

phenotype$phenotype_tumor = grepl("Tumor", phenotype$sra.sample_attributes)
mod_SRP212704 = model.matrix(~phenotype_tumor, data = phenotype)

geneCounts_SRP212704 = assays(rse_gene_SRP212704)$raw_counts
geneCounts_SRP212704 = geneCounts_SRP212704 / colSums(geneCounts_SRP212704)

geneCounts_2 = cbind(geneCounts_SRP118614, geneCounts_SRP212704)
mod_2 = rbind(mod_SRP118614, mod_SRP212704)

numZeros = apply(geneCounts_2, 1, function(x) sum(x == 0))
geneCounts_2 = geneCounts_2[numZeros == 0 ,]

n_sv[2] = num.sv(geneCounts_SRP212704, mod_SRP212704, method = "be")
n_sv_agg[2] = num.sv(geneCounts_2, mod_2, method = "be")


#3
#Overall design: We sequenced the transcriptome (polyA+) of 20 prostate cancer 
#tumors and 10 matched normal tissues using Illumina GAII platform. 
#Then we used bioinformatic approaches to identify prostate cancer specific 
#aberrations which include gene fusion, alternative splicing, somatic mutation
proj_info <- subset(
  human_projects,
  project == "SRP002628"
)
rse_gene_SRP002628 <- create_rse(proj_info)
phenotype = data.frame(colData(rse_gene_SRP002628))
genotype = data.frame(rowData(rse_gene_SRP002628))

phenotype$phenotype_tumor = grepl("Prostate cancer tissue", phenotype$sra.sample_attributes)
mod_SRP002628 = model.matrix(~phenotype_tumor, data = phenotype)

geneCounts_SRP002628 = assays(rse_gene_SRP002628)$raw_counts
geneCounts_SRP002628 = geneCounts_SRP002628 / colSums(geneCounts_SRP002628)

geneCounts_3 = cbind(geneCounts_SRP118614, geneCounts_SRP212704, geneCounts_SRP002628)
mod_3 = rbind(mod_SRP118614, mod_SRP212704, mod_SRP002628)

numZeros = apply(geneCounts_3, 1, function(x) sum(x == 0))
geneCounts_3 = geneCounts_3[numZeros == 0 ,]

n_sv[3] = num.sv(geneCounts_SRP002628, mod_SRP002628, method = "be")
n_sv_agg[3] = num.sv(geneCounts_3, mod_3, method = "be")

#try plotting the SV variables
svobj = sva(geneCounts_3, mod_3, mod_3[, 1], n.sv = n_sv_agg[3])
annotate = c(rep("study1", ncol(geneCounts_SRP118614)), 
             rep("study2", ncol(geneCounts_SRP212704)), 
             rep("study3", ncol(geneCounts_SRP002628)))
qplot(svobj$sv[,1], svobj$sv[,2], color = annotate)

#4
#To evaluate TERC expression in an RNAseq dataset derived from 25 fresh 
#frozen prostate cancer samples and matched benign tissues
proj_info <- subset(
  human_projects,
  project == "SRP114904"
)
rse_gene_SRP114904 <- create_rse(proj_info)
phenotype = data.frame(colData(rse_gene_SRP114904))
genotype = data.frame(rowData(rse_gene_SRP114904))

phenotype$phenotype_tumor = grepl("Prostate Cancer", phenotype$sra.sample_attributes)
mod_SRP114904 = model.matrix(~phenotype_tumor, data = phenotype)


geneCounts_SRP114904 = assays(rse_gene_SRP114904)$raw_counts
geneCounts_SRP114904 = geneCounts_SRP114904 / colSums(geneCounts_SRP114904)

geneCounts_4 = cbind(geneCounts_SRP118614, geneCounts_SRP212704, geneCounts_SRP002628, geneCounts_SRP114904)
mod_4 = rbind(mod_SRP118614, mod_SRP212704, mod_SRP002628, mod_SRP114904)

numZeros = apply(geneCounts_4, 1, function(x) sum(x == 0))
geneCounts_4 = geneCounts_4[numZeros == 0 ,]

n_sv[4] = num.sv(geneCounts_SRP114904, mod_SRP114904, method = "be")
n_sv_agg[4] = num.sv(geneCounts_4, mod_4, method = "be")


#We then consider the following:  
#Is there a difference between study 1 and study 2 just using control samples with SVs? 
#Run DE: Do adj. p-values look uniform?
mod_2controls = c(rep(0, ncol(geneCounts_SRP118614)), 
                  rep(1, ncol(geneCounts_SRP212704)))
controls = which(mod_2[, 2] == 0)
geneCounts_2controls = geneCounts_2[, controls]
mod_2controls = mod_2controls[controls]
mod_2controls = model.matrix(~ mod_2controls)

n.sv = num.sv(geneCounts_2controls, mod_2controls, method = "be")
svobj = sva(geneCounts_2controls, mod_2controls, mod_2controls[, 1], n.sv = n.sv)
mod_2controls_sv = cbind(mod_2controls, svobj$sv)
#run limma for DE
fit = lmFit(geneCounts_2controls, mod_2controls_sv)
eb = eBayes(fit)
result = topTable(eb, adjust = "BH", number = nrow(geneCounts_2controls))

plot(1:nrow(result)/(nrow(result)+1), sort(result$adj.P.Val))
abline(a = 0, b = 1)

