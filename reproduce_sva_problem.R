library(recount3)
#from https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
#filter for humans, sample size >= 50. 
studyMetdata = read.csv("recount3_selection_2022-08-09 16_01_21.csv") 
#get rid of single cell experiments. 
studyMetdata = studyMetdata[-grep("SINGLE", toupper(studyMetdata$study_abstract)) ,]
human_projects <- available_projects()


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
geneCounts_SRP118614 = geneCounts_SRP118614 / genotype$bp_length #normalize by gene length
#numZeros = apply(geneCounts, 1, function(x) sum(x == 0))
#geneCounts = geneCounts[numZeros == 0 ,]

n.sv = num.sv(geneCounts_SRP118614, mod_SRP118614, method = "be")


#add more!

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
geneCounts_SRP212704 = geneCounts_SRP212704 / genotype$bp_length

geneCounts_2 = cbind(geneCounts_SRP118614, geneCounts_SRP212704)
mod_2 = rbind(mod_SRP118614, mod_SRP212704)

numZeros = apply(geneCounts_2, 1, function(x) sum(x == 0))
geneCounts_2 = geneCounts_2[numZeros == 0 ,]

n.sv = num.sv(geneCounts_2, mod_2, method = "be")
svobj = sva(geneCounts_2, mod_2, mod_2[, 1], n.sv = n.sv)

#3
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
geneCounts_SRP002628 = geneCounts_SRP002628 / genotype$bp_length

geneCounts_3 = cbind(geneCounts_SRP118614, geneCounts_SRP212704, geneCounts_SRP002628)
mod_3 = rbind(mod_SRP118614, mod_SRP212704, mod_SRP002628)

numZeros = apply(geneCounts_3, 1, function(x) sum(x == 0))
geneCounts_3 = geneCounts_3[numZeros == 0 ,]

n.sv = num.sv(geneCounts_3, mod_3, method = "be")
svobj = sva(geneCounts_3, mod_3, mod_3[, 1], n.sv = n.sv)
qplot(1:154, svobj$sv)

#4
proj_info <- subset(
  human_projects,
  project == "SRP151104"
)
rse_gene_SRP151104 <- create_rse(proj_info)
phenotype = data.frame(colData(rse_gene_SRP151104))
genotype = data.frame(rowData(rse_gene_SRP151104))
mod_SRP151104 = matrix(rep(1, nrow(phenotype) * 2), nrow = nrow(phenotype))
geneCounts_SRP151104 = assays(rse_gene_SRP151104)$raw_counts
geneCounts_SRP151104 = geneCounts_SRP151104 / genotype$bp_length

geneCounts_4 = cbind(geneCounts_SRP118614, geneCounts_SRP212704, geneCounts_SRP002628, geneCounts_SRP151104)
mod_4 = rbind(mod_SRP118614, mod_SRP212704, mod_SRP002628, mod_SRP151104)

numZeros = apply(geneCounts_4, 1, function(x) sum(x == 0))
geneCounts_4 = geneCounts_4[numZeros == 0 ,]

n.sv = num.sv(geneCounts_4, mod_4, method = "be")
svobj = sva(geneCounts_4, mod_4, mod_4[, 1], n.sv = n.sv)
qplot(1:187, svobj$sv)
