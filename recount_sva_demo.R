#http://rna.recount.bio/docs/quick-access.html

#BiocManager::install("recount3")
library("recount3")
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
rse_gene_SRP009615 <- create_rse(proj_info)

rse_gene_SRP009615

colData(rse_gene_SRP009615)
rowData(rse_gene_SRP009615)
assays(rse_gene_SRP009615)$raw_counts

#snapcount: not sure if this is userful. 
BiocManager::install("snapcount")
library("snapcount")
qb <- QueryBuilder(compilation = "gtex", regions = "chr1:1-100000")
qb <- set_row_filters(qb, samples_count >= 20)
query = query_gene(qb)
colData(query)
rowData(query)



#sva
#sva tutorial and 
#http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html

BiocManager::install("sva")

#sample data
#BiocManager::install("bladderbatch")
library(bladderbatch)
data(bladderdata)

#BiocManager::install("pamr")
library(pamr)
library(limma)
library(sva)

pheno = pData(bladderEset)
edata  = exprs(bladderEset)
mod = model.matrix(~as.factor(cancer), data = pheno)
mod0 = model.matrix(~1, data = pheno)

#number of surrogate variables to be estimated
n.sv = num.sv(edata, mod, method = "leek")

#est surrogate variables
svobj = sva(edata, mod, mod0, n.sv = n.sv)

#examine differential gene expression, with out surrogate variable analysis
pValues = f.pvalue(edata, mod, mod0)
qValues = p.adjust(pValues, method = "BH")


#now with surrogate variables
modSv = cbind(mod, svobj$sv)
mod0Sv = cbind(mod0, svobj$sv)
pValuesSv = f.pvalue(edata, modSv, mod0Sv)
qValuesSv = p.adjust(pValuesSv, method = "BH")




