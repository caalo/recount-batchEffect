library(recount3)
library(snapcount)
library(tidyverse)
library(ggforce)
library(DESeq2)
theme_set(theme_bw())

genes = c("CD99", "ACTB", "TBP", "GAPDH", "TFRC", "KRAS", "TP53", "EGFR", "MTOR", "POLE4")

analysis = data.frame(gene = character(),
                      expression = numeric(),
                      coverage = numeric(),
                      tissue = character())
for(g in genes) {
  cat(g, "\n")
  sb <- QueryBuilder(compilation="gtex", regions=g)
  gene.query <- query_gene(sb)
  gene.phenotype = data.frame(colData(gene.query))
  gene.genotype = data.frame(rowData(gene.query))
  gene.expression = assays(gene.query)$counts
  gene.length = as.numeric(strsplit(gene.genotype$gene_id.gene_name.gene_type.bp_length, ":")[[1]][4])
  
  
  gene_analysis = data.frame(gene = g,
                            expression = gene.expression[1 ,] / gene.length, #divide expression by number of bases in gene. 
                            coverage = gene.phenotype$bases / (3.2 * 10^9), #not right
                            tissue = gene.phenotype$Body_Site)
  
  analysis = rbind(analysis, gene_analysis)
}

pdf("gene_normalization.pdf", width = 8, height = 6)
for(page in 1:9) {
  plot_df = analysis %>% filter(gene == "KRAS")
  g = ggplot(plot_df, aes(x = log10(coverage), y = log10(expression), color = tissue)) +
        geom_point(alpha = .5, show.legend = FALSE) + 
        geom_smooth(method = "lm", show.legend = FALSE) + ggtitle("KRAS") +
        facet_wrap_paginate(~tissue, ncol = 3, nrow = 2, page = page)
  print(g)
}
dev.off()



coefs_analysis = expand.grid(tissue = unique(analysis$tissue),
                             gene = genes)
coefs_analysis$coef = NA
coefs_analysis$n_samples = NA
coefs_analysis$expr_level = NA

for(i in 1:nrow(coefs_analysis)) {
  t = coefs_analysis$tissue[i]
  g = coefs_analysis$gene[i]
  tissue_analysis = analysis %>% filter(tissue == t & gene == g)
  tissue_analysis$logExpression = log10(tissue_analysis$expression)
  tissue_analysis$logCoverage = log10(tissue_analysis$coverage)
  model = lm(logExpression ~ logCoverage, tissue_analysis)
  coefs_analysis$coef[i] = model$coefficients["logCoverage"]
  coefs_analysis$n_samples[i] = nrow(tissue_analysis)
  coefs_analysis$expr_level[i] = predict.lm(model, data.frame(logCoverage = c(.25)))
}

coefs_analysis = coefs_analysis %>% filter(n_samples >= 50)
pdf("gene_normalization_coefs.pdf", width=8, height=8)
ggplot(coefs_analysis, aes(coef)) + geom_density() + facet_wrap(~gene) 

coefs_analysis$expr_bin = cut(coefs_analysis$expr_level, breaks = quantile(coefs_analysis$expr_level, seq(0, 1, .2)))
coefs_analysis$expr_bin2 = cut(coefs_analysis$expr_level, breaks = 5)

ggplot(coefs_analysis, aes(coef)) + geom_density() + facet_wrap(~expr_bin) + ggtitle("Bin by expression at .25 log10(coverage)", subtitle = "5 quantiles")
ggplot(coefs_analysis, aes(coef)) + geom_density() + facet_wrap(~expr_bin2) + ggtitle("Bin by expression at .25 log10(coverage)", subtitle = "5 evenly spaced cutoffs")

dev.off()
#global scale factor?

