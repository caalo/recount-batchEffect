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

#BiocManager::install("sva")

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

#sva by hand: (irwsva.build function)

edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  
  n <- length(p)
  transf <- match.arg(transf)
  
  if(transf=="probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1-eps)
    x <- qnorm(p)
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0*dnorm(x)/y
  }
  
  if(transf=="logit") {
    x <- log((p+eps)/(1-p+eps))
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1+exp(x))^2
    lfdr <- pi0 * dx/y
  }
  
  if(trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if(monotone) {	
    lfdr <- lfdr[order(p)]
    lfdr <- mono(lfdr)
    lfdr <- lfdr[rank(p)]
  }
  
  return(lfdr)
}

mono <- function(lfdr){
  .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
}

dat = edata
n <- ncol(dat)
m <- nrow(dat)
if(is.null(mod0)){mod0 <- mod[,1]}
Id <- diag(n)
resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))  
uu <- eigen(t(resid)%*%resid)
vv <- uu$vectors
ndf <- n - dim(mod)[2]

pprob <- rep(1,m)
one <- rep(1,n)
Id <- diag(n)
df1 <- dim(mod)[2] + n.sv
df0 <- dim(mod0)[2]  + n.sv

i = 1 #out of B iterations
#to find P(b_i = 0 | ...) 
mod.b <- cbind(mod,uu$vectors[,1:n.sv])
mod0.b <- cbind(mod0,uu$vectors[,1:n.sv])
ptmp <- f.pvalue(dat,mod.b,mod0.b)
pprob.b <- (1-edge.lfdr(ptmp))





#examine differential gene expression, with out surrogate variable analysis
pValues = f.pvalue(edata, mod, mod0)
qValues = p.adjust(pValues, method = "BH")


#now with surrogate variables
modSv = cbind(mod, svobj$sv)
mod0Sv = cbind(mod0, svobj$sv)
pValuesSv = f.pvalue(edata, modSv, mod0Sv)
qValuesSv = p.adjust(pValuesSv, method = "BH")



