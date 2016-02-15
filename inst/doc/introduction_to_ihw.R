## ---- message=FALSE, warning=FALSE---------------------------------------
library("methods")
library("airway")
library("DESeq2")
data("airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)

## ------------------------------------------------------------------------
colnames(res)

## ---- message=FALSE, warning=FALSE---------------------------------------
library("IHW")
ihw_res <- ihw(res$pvalue, res$baseMean, alpha = 0.1)

## ------------------------------------------------------------------------
rejections(ihw_res)

## ------------------------------------------------------------------------
head(adj_pvalues(ihw_res))
sum(adj_pvalues(ihw_res) <= 0.1, na.rm = TRUE) == rejections(ihw_res)

## ------------------------------------------------------------------------
padj_bh <- p.adjust(res$pvalue, method = "BH")
sum(padj_bh <= 0.1, na.rm = TRUE)

## ------------------------------------------------------------------------
head(weights(ihw_res))

## ------------------------------------------------------------------------
weights(ihw_res, levels_only=TRUE)

## ------------------------------------------------------------------------
plot_ihw(ihw_res)

## ------------------------------------------------------------------------
colnames(as.data.frame(ihw_res))

## ------------------------------------------------------------------------
set.seed(1)
X   <- runif(100000, min=0, max=2.5) #covariate
H   <- rbinom(100000,1,0.1)            #hypothesis true or false
Z   <- rnorm(100000, H*X)              #Z-score
pvalue <- 1-pnorm(Z)                  #pvalue
sim <- data.frame(X=X, H=H, Z=Z, pvalue=pvalue)

## ------------------------------------------------------------------------
sum(p.adjust(sim$pvalue, method="BH") <= 0.1)

## ------------------------------------------------------------------------
filter_threshold <- 0.1
pvalue_filt <- pvalue[pvalue <= filter_threshold]

## ------------------------------------------------------------------------
sum(p.adjust(pvalue_filt, method="BH", n=length(pvalue)) <= 0.1)  

## ------------------------------------------------------------------------
nbins <- 20
sim$group <- groups_by_filter(sim$X, nbins)
m_groups <- table(sim$group)

## ------------------------------------------------------------------------
sim_filtered <- subset(sim, sim$pvalue <= filter_threshold) 
ihw_filt <- ihw(sim_filtered$pvalue, sim_filtered$group, .1, m_groups = m_groups)
rejections(ihw_filt)

