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

