## ------------------------------------------------------------------------
library("DESeq2")

## ------------------------------------------------------------------------
library("airway")
data("airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
dds <- DESeq(dds)
dds <- results(dds)
dds <-subset(dds, !is.na(pvalue))

## ------------------------------------------------------------------------
library("ddhw")
res_milp <- ddhw(dds$pvalue,dds$baseMean, 10, 0.1)
weights(res_milp)
sum(adj_pvalues(res_milp) < 0.1)

## ------------------------------------------------------------------------
res_subplex <- ddhw(dds$pvalue,dds$baseMean, 10, 0.1, optim_method="subplex")
weights(res_subplex)
sum(adj_pvalues(res_subplex) < 0.1)

## ------------------------------------------------------------------------
sum(p.adjust(dds$pvalue, method="BH") < 0.1)

