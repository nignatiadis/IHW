## ---- message=FALSE, warning=FALSE---------------------------------------
library("ggplot2")
library("methods")
library("airway")
library("DESeq2")
data("airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
dds <- DESeq(dds)
de_res <- as.data.frame(results(dds))

## ------------------------------------------------------------------------
colnames(de_res)

## ---- message=FALSE, warning=FALSE---------------------------------------
library("IHW")
ihw_res <- ihw(pvalue ~ baseMean,  data=de_res, alpha = 0.1)

## ------------------------------------------------------------------------
rejections(ihw_res)

## ------------------------------------------------------------------------
head(adj_pvalues(ihw_res))
sum(adj_pvalues(ihw_res) <= 0.1, na.rm = TRUE) == rejections(ihw_res)

## ------------------------------------------------------------------------
padj_bh <- p.adjust(de_res$pvalue, method = "BH")
sum(padj_bh <= 0.1, na.rm = TRUE)

## ------------------------------------------------------------------------
head(weights(ihw_res))

## ------------------------------------------------------------------------
weights(ihw_res, levels_only=TRUE)

## ----plot_ihw_res, fig.width = 5, fig.height = 3-------------------------
plot(ihw_res)

## ------------------------------------------------------------------------
ihw_res_df <- as.data.frame(ihw_res)
colnames(ihw_res_df)

## ----eval=FALSE----------------------------------------------------------
#  ihw_bonf <- ihw(pvalue ~ baseMean, data=de_res, alpha = 0.1, adjustment_type = "bonferroni")

## ----rkpvslogp, fig.width = 6, fig.height = 3----------------------------
de_res <- na.omit(de_res)
de_res$geneid <- as.numeric(as.numeric(gsub("ENSG[+]*", "", rownames(de_res))))

# set up data frame for plotting
df <- rbind(data.frame(pvalue = de_res$pvalue, covariate = rank(de_res$baseMean)/nrow(de_res), 
                       covariate_type="base mean"),
            data.frame(pvalue = de_res$pvalue, covariate = rank(de_res$geneid)/nrow(de_res), 
                       covariate_type="gene id"))

ggplot(df, aes(x=covariate, y = -log10(pvalue))) +
                         geom_hex(bins = 100) + 
                         facet_grid( . ~ covariate_type)

## ----pvalhistogram, fig.width = 3, fig.height = 3------------------------
ggplot(de_res, aes(x=pvalue)) + geom_histogram(binwidth=0.025, boundary = 0)

## ----goodhist, fig.width = 8, fig.height = 5-----------------------------
de_res$baseMean_group <- groups_by_filter(de_res$baseMean,8)

ggplot(de_res, aes(x=pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMean_group, nrow=2)

## ----goodecdf, fig.width = 5, fig.height = 3-----------------------------
ggplot(de_res, aes(x = pvalue, col = baseMean_group)) + stat_ecdf(geom = "step") 

## ----badhist, fig.width = 8, fig.height = 5------------------------------
de_res$lfc_group <- groups_by_filter(abs(de_res$log2FoldChange),8)

ggplot(de_res, aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ lfc_group, nrow=2)

## ----badecdf, fig.width = 5, fig.height = 3------------------------------
ggplot(de_res, aes(x = pvalue, col = lfc_group)) + stat_ecdf(geom = "step") 

## ------------------------------------------------------------------------
set.seed(1)
X   <- runif(100000, min=0, max=2.5)   # covariate
H   <- rbinom(100000,1,0.1)            # hypothesis true or false
Z   <- rnorm(100000, H*X)              # Z-score
pvalue <- 1-pnorm(Z)                   # pvalue
sim <- data.frame(X=X, H=H, Z=Z, pvalue=pvalue)

## ------------------------------------------------------------------------
padj <- p.adjust(sim$pvalue, method="BH")
sum( padj <= 0.1)

## ------------------------------------------------------------------------
filter_threshold <- 0.1
selected <- which(pvalue <= filter_threshold)
pvalue_filt <- pvalue[selected]

## ----compareselected, fig.width=3, fig.height=3--------------------------
padj_filt <- p.adjust(pvalue_filt, method = "BH", n = length(pvalue))
qplot(padj[selected], padj_filt)
sum(padj_filt <= 0.1)  

## ----justcheck, echo=FALSE-----------------------------------------------
stopifnot(max(abs(padj[selected]- padj_filt)) <= 0.001)

## ------------------------------------------------------------------------
nbins <- 20
sim$group <- groups_by_filter(sim$X, nbins)
m_groups <- table(sim$group)

## ------------------------------------------------------------------------
sim_filtered <- subset(sim, sim$pvalue <= filter_threshold) 
ihw_filt <- ihw(pvalue ~ group, data=sim_filtered,
                alpha = .1, m_groups = m_groups)
rejections(ihw_filt)

