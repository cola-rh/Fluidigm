
## Hierarchical consensus partitioning analysis on dataset 'Fluidigm'

Runnable code:

```r
library(cola)
library(scRNAseq)
data = ReprocessedFluidigmData()

anno = colData(data)[, c("Biological_Condition", "Coverage_Type", "Cluster1", "Cluster2")]
anno = as.data.frame(anno)

mat = assays(data)$rsem_tpm[, anno$Coverage_Type == "High"]
count = assays(data)$rsem_counts[, anno$Coverage_Type == "High"]
l = apply(count, 1, function(x) sum(x > 0)/length(x) > 0.1)
mat = log2(mat[l, ] + 1)
mat = adjust_matrix(mat)

anno = anno[anno$Coverage_Type == "High", ]

rh = hierarchical_partition(mat, cores = 4, anno = anno)
saveRDS(rh, file = "Fluidigm_cola_rh.rds")

cola_report(rh, output = "Fluidigm_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'Fluidigm'")
```

[The HTML report is here.](https://cola-rh.github.io/Fluidigm/Fluidigm_cola_rh_report/cola_hc.html) (generated by __cola__ version 2.0.0, R 4.1.0.)

