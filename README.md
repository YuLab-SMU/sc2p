# SC2P: identify cell type specific pathways

> SC (single cell) + CP (clusterProfiler) => SCCP => SC2P (single cell to pathway)

1. score cell states by pathways (e.g., GSVA)
2. identify cell type significant pathways 

```
> eres <- scKEGG.Seurat(pbmc)
'select()' returned 1:1 mapping between keys and columns
Estimating GSVA scores for 356 gene sets.
Estimating ECDFs with Gaussian kernels
  |======================================================================| 100%
> head(eres, 3)
        ID     logFC   AveExpr         t       P.Value     adj.P.Val         B
1 hsa00290 0.5263822 0.5108176 181.02353  0.000000e+00  0.000000e+00 3412.9098
2 hsa00290 0.5656179 0.5108176  26.81641 2.885105e-140 1.027097e-137  309.6167
3 hsa00290 0.5369375 0.5108176 130.95099  0.000000e+00  0.000000e+00 2646.4178
        ident                                 Description
1 Naive CD4 T Valine, leucine and isoleucine biosynthesis
2    Platelet Valine, leucine and isoleucine biosynthesis
3           B Valine, leucine and isoleucine biosynthesis
```

