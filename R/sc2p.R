
pathway_score <- function(mat, gs, method, ...) {
    if (method == 'gsva') {
        library(GSVA)
        gsva_par <- gsvaParam(mat, gs)
        res <- gsva(gsva_par) # todo: BiocParallel 
    }

    return(res)    
}

pathway_sig <- function(score, cls) {
    library(limma)
    library(Biobase)
    pset <- ExpressionSet(assayData = score)
    pData(pset) <- cls
    design <- model.matrix(~0+ident, pset)
    plm <- lmFit(pset, design=design)
    eb <- eBayes(plm)
    ids <- levels(cls$ident)
    
    res <- lapply(seq_along(ids), function(i) {
        id <- sprintf("ident%s", ids[i])
        y <- topTable(eb, coef = id, number=5) # <--- needs to change
        y$ident <- ids[i]
        return(y)
    }) |> yulab.utils::rbindlist()

    # res$ID <- rownames(res) # limma output may add suffix in the rownames
    
    res$ID <- sub("^(\\D+\\d{5}).*", "\\1", rownames(res)) # for KEGG only

    return(res)
}


sc2p.Seurat <- function(x, gset, TERM2NAME) {
    y <- x[[Assays(x)]]$counts
    cell_score <- pathway_score(y, gset, method = 'gsva')
    cls <- FetchData(x, vars = 'ident')
    res <- pathway_sig(cell_score, cls)
    if (missing(TERM2NAME)) return(res)
    names(TERM2NAME) <- c("ID", "Description")
    merge(res, TERM2NAME,by="ID")
}


scKEGG.Seurat <- function(x) {
    ## gene set
    kk <- clusterProfiler:::download_KEGG('hsa')
    gs <- kk$KEGGPATHID2EXTID
    gene <- clusterProfiler::bitr(gs$to, 'ENTREZID', 'SYMBOL', 'org.Hs.eg.db')
    gs2 <- merge(gs, gene, by.x = 'to', by.y='ENTREZID')
    gset <- split(gs2$SYMBOL, gs2$from)

    sc2p.Seurat(x, gset, TERM2NAME=kk$KEGGPATHID2NAME)
}

