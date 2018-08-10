#' Perform KEEG Pathways Enrichment
#'
#' @param Tool
#' @param ID
#' @param Width,Height,Res,Unit,image.format
#' @param env
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#' @inheritParams GOnto
#'
#' @return Enriched terms.
#' @export
#'
#' @import pathview
#'
#' @examples
#' \dontrun{
#' KEEG_ENRICH(Tool = "edgeR", env = "env name without quotes")
#' }
KEEG_ENRICH <- function(Tool = "edgeR",
                        ID = "GeneID",
                        Width = 8,
                        Height = 4,
                        Res = 300,
                        Unit = "in",
                        image.format = "png",
                        env){

    # library(pathview, gage, gageData)

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (missing("workDir")){
        workDir <- string_vars[["envir_link"]]$workDir
    }

    # assign("PATH", file.path(workDir, "GDCRtools", toupper(string_vars[["envir_link"]]$tumor), "Analyses"),
    #        envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    dir.create(path = paste0(PATH, "/KEGG_GAGE_", tolower(FinalData)), showWarnings = FALSE)
    dir.create(path = paste0(PATH, "/KEGG_GAGE_", tolower(FinalData), "/aux_files"), showWarnings = FALSE)
    DIR <- paste0(PATH, "/KEGG_GAGE_", tolower(FinalData))
    #http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html

    #kegg.sets.hs list of 229 elements with gene Entrez IDs for a single KEGG pathway
    # e.g. “Global Map” and “Human Diseases”
    kegg.sets.hs <- GDCRtools::kegg.sets.hs
    #for only sinaling and metabolic pathways
    sigmet.idx.hs <- GDCRtools::sigmet.idx.hs
    # gene sets of sinaling and metabolic pathways only
    kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
    # head(kegg.sets.hs, 3)

    File <- paste("resultadosDE", Tool, sep = ".")
    # File2 <- paste("Results.Completed", Tool, sep = ".")

    #input DE genes
    resultadosDE <- get(File, envir = string_vars[["envir_link"]])[[pairName]]
    # Results.Completed <- get(File, envir = string_vars[["envir_link"]])[[pairName]]

    # generate annotation table
    annotation_table <- GDCRtools::annotation_table

    if (tolower(dataBase) == "gdc") {
        rownames(annotation_table) <- annotation_table$ensembl
        resultadosDE$GeneID <- annotation_table[resultadosDE$ensembl, "GeneID"]
    }

    foldchanges <- resultadosDE[, "log2FC", drop = FALSE]
    rownames(foldchanges) <- resultadosDE$GeneID
    # names(foldchanges) <- resultadosDE$GeneID

    keggres <- gage::gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)

    # Get the pathways
    # library(magrittr)
    keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater)
    keggrespathways <- na.exclude(keggrespathways)
    keggrespathways <- as.character(keggrespathways$id)

    # Get the IDs.
    keggresids <- unname(sapply(keggrespathways, function(w){
        paste0(unlist(strsplit(w, " "))[1])}))

    # get user wdir
    wdir <- getwd()
    setwd(DIR)

    DIR2 <- file.path(DIR, "aux_files")

    # plot multiple pathways (plots saved to disk and returns a throwaway list object)
    tmp <- sapply(keggresids, function(pid) pathview::pathview(gene.data = foldchanges,
                                                               pathway.id = pid,
                                                               species = "hsa",
                                                               kegg.dir = DIR2))
    # return to user wdir
    setwd(wdir)

    #see Plotting fold changes in genomic space in http://www.bioconductor.org/help/workflows/rnaseqGene/
    # #for ids
    # File$name <-    mapIds(org.Hs.eg.db,
    #                     keys=row.names(File),
    #                     column="GENENAME",
    #                     column="ENTREZID",
    #                     column="SYMBOL",
    #                     keytype="ENSEMBL",
    #                     multiVals="first")

    gc()
    message("Done!\n")
}
