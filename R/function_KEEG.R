#' Perform KEEG Pathways Enrichment
#'
#' @param Tool
#' @param FDR_cutoff
#' @param ID
#' @param pairName
#' @param Width,Height,Res,Unit,image_format
#' @param env
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#' @inheritParams GOnto
#'
#' @return Enriched terms.
#' @export
#'
#' @examples
#' \dontrun{
#' KEEG_ENRICH(Tool = "edgeR", env = "env name without quotes")
#' }
KEEG_ENRICH <- function(Tool = "edgeR",
                        FDR_cutoff = 0.05,
                        ID = "GeneID",
                        pairName = "G2_over_G1",
                        Width = 8,
                        Height = 4,
                        Res = 300,
                        Unit = "in",
                        image_format = "png",
                        env){

    # library(pathview, gage, gageData)

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))
    groupGen <- string_vars[["envir_link"]]$groupGen
    Name <- string_vars[["envir_link"]]$Name

    # assign("PATH", file.path(workDir, "GDCtools", toupper(string_vars[["envir_link"]]$tumor), "Analyses"),
    #        envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    if (grepl("crosstable", tolower(Tool))) {
        if (tolower(Tool) == "crosstable.deseq2") {
            DIR <- paste0(PATH, "/CrossData_deseq2")
        } else if (tolower(Tool) == "crosstable.edger") {
            DIR <- paste0(PATH, "/CrossData_edger")
        } else if (tolower(Tool) == "crosstable.ebseq") {
            DIR <- paste0(PATH, "/CrossData_ebseq")
        }
        dir.create(file.path(DIR, paste0("Ontology_Results", tolower(groupGen))), showWarnings = FALSE)
        DIR <- file.path(DIR, paste0("Ontology_Results", tolower(groupGen)))

        File <- "resultadosDE_crossed"

    } else {
        File <- paste("resultadosDE", Tool, sep = ".")

        dir.create(paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_",
                          Tool, "_", toupper(Name)), showWarnings = FALSE)
        DIR <- paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_", Tool, "_", toupper(Name))

    }

    DIR <- paste0(DIR, "/KEGG_GAGE_", tolower(Tool))

    dir.create(DIR, showWarnings = FALSE)
    dir.create(file.path(DIR, "aux_files"), showWarnings = FALSE)


    #http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html

    #kegg.sets.hs list of 229 elements with gene Entrez IDs for a single KEGG pathway
    # e.g. “Global Map” and “Human Diseases”
    if (!exists("kegg.sets.hs", envir = get(envir_link))) {
        data(kegg.sets.hs, package = "gageData", envir = get(envir_link))
        data(sigmet.idx.hs, package = "gageData", envir = get(envir_link))
    }
    kegg.sets.hs <- string_vars[["envir_link"]]$kegg.sets.hs
    #for only sinaling and metabolic pathways
    sigmet.idx.hs <- string_vars[["envir_link"]]$sigmet.idx.hs
    # gene sets of sinaling and metabolic pathways only
    kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
    # head(kegg.sets.hs, 3)

    #input DE genes
    resultadosDE <- get(File, envir = string_vars[["envir_link"]])[[pairName]]
    # Results_Completed <- get(File, envir = string_vars[["envir_link"]])[[pairName]]

    # generate annotation table
    annotation_table <- GDCtools::annotation_table

    dataBase <- string_vars[["envir_link"]]$dataBase

    if (tolower(dataBase) == "gdc") {
        rownames(annotation_table) <- annotation_table$ensembl
        resultadosDE$GeneID <- annotation_table[resultadosDE$ensembl, "GeneID"]
    }

    # foldchanges <- resultadosDE[, "log2FC", drop = FALSE]
    # rownames(foldchanges) <- resultadosDE$GeneID
    foldchanges <- resultadosDE[, "log2FC", drop = TRUE]
    names(foldchanges) <- resultadosDE$GeneID

    keggres <- gage::gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)

    # Get the pathways
    # library(magrittr)
    keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater)
    keggrespathways <- na.exclude(keggrespathways)
    keggrespathways$id<- as.character(keggrespathways$id)
    keggrespathways <- keggrespathways[keggrespathways$q.val < FDR_cutoff, ]

    if (nrow(keggrespathways) == 0) {
        stop(message("There is not any kegg pathway with FDR < ", FDR_cutoff, "\n"))
    }

    # Get the IDs.
    keggresids <- unname(sapply(keggrespathways, function(w){
        paste0(unlist(strsplit(w, " "))[1])}))

    # get user wdir
    wdir <- getwd()
    setwd(DIR)

    DIR2 <- file.path(DIR, "aux_files")

    data(bods, package = "pathview")
    data(korg, package = "pathview")

    # importing needed funtions
    columns <- AnnotationDbi::columns
    select <- AnnotationDbi::select

    # plot multiple pathways (plots saved to disk and returns a throwaway list object)
    tmp <- sapply(keggresids, function(pid) pathview::pathview(gene.data = foldchanges,
                                                               pathway.id = pid,
                                                               species = "hsa",
                                                               kegg.dir = DIR2))
    remove(bods, korg, envir = .GlobalEnv)

    # return to user wdir
    setwd(wdir)


    gc()
    message("Done!\n")
}
