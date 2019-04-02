#' Gene Set Enrichment Analysis
#'
#' @param FDR_cutoff
#' @param Width,Height,Res,Unit,image_format
#' @param Tool
#' @param ID
#' @param pairName
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
#' GSEA(Tool = "edgeR", env = "env name without quotes")
#' }
GSEA <- function(FDR_cutoff = 0.05,
                 Width = 10,
                 Height = 3,
                 Res = 500,
                 Unit = "in",
                 image_format = "png",
                 Tool = "edgeR",
                 ID = "GeneID",
                 pairName = "G2_over_G1",
                 env){

    message("Performing Gene Set Enrichment Analysis - over Reactome...\n")
    # Looks for enrichment by fold direction
    #List of input id types

    # BiocParallel::register(BiocParallel::SerialParam())

    BiocParallel::register(BiocParallel::SnowParam(1))

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    if (missing(Tool)){Tool <- string_vars[["envir_link"]]$Tool}

    Name <- string_vars[["envir_link"]]$Name
    # dataBase <- string_vars[["envir_link"]]$dataBase
    groupGen <- string_vars[["envir_link"]]$groupGen

    if (grepl("crosstable", tolower(Tool))) {
        TCGAExpression <- string_vars[["envir_link"]]$Results_Completed_crossed
    } else {
        TCGAExpression <- eval(parse(text= paste0("string_vars[['envir_link']]$Results_Completed.",
                                                  Tool)))
    }

    TCGAExpression <- TCGAExpression[[pairName]]

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
    } else {
        dir.create(paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_",
                          Tool, "_", toupper(Name)), showWarnings = FALSE)
        DIR <- paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_", Tool, "_", toupper(Name))
    }

    dir.create(file.path(DIR, "GSEA_Output"), showWarnings = FALSE)

    suppressPackageStartupMessages(input.types <- as.matrix(x = AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)))
    write.table(input.types, paste0(DIR,
                                  "/GSEA_Output/input_types",
                                  ID, "_", pairName, ".txt"), row.names = FALSE)

    # GeneSetVector <- string_vars[["envir_link"]]$TCGAExpression$FC

    GeneSetVector <- TCGAExpression$FC

    if (ID == "GeneID"){
        # It demands FOld change in one vector decreasing ordered
        names(GeneSetVector) <- TCGAExpression$GeneID
        GeneSetVector <- GeneSetVector[order(GeneSetVector, decreasing = TRUE)]

    } else if(ID == "GeneSymbol"){
        # It demands FOld change in one vector decreasing ordered
        names(GeneSetVector) <- TCGAExpression$GeneSymbol
        GeneSetVector <- GeneSetVector[order(GeneSetVector, decreasing = TRUE)]

        # Again, change the names for ID using HUGO names
        ids <- clusterProfiler::bitr(names(GeneSetVector), fromType = "SYMBOL",
                                     toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
        GeneSetVector <- GeneSetVector[ids$ALIAS]
        names(GeneSetVector) <- ids$ENTREZID

    } else if (ID == "Ensembl"){
        names(GeneSetVector) <- rownames(TCGAExpression)
        GeneSetVector <- GeneSetVector[order(GeneSetVector, decreasing = TRUE)]

        # Again, change the names for ID using HUGO names
        ids <- clusterProfiler::bitr(names(GeneSetVector), fromType = 'ENSEMBL',
                                                      toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
        GeneSetVector <- GeneSetVector[ids$ENSEMBL]
        names(GeneSetVector) <- ids$ENTREZID
    }  else if (tolower(ID) == "refgene") {
        names(GeneSetVector) <- TCGAExpression$Ensembl
        GeneSetVector <- GeneSetVector[order(GeneSetVector, decreasing = TRUE)]

        # Again, change the names for ID using HUGO names
        ids <- clusterProfiler::bitr(names(GeneSetVector), fromType = 'REFSEQ',
                                     toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
        GeneSetVector <- GeneSetVector[ids$ENSEMBL]
        names(GeneSetVector) <- ids$ENTREZID
    }

    # Perform GSEA
    GSEA_REACTOME <- ReactomePA::gsePathway(GeneSetVector,
                                pvalueCutoff = FDR_cutoff, nPerm = 10000,
                                pAdjustMethod = "BH", verbose = TRUE)

    #Get summary
    GSEA_REACTOME_summary_fdr <- as.data.frame(GSEA_REACTOME)

    # Write
    write.csv(GSEA_REACTOME_summary_fdr, file = paste0(DIR,
                                                   "/GSEA_Output/REAC.GSEA.enrichment",
                                                   ID, "_", pairName, ".csv"), row.names = FALSE)
    assign("GeneSetVector", GeneSetVector, envir = get(envir_link))
    # GSEA_REACTOME_summary_fdr <- GSEA_REACTOME.summary[which(GSEA_REACTOME.summary$p.adjust < FDR_cutoff), ]


    # REACTOME-GSEA Plot
    if (nrow(GSEA_REACTOME_summary_fdr) != 0) {

        if (nrow(GSEA_REACTOME_summary_fdr) >= 15){
            PlotNowREACT <- GSEA_REACTOME_summary_fdr[1:15, ]
        } else{
            PlotNowREACT <- GSEA_REACTOME_summary_fdr
        }

        PlotNowREACT[, 2] <- gsub(" of", "", PlotNowREACT[, 2])
        # large <- lengths(strsplit(PlotNowREACT[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowREACT[, 2], " ")) > 7
        PlotNowREACT[large, 2] <- unname(sapply(PlotNowREACT[large, 2],
                                             function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                               collapse = " ")}))

        # is it better to use the manually adjusted?
        log_10_GSEA <- -log10(as.numeric(PlotNowREACT[, "p.adjust"]))
        new_inf <- log_10_GSEA[order(log_10_GSEA, decreasing = TRUE)]
        log_10_GSEA[log_10_GSEA == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

        longest_word <- max(stringr::str_count(PlotNowREACT$Description))


        Gene_Ratio <- unname(sapply(PlotNowREACT$leading_edge,
                                    function(w){unlist(strsplit(w, ", "))[1]}))
        Gene_Ratio <- round(as.numeric(gsub("tags=|%", "", Gene_Ratio))/100, 2)


        if (longest_word > 50) {
            longest_Width <- Width * 1.5
        } else {
            longest_Width <- Width
        }

        if (tolower(image_format) == "png") {
            png(filename = paste0(DIR,
                              "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", ID,
                              "_", pairName, ".png"),
                width = longest_Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = paste0(DIR,
                              "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", ID,
                              "_", pairName, ".svg"),
                width = longest_Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }

        Count <- ceiling(PlotNowREACT$setSize * Gene_Ratio)


        p <- ggplot2::ggplot(PlotNowREACT, ggplot2::aes(x = log_10_GSEA,
                                                     y = forcats::fct_reorder(Description, log_10_GSEA))) +
            ggplot2::geom_point(ggplot2::aes(size = Count, color = Gene_Ratio)) +
            ggplot2::scale_colour_gradient(limits=c(0, 1), low="red", high = "blue") +
            ggplot2::labs(y = "", x = "-log(FDR)") +
            ggplot2::theme_bw(base_size = 10) +
            ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold",
                                                                size = 16),
                           axis.text = ggplot2::element_text(face = "bold",
                                                             color = "#011600", size = 12),
                           title = ggplot2::element_text(face = "bold",
                                                         size = 18),
                           plot.title = ggplot2::element_text(hjust = 0.5))
        print(p)
        dev.off()

    } else {
        message(cat("There is nothing to show in Gene Set Enrichment Analysis - over Reactome"))
    }
    gc()
    message("Done!\n")
}
