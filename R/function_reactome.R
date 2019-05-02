#' DESEASE-ONTOLOGY and REACTOME ENRICHMENT
#'
#' @param p_cutoff
#' @param FDR_cutoff
#' @param Width,Height,Res,Unit,image_format
#' @param Tool
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
#' DO_React_enrich(Tool = "edgeR", env = "env name without quotes")
#' }
DO_React_enrich <- function(p_cutoff = 0.05,
                           FDR_cutoff = 0.05,
                           Width = 8,
                           Height = 4,
                           Res = 300,
                           Unit = "in",
                           image_format = "png",
                           Tool = "edgeR",
                           pairName = "G2_over_G1",
                           env){

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
    # pairName <- string_vars[["envir_link"]]$pairName
    dataBase <- string_vars[["envir_link"]]$dataBase
    groupGen <- string_vars[["envir_link"]]$groupGen

    if (grepl("crosstable", tolower(Tool))) {
        TCGAExpression <- string_vars[["envir_link"]]$Results_Completed_crossed
    } else {
        TCGAExpression <- eval(parse(text= paste0("string_vars[['envir_link']]$Results_Completed.",
                                                  Tool)))
    }
    TCGAExpression <- TCGAExpression[[pairName]]

    if (grepl("crosstable", tolower(Tool))) {
        DEGenes_all <- string_vars[["envir_link"]]$resultadosDE_crossed[[pairName]]
    } else {
        DEGenes_all <- get(paste("resultadosDE", Tool, sep = "."),
                           envir = string_vars[["envir_link"]])[[pairName]]
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
    } else {
        dir.create(paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_",
                          Tool, "_", toupper(Name)), showWarnings = FALSE)
        DIR <- paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_", Tool, "_", toupper(Name))
    }

    dir.create(file.path(DIR, "REACTOME_Output"), showWarnings = FALSE)
    # dir.create(file.path(DIR, "REACTOME_Output", "PathwayMaps"),
    #            showWarnings = FALSE)
    dir.create(file.path(DIR, "DO_Output"), showWarnings = FALSE)

    # "GeneID" equal to "entrez gene id"
    if (tolower(dataBase) == "gdc"){
        geneIDList_DE <- clusterProfiler::bitr(DEGenes_all$ensembl, fromType = "ENSEMBL",
                                               toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
        geneIDList_Universe <- clusterProfiler::bitr(TCGAExpression$ensembl, fromType="ENSEMBL",
                                                     toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
    } else {
        geneIDList_DE <- DEGenes_all$GeneID
        geneIDList_Universe <- TCGAExpression$GeneID
    }

    # DO ####

    # Perform DO Enrichment analysis
    message("\nPerforming DO enrichment...\n")
    DO_Enrichment <- DOSE::enrichDO(gene = geneIDList_DE,
                              ont = "DO",
                              pvalueCutoff  = FDR_cutoff,
                              pAdjustMethod = "BH",
                              universe = geneIDList_Universe,
                              minGSSize = 2, #at least a pair
                              qvalueCutoff = 1,
                              readable = TRUE)

    # Makes summary
    DO_Enriched_summary <- as.data.frame(DO_Enrichment)

    # Write
    write.csv(DO_Enriched_summary,
              file=paste0(DIR,
                          "/DO_Output/DO_enrichment_", pairName, ".csv"),
              row.names = FALSE)

    #Get fdr <0.05
    DO_Enriched_fdr <- DO_Enriched_summary[which(DO_Enriched_summary$p.adjust < FDR_cutoff), ]

    # DO Plot
    message("Plotting DO enrichment...\n")

    if (nrow(DO_Enriched_fdr) != 0){

        if (nrow(DO_Enriched_fdr) >= 10){
            PlotNowDO <- DO_Enriched_fdr[1:10, ]
        } else{
            PlotNowDO <- DO_Enriched_fdr
        }

        PlotNowDO[, 2] <- gsub(" of", "", PlotNowDO[, 2])
        # large <- lengths(strsplit(PlotNowREACT[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowDO[, 2], " ")) > 7
        PlotNowDO[large, 2] <- unname(sapply(PlotNowDO[large, 2],
                                                function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                                  collapse = " ")}))

        log_10_DO <- -log10(as.numeric(PlotNowDO[, "p.adjust"]))
        new_inf <- log_10_DO[order(log_10_DO, decreasing = TRUE)]
        log_10_DO[log_10_DO == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

        longest_word <- max(stringr::str_count(PlotNowDO$Description))

        Gene_Ratio1 <- unname(sapply(PlotNowDO$GeneRatio,
                                    function(w){unlist(strsplit(w, "/"))[1]}))
        Gene_Ratio2 <- unname(sapply(PlotNowDO$GeneRatio,
                                     function(w){unlist(strsplit(w, "/"))[2]}))

        Gene_Ratio <- round(as.numeric(Gene_Ratio1)/as.numeric(Gene_Ratio2), 2)

        if (longest_word > 50) {
            longest_Width <- Width * 1.5
        } else {
            longest_Width <- Width
        }

        if (tolower(image_format) == "png") {
            png(file = paste0(DIR,
                              "/DO_Output/DO_EnrichPlot_10first_", pairName, ".png"),
                width = longest_Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(file = paste0(DIR,
                              "/DO_Output/DO_EnrichPlot_10first_", pairName, ".svg"),
                width = longest_Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }

        p <- ggplot2::ggplot(PlotNowDO, ggplot2::aes(x = log_10_DO,
                                                     y = forcats::fct_reorder(Description, log_10_DO))) +
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
         message(cat("There is nothing to show in Disease Ontology\n"))
    }

    # REACTOME ####
    # perform enrichment in REACTOME
    message("Performing reactome enrichment...\n")
    REACTOME_Enriched <- ReactomePA::enrichPathway(geneIDList_DE,
                                       pvalueCutoff = 1,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 1,
                                       universe = geneIDList_Universe,
                                       minGSSize = 2,
                                       readable = TRUE)
    # Makes summary
    REACTOME_Enriched_summary <- as.data.frame(REACTOME_Enriched)

    # Write
    write.csv(REACTOME_Enriched_summary,
              file=paste0(DIR,
                         "/REACTOME_Output/REAC_enrichment_", pairName, ".csv"), row.names = FALSE)

    #Get fdr <0.05
    REACTOME_Enriched_fdr <- REACTOME_Enriched_summary[which(REACTOME_Enriched_summary$p.adjust < FDR_cutoff), ]
    REACTOME_Enriched_fdr <- REACTOME_Enriched_fdr[order(REACTOME_Enriched_fdr$p.adjust, decreasing = FALSE), ]

    # REACTOME Plot
    message("Plotting Reactome enrichment\n")
    if (nrow(REACTOME_Enriched_fdr) != 0) {

        if (nrow(REACTOME_Enriched_fdr) >= 10){
            PlotNowREACT <- REACTOME_Enriched_fdr[1:10, ]
        } else{
            PlotNowREACT <- REACTOME_Enriched_fdr
        }

        PlotNowREACT[, 2] <- gsub(" of", "", PlotNowREACT[, 2])
        # large <- lengths(strsplit(PlotNowREACT[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowREACT[, 2], " ")) > 7
        PlotNowREACT[large, 2] <- unname(sapply(PlotNowREACT[large, 2],
                                                function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                                  collapse = " ")}))

        log_10_REACT <- -log10(as.numeric(PlotNowREACT[, "p.adjust"]))
        new_inf <- log_10_REACT[order(log_10_REACT, decreasing = TRUE)]
        log_10_REACT[log_10_REACT == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

        longest_word <- max(stringr::str_count(PlotNowREACT$Description))

        Gene_Ratio1 <- unname(sapply(PlotNowREACT$GeneRatio,
                                     function(w){unlist(strsplit(w, "/"))[1]}))
        Gene_Ratio2 <- unname(sapply(PlotNowREACT$GeneRatio,
                                     function(w){unlist(strsplit(w, "/"))[2]}))

        Gene_Ratio <- round(as.numeric(Gene_Ratio1)/as.numeric(Gene_Ratio2), 2)

        if (longest_word > 50) {
            longest_Width <- Width * 1.5
        } else {
            longest_Width <- Width
        }

        if (tolower(image_format) == "png") {
            png(file = paste0(DIR,
                              "/REACTOME_Output/REACT_EnrichPlot_10first_", pairName, ".png"),
                width = longest_Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(file = paste0(DIR,
                              "/REACTOME_Output/REACT_EnrichPlot_10first_", pairName, ".svg"),
                width = longest_Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }

        p <- ggplot2::ggplot(PlotNowREACT, ggplot2::aes(x = log_10_REACT,
                                                     y = forcats::fct_reorder(Description, log_10_REACT))) +
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
        message(cat("There is nothing to show from Reactome Enrichment"))
    }

    # Write
    write.csv(REACTOME_Enriched_summary,
              file=paste0(DIR,
                          "/REACTOME_Output/REAC_enrichment_under_cutoff_", pairName, ".csv"),
              row.names = FALSE)

    gc()
    message("Done!\n")
}
