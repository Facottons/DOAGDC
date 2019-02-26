#' DESEASE-ONTOLOGY and REACTOME ENRICHMENT
#'
#' @param p.cutoff
#' @param FDR.cutoff
#' @param Width,Height,Res,Unit,image.format
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
#' DO_REAC_ENRICH (Tool = "edgeR", env = "env name without quotes")
#' }
DO_REAC_ENRICH <- function(p.cutoff = 0.05,
                           FDR.cutoff = 0.05,
                           Width = 8,
                           Height = 4,
                           Res = 300,
                           Unit = "in",
                           image.format = "png",
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
        TCGAExpression <- string_vars[["envir_link"]]$Results.Completed.crossed
    } else {
        TCGAExpression <- eval(parse(text= paste0("string_vars[['envir_link']]$Results.Completed.",
                                                  Tool)))
    }
    TCGAExpression <- TCGAExpression[[pairName]]

    if (grepl("crosstable", tolower(Tool))) {
        DEGenes.all <- string_vars[["envir_link"]]$resultadosDE.crossed[[pairName]]
    } else {
        DEGenes.all <- get(paste("resultadosDE", Tool, sep = "."),
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
        geneIDList.DE <- clusterProfiler::bitr(DEGenes.all$ensembl, fromType = "ENSEMBL",
                                               toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
        geneIDList.Universe <- clusterProfiler::bitr(TCGAExpression$ensembl, fromType="ENSEMBL",
                                                     toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
    } else {
        geneIDList.DE <- DEGenes.all$GeneID
        geneIDList.Universe <- TCGAExpression$GeneID
    }

    # DO ####

    # Perform DO Enrichment analysis
    message("\nPerforming DO enrichment...\n")
    DO.Enrichment <- DOSE::enrichDO(gene = geneIDList.DE,
                              ont = "DO",
                              pvalueCutoff  = FDR.cutoff,
                              pAdjustMethod = "BH",
                              universe = geneIDList.Universe,
                              minGSSize = 2, #at least a pair
                              qvalueCutoff = 1,
                              readable = TRUE)

    # Makes summary
    DO.Enriched.summary <- as.data.frame(DO.Enrichment)

    # Write
    write.csv(DO.Enriched.summary,
              file=paste0(DIR,
                          "/DO_Output/DO_enrichment_", pairName, ".csv"),
              row.names = FALSE)

    #Get fdr <0.05
    DO.Enriched.fdr <- DO.Enriched.summary[which(DO.Enriched.summary$p.adjust < FDR.cutoff), ]

    # DO Plot
    message("Plotting DO enrichment...\n")

    if (nrow(DO.Enriched.fdr) != 0){

        if (nrow(DO.Enriched.fdr) >= 10){
            PlotNowDO <- DO.Enriched.fdr[1:10, ]
        } else{
            PlotNowDO <- DO.Enriched.fdr
        }

        PlotNowDO[, 2] <- gsub(" of", "", PlotNowDO[, 2])
        # large <- lengths(strsplit(PlotNowREACT[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowDO[, 2], " ")) > 7
        PlotNowDO[large, 2] <- unname(sapply(PlotNowDO[large, 2],
                                                function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                                  collapse = " ")}))

        log.10.DO <- -log10(as.numeric(PlotNowDO[, "p.adjust"]))
        new.inf <- log.10.DO[order(log.10.DO, decreasing = TRUE)]
        log.10.DO[log.10.DO == "Inf"] <- (new.inf[new.inf != "Inf"][1]+1)

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

        if (tolower(image.format) == "png") {
            png(file = paste0(DIR,
                              "/DO_Output/DO_EnrichPlot_10first_", pairName, ".png"),
                width = longest_Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(file = paste0(DIR,
                              "/DO_Output/DO_EnrichPlot_10first_", pairName, ".svg"),
                width = longest_Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }

        p <- ggplot2::ggplot(PlotNowDO, ggplot2::aes(x = log.10.DO,
                                                     y = forcats::fct_reorder(Description, log.10.DO))) +
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

        # if (tolower(image.format) == "png") {
        #     png(file = paste0(DIR,
        #                       "/DO_Output/DO_EnrichPlot_10first_", pairName, "_2.png"),
        #         width = longest_Width, height = Height, res = Res, units = Unit)
        # } else if (tolower(image.format) == "svg") {
        #     svg(file = paste0(DIR,
        #                       "/DO_Output/DO_EnrichPlot_10first_", pairName, "_2.svg"),
        #         width = longest_Width, height = Height, onefile = TRUE)
        # } else {
        #     stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        # }
        # # par(oma = c(0,10,0,0),lwd = 2.5)
        # par(mar = c(4.5,18,2,2),lwd = 2.5)
        #
        # barplot(log.10.DO,
        #         names=PlotNowDO[,"Description"], las=2, horiz=TRUE, xlab = "-log(FDR)",
        #         main="Disease Ontology", lwd=2, cex.lab=1, cex.axis=2, cex.main=1.4,
        #         axes=FALSE, col = RColorBrewer::brewer.pal(8,"Set1")[5],
        #         border="white",
        #         cex.names=0.8, space=0.001)
        #
        # if (-log10(min(as.numeric(PlotNowDO[, "p.adjust"]))) == "Inf"){
        #     abline(v = 1:floor(0), col = "white", lwd = 4)
        # } else {
        #     abline(v=1:floor(-log10(min(as.numeric(PlotNowDO[, "p.adjust"]), na.rm = TRUE))),
        #            col="white", lwd=4)
        # }
        #
        # axis(side = 1, 0:(max(log.10.DO)+2),
        #      lwd = 1.5, cex.axis = 1.2)
        #
        # abline(v = -log10(FDR.cutoff), lwd = 3)
        #
        # par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 1, 0), new=TRUE)
        # plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        #
        # legend("topright", "-log(FDR.cutoff)", lty = 1, lwd = 4,
        #        bty = "n", cex = 0.5)
        # dev.off()

    } else {
         message(cat("There is nothing to show in Disease Ontology\n"))
    }

    # REACTOME ####
    # perform enrichment in REACTOME
    message("Performing reactome enrichment...\n")
    REACTOME.Enriched <- ReactomePA::enrichPathway(geneIDList.DE,
                                       pvalueCutoff = 1,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 1,
                                       universe = geneIDList.Universe,
                                       minGSSize = 2,
                                       readable = TRUE)
    # Makes summary
    REACTOME.Enriched.summary <- as.data.frame(REACTOME.Enriched)

    # Write
    write.csv(REACTOME.Enriched.summary,
              file=paste0(DIR,
                         "/REACTOME_Output/REAC_enrichment_", pairName, ".csv"), row.names = FALSE)

    #Get fdr <0.05
    REACTOME.Enriched.fdr <- REACTOME.Enriched.summary[which(REACTOME.Enriched.summary$p.adjust < FDR.cutoff), ]
    REACTOME.Enriched.fdr <- REACTOME.Enriched.fdr[order(REACTOME.Enriched.fdr$p.adjust, decreasing = FALSE), ]

    # REACTOME Plot
    message("Plotting Reactome enrichment\n")
    if (nrow(REACTOME.Enriched.fdr) != 0) {

        if (nrow(REACTOME.Enriched.fdr) >= 10){
            PlotNowREACT <- REACTOME.Enriched.fdr[1:10, ]
        } else{
            PlotNowREACT <- REACTOME.Enriched.fdr
        }

        PlotNowREACT[, 2] <- gsub(" of", "", PlotNowREACT[, 2])
        # large <- lengths(strsplit(PlotNowREACT[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowREACT[, 2], " ")) > 7
        PlotNowREACT[large, 2] <- unname(sapply(PlotNowREACT[large, 2],
                                                function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                                  collapse = " ")}))

        log.10.REACT <- -log10(as.numeric(PlotNowREACT[, "p.adjust"]))
        new.inf <- log.10.REACT[order(log.10.REACT, decreasing = TRUE)]
        log.10.REACT[log.10.REACT == "Inf"] <- (new.inf[new.inf != "Inf"][1]+1)

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

        if (tolower(image.format) == "png") {
            png(file = paste0(DIR,
                              "/REACTOME_Output/REACT_EnrichPlot_10first_", pairName, ".png"),
                width = longest_Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(file = paste0(DIR,
                              "/REACTOME_Output/REACT_EnrichPlot_10first_", pairName, ".svg"),
                width = longest_Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }

        p <- ggplot2::ggplot(PlotNowREACT, ggplot2::aes(x = log.10.REACT,
                                                     y = forcats::fct_reorder(Description, log.10.REACT))) +
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

        # if (tolower(image.format) == "png") {
        #     png(file = paste0(DIR,
        #                       "/REACTOME_Output/REACT_EnrichPlot_10first_", pairName, "_2.png"),
        #         width = longest_Width, height = Height, res = Res, units = Unit)
        # } else if (tolower(image.format) == "svg") {
        #     svg(file = paste0(DIR,
        #                       "/REACTOME_Output/REACT_EnrichPlot_10first_", pairName, "_2.svg"),
        #         width = longest_Width, height = Height, onefile = TRUE)
        # } else {
        #     stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        # }
        # # par(oma = c(0,10,0,0),lwd = 2.5)
        # par(mar = c(4.5,20,2,2),lwd = 2.5)
        # barplot(log.10.REACT,
        #         names = PlotNowREACT[, "Description"], las = 2, horiz = TRUE,
        #         xlab = "-log(FDR)", main = "REACTOME", lwd = 2, cex.lab = 1.1,
        #         cex.axis = 1, cex.main = 1.5, axes = FALSE,
        #         col = RColorBrewer::brewer.pal(8, "Set1")[4], border = "white",
        #         cex.names = 0.8, space = 0.001)
        #
        # if (-log10(min(as.numeric(PlotNowREACT[, "p.adjust"]))) == "Inf"){
        #     abline(v = 1:floor(0), col = "white", lwd = 4)
        # } else {
        #     abline(v = 1:floor(-log10(min(PlotNowREACT[, "p.adjust"]))),
        #            col = "white", lwd = 4)
        # }
        #
        # axis(side = 1, lwd = 1.5, cex.axis = 1.1,
        #      0:(max(log.10.REACT)+2))
        #
        # abline(v = -log10(FDR.cutoff), lwd = 3)
        #
        # par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 1, 0), new=TRUE)
        # plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        #
        # legend("topright", "-log(FDR.cutoff)", lty = 1, lwd = 4,
        #        bty = "n", cex = 0.5)
        #
        # dev.off()

        # # Plot Specific Pathways with fold-changes
        # for(kiwi in 1:length(REACTOME.Enriched.fdr[, "Description"])){
        #     x <- gsub(" ", "_", REACTOME.Enriched.fdr[kiwi, "Description"])
        #     x <- gsub("[(]", "_", x)
        #     x <- gsub("[)]", "_", x)
        #     x <- gsub("[/]", "_", x)
        #     x <- gsub("[']", "", x)
        #     x <- gsub("[+]", "plus", x)
        #     x <- gsub("[,]", "_", x)
        #     if (kiwi != 3) {
        #         png(file = paste0("./Ontology_Results_", Tool, "_",
        #                           toupper(Name), "/REACTOME_Output/PathwayMaps/React_",
        #                           x, "_", ID, ".png"),
        #             width = Width, heigh = 4, units = Unit, res = Res)
        #         par(mar = c(2,2,2,6),lwd = 2.5)
        #         tryCatch(viewPathway(pathName = REACTOME.Enriched.fdr[kiwi, "Description"],
        #                     readable = TRUE, organism = "human", vertex.label.font = 1,
        #                     foldChange = GeneSetVector[unique(names(GeneSetVector))]),
        #                  error = function(e) e)
        #         dev.off()
        #     }
        # }
    } else {
        message(cat("There is nothing to show from Reactome Enrichment"))
    }

    # Write
    write.csv(REACTOME.Enriched.summary,
              file=paste0(DIR,
                          "/REACTOME_Output/REAC_enrichment_under_cutoff_", pairName, ".csv"),
              row.names = FALSE)

    #
    #
    # Why Viewpathway does not work with:
    #   "Metabolism of steroid hormones and vitamin D"
    #       REACTOME.Enriched.fdr[3, "Description"]
    #           from TCGT
    #           Neurotoxicity of clostridium toxins mais esse[16, "Description"]
    #
    #

    # perform enrichment in REACTOME with cutoffs
    # REACTOME.Enriched.2 <- ReactomePA::enrichPathway(geneIDList.DE,
    #                                      organism = "human",
    #                                      pvalueCutoff = p.cutoff,
    #                                      pAdjustMethod = "BH",
    #                                      qvalueCutoff = 1,
    #                                      universe = geneIDList.Universe,
    #                                      minGSSize = 2,
    #                                      readable = TRUE)

    # # Plot Enrichment Map With Names
    # png(file = paste0(DIR,
    #                   "/REACTOME_Output/REACT_EnrichMap2_", ID, ".png"),
    #     #width = 6000, heigh = 6000, res = Res)
    #     width = Width, heigh = 12, units = Unit, res = Res)
    # par(mar = c(2,8,2,8), lwd = 2.5)
    # suppressWarnings(enrichMap(x = REACTOME.Enriched.2,
    #           layout = layout.kamada.kawai,
    #           vertex.label.cex = 0.5,
    #           vertex.label.font = 0.5,
    #           #vertex.label.color = "#666666",
    #           foldChange = NULL, fixed = FALSE,
    #           col.bin = 10))
    # dev.off()

    #
    # # Plot Enrichment Map With Names on graph plots
    # png(file = paste0(DIR,
    #                   "/REACTOME_Output/REACT_EnrichMap_", ID, ".png"),
    #     #width = 6000, heigh = 6000, res = Res)
    #     width = Width, heigh = 6, units = Unit, res = Res)
    # par(mar = c(2,8,2,8), lwd = 2.5)
    # suppressWarnings(enrichMap(x = REACTOME.Enriched.2,
    #                            layout = layout.kamada.kawai,
    #                            vertex.label.cex = 0.5,
    #                            vertex.label.font = 0.5,
    #                            #vertex.label.color = "#666666",
    #                            foldChange = NULL,
    #                            col.bin = 10))
    # dev.off()
    #
    # # Plot the cnetplot for the First Cat
    # png(file = paste0("./Ontology_Results_", Tool, "_",
    #                   toupper(Name), "/REACTOME_Output/REACT_cnetplot_", ID, ".png"),
    #     width = Width, heigh = 8, units = Unit, res = Res)
    # suppressWarnings(cnetplot(REACTOME.Enriched, categorySize = "Count",
    #                           foldChange = NULL,
    #          showCategory = 1))
    # dev.off()

    gc()
    message("Done!\n")
}
