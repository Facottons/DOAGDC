#' Gene Set Enrichment Analysis
#'
#' @param FDR.cutoff
#' @param Width,Height,Res,Unit,image.format
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
GSEA <- function(FDR.cutoff = 0.05,
                 Width = 10,
                 Height = 3,
                 Res = 500,
                 Unit = "in",
                 image.format = "png",
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
        TCGAExpression <- string_vars[["envir_link"]]$Results.Completed.crossed
    } else {
        TCGAExpression <- eval(parse(text= paste0("string_vars[['envir_link']]$Results.Completed.",
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
    GSEA.REACTOME <- ReactomePA::gsePathway(GeneSetVector,
                                pvalueCutoff = FDR.cutoff, nPerm = 10000,
                                pAdjustMethod = "BH", verbose = TRUE)

    #Get summary
    GSEA.REACTOME.summary.fdr <- as.data.frame(GSEA.REACTOME)

    # Write
    write.csv(GSEA.REACTOME.summary.fdr, file = paste0(DIR,
                                                   "/GSEA_Output/REAC.GSEA.enrichment",
                                                   ID, "_", pairName, ".csv"), row.names = FALSE)
    assign("GeneSetVector", GeneSetVector, envir = get(envir_link))
    # GSEA.REACTOME.summary.fdr <- GSEA.REACTOME.summary[which(GSEA.REACTOME.summary$p.adjust < FDR.cutoff), ]


    # REACTOME-GSEA Plot
    if (nrow(GSEA.REACTOME.summary.fdr) != 0) {

        if (nrow(GSEA.REACTOME.summary.fdr) >= 15){
            PlotNowREACT <- GSEA.REACTOME.summary.fdr[1:15, ]
        } else{
            PlotNowREACT <- GSEA.REACTOME.summary.fdr
        }

        PlotNowREACT[, 2] <- gsub(" of", "", PlotNowREACT[, 2])
        # large <- lengths(strsplit(PlotNowREACT[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowREACT[, 2], " ")) > 7
        PlotNowREACT[large, 2] <- unname(sapply(PlotNowREACT[large, 2],
                                             function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                               collapse = " ")}))

        # is it better to use the manually adjusted?
        log.10.GSEA <- -log10(as.numeric(PlotNowREACT[, "p.adjust"]))
        new.inf <- log.10.GSEA[order(log.10.GSEA, decreasing = TRUE)]
        log.10.GSEA[log.10.GSEA == "Inf"] <- (new.inf[new.inf != "Inf"][1]+1)

        longest_word <- max(stringr::str_count(PlotNowREACT$Description))


        Gene_Ratio <- unname(sapply(PlotNowREACT$leading_edge,
                                    function(w){unlist(strsplit(w, ", "))[1]}))
        Gene_Ratio <- round(as.numeric(gsub("tags=|%", "", Gene_Ratio))/100, 2)


        if (longest_word > 50) {
            longest_Width <- Width * 1.5
        } else {
            longest_Width <- Width
        }

        if (tolower(image.format) == "png") {
            png(filename = paste0(DIR,
                              "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", ID,
                              "_", pairName, ".png"),
                width = longest_Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = paste0(DIR,
                              "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", ID,
                              "_", pairName, ".svg"),
                width = longest_Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }

        Count <- ceiling(PlotNowREACT$setSize * Gene_Ratio)


        p <- ggplot2::ggplot(PlotNowREACT, ggplot2::aes(x = log.10.GSEA,
                                                     y = forcats::fct_reorder(Description, log.10.GSEA))) +
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
        #     png(filename = paste0(DIR,
        #                           "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", ID,
        #                           "_", pairName, "_2.png"),
        #         width = longest_Width, height = Height, res = Res, units = Unit)
        # } else if (tolower(image.format) == "svg") {
        #     svg(filename = paste0(DIR,
        #                           "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", ID,
        #                           "_", pairName, "_2.svg"),
        #         width = longest_Width, height = Height, onefile = TRUE)
        # } else {
        #     stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        # }
        # par(mar = c(4.5,10,2,2),lwd = 2.5)
        # barplot(log.10.GSEA, horiz = TRUE,
        #         xlab = "", main = "REACTOME-GSEA", lwd = 2, cex.lab = 1.2,
        #         cex.axis = 0.7, cex.main = 1.1, axes = FALSE,
        #         col = RColorBrewer::brewer.pal(8, "Set1")[5], border = "white",
        #         cex.names = 0.8, space = 0.001)
        # title(xlab = "-log(FDR)", line = 2.5, cex.lab = 0.7, family = "Calibri Light")
        # if (-log10(min(as.numeric(PlotNowREACT[, "p.adjust"]))) == "Inf"){
        #     abline(v = 1:floor(0), col = "white", lwd = 4)
        # } else {
        #     abline(v = 1:floor(-log10(min(as.numeric(PlotNowREACT[, "p.adjust"]), na.rm = TRUE))),
        #            col = "white", lwd = 4)
        # }
        #
        # axis(side = 1, 0:(max(log.10.GSEA)+2),
        #      lwd = 1.5, cex.axis = 0.9)
        #
        # axis(2, at=c(0.5:9.5), lwd = 0, line = -1, cex.axis = 0.6,
        #      labels = PlotNowREACT[, "Description"], las = 2)
        #
        # abline(v = -log10(FDR.cutoff), lwd = 3)
        #
        # par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 1, 0), new=TRUE)
        # plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        #
        # legend("topright", "-log(FDR.cutoff) ", lty = 1, lwd = 2,
        #        bty = "n", cex = 0.3)
        #
        # dev.off()
    } else {
        message(cat("There is nothing to show in Gene Set Enrichment Analysis - over Reactome"))
    }
    gc()
    message("Done!\n")
}
