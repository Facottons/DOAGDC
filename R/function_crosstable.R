#' Venn diagram of differential expression genes list
#'
#' @param FinalData A character string indicating the name of which diferential
#'   expression package should be used to get the statistical values in the
#'   final list. The default is "EBSeq".
#' @param n_pack A numerical value indicating the number of expression analysis
#'   to be used in venn diagram. It is expected the number 2 or 3. The default
#'   is 3.
#' @param packageNames A character vector indicating the names of at least two
#'   diferential expression packages used in previous steps: "DESeq2", "edgeR",
#'   "DESeq2, or "All".
#' @param pairName
#' @param Width,Height,Res,Unit
#' @param Colors A character vector indicating the colors to be used in the venn
#'   diagram. The default is c('green', 'blue', "red").
#' @param VennDiagram_imagetype A character string indicating the image_format
#'   (e.g. "tiff", "png" or "svg"). The default is "png".
#' @param workDir
#' @param env
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#'
#' @return a list of differentially expressed genes in common between two or
#'   three differential expression analysis packages.
#'
#' @importFrom VennDiagram venn.diagram
#' @export
#'
#' @examples
#' \dontrun{
#' venn_diagram("DESeq2", 2, packageNames = c("edgeR", "DESeq2"), env = "env name without quotes")
#' }
venn_diagram <- function(FinalData, n_pack = 3,
                      packageNames,
                      pairName = "G2_over_G1",
                      Width = 2000,
                      Height = 2000,
                      Res = 300,
                      Unit = "px",
                      Colors = c('green', 'blue', "red"),
                      VennDiagram_imagetype = "png",
                      workDir,
                      env) {


    # local functions ####

    volcano <- function(Results_Completed_local){


        # START VOLCANO PLOT
        Results_Completed_local$Colour = rgb(100, 100, 100, 50, maxColorValue = 255)

        # Set new column values to appropriate colours
        Results_Completed_local$Colour[Results_Completed_local$log2FC >= log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff] <- rgb(222, 22, 22, 50, maxColorValue = 255)
        Results_Completed_local$Colour[Results_Completed_local$log2FC < -log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff] <- rgb(56, 50, 237, 50, maxColorValue = 255)

        ####Volcano Plot
        axislimits_x <- ceiling(max(c(-min(Results_Completed_local$log2FC, na.rm = TRUE) - 1,
                                      max(Results_Completed_local$log2FC, na.rm = TRUE) + 1)))

        log_10_FDR <- -log10(Results_Completed_local$FDR)
        new_inf <- log_10_FDR[order(log_10_FDR, decreasing = TRUE)]
        log_10_FDR[log_10_FDR == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        axislimits_y <- ceiling(max(log_10_FDR, na.rm = TRUE)) + 1

        message("Start volcano plot...")
        #Volcano Plot
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, paste0("VolcanoPlot_Basic_",
                                                 comb_name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, paste0("VolcanoPlot_Basic_",
                                                 comb_name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(4,6,3,2), mgp = c(2,.7,0), tck = -0.01)
        plot(Results_Completed_local$log2FC, log_10_FDR, axes = FALSE,
             xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
             xlab = expression('log'[2]*'(FC)'),
             ylab = "",
             # main = "Volcano Plot",
             cex.lab = 1.5, cex.main = 2,
             cex.sub = 2,
             pch = 16, col = Results_Completed_local$Colour, cex = 2.5, las = 1)
        title(ylab = "-log(FDR)", line = 4, cex.lab = 1.5, family = "Calibri Light")
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5, las = 1)
        abline(v = log2(FC_cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(v = -log2(FC_cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(h = -log10(FDR_cutoff), col = "black", lwd = 4, lty = 3)
        text(x = -axislimits_x + 0.3, y = axislimits_y/10, labels =
                 length(Results_Completed_local$Colour[Results_Completed_local$log2FC <= -log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff]),
             cex = 1, col = "blue")
        text(x = axislimits_x - 0.4, y = axislimits_y/10, labels =
                 length(Results_Completed_local$Colour[Results_Completed_local$log2FC >= log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff]),
             cex = 1, col = "red")
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

        legend("topright", border = FALSE, bty = "n",
               legend = c(paste0("-log(", FDR_cutoff, ") = ", round(-log10(FDR_cutoff), 2)),
                        paste0("\u00B1", " log", "\u2082","(", FC_cutoff, ") = ",log2(FC_cutoff)),
                        "UP", "DOWN"), pt.cex = c(0, 0, 1.8, 1.8),
               lty = c(3, 6, 0, 0), pch = c(-1, -1, 16, 16), cex = c(0.8, 0.8, 0.8, 0.8), lwd = c(2, 2, 5, 5),
               col = c("Black", "Black", "red3", "blue3"))
        dev.off()
    }


    # code ####
    if (missing(env)) {stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))) {
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    dir.create(path = paste0(PATH, "/CrossData_", tolower(FinalData)), showWarnings = FALSE)
    DIR <- paste0(PATH, "/CrossData_", tolower(FinalData))

    if (is.null(FinalData)) {
        stop("Please insert which DE method is going to be used in result!! (EBSEq, DESeq2 or edgeR?)")
    }

    if (tolower(FinalData) == "ebseq") {
        DEcross <- string_vars[["envir_link"]]$resultadosDE_EBSeq[[pairName]]
        complete_cross <- string_vars[["envir_link"]]$Results_Completed_EBSeq[[pairName]]
    } else if (tolower(FinalData) == "edger") {
        DEcross <- string_vars[["envir_link"]]$resultadosDE_edgeR[[pairName]]
        complete_cross <- string_vars[["envir_link"]]$Results_Completed_edgeR[[pairName]]
    } else if (tolower(FinalData) == "deseq2") {
        DEcross <- string_vars[["envir_link"]]$resultadosDE_DESeq2[[pairName]]
        complete_cross <- string_vars[["envir_link"]]$Results_Completed_DESeq2[[pairName]]
    }

    if (tolower(packageNames) == "all") {
        packageNames <- c("edgeR", "DESeq2", "EBSeq")
    }

    if (n_pack == 3) {
        p1 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE_",
                     packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE_",
                     packageNames[tolower(packageNames) != tolower(FinalData)][1], "[[pairName]]")
        p3 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE_",
                     packageNames[tolower(packageNames) != tolower(FinalData)][2], "[[pairName]]")
        p1completed <- paste0("string_vars[['envir_link']]",
                              "$Results_Completed_",
                              packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2completed <- paste0("string_vars[['envir_link']]",
                              "$Results_Completed_",
                              packageNames[tolower(packageNames) != tolower(FinalData)][1], "[[pairName]]")
        p3completed <- paste0("string_vars[['envir_link']]",
                              "$Results_Completed_",
                              packageNames[tolower(packageNames) != tolower(FinalData)][2], "[[pairName]]")
        DEcross <- DEcross[which(rownames(eval(parse(text = p1))) %in% rownames(eval(parse(text = p2)))), ]
        DEcross <- DEcross[which(rownames(DEcross) %in% rownames(eval(parse(text = p3)))), ]
        complete_cross <- complete_cross[which(rownames(eval(parse(text = p1completed))) %in% rownames(eval(parse(text = p2completed)))), ]
        complete_cross <- complete_cross[which(rownames(complete_cross) %in% rownames(eval(parse(text = p3completed)))), ]
        write.csv(x = DEcross, file = file.path(DIR,
                                               "DE_crossData.csv"))
        write.csv(x = complete_cross, file = file.path(DIR,
                                                      "complete_crossData.csv"))

        tmp <- paste0(tolower(packageNames), c(rep("_Vs_", 2), ""), collapse = "")

        #ploting venn diagram
        VennDiagram::venn.diagram(
            x = list(rownames(eval(parse(text = p1))), rownames(eval(parse(text = p2))),
                     rownames(eval(parse(text = p3)))),
            category.names = c(packageNames[tolower(packageNames) == tolower(FinalData)],
                               packageNames[tolower(packageNames) != tolower(FinalData)][1],
                               packageNames[tolower(packageNames) != tolower(FinalData)][2]),
            filename = file.path(DIR, paste0("venn_diagramm_",
                                             pairName, "_",
                                             tmp,".", tolower(VennDiagram_imagetype))),
            resolution = Res,
            output = TRUE, imagetype = VennDiagram_imagetype, height = Height, width = Width,
            units = Unit,
            compression = "lzw",
            lwd = 2, lty = 'blank',
            fill = Colors,
            cex = 1.2, print.mode = c("raw", "percent"), sigdigs = 2,
            fontface = "bold", fontfamily = "sans",
            cat.cex = 1.6, cat.fontface = "bold", cat.default.pos = "outer",
            cat.pos = c(-40, 40, 135), cat.dist = c(0.07, 0.07, 0.07),
            cat.fontfamily = "sans")
    }
    else if (n_pack == 2) {
        p1 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE_",
                     packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE_",
                     packageNames[tolower(packageNames) != tolower(FinalData)][1], "[[pairName]]")
        p1completed <- paste0("string_vars[['envir_link']]",
                              "$Results_Completed_",
                              packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2completed <- paste0("string_vars[['envir_link']]",
                              "$Results_Completed_",
                              packageNames[tolower(packageNames) != tolower(FinalData)], "[[pairName]]")
        DEcross <- DEcross[which(rownames(eval(parse(text = p1))) %in% rownames(eval(parse(text = p2)))), ]
        complete_cross <- complete_cross[which(rownames(eval(parse(text = p1completed))) %in% rownames(eval(parse(text = p2completed)))), ]
        write.csv(x = DEcross, file = file.path(DIR,
                                               "DE_crossData.csv"))
        write.csv(x = complete_cross, file = file.path(DIR,
                                                      "complete_crossData.csv"))

        tmp <- paste0(tolower(packageNames), c("_Vs_", ""), collapse = "")

        #ploting venn diagram
        VennDiagram::venn.diagram(
            x = list(rownames(eval(parse(text = p1))), rownames(eval(parse(text = p2)))),
            category.names = c(packageNames[tolower(packageNames) == tolower(FinalData)],
                               packageNames[tolower(packageNames) != tolower(FinalData)]),
            filename = file.path(DIR, paste0("venn_diagramm_",
                                             pairName, "_",
                                             tmp,".", tolower(VennDiagram_imagetype))),
            resolution = Res,
            output = TRUE, imagetype = VennDiagram_imagetype, height = Height, width = Width,
            units = Unit,
            compression = "lzw",
            lwd = 2, lty = 'blank',
            fill = Colors,
            cex = 1.2, print.mode = c("raw", "percent"), sigdigs = 2,
            fontface = "bold", fontfamily = "sans",
            cat.cex = 1.6, cat.fontface = "bold", cat.default.pos = "outer",
            cat.pos = c(-35, 30), cat.dist = c(0.055, -0.05),
            cat.fontfamily = "sans")
    }

    tested1 <- vector("list", 1)
    names(tested1) <- pairName
    tested1[[pairName]] <- DEcross

    tested2 <- vector("list", 1)
    names(tested2) <- pairName
    tested2[[pairName]] <- complete_cross

    suppressWarnings(volcano(complete_cross))

    assign("FinalData", FinalData, envir = get(envir_link))
    assign("Results_Completed_crossed", tested2, envir = get(envir_link))
    assign("resultadosDE_crossed", tested1, envir = get(envir_link))
    if (tolower(FinalData) == "ebseq") {
        string_vars[["envir_link"]]$Tool <- "CrossTable_EBSeq"
    } else if (tolower(FinalData) == "edger") {
        string_vars[["envir_link"]]$Tool <- "CrossTable_edgeR"
    } else if (tolower(FinalData) == "deseq2") {
        string_vars[["envir_link"]]$Tool <- "CrossTable_DESeq2"
    }

    message("Done!\n")
}
