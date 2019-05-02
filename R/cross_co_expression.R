#' Cross co-expression gene list against differential expression genes list
#'
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
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' cross_co_expression(pairName = "G2_over_G1", env = "env name without quotes")
#' }
cross_co_expression <- function(pairName = "G2_over_G1",
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


    dir.create(path = paste0(PATH, "/CrossData_", tolower(FinalData), "/co_expression"), showWarnings = FALSE)
    DIR <- paste0(PATH, "/CrossData_", tolower(FinalData), "/co_expression")

    if (is.null(FinalData)) {
        stop("You must run the 'venn_diagram' fucntion first!")
    }

    # FinalData from venn
    DEcross <- string_vars[["envir_link"]]$resultadosDE_crossed
    # coexp_genes from co_expression
    coexp_genes <- string_vars[["envir_link"]]$coexp_genes

    # save
    DEcross <- DEcross[intersect(rownames(DEcross), coexp_genes), ]
    write.csv(x = DEcross, file = file.path(DIR,
                                            "DE_Co_exp_crossData.csv"))

    #ploting venn diagram
    VennDiagram::venn.diagram(
        x = list(rownames(DEcross), coexp_genes),
        category.names = c(paste0("Venn (", string_vars[["envir_link"]]$FinalData, ")"), "Co_expression"),
        filename = file.path(DIR, paste0("Co_exp_venn_diagramm_",
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


    suppressWarnings(volcano(DEcross))

    tested1 <- vector("list", 1)
    names(tested1) <- pairName
    tested1[[pairName]] <- DEcross

    assign("resultadosDE_crossed_Co", tested1, envir = get(envir_link))
    string_vars[["envir_link"]]$Tool <- "CrossTable_Co_exp"


    message("Done!\n")
}
