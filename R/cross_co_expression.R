#' Cross co-expression gene list against differential expression genes list
#'
#' @param pair_name
#' @param final_data A character string indicating the name of which diferential
#'    expression package should be used to get the statistical values in the
#'    final list. The default is "EBSeq".
#' @param width,height,res,unit
#' @param colors A character vector indicating the colors to be used in the
#'    venn diagram. The default is c('green', 'blue', "red").
#' @param venn_diagram_imagetype A character string indicating the image_format
#'    (e.g. "tiff", "png" or "svg"). The default is "png".
#' @param work_dir
#' @param fc_cutoff
#' @param fdr_cutoff
#' @param env
#' @inheritParams dea_ebseq
#' @inheritParams groups_identification_mclust
#' @inheritParams concatenate_exon
#'
#' @return a list of genes co-expressed and differentially expressed genes
#'    inside the determined environment name for further use.
#'
#' @export
#'
#' @importFrom VennDiagram venn.diagram
cross_co_expression <- function(pair_name = "G2_over_G1",
                                final_data,
                                width = 2000,
                                height = 2000,
                                res = 300,
                                unit = "px",
                                colors = c("green", "blue", "red"),
                                venn_diagram_imagetype = "png",
                                work_dir,
                                fc_cutoff = 2,
                                fdr_cutoff = 0.05,
                                env) {

    # local functions ####
    volcano <- function(results_completed_local) {

        # START VOLCANO PLOT
        results_completed_local$colour <- hsv(0, 0, .39, .2)

        # Set new column values to appropriate colours
        fc_up <- results_completed_local$log2FC >= log2(fc_cutoff)
        fdr_ok <- results_completed_local$FDR <= fdr_cutoff
        results_completed_local$colour[fc_up & fdr_ok] <- hsv(0, .9, .87, .5)
        fc_down <- results_completed_local$log2FC < -log2(fc_cutoff)
        tmp <- hsv(.67, .79, .93, .2)
        results_completed_local$colour[fc_down & fdr_ok] <- tmp

        #### Volcano Plot
        axislimits_x <- ceiling(max(c(
            -min(results_completed_local$log2FC,
                na.rm = TRUE
            ) - 1,
            max(results_completed_local$log2FC,
                na.rm = TRUE
            ) + 1
        )))

        log_10_fdr <- -log10(results_completed_local$FDR)
        new_inf <- log_10_fdr[order(log_10_fdr, decreasing = TRUE)]
        log_10_fdr[log_10_fdr == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        axislimits_y <- ceiling(max(log_10_fdr, na.rm = TRUE)) + 1

        message("Start volcano plot...")
        # Volcano Plot
        if (tolower(venn_diagram_imagetype) == "png") {
            png(
                filename = file.path(dir, paste0(
                    "VolcanoPlot_Basic_",
                    ".png"
                )),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(venn_diagram_imagetype) == "svg") {
            svg(
                filename = file.path(dir, paste0(
                    "VolcanoPlot_Basic_",
                    ".svg"
                )),
                width = width, height = height, onefile = TRUE
            )
        }
        par(mar = c(4, 6, 3, 2), mgp = c(2, .7, 0), tck = -0.01)
        plot(results_completed_local$log2FC, log_10_fdr,
            axes = FALSE,
            xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
            xlab = expression("log"[2] * "(FC)"),
            ylab = "",
            # main = "Volcano Plot",
            cex.lab = 1.5, cex.main = 2,
            cex.sub = 2,
            pch = 16, col = results_completed_local$colour,
            cex = 2.5, las = 1
        )
        title(
            ylab = "-log(FDR)", line = 4, cex.lab = 1.5,
            family = "Calibri Light"
        )
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5, las = 1)
        abline(
            v = log2(fc_cutoff), col = "black", lty = 6, cex = 0.8,
            lwd = 4
        )
        abline(
            v = -log2(fc_cutoff), col = "black", lty = 6, cex = 0.8,
            lwd = 4
        )
        abline(h = -log10(fdr_cutoff), col = "black", lwd = 4, lty = 3)
        text(
            x = -axislimits_x + 0.3,
            y = axislimits_y / 10,
            labels = length(results_completed_local$colour[fc_down & fdr_ok]),
            cex = 1, col = "blue"
        )
        text(
            x = axislimits_x - 0.4,
            y = axislimits_y / 10,
            labels = length(results_completed_local$colour[fc_up & fdr_ok]),
            cex = 1, col = "red"
        )
        par(
            fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
            mar = c(0, 0, 0, 0), new = TRUE
        )
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

        legend("topright",
            border = FALSE, bty = "n",
            legend = c(
                paste0(
                    "-log(", fdr_cutoff, ") = ",
                    round(-log10(fdr_cutoff), 2)
                ),
                paste0(
                    "\u00B1", " log", "\u2082", "(",
                    fc_cutoff, ") = ", log2(fc_cutoff)
                ),
                "UP", "DOWN"
            ), pt.cex = c(0, 0, 1.8, 1.8),
            lty = c(3, 6, 0, 0), pch = c(-1, -1, 16, 16),
            cex = c(0.8, 0.8, 0.8, 0.8), lwd = c(2, 2, 5, 5),
            col = c("Black", "Black", "red3", "blue3")
        )
        dev.off()
    }


    # code ####
    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please insert",
            " the 'env' name and try again!"
        ))
    }

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("name_e", envir = get(envir_link))) {
        path <- file.path(
            string_vars[["envir_link"]]$path,
            string_vars[["envir_link"]]$name_e
        )
    } else {
        path <- string_vars[["envir_link"]]$path
    }


    dir.create(path = paste0(
        path, "/CrossData_", tolower(final_data),
        "/co_expression"
    ), showWarnings = FALSE)
    dir <- paste0(path, "/CrossData_", tolower(final_data), "/co_expression")

    if (is.null(final_data)) {
        stop("You must run the 'venn_diagram' fucntion first!")
    }

    # final_data from venn
    de_cross <- string_vars[["envir_link"]]$resultados_de_crossed
    # coexp_genes from co_expression
    coexp_genes <- string_vars[["envir_link"]]$coexp_genes

    # save
    de_cross <- de_cross[intersect(rownames(de_cross), coexp_genes), ]
    write.csv(x = de_cross, file = file.path(
        dir,
        "DE_Co_exp_crossData.csv"
    ))

    # ploting venn diagram
    VennDiagram::venn.diagram(
        x = list(rownames(de_cross), coexp_genes),
        category.names = c(
            paste0(
                "Venn (",
                string_vars[["envir_link"]]$final_data, ")"
            ),
            "Co_expression"
        ),
        filename = file.path(dir, paste0(
            "Co_exp_venn_diagramm_",
            pair_name, "_",
            tmp, ".",
            tolower(venn_diagram_imagetype)
        )),
        resolution = res,
        output = TRUE, imagetype = venn_diagram_imagetype,
        height = height, width = width,
        units = unit,
        compression = "lzw",
        lwd = 2, lty = "blank",
        fill = colors,
        cex = 1.2, print.mode = c("raw", "percent"), sigdigs = 2,
        fontface = "bold", fontfamily = "sans",
        cat.cex = 1.6, cat.fontface = "bold", cat.default.pos = "outer",
        cat.pos = c(-35, 30), cat.dist = c(0.055, -0.05),
        cat.fontfamily = "sans"
    )


    suppressWarnings(volcano(de_cross))

    tested1 <- vector("list", 1)
    names(tested1) <- pair_name
    tested1[[pair_name]] <- de_cross

    assign("resultados_de_crossed_Co", tested1, envir = get(envir_link))
    string_vars[["envir_link"]]$tool <- "CrossTable_Co_exp"


    message("Done!\n")
}
