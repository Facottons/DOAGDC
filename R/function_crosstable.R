#' Venn diagram of differential expression genes list
#'
#' @param final_data A character string indicating the name of which diferential
#'    expression package should be used to get the statistical values in the
#'    final list. The default is "EBSeq".
#' @param n_pack A numerical value indicating the number of expression analysis
#'    to be used in venn diagram. It is expected the number 2 or 3. The default
#'    is 3.
#' @param package_names A character vector indicating the names of at least two
#'    diferential expression packages used in previous steps: "DESeq2", "edgeR",
#'    "DESeq2, or "All".
#' @param pair_name
#' @param width,height,res,unit
#' @param colors A character vector indicating the colors to be used in the
#'    venn diagram. The default is c('green', 'blue', "red").
#' @param venn_diagram_imagetype A character string indicating the
#'    venn_diagram_imagetype (e.g. "tiff", "png" or "svg").
#'    The default is "png".
#' @param work_dir
#' @param fc_cutoff
#' @param fdr_cutoff
#' @param env
#' @inheritParams dea_ebseq
#' @inheritParams groups_identification_mclust
#' @inheritParams concatenate_exon
#'
#' @return a list of differentially expressed genes in common between two or
#'    three differential expression analysis packages.
#'
#' @importFrom VennDiagram venn.diagram
#' @export
#'
#' @examples
#' library(DOAGDC)
#'
#' # data already downloaded using the 'download_gdc' function
#' concatenate_expression("gene",
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
#'
#' # separating gene HIF3A expression data patients in two groups
#' groups_identification_mclust("gene", 2,
#'     name = "HIF3A",
#'     modelName = "E",
#'     env = CHOL_LEGACY_gene_tumor_data,
#'     tumor = "CHOL"
#' )
#'
#' # load not normalized data
#' concatenate_expression("gene",
#'     normalization = FALSE,
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     env = CHOL_LEGACY_gene_tumor_data,
#'     work_dir = "~/Desktop"
#' )
#'
#' # start DE analysis
#' dea_edger(
#'     name = "HIF3A",
#'     group_gen = "mclust",
#'     env = CHOL_LEGACY_gene_tumor_data
#' )
#'
#' dea_ebseq(
#'     pair_name = "G2_over_G1",
#'     rounds = 2,
#'     name = "HIF3A",
#'     group_gen = "mclust",
#'     env = CHOL_LEGACY_gene_tumor_data
#' )
#'
#' # run the Venn diagram
#' venn_diagram(
#'     final_data = "edgeR",
#'     n_pack = 2,
#'     package_names = c("EBSeq", "edgeR"),
#'     env = CHOL_LEGACY_gene_tumor_data
#' )
venn_diagram <- function(final_data,
                         n_pack = 3,
                         package_names,
                         pair_name = "G2_over_G1",
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
        fdr_ok <- results_completed_local$fdr <= fdr_cutoff
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

        log_10_fdr <- -log10(results_completed_local$fdr)
        new_inf <- log_10_fdr[order(log_10_fdr, decreasing = TRUE)]
        log_10_fdr[log_10_fdr == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        axislimits_y <- ceiling(max(log_10_fdr, na.rm = TRUE)) + 1

        message("Start volcano plot...")
        # Volcano Plot
        if (tolower(venn_diagram_imagetype) == "png") {
            png(
                filename = file.path(dir, paste0(
                    "VolcanoPlot_Basic_",
                    comb_name, ".png"
                )),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(venn_diagram_imagetype) == "svg") {
            svg(
                filename = file.path(dir, paste0(
                    "VolcanoPlot_Basic_",
                    comb_name, ".svg"
                )),
                width = width, height = height, onefile = TRUE
            )
        } else {
            tmp <- "Insert a valid venn_diagram_imagetype! ('png' or 'svg')"
            stop(message(tmp))
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
            x = -axislimits_x + 0.3, y = axislimits_y / 10, labels =
                length(results_completed_local$colour[fc_down & fdr_ok]),
            cex = 1, col = "blue"
        )
        text(
            x = axislimits_x - 0.4, y = axislimits_y / 10, labels =
                length(results_completed_local$colour[fc_up & fdr_ok]),
            cex = 1, col = "red"
        )
        par(
            fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
            new = TRUE
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

    comb_name <- paste0(tolower(package_names), collapse = "_")

    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please",
            " insert the 'env' name and try again!"
        ))
    }

    ev <- deparse(substitute(env))
    sv <- list(ev = get(ev))

    if (exists("name_e", envir = get(ev))) {
        path <- file.path(
            sv[["ev"]]$path,
            sv[["ev"]]$name_e
        )
    } else {
        path <- sv[["ev"]]$path
    }

    dir.create(
        path = paste0(path, "/CrossData_", tolower(final_data)),
        showWarnings = FALSE
    )
    dir <- paste0(path, "/CrossData_", tolower(final_data))

    if (is.null(final_data)) {
        stop(
            "Please insert which DE method is going",
            " to be used in result!! (EBSEq, DESeq2 or edgeR?)"
        )
    }

    if (tolower(final_data) == "ebseq") {
        de_cross <- sv[["ev"]]$resultados_de_ebseq[[pair_name]]
        complete_cross <- sv[["ev"]]$results_completed_ebseq[[pair_name]]
    } else if (tolower(final_data) == "edger") {
        de_cross <- sv[["ev"]]$resultados_de_edger[[pair_name]]
        complete_cross <- sv[["ev"]]$results_completed_edger[[pair_name]]
    } else if (tolower(final_data) == "deseq2") {
        de_cross <- sv[["ev"]]$resultados_de_deseq2[[pair_name]]
        complete_cross <- sv[["ev"]]$results_completed_deseq2[[pair_name]]
    }

    if ("all" %in% tolower(package_names)) {
        package_names <- c("edger", "deseq2", "ebseq")
    } else {
        package_names <- tolower(package_names)
    }

    if (n_pack == 3) {
        p1 <- paste0(
            "sv[['ev']]",
            "$resultados_de_",
            package_names[tolower(package_names) == tolower(final_data)],
            "[[pair_name]]"
        )
        p2 <- paste0(
            "sv[['ev']]",
            "$resultados_de_",
            package_names[tolower(package_names) != tolower(final_data)][1],
            "[[pair_name]]"
        )
        p3 <- paste0(
            "sv[['ev']]",
            "$resultados_de_",
            package_names[tolower(package_names) != tolower(final_data)][2],
            "[[pair_name]]"
        )
        p1completed <- paste0(
            "sv[['ev']]",
            "$results_completed_",
            package_names[tolower(package_names) == tolower(final_data)],
            "[[pair_name]]"
        )
        p2completed <- paste0(
            "sv[['ev']]",
            "$results_completed_",
            package_names[tolower(package_names) != tolower(final_data)][1],
            "[[pair_name]]"
        )
        p3completed <- paste0(
            "sv[['ev']]",
            "$results_completed_",
            package_names[tolower(package_names) != tolower(final_data)][2],
            "[[pair_name]]"
        )

        tmp1 <- rownames(eval(parse(text = p1)))
        tmp2 <- rownames(eval(parse(text = p2)))
        tmp3 <- rownames(eval(parse(text = p3)))

        de_cross <- de_cross[which(tmp1 %in% tmp2), ]
        de_cross <- de_cross[which(rownames(de_cross) %in% tmp3), ]
        tmp4 <- rownames(eval(parse(text = p1completed)))
        tmp5 <- rownames(eval(parse(text = p2completed)))
        complete_cross <- complete_cross[which(tmp4 %in% tmp5), ]
        tmp4 <- rownames(complete_cross)
        tmp5 <- rownames(eval(parse(text = p3completed)))
        complete_cross <- complete_cross[which(tmp4 %in% tmp5), ]
        write.csv(x = de_cross, file = file.path(
            dir,
            "DE_crossData.csv"
        ))
        write.csv(x = complete_cross, file = file.path(
            dir,
            "complete_crossData.csv"
        ))

        tmp <- paste0(tolower(package_names), c(rep("_Vs_", 2), ""),
            collapse = ""
        )

        # ploting venn diagram
        VennDiagram::venn.diagram(
            x = list(tmp1, tmp2, tmp3),
            category.names = c(
                package_names[tolower(package_names) == tolower(final_data)],
                package_names[tolower(package_names) != tolower(final_data)][1],
                package_names[tolower(package_names) != tolower(final_data)][2]
            ),
            filename = file.path(dir, paste0(
                "venn_diagramm_",
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
            cat.pos = c(-40, 40, 135), cat.dist = c(0.07, 0.07, 0.07),
            cat.fontfamily = "sans"
        )
    }
    else if (n_pack == 2) {
        p1 <- paste0(
            "sv[['ev']]",
            "$resultados_de_",
            package_names[tolower(package_names) == tolower(final_data)],
            "[[pair_name]]"
        )
        p2 <- paste0(
            "sv[['ev']]",
            "$resultados_de_",
            package_names[tolower(package_names) != tolower(final_data)][1],
            "[[pair_name]]"
        )
        p1completed <- paste0(
            "sv[['ev']]",
            "$results_completed_",
            package_names[tolower(package_names) == tolower(final_data)],
            "[[pair_name]]"
        )
        p2completed <- paste0(
            "sv[['ev']]",
            "$results_completed_",
            package_names[tolower(package_names) != tolower(final_data)],
            "[[pair_name]]"
        )
        tmp1 <- rownames(eval(parse(text = p1)))
        tmp2 <- rownames(eval(parse(text = p2)))
        de_cross <- de_cross[which(tmp1 %in% tmp2), ]
        tmp1 <- rownames(eval(parse(text = p1completed)))
        tmp2 <- rownames(eval(parse(text = p2completed)))
        complete_cross <- complete_cross[which(tmp1 %in% tmp2), ]
        write.csv(x = de_cross, file = file.path(
            dir,
            "DE_crossData.csv"
        ))
        write.csv(x = complete_cross, file = file.path(
            dir,
            "complete_crossData.csv"
        ))

        tmp <- paste0(tolower(package_names), c("_Vs_", ""), collapse = "")

        # ploting venn diagram
        VennDiagram::venn.diagram(
            x = list(rownames(eval(parse(text = p1))),
                        rownames(eval(parse(text = p2)))),
            category.names = c(
                package_names[tolower(package_names) == tolower(final_data)],
                package_names[tolower(package_names) != tolower(final_data)]
            ),
            filename = file.path(dir, paste0(
                "venn_diagramm_",
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
            fill = colors[1:2],
            cex = 1.2, print.mode = c("raw", "percent"), sigdigs = 2,
            fontface = "bold", fontfamily = "sans",
            cat.cex = 1.6, cat.fontface = "bold", cat.default.pos = "outer",
            cat.pos = c(-35, 30), cat.dist = c(0.055, -0.05),
            cat.fontfamily = "sans"
        )
    }

    tested1 <- vector("list", 1)
    names(tested1) <- pair_name
    tested1[[pair_name]] <- de_cross

    tested2 <- vector("list", 1)
    names(tested2) <- pair_name
    tested2[[pair_name]] <- complete_cross

    suppressWarnings(volcano(complete_cross))

    assign("final_data", final_data, envir = get(ev))
    assign("results_completed_crossed", tested2, envir = get(ev))
    assign("resultados_de_crossed", tested1, envir = get(ev))
    if (tolower(final_data) == "ebseq") {
        sv[["ev"]]$tool <- "crosstable_ebseq"
    } else if (tolower(final_data) == "edger") {
        sv[["ev"]]$tool <- "crosstable_edger"
    } else if (tolower(final_data) == "deseq2") {
        sv[["ev"]]$tool <- "crosstable_deseq2"
    }

    message("Done!\n")
}
