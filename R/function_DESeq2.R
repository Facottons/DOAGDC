#' Run DESeq2 gene Differential Expression Analysis (DEA).
#'
#' @param name
#' @param core_number A numeric value indicating how many CPU cores should be
#'    used in the analysis. The default value is 2.
#' @param test A character string indicating which test should be used:
#'    \code{"LRT", "wald" or "Default Test"}. The default is \code{"Default
#'    Test"}.
#' @param clinical_pair
#' @param group_gen
#' @param fc_cutoff
#' @param work_dir
#' @param tumor
#' @param fdr_cutoff
#' @param width,height,res,unit,image_format
#' @param env
#' @param cooks_cutoff Cooks distance remove outliers from the analysis; More
#'    details in DESeq2 \link{results} page. The default is \code{FALSE}.
#' @inheritParams download_gdc
#' @inheritParams dea_ebseq
#' @inheritParams concatenate_exon
#' @inheritParams groups_identification_mclust
#'
#' @return A matrix with DE genes in row and statistical values in columns.
#' @export
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel register
#' @importFrom BiocParallel SnowParam
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
#' dea_deseq2(
#'     name = "HIF3A",
#'     test = "LRT",
#'     env = CHOL_LEGACY_gene_tumor_data,
#'     group_gen = "mclust"
#' )
dea_deseq2 <- function(name,
                       core_number = 2,
                       test = "Default Test",
                       group_gen,
                       clinical_pair,
                       fc_cutoff = 2,
                       work_dir,
                       tumor,
                       fdr_cutoff = 0.05,
                       width = 2000,
                       height = 1500,
                       res = 300,
                       unit = "px",
                       image_format = "png",
                       env,
                       cooks_cutoff = FALSE) {

    # DESeq2 normalization does not account for gene length, and there are
    # sound reasons for making that choice when using the data for statistical
    # hypothesis testing. Visualization based on regularized log transformation
    # should not be biased based on gene length.  However, gene expression in
    # RNA-seq does have a gene length bias "built-in"; this is a function of
    # the "count" nature of RNA-seq and not due to any software processing of
    # the data.


    # #verifying if the package is already installed

    # local functions ####
    deseq_plots <- function(dds) {
        if (tolower(image_format) == "png") {
            png(
                filename = paste0(dir, "/plotDispEsts.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(dir, "/plotDispEsts.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        DESeq2::plotDispEsts(dds,
            main = paste0("Test = ", toupper(test)),
            las = 1, ylab = "Dispersion",
            xlab = "Mean of normalized counts"
        )
        dev.off()

        if (tolower(image_format) == "png") {
            png(
                filename = paste0(
                    dir,
                    "/deseq2_MAplot.png"
                ),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(
                    dir,
                    "/deseq2_MAplot.svg"
                ),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        DESeq2::plotMA(dds,
            ylim = c(-log2(fc_cutoff) - 1, log2(fc_cutoff) + 1),
            main = "DESeq2",
            las = 1, ylab = "Log2FC",
            xlab = "Mean of normalized counts"
        )
        abline(h = log2(fc_cutoff), col = "dodgerblue", lwd = 2)
        abline(h = -log2(fc_cutoff), col = "dodgerblue", lwd = 2)
        dev.off()

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "transformation_effect.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "transformation_effect.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
            stop(message(tmp))
        }
        par(mar = c(4, 4, 1, 2))
        px <- DESeq2::counts(dds)[, 1] / DESeq2::sizeFactors(dds)[1]
        ord <- order(px)
        ord <- ord[px[ord] < 150]
        ord <- ord[seq(1, length(ord), length = 50)]
        vstcol <- c("blue", "black")
        matplot(px[ord], cbind(
            SummarizedExperiment::assay(var_trans)[, 1],
            log2(px)
        )[ord, ],
        type = "l", lty = 1,
        col = vstcol, xlab = "n", ylab = "f(n)", las = 1
        )
        legend("bottomright",
            legend = c(
                expression("variance stabilizing transformation"),
                expression(log[2](n / s[1]))
            ), fill = vstcol,
            cex = 0.8
        )
        dev.off()
    }

    volcano <- function(results_completed, pairs) {
        if (pairs > 0) {
            comb_name <- names(results_completed[pairs])
            results_compl_local <- results_completed[[pairs]]
        } else {
            comb_name <- "G2_over_G1"
            results_compl_local <- results_completed
        }

        # START VOLCANO PLOT
        results_compl_local$colour <- hsv(0, 0, .39, 0.5)

        # Set new column values to appropriate colours
        tmp <- results_compl_local$log2FC >= log2(fc_cutoff)
        tmp1 <- results_compl_local$FDR <= fdr_cutoff
        results_compl_local$colour[tmp & tmp1] <- hsv(0, .9, .87, .5)
        tmp <- results_compl_local$log2FC < -log2(fc_cutoff)
        results_compl_local$colour[tmp & tmp1] <- hsv(.67, .79, .93, .5)

        #### Volcano Plot
        axislimits_x <- ceiling(max(c(
            -min(results_compl_local$log2FC, na.rm = TRUE) - 1,
            max(results_compl_local$log2FC, na.rm = TRUE) + 1
        )))

        log_10_fdr <- -log10(results_compl_local$FDR)
        new_inf <- log_10_fdr[order(log_10_fdr, decreasing = TRUE)]
        log_10_fdr[log_10_fdr == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        axislimits_y <- ceiling(max(log_10_fdr, na.rm = TRUE)) + 1

        message("Start volcano plot...")
        # Volcano Plot
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, paste0(
                    "VolcanoPlot_Basic_",
                    comb_name, ".png"
                )),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, paste0(
                    "VolcanoPlot_Basic_",
                    comb_name, ".svg"
                )),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid ",
                "image_format! ('png' or 'svg')"
            ))
        }
        par(mar = c(4, 6, 3, 2), mgp = c(2, .7, 0), tck = -0.01)
        plot(results_compl_local$log2FC, log_10_fdr,
            axes = FALSE,
            xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
            xlab = expression("log"[2] * "(fc)"),
            ylab = "",
            cex.lab = 1.5, cex.main = 2,
            cex.sub = 2,
            pch = 16, col = results_compl_local$colour, cex = 2.5, las = 1
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
        tmp1 <- results_compl_local$log2FC <= -log2(fc_cutoff)
        tmp2 <- results_compl_local$FDR <= fdr_cutoff
        text(
            x = -axislimits_x + 0.3, y = axislimits_y / 10,
            labels = length(results_compl_local$colour[tmp1 & tmp2]),
            cex = 1, col = "blue"
        )

        tmp1 <- results_compl_local$log2FC >= log2(fc_cutoff)

        text(
            x = axislimits_x - 0.4, y = axislimits_y / 10,
            labels = length(results_compl_local$colour[tmp1 & tmp2]),
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
                    "\u00B1", " log", "\u2082", "(", fc_cutoff,
                    ") = ", log2(fc_cutoff)
                ),
                "UP", "DOWN"
            ), pt.cex = c(0, 0, 1.8, 1.8),
            lty = c(3, 6, 0, 0), pch = c(-1, -1, 16, 16),
            cex = c(0.8, 0.8, 0.8, 0.8), lwd = c(2, 2, 5, 5),
            col = c("Black", "Black", "red3", "blue3")
        )
        dev.off()
    }

    # Code ####

    BiocParallel::register(BiocParallel::SnowParam(core_number))

    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please insert",
            " the 'env' name and try again!"
        ))
    }

    ev <- deparse(substitute(env))
    sv <- list(ev = get(ev))

    name <- gsub("-", "_", name)

    work_dir <- ifelse(missing("work_dir"), sv[["ev"]]$work_dir, work_dir)

    data_base <- sv[["ev"]]$data_base

    assign("path", file.path(
        work_dir, "DOAGDC",
        toupper(sv[["ev"]]$tumor),
        "Analyses"
    ), envir = get(ev))

    if (missing(group_gen)) {
        stop(message(
            "Please insert a valid 'group_gen' value!!",
            " ('mlcust', 'coxHR' or 'clinical')"
        ))
    }

    assign("group_gen", group_gen, envir = get(ev))

    path <- ifelse(
        exists("name_e", envir = get(ev)),
        file.path(
            sv[["ev"]]$path,
            sv[["ev"]]$name_e
        ), sv[["ev"]]$path
    )

    dir.create(path, showWarnings = FALSE)

    dir.create(paste0(
        path, "/DESeq2_Results.", tolower(group_gen), "_",
        toupper(name)
    ),
    showWarnings = FALSE
    )
    dir <- paste0(
        path, "/DESeq2_Results.", tolower(group_gen), "_",
        toupper(name)
    )

    # Cooks distance remove outliers from the analysis, it looks to see how
    # much each sample contributes to a genes overall value fold change, with
    # samples that cause extreme effects removed.

    # Preparing groups data
    if (tolower(group_gen) == "mclust") {
        grupos_deseq2 <- as.data.frame(sv[["ev"]]$groups)
        grupos_deseq2[, "Selected_classification"] <- factor(
            grupos_deseq2[, "Selected_classification"]
        )
        colnames(grupos_deseq2) <- c("condition", "type")
        grupos_deseq2[, 2] <- c("paired-end")
    } else if (tolower(group_gen) == "clinical") {
        grupos_deseq2 <- as.data.frame(
            sv[["ev"]]$clinical_groups_clinical[[clinical_pair]]
        )
        grupos_deseq2[, "Selected_classification"] <- factor(
            grupos_deseq2[, "Selected_classification"]
        )
        colnames(grupos_deseq2) <- c("condition", "type")
        grupos_deseq2[, 2] <- c("paired-end")
    } else if (tolower(group_gen) == "coxhr") {
        grupos_deseq2 <- as.data.frame(sv[["ev"]]$clinical_groups)
        grupos_deseq2$classification <- gsub(
            "low", "1",
            grupos_deseq2$classification
        )
        grupos_deseq2$classification <- gsub(
            "high", "2",
            grupos_deseq2$classification
        )
        grupos_deseq2[, "classification"] <- factor(
            grupos_deseq2[, "classification"]
        )
        colnames(grupos_deseq2) <- c("type", "condition")
        grupos_deseq2[, 1] <- c("paired-end")
        grupos_deseq2 <- grupos_deseq2[, c(2, 1)]
    }


    # check patient in common
    tmp1 <- rownames(grupos_deseq2)
    tmp2 <- colnames(sv[["ev"]]$gene_tumor_not_normalized)

    grupos_deseq2 <- grupos_deseq2[tmp1 %in% tmp2, ]

    assign("cond_heatmap", grupos_deseq2[, 1], envir = get(ev))

    # selecting specifics patients
    completed_matrix <- sv[["ev"]]$gene_tumor_not_normalized[, rownames(
        grupos_deseq2
    )]

    message("Filtering data ...")
    # remove residual data exclude rows with all zeroes (genes with no counts)
    filter_rows <- rowSums(completed_matrix) >= 10
    completed_matrix <- completed_matrix[filter_rows, ]
    message("It was filtered = ", as.numeric(table(filter_rows)[1]), " genes!")

    completed_matrix <- round(completed_matrix)
    storage.mode(completed_matrix) <- "integer"

    data_base <- sv[["ev"]]$data_base
    # Importing data into DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = completed_matrix,
        colData = grupos_deseq2,
        design = ~condition
    )

    # separate conditoins for more than 1 pair
    group2_number <- max(as.numeric(levels(grupos_deseq2[, 1])))
    combinations <- combn(1:group2_number, 2)
    tested <- vector("list", group2_number)
    resultados_de <- vector("list", group2_number)
    results_completed <- vector("list", group2_number)
    combinations_names <- apply(combinations, 2, function(x) {
        paste0("G", x[2], "_over_", "G", x[1])
    })
    names(tested) <- combinations_names
    names(resultados_de) <- combinations_names
    names(results_completed) <- combinations_names
    if (group2_number > 2) {
        message("Performing differential expressed analysis\n")
        for (pairs in seq_len(ncol(combinations))) {
            if (tolower(test) == "lrt") {
                dds <- DESeq2::DESeq(dds,
                    test = "LRT",
                    reduced = ~1, parallel = TRUE
                )
            } else if (tolower(test) == "wald") {
                dds <- DESeq2::DESeq(dds, test = "Wald", parallel = TRUE)
            } else {
                dds <- DESeq2::estimateSizeFactors(dds)
                dds <- DESeq2::estimateDispersions(dds)
                dds <- DESeq2::nbinomWaldTest(dds)
            }

            normalized_expression <- DESeq2::counts(dds, normalized = TRUE)
            var_trans <- DESeq2::varianceStabilizingTransformation(dds,
                blind = TRUE
            )

            deseq_plots(dds)


            res <- DESeq2::results(dds,
                parallel = TRUE, addMLE = FALSE,
                cooks_cutoff = cooks_cutoff,
                contrast = c(
                    "condition",
                    as.character(combinations[2, pairs]),
                    as.character(combinations[1, pairs])
                )
            )
            # , alpha = fdr_cutoff) alfa is the FDR

            if (tolower(image_format) == "png") {
                png(
                    filename = file.path(dir, "MA_plot1.png"),
                    width = width, height = height, res = res, units = unit
                )
            } else if (tolower(image_format) == "svg") {
                svg(
                    filename = file.path(dir, "MA_plot1.svg"),
                    width = width, height = height, onefile = TRUE
                )
            } else {
                tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
                stop(message(tmp))
            }
            DESeq2::plotMA(res,
                main = "DESeq2", ylim = c(
                    -log2(fc_cutoff) - 1,
                    log2(fc_cutoff) + 1
                ),
                las = 1,
                ylab = "Log2FC",
                xlab = "Mean of normalized counts"
            )
            abline(h = log2(fc_cutoff), col = "dodgerblue", lwd = 2)
            abline(h = -log2(fc_cutoff), col = "dodgerblue", lwd = 2)
            dev.off()

            results_compl_local <- as.data.frame(res)

            if (tolower(data_base) == "legacy") {
                gene_symbol <- strsplit(row.names(results_compl_local),
                    split = "\\|"
                )
                gene_symbol <- as.data.frame(gene_symbol)
                gene_symbol <- t(gene_symbol)

                results_compl_local$gene_symbol <- gene_symbol[, 1]
                results_compl_local$geneid <- gene_symbol[, 2]

                colnames(results_compl_local)[c(2, 5, 6)] <- c(
                    "log2FC",
                    "Pvalue", "FDR"
                )

                results_compl_local$fc <- 2** (results_compl_local$log2FC)

                results_compl_local <- results_compl_local[, c(
                    7, 8, 5,
                    6, 9, 2, 1, 3, 4
                )]
            } else {
                colnames(results_compl_local)[c(2, 5, 6)] <- c(
                    "log2FC",
                    "Pvalue", "FDR"
                )

                results_compl_local$fc <- 2** (results_compl_local$log2FC)

                results_compl_local <- results_compl_local[, c(
                    5, 6, 7,
                    2, 1, 3, 4
                )]
            }

            results_completed[[pairs]] <- results_compl_local

            # filtering threashold
            tmp <- results_compl_local$FDR < fdr_cutoff
            res_filtered <- results_compl_local[tmp, ]
            tmp <- abs(res_filtered$log2FC) > log2(fc_cutoff)
            res_filtered <- res_filtered[tmp, ]
            res_filtered <- na.exclude(res_filtered)
            resultados_de[[pairs]] <- as.data.frame(res_filtered)

            write.csv(resultados_de[[pairs]], file = paste0(dir, "/ResultsDE_",
                                                combinations_names, ".csv"))
            write.csv(x = normalized_expression, file = paste0(
                dir,
                "/Normalized_Expression_", tolower(test), ".csv"
            ))

            suppressWarnings(volcano(results_completed, pairs))
        }
    } else {
        message("Performing differential expressed analysis\n")
        if (tolower(test) == "lrt" || tolower(test) == "wald") {
            dds <- DESeq2::DESeq(dds,
                test = toupper(test),
                reduced = ~1, parallel = TRUE
            )
        } else {
            dds <- DESeq2::estimateSizeFactors(dds)
            dds <- DESeq2::estimateDispersions(dds)
            dds <- DESeq2::nbinomWaldTest(dds)
        }

        # save normalized counts
        normalized_expression <- DESeq2::counts(dds, normalized = TRUE)

        # rlogTransformation() stabilizes the variance across the range of
        # counts, so genes have a nearly equal effect on the distances and in
        # the PCA plot for example. the rlog is not biased towards long genes
        # (high count genes). varianceStabilizingTransformation is a faster
        # choice
        var_trans <- DESeq2::varianceStabilizingTransformation(dds,
            blind = TRUE
        )

        deseq_plots(dds)

        res <- DESeq2::results(dds,
            parallel = TRUE, addMLE = FALSE,
            cooks_cutoff = cooks_cutoff,
            contrast = c("condition", "2", "1")
        )

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "MA_plot1.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "MA_plot1.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
            stop(message(tmp))
        }
        # MA plot. Points which fall out of the window are plotted as open
        # triangles pointing either up or down
        DESeq2::plotMA(res,
            main = "DESeq2", ylim = c(
                -log2(fc_cutoff) - 1,
                log2(fc_cutoff) + 1
            ),
            las = 1, ylab = "Log2FC",
            xlab = "Mean of normalized counts"
        )
        abline(h = log2(fc_cutoff), col = "dodgerblue", lwd = 2)
        abline(h = -log2(fc_cutoff), col = "dodgerblue", lwd = 2)
        dev.off()

        results_compl_local <- as.data.frame(res)

        if (tolower(data_base) == "legacy") {
            gene_symbol <- strsplit(row.names(results_compl_local),
                split = "\\|"
            )
            gene_symbol <- as.data.frame(gene_symbol)
            gene_symbol <- t(gene_symbol)

            results_compl_local$gene_symbol <- gene_symbol[, 1]
            results_compl_local$geneid <- gene_symbol[, 2]

            colnames(results_compl_local)[c(2, 5, 6)] <- c(
                "log2FC",
                "Pvalue", "FDR"
            )

            results_compl_local$fc <- 2 ** (results_compl_local$log2FC)
            results_compl_local$fc <- ifelse(results_compl_local$fc < 1,
                (-1 / results_compl_local$fc),
                results_compl_local$fc
            )

            results_compl_local <- results_compl_local[, c(
                7, 8, 5, 6,
                9, 2, 1, 3, 4
            )]
        } else {
            colnames(results_compl_local)[c(2, 5, 6)] <- c(
                "log2FC",
                "Pvalue", "FDR"
            )

            results_compl_local$fc <- 2 ** (results_compl_local$log2FC)

            results_compl_local$fc <- ifelse(results_compl_local$fc < 1,
                (-1 / results_compl_local$fc),
                results_compl_local$fc
            )

            results_compl_local <- results_compl_local[, c(
                5, 6, 7, 2,
                1, 3, 4
            )]
        }

        results_completed[[1]] <- results_compl_local

        # filtering threashold
        tmp <- results_compl_local$FDR < fdr_cutoff
        res_filtered <- results_compl_local[tmp, ]
        tmp <- abs(res_filtered$log2FC) > log2(fc_cutoff)
        res_filtered <- res_filtered[tmp, ]
        res_filtered <- na.exclude(res_filtered)
        resultados_de[[1]] <- as.data.frame(res_filtered)

        write.csv(resultados_de[[1]], file = paste0(
            dir, "/ResultsDE_",
            combinations_names, ".csv"
        ))
        write.csv(x = normalized_expression, file = paste0(
            dir,
            "/Normalized_Expression_",
            tolower(test), ".csv"
        ))

        suppressWarnings(volcano(results_completed[[1]], 0))
    }

    assign("normalized_expression_deseq2",
        normalized_expression,
        envir = get(ev)
    )
    assign("results_compl_deseq2",
        results_completed,
        envir = get(ev)
    )
    assign("resultados_de_deseq2", resultados_de, envir = get(ev))

    assign("tool", "deseq2", envir = get(ev))

    message("Done!")
}
