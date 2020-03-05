#' Run EBSeq gene Differential Expression Analysis (DEA).
#'
#' @param name
#' @param work_dir
#' @param env
#' @param tumor
#' @param group_gen A character string of the groups generation function:
#'    \itemize{\item{"mclust" - }{\code{groups_identification_mclust()};}
#'    \item{"CoxHR" - }{\code{groups_identification_coxHR()};} \item{"clinical"
#'    - }{\code{groups_identification_clinical()}.} }
#' @param clinical_pair A character string containing one of the group pairs
#'    selected after statistical analysis runned in \code{clinical_terms()}
#'    function.
#' @param pair_name A character string indicating which condition name should be
#'    used. When there are only two groups the default is \code{"G2_over_G1"}.
#' @param rounds Numerical value indicating the number of iterations. It is
#'    recommended to check the Alpha and Beta convergence plots in output and
#'    adjust this value until the hyper-parameter estimations converged. The
#'    default is \code{7}.
#' @param norm_type A character string indicating which EBSeq normalization
#'    factors type should be used in the analysis "QuantileNorm" or
#'    "MedianNorm". 'The default is \code{"all"}.
#' @param ebtest_qtrm,ebtest_qtrm_cut Numerical value. It is removed from the
#'    analysis genes with ebtest_qtrm th quantile < = ebtest_qtrm_cut. More
#'    details in EBSeq \link{EBTest} page. The default is \code{ebtest_qtrm =
#'    0.75} and \code{ebtest_qtrm_cut = 10}.
#' @param p_cutoff Numerical value indicating the maximum value for p-values.
#'    The default is \code{0.05}.
#' @param fdr_cutoff Numerical value indicating the maximum value for FDR
#'    values. The default is \code{0.05}.
#' @param fc_cutoff Numerical value indicating the maximum value for Fold
#'    Change '(FC). The default is \code{2}.
#' @param width,height,res,unit,image_format
#' @param bullard_quantile Numerical value indicating the quantile for the
#'    Bullard's normalization. The default is \code{0.75}.
#' @inheritParams concatenate_exon
#' @inheritParams download_gdc
#' @inheritParams groups_identification_mclust
#'
#' @return A matrix with DE genes in row and statistical values in columns.
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
#' # considering concatenate_expression and groups_identification already runned
#' dea_ebseq(
#'     pair_name = "G2_over_G1",
#'     rounds = 7,
#'     name = "HIF3A",
#'     env = CHOL_LEGACY_gene_tumor_data
#' )
dea_ebseq <- function(name, work_dir, env, tumor,
                      group_gen,
                      clinical_pair,
                      pair_name = "G2_over_G1",
                      rounds = 7,
                      norm_type = "QuantileNorm",
                      ebtest_qtrm = 0.75,
                      ebtest_qtrm_cut = 10,
                      p_cutoff = 0.05,
                      fdr_cutoff = 0.05,
                      fc_cutoff = 2,
                      width = 2000,
                      height = 1500,
                      res = 300,
                      unit = "px",
                      image_format = "png",
                      bullard_quantile = 0.75) {

    # local functions ####
    ebseq_plots <- function(test, pairs) {
        comb_name <- ifelse(pairs > 1, names(tested)[pairs], pair_name)

        # QQP
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "QQplot.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "QQplot.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 2))
        EBSeq::QQP(test)
        dev.off()
        # DenNHist
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "DenNHistplot.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "DenNHistplot.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 2))
        EBSeq::DenNHist(test)
        dev.off()

        ###### PlotPostVsRawFC plot
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "PostVsRawFCplot.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "PostVsRawFCplot.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        EBSeq::PlotPostVsRawFC(test, fc)
        dev.off()

        # Get the iters and generate a table    iters = convergência??
        iters <- cbind(test$Alpha, test$Beta, test$P)
        rownames(iters) <- rownames(test$Alpha)
        colnames(iters) <- c("Alpha", "Beta_Ng1", "P")
        write.csv(iters, file = paste0(
            dir, "/iters.results_",
            comb_name, ".csv"
        ))

        # Print the plot
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "Convergence.plot.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "Convergence.plot.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 3))
        plot(test$Alpha,
            frame = FALSE, xlab = "Rounds",
            ylab = "Alpha", cex = 1.1, col = "#FECB92", pch = 16,
            cex.axis = 1.5,
            cex.lab = 1.5
        )
        lines(test$Alpha, col = "#FDB462", lwd = 2)
        plot(test$Beta,
            frame = FALSE, xlab = "Rounds",
            ylab = "Beta", cex = 1.1, col = "#FECB92", pch = 16,
            cex.axis = 1.5,
            cex.lab = 1.5
        )
        lines(test$Beta, col = "#FDB462", lwd = 2)
        plot(test$P,
            frame = FALSE, xlab = "Rounds",
            ylab = "P", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5,
            cex.lab = 1.5
        )
        lines(test$P, col = "#FDB462", lwd = 2)
        dev.off()

        # Generate table with values of PPEE and PPDE
        ebseq_results <- as.data.frame(EBSeq::GetPPMat(test))
        if (tolower(data_base) == "legacy") {
            gene_symbol <- strsplit(row.names(ebseq_results), split = "\\|")
            gene_symbol <- as.data.frame(gene_symbol)
            gene_symbol <- t(gene_symbol)
            ebseq_results[, 2] <- gene_symbol[, 1]
            colnames(ebseq_results)[1:2] <- c("FDR", "gene_symbol")
            ebseq_results$geneid <- gene_symbol[, 2]
            ebseq_results$post_fc <- fc$PostFC
            ebseq_results$log2_fc <- log2(ebseq_results$post_fc)

            ebseq_results$fc <- ifelse(ebseq_results$post_fc < 1,
                -1 / ebseq_results$post_fc,
                ebseq_results$post_fc
            )

            ebseq_results <- ebseq_results[, c(2, 3, 1, 6, 5, 4)]
        } else {
            colnames(ebseq_results)[1] <- c("FDR")
            ebseq_results$post_fc <- fc$PostFC
            ebseq_results$log2_fc <- log2(ebseq_results$post_fc)

            ebseq_results$fc <- ifelse(ebseq_results$post_fc < 1,
                -1 / ebseq_results$post_fc,
                ebseq_results$post_fc
            )
            ebseq_results <- ebseq_results[, c(1, 4, 3)]
        }

        if (exists("results_compl_ebseq", envir = get(ev))) {
            results_completed <- sv[["ev"]]$results_compl_ebseq
            resultados_de <- sv[["ev"]]$resultados_de_EBSeq
        }

        results_completed[[comb_name]] <- ebseq_results
        assign("results_compl_ebseq", results_completed,
            envir = get(ev)
        )

        # for DE
        resultados_de <- ebseq_results[abs(ebseq_results$fc) > fc_cutoff, ]
        resultados_de <- resultados_de[resultados_de$FDR < fdr_cutoff, ]
        tmp <- order(resultados_de$FDR, -resultados_de$fc)
        resultados_de <- resultados_de[tmp, ]
        write.csv(resultados_de, file = paste0(
            dir, "/ResultsDE_",
            comb_name, ".csv"
        ))
        resultados_de[[comb_name]] <- resultados_de
        assign("resultados_de_ebseq", resultados_de, envir = get(ev))
    }

    volcano <- function(results_completed, pairs) {
        comb_name <- pairs
        ebseq_results <- results_completed[[pairs]]

        # START VOLCANO PLOT
        ebseq_results$colour <- hsv(0, 0, .39, 0.5)

        # Set new column values to appropriate colours
        tmp <- ebseq_results$FDR <= fdr_cutoff
        tmp1 <- ebseq_results$log2_fc >= log2(fc_cutoff)
        ebseq_results$colour[tmp1 & tmp] <- hsv(0, .9, .87, .5)
        tmp1 <- ebseq_results$log2_fc <= -log2(fc_cutoff)
        ebseq_results$colour[tmp1 & tmp] <- hsv(.67, .79, .93, .5)

        #### Volcano Plot
        axislimits_x <- ceiling(max(c(
            -min(ebseq_results$log2_fc, na.rm = TRUE) - 1,
            max(ebseq_results$log2_fc, na.rm = TRUE) + 1
        )))

        log_10_fdr <- -log10(ebseq_results$FDR)
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
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(4, 6, 3, 2), mgp = c(2, .7, 0), tck = -0.01)
        plot(ebseq_results$log2_fc, log_10_fdr,
            axes = FALSE,
            xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
            xlab = expression("log"[2] * "(FC)"),
            ylab = "",
            cex.lab = 1.5, cex.main = 2,
            cex.sub = 2,
            pch = 16, col = ebseq_results$colour, cex = 2.5, las = 1
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
        tmp <- ebseq_results$FDR <= fdr_cutoff
        tmp1 <- ebseq_results$log2_fc <= -log2(fc_cutoff)
        text(
            x = -axislimits_x + 0.3,
            y = axislimits_y / 10,
            labels = length(ebseq_results$colour[tmp1 & tmp]),
            cex = 1, col = "blue"
        )
        tmp1 <- ebseq_results$log2_fc >= log2(fc_cutoff)
        text(
            x = axislimits_x - 0.4,
            y = axislimits_y / 10,
            labels = length(ebseq_results$colour[tmp1 & tmp]),
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

    sw <- function(x) {
        suppressWarnings(x)
    }

    # Code ####
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
    ),
    envir = get(ev)
    )

    assign("group_gen", group_gen, envir = get(ev))

    path <- ifelse(
        exists("name_e", envir = get(ev)),
        file.path(
            sv[["ev"]]$path,
            sv[["ev"]]$name_e
        ), sv[["ev"]]$path
    )

    dir.create(path, showWarnings = FALSE)

    # creating the dir to outputs
    dir.create(paste0(
        path, "/EBSeq_Results.", tolower(group_gen), "_",
        toupper(name)
    ),
    showWarnings = FALSE
    )
    dir <- paste0(
        path, "/EBSeq_Results.", tolower(group_gen), "_",
        toupper(name)
    )

    # level[1] over level[2]
    levels <- gsub("[Gg]", "", unlist(strsplit(pair_name, "_over_")),
        perl = TRUE
    )

    # For groups
    if (tolower(group_gen) == "mclust") {
        grupos_ebseq <- as.data.frame(sv[["ev"]]$groups)
        group2_number <- max(grupos_ebseq[, "Selected_classification"])
        grupos_ebseq[, "Selected_classification"] <- factor(
            grupos_ebseq[, "Selected_classification"],
            levels = levels
        )
        colnames(grupos_ebseq) <- c("condition", "type")
        grupos_ebseq[, 2] <- c("paired-end")
        grupos_ebseq <- na.omit(grupos_ebseq)
    } else if (tolower(group_gen) == "clinical") {
        grupos_ebseq <- as.data.frame(
            sv[["ev"]]$clinical_groups_clinical[[clinical_pair]]
        )
        group2_number <- max(grupos_ebseq[, "Selected_classification"])

        grupos_ebseq[, "Selected_classification"] <- factor(
            grupos_ebseq[, "Selected_classification"]
        )

        colnames(grupos_ebseq) <- c("condition", "type")
        grupos_ebseq[, 2] <- c("paired-end")
    } else if (tolower(group_gen) == "coxhr") {
        grupos_ebseq <- as.data.frame(sv[["ev"]]$clinical_groups)
        grupos_ebseq$classification <- gsub(
            "low", "1",
            grupos_ebseq$classification
        )
        grupos_ebseq$classification <- gsub(
            "high", "2",
            grupos_ebseq$classification
        )
        group2_number <- max(as.numeric(grupos_ebseq[, "classification"]))
        grupos_ebseq[, "classification"] <- factor(
            grupos_ebseq[, "classification"]
        )
        colnames(grupos_ebseq) <- c("type", "condition")
        grupos_ebseq[, 1] <- c("paired-end")
        grupos_ebseq <- grupos_ebseq[, c(2, 1)]
    } else {
        tmp <- paste0(
            "Please insert a valid 'group_gen' value!!",
            " ('mclust', 'coxHR' or 'clinical')"
        )
        stop(message(tmp))
    }

    # check patient in common
    tmp1 <- colnames(sv[["ev"]]$gene_tumor_not_normalized)
    tmp <- rownames(grupos_ebseq) %in% tmp1

    grupos_ebseq <- grupos_ebseq[tmp, ]

    assign("cond_heatmap", grupos_ebseq[, 1], envir = get(ev))

    # selecting specifics patients
    completed_matrix <- sv[["ev"]]$gene_tumor_not_normalized[, rownames(
        grupos_ebseq
    )]

    if (tolower(norm_type) == "quantilenorm") {
        sizes <- EBSeq::QuantileNorm(completed_matrix, bullard_quantile)
    } else if (tolower(norm_type) == "mediannorm") {
        sizes <- EBSeq::MedianNorm(completed_matrix)
    } else {
        stop(message("Please insert a valid 'norm_type' argument!!"))
    }

    if (!exists("results_compl_ebseq", envir = get(ev))) {
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
    }

    message("Starting EBSeq analysis...")
    eb_out <- EBSeq::EBTest(
        Data = completed_matrix,
        Conditions = grupos_ebseq[, "condition"],
        # filtrando ruído
        Qtrm = ebtest_qtrm,
        QtrmCut = ebtest_qtrm_cut,
        sizeFactors = sizes,
        maxround = rounds
    )
    # generating Fold Change
    fc <- EBSeq::PostFC(eb_out)

    ebseq_plots(eb_out, 0)
    sw(volcano(sv[["ev"]]$results_compl_ebseq, pair_name))

    # A disadvantage to and serious risk of using fold change (for DE) is that
    # it is biased [2] and may miss differentially expressed genes with large
    # differences (B-A) but small ratios (A/B), leading to a high miss rate at
    # high intensities. source: wiki

    # Generatin' normalized data matrix
    normalized_expression <- EBSeq::GetNormalizedMat(completed_matrix, sizes)
    write.csv(normalized_expression, file = paste0(
        dir,
        "/Normalized_Expression.csv"
    ))

    assign("normalized_expression_ebseq", normalized_expression,
        envir = get(ev)
    )
    assign("tool", "ebseq", envir = get(ev))

    message("Done!")
}
