#' Run edgeR gene Differential Expression Analysis (DEA).
#'
#' @param name
#' @param method A character string indicating which method should be used:
#'    \code{"exacttest"} or \code{"glmlrt"}. The default is \code{"exacttest"}.
#' @param clinical_pair
#' @param group_gen
#' @param fc_cutoff
#' @param work_dir
#' @param env
#' @param fdr_cutoff
#' @param width,height,res,unit,image_format
#' @inheritParams concatenate_exon
#' @inheritParams dea_ebseq
#' @inheritParams download_gdc
#' @inheritParams concatenate_exon
#' @inheritParams groups_identification_mclust
#'
#' @return A matrix with DE genes in row and statistical values in columns.
#'
#' @import edgeR
#' @export
#'
#' @examples
#' library(DOAGDC)
#'
#' # data already downloaded using the 'download_gdc' function
#' concatenate_expression("gene",
#'    name = "HIF3A",
#'    data_base = "legacy",
#'    tumor = "CHOL",
#'    work_dir = "~/Desktop"
#' )
#'
#' # separating gene HIF3A expression data patients in two groups
#' groups_identification_mclust("gene", 2,
#'    name = "HIF3A",
#'    modelName = "E",
#'    env = CHOL_LEGACY_gene_tumor_data,
#'    tumor = "CHOL"
#' )
#'
#' # load not normalized data
#' concatenate_expression("gene",
#'    normalization = FALSE,
#'    name = "HIF3A",
#'    data_base = "legacy",
#'    tumor = "CHOL",
#'    env = CHOL_LEGACY_gene_tumor_data,
#'    work_dir = "~/Desktop"
#' )
#'
#' # start DE analysis
#' dea_edger(
#'    name = "HIF3A",
#'    group_gen = "mclust",
#'    env = CHOL_LEGACY_gene_tumor_data
#' )
dea_edger <- function(name,
                    method = "exacttest",
                    clinical_pair,
                    group_gen,
                    fc_cutoff = 2,
                    work_dir, env,
                    fdr_cutoff = 0.05,
                    width = 2000,
                    height = 1500,
                    res = 300,
                    unit = "px",
                    image_format = "png") {

    # local functions ####
    sw <- function(x) {
        suppressWarnings(x)
    }

    exact_groups_fix <- function(tested, pairs) {
        if (pairs > 0) {
            comb_name <- names(tested)[pairs]
            tested_local <- tested[[pairs]]
        } else {
            comb_name <- "G2_over_G1"
            tested_local <- tested
        }

        table_de <- edgeR::topTags(tested_local,
            n = nrow(tested_local$table),
            adjust.method = "BH", sort.by = "PValue"
        )$table

        table_de$fc <- 2**table_de[, "logFC"]
        table_de <- table_de[, c(5, 1, 2, 3, 4)]
        colnames(table_de)[c(2, 5)] <- c("log2FC", "fdr")

        if (tolower(data_base) == "legacy") {
            gene_symbol <- strsplit(row.names(table_de), split = "\\|")
            gene_symbol <- as.data.frame(gene_symbol)
            gene_symbol <- t(gene_symbol)

            table_de[, "gene_symbol"] <- gene_symbol[, 1]
            table_de[, "geneid"] <- gene_symbol[, 2]
            table_de <- table_de[, c(6, 7, 4, 5, 1, 2, 3)]
        } else {
            table_de <- table_de[, c(4, 5, 1, 2, 3)]
        }

        results_completed_local <- table_de
        table_de <- table_de[table_de$fdr < fdr_cutoff, ]
        table_de <- table_de[abs(table_de$log2FC) > log2(fc_cutoff), ]

        write.csv(
            x = table_de, file = paste0(
                dir, "/ResultsDE_exactTest_",
                comb_name, ".csv"
            ),
            row.names = TRUE
        )

        if (pairs > 0) {
            sv[["ev"]]$resultados_de_edger[[pairs]] <- table_de
            sv[["ev"]]$results_compl_edger[[pairs]] <- results_completed_local
        } else {
            results_completed <- vector("list", 1)
            resultados_de <- vector("list", 1)
            names(results_completed) <- comb_name
            names(resultados_de) <- comb_name

            results_completed[[1]] <- results_completed_local
            resultados_de[[1]] <- table_de
            assign("results_compl_edger", results_completed,
                envir = get(ev)
            )
            assign("resultados_de_edger", resultados_de,
                envir = get(ev)
            )
        }

        # plots
        message("Plotting...")
        de <- edgeR::decideTestsDGE(tested_local)
        detags <- rownames(dge)[as.logical(de)]
        assign(paste0("detags_", comb_name), detags, envir = get(ev))

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, paste0(
                    "exactTest_",
                    comb_name, ".png"
                )),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, paste0(
                    "exactTest_",
                    comb_name, ".svg"
                )),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format!",
                " ('png' or 'svg')"
            ))
        }
        par(mar = c(5.1, 4.1, 2, 2.1))
        edgeR::plotSmear(tested_local,
            de.tags = detags, las = 1,
            ylab = "log2FC"
        )
        abline(h = c(-log2(fc_cutoff), log2(fc_cutoff)), col = "blue", cex = 2)
        legend("topright",
            border = FALSE, bty = "n",
            legend = c(
                "p.adj > fdr_cutoff", "p.adj < fdr_cutoff",
                "log2FC_cutoff"
            ),
            lty = c(0, 0, 1), pch = c(20, 20, -1), cex = 0.8, lwd = c(1, 1, 2),
            col = c("Black", "Red", "blue")
        )
        dev.off()
    }

    glmlrt_groups_fix <- function(tested, pairs) {
        if (pairs > 0) {
            comb_name <- names(tested)[pairs]
            tested_local <- tested[[pairs]]
        } else {
            comb_name <- "G2_over_G1"
            tested_local <- tested
        }

        table_de <- edgeR::topTags(tested_local,
            n = nrow(tested_local$table),
            adjust.method = "BH", sort.by = "PValue"
        )$table
        colnames(table_de)[c(1, 2, 5)] <- c("log2FC", "log2CPM", "fdr")
        table_de[, "fc"] <- 2**table_de[, "log2FC"]
        table_de <- table_de[order(rownames(table_de)), ]

        if (tolower(data_base) == "legacy") {
            gene_symbol <- strsplit(row.names(table_de), split = "\\|")
            gene_symbol <- as.data.frame(gene_symbol)
            gene_symbol <- t(gene_symbol)

            table_de[, "gene_symbol"] <- gene_symbol[, 1]
            table_de[, "geneid"] <- gene_symbol[, 2]
            table_de <- table_de[, c(7, 8, 4, 5, 6, 1, 3, 2)]
        } else {
            table_de <- table_de[, c(4, 5, 1, 2, 3)]
        }

        results_completed_local <- table_de
        table_de <- table_de[table_de$fdr < fdr_cutoff, ]
        table_de <- table_de[abs(table_de$log2FC) > log2(fc_cutoff), ]
        write.csv(
            x = table_de, file = paste0(
                dir, "/ResultsDE_glmLRT_",
                comb_name, ".csv"
            ),
            row.names = FALSE
        )

        if (pairs > 0) {
            sv[["ev"]]$resultados_de_edger[[pairs]] <- table_de
            sv[["ev"]]$results_compl_edger[[pairs]] <- results_completed_local
        } else {
            results_completed <- vector("list", 1)
            resultados_de <- vector("list", 1)
            names(results_completed) <- comb_name
            names(resultados_de) <- comb_name

            results_completed[[1]] <- results_completed_local
            resultados_de[[1]] <- table_de
            assign("results_compl_edger", results_completed,
                envir = get(ev)
            )
            assign("resultados_de_edger", resultados_de,
                envir = get(ev)
            )
        }

        # plots
        message("Plotting...")
        de <- edgeR::decideTestsDGE(tested_local)
        detags <- rownames(dge)[as.logical(de)]
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, paste0(
                    "glmLRT_",
                    comb_name, ".png"
                )),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, paste0(
                    "glmLRT_",
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
        edgeR::plotSmear(tested_local,
            de.tags = detags, las = 1,
            ylab = "log2FC"
        )
        abline(h = c(-log2(fc_cutoff), log2(fc_cutoff)), col = "blue", cex = 2)
        legend("topright",
            border = FALSE, bty = "n",
            legend = c(
                "p.adj > fdr_cutoff", "p.adj < fdr_cutoff",
                "logfc_cutoff"
            ),
            lty = c(0, 0, 1), pch = c(20, 20, -1), cex = 0.8, lwd = c(1, 1, 2),
            col = c("Black", "Red", "blue")
        )
        dev.off()
    }

    volcano <- function(results_completed, pairs) {
        if (pairs > 0) {
            comb_name <- names(results_completed)[pairs]
            results_completed_local <- results_completed[[pairs]]
        } else {
            comb_name <- "G2_over_G1"
            results_completed_local <- results_completed[[1]]
        }

        # START VOLCANO PLOT
        results_completed_local$colour <- hsv(0, 0, .39, 0.5)

        # Set new column values to appropriate colours
        tmp1 <- results_completed_local$log2FC >= log2(fc_cutoff)
        tmp2 <- results_completed_local$fdr <= fdr_cutoff
        results_completed_local$colour[tmp1 & tmp2] <- hsv(0, .9, .87, .5)

        tmp1 <- results_completed_local$log2FC < -log2(fc_cutoff)
        results_completed_local$colour[tmp1 & tmp2] <- hsv(.67, .79, .93, .5)

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
            tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
            stop(message(tmp))
        }
        par(mar = c(4, 6, 3, 2), mgp = c(2, .7, 0), tck = -0.01)
        plot(results_completed_local$log2FC, log_10_fdr,
            axes = FALSE,
            xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
            xlab = expression("log"[2] * "(FC)"),
            ylab = "",
            cex.lab = 1.5, cex.main = 2,
            cex.sub = 2,
            pch = 16, col = results_completed_local$colour, cex = 2.5, las = 1
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
                length(results_completed_local$colour[tmp1 & tmp2]),
            cex = 1, col = "blue"
        )

        tmp1 <- results_completed_local$log2FC >= log2(fc_cutoff)
        text(
            x = axislimits_x - 0.4, y = axislimits_y / 10,
            labels = length(results_completed_local$colour[tmp1 & tmp2]),
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
                    "\u00B1", " log", "\u2082", "(", fc_cutoff, ") = ",
                    log2(fc_cutoff)
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
    if (missing(env)) {
        tmp <- paste0(
            "The 'env' argument is missing, please insert",
            "the 'env' name and try again!"
        )
        stop(message(tmp))
    }

    ev <- deparse(substitute(env))

    sv <- list(ev = get(ev))

    name <- gsub("-", "_", name)

    if (missing("work_dir")) {
        work_dir <- sv[["ev"]]$work_dir
    }

    data_base <- sv[["ev"]]$data_base

    assign("path", file.path(
        work_dir, "DOAGDC",
        toupper(sv[["ev"]]$tumor), "Analyses"
    ),
    envir = get(ev)
    )

    assign("group_gen", group_gen, envir = get(ev))

    if (exists("name_e", envir = get(ev))) {
        path <- file.path(
            sv[["ev"]]$path,
            sv[["ev"]]$name_e
        )
        dir.create(path, showWarnings = FALSE)
    } else {
        path <- sv[["ev"]]$path
    }

    # creating the dir to outputs
    dir.create(paste0(
        path, "/edgeR_Results.",
        tolower(group_gen), "_", toupper(name)
    ),
    showWarnings = FALSE
    )
    dir <- paste0(
        path, "/edgeR_Results.", tolower(group_gen), "_",
        toupper(name)
    )

    # for the groups
    if (tolower(group_gen) == "mclust") {
        # groups from mclust
        grupos_edger <- as.data.frame(sv[["ev"]]$groups)
        grupos_edger[, "Selected_classification"] <- factor(
            grupos_edger[, "Selected_classification"]
        )
        colnames(grupos_edger) <- c("conditions", "type")
        grupos_edger[, 2] <- c("paired-end")
    } else if (tolower(group_gen) == "clinical") {
        stop()

        grupos_edger <- as.data.frame(
            sv[["ev"]]$clinical_groups_clinical[[clinical_pair]]
        )
        grupos_edger[, "Selected_classification"] <- factor(
            grupos_edger[, "Selected_classification"]
        )
        colnames(grupos_edger) <- c("conditions", "type")
        grupos_edger[, 2] <- c("paired-end")
    } else if (tolower(group_gen) == "coxhr") {
        grupos_edger <- as.data.frame(sv[["ev"]]$clinical_groups)
        grupos_edger$classification <- gsub(
            "low", "1",
            grupos_edger$classification
        )
        grupos_edger$classification <- gsub(
            "high", "2",
            grupos_edger$classification
        )
        grupos_edger[, "classification"] <- factor(
            grupos_edger[, "classification"]
        )
        colnames(grupos_edger) <- c("type", "conditions")
        grupos_edger[, 1] <- c("paired-end")
        grupos_edger <- grupos_edger[, c(2, 1)]
    } else {
        stop(message(
            "Please insert a valid 'group_gen' value!!",
            " ('mlcust', coxHR or 'clinical')"
        ))
    }

    # check patient in common
    tmp1 <- rownames(grupos_edger)
    tmp2 <- colnames(sv[["ev"]]$gene_tumor_not_normalized)

    grupos_edger <- grupos_edger[tmp1 %in% tmp2, ]

    assign("cond_heatmap", grupos_edger[, 1], envir = get(ev))

    # selecting specifics patients
    completed_matrix <- sv[["ev"]]$gene_tumor_not_normalized[, rownames(
        grupos_edger
    )]

    message("Filtering data ...")
    # remove residual data exclude rows with less than 10
    filter_rows <- rowSums(completed_matrix) >= 10
    completed_matrix <- completed_matrix[filter_rows, ]
    message("It was filtered = ", as.numeric(table(filter_rows)[1]), " genes!")

    # round numbers to whole counts
    completed_matrix <- round(completed_matrix, 0)

    # Cria o objeto DGEList
    dge <- edgeR::DGEList(
        counts = completed_matrix,
        group = grupos_edger[, 1]
    )

    # Tdge$countso apply TMM normalization, it is convenient to create a DGEList
    # Robinson, M.D. and Oshlack, A. (2010). A scaling normalization method for
    # differential expres- sion analysis of RNA-seq data. Genome Biology 11,
    # R25. The calcNormFactors function normalizes for RNA composition by
    # finding a set of scaling factors for the library sizes that minimize the
    # log-fold changes between the samples for most genes. The default method
    # for computing these scale factors uses a trimmed mean of Mvalues (TMM)
    # between each pair of samples [26]. We call the product of the original
    # library size and the scaling factor the effective library size. The
    # effective library size replaces the original library size in all downsteam
    # analyses. TMM normalization is applied to this dataset to account for
    # compositional difference between the libraries.

    # binomial models are fitted and dispersion estimates are obtained, we can
    # proceed with testing procedures for determining differential expression
    # using the exact test. The GLM likelihood ratio test is based on the idea
    # of fitting negative binomial GLMs with the Cox-Reid dispersion estimates.
    design <- model.matrix(~ grupos_edger[, 1])
    colnames(design) <- paste0("G", seq(1, max(as.numeric(grupos_edger[, 1]))))
    rownames(design) <- colnames(dge$counts)

    # PS: CPM = counts-per-million The logCPM values can optionally be converted
    # to RPKM or FPKM by subtracting log2 of gene length, see rpkm()
    if (tolower(method) == "exacttest") {
        dge <- edgeR::calcNormFactors(dge)

        normalized_expression_tmm <- dge$counts

        # For general experiments (with multiple factors), edgeR uses the
        # Cox-Reid profile-adjusted likelihood (CR) method in estimating
        # dispersions To estimate common dispersion, trended dispersions and
        # tagwise dispersions in one run:
        message("Estimating common dispersion...")
        dge <- edgeR::estimateCommonDisp(dge)
        group2_number <- max(as.numeric(levels(grupos_edger[, 1])))
        if (group2_number > 2) {
            message(
                "There are more than two group ",
                "combinations, this may take a while...\n"
            )

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

            assign("results_compl_edger", results_completed,
                envir = get(ev)
            )
            assign("resultados_de_edger", resultados_de,
                envir = get(ev)
            )

            count <- 0
            pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
            for (pairs in seq(1, ncol(combinations))) {
                count <- count + 1
                tested[[pairs]] <- edgeR::exactTest(dge,
                    pair = combinations[, pairs]
                )
                exact_groups_fix(tested, pairs)
                sw(volcano(
                    sv[["ev"]]$results_compl_edger,
                    pairs
                ))

                setTxtProgressBar(pb, count)
            }
            close(pb)
        } else {
            tested <- edgeR::exactTest(dge)
            exact_groups_fix(tested, 0)
            sw(volcano(sv[["ev"]]$results_compl_edger, 0))
        }
    } else if (tolower(method) == "glmlrt") {
        dge <- edgeR::calcNormFactors(dge)

        normalized_expression_tmm <- dge$counts

        message("Estimating common dispersion...")
        dge <- edgeR::estimateGLMCommonDisp(dge, design)
        adge_list <- edgeR::estimateGLMTagwiseDisp(dge, design)
        aglm_fit <- edgeR::glmFit(adge_list, design,
            dispersion = adge_list$tagwise.dispersion,
            prior.count.total = 0
        )
        # To perform likelihood ratio tests: fits the negative binomial GLM for
        # each tag and produces an object of class DGEGLM with some new
        # components
        group2_number <- max(as.numeric(levels(grupos_edger[, 1])))
        if (group2_number > 2) {
            message(
                "There are more than two group combinations, ",
                "this may take a while...\n"
            )

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

            assign("results_compl_edger", results_completed,
                envir = get(ev)
            )
            assign("resultados_de_edger", resultados_de,
                envir = get(ev)
            )

            count <- 0
            pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
            for (pairs in seq(1, ncol(combinations))) {
                if (combinations[1, pairs] == 1) {
                    tested[[pairs]] <- edgeR::glmLRT(aglm_fit, coef = pairs + 1)
                } else {
                    count <- count + 1
                    vec <- numeric(ncol(combinations))
                    vec[combinations[1, pairs]] <- -1
                    vec[combinations[2, pairs]] <- 1

                    tested[[pairs]] <- edgeR::glmLRT(aglm_fit, contrast = vec)
                }

                glmlrt_groups_fix(tested, pairs)
                sw(volcano(
                    sv[["ev"]]$results_compl_edger,
                    pairs
                ))

                setTxtProgressBar(pb, count)
            }
            close(pb)
        } else {
            tested <- edgeR::glmLRT(aglm_fit, coef = 2)
            glmlrt_groups_fix(tested, 0)
            sw(volcano(sv[["ev"]]$results_compl_edger, 0))
        }
    } else {
        stop("Please insert a valid method!!! ('exactTest' or 'glmLRT')")
    }

    write.csv(
        x = normalized_expression_tmm,
        file = file.path(dir, "TMM_Normalized_Expression.csv")
    )
    assign("normalized_expression_edger", normalized_expression_tmm,
        envir = get(ev)
    )

    assign("tool", "edger", envir = get(ev))

    message("Done!\n")
}
