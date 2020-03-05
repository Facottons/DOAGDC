#' Separate patients in groups
#'
#' \code{groups_identification_mclust} is a function designed to separate
#' patients in groups, powered by mclust.
#'
#' @param data_type
#' @param group_number Numerical value indicating how many groups should be
#'    generated.
#' @param modelName A character string indicating which mclust model
#'    name will be
#'    used. For more details please check \link{mclustModelNames} help file.
#' @param uncertainty_cutoff Numerical value indicating which uncertainty value
#'    for the separation should be tolerated. Patients over this threshold will
#'    be removed from the analysis.
#' @param rerun_plots Logical value where TRUE indicate that the function
#'    should
#'    run the step of group generation using the \code{uncertainty_cutoff}
#'    parameter for filtering the data. The default is \code{FALSE}.
#' @param n_breaks Numerical value giving the number of cells for
#'    the \code{hist} bars. As default \code{n_breaks = 55}.
#' @param width,height,res,unit Graphical parameters. See \link{par} for more
#'    details. As default \code{width = 2000, height = 1500, res =
#'    300 and unit = "px"}.
#' @param image_format A character string indicating which image_format will be
#'    used. It could be "png" or "svg". The only unit available in "svg" is
#'    inches ('in'). The default is "png".
#' @param save_data
#' @param env
#' @param tumor
#' @param data_base
#' @param work_dir
#' @param name
#' @inheritParams download_gdc
#' @inheritParams concatenate_exon
#'
#' @return the groups generated after using the mclust analysis.
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales alpha
#' @importFrom scales percent
#' @importFrom mclust cdensV
#'
#' @examples
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
groups_identification_mclust <- function(data_type,
                                        group_number = "",
                                        modelName = NULL,
                                        uncertainty_cutoff = 0.05,
                                        rerun_plots = FALSE,
                                        n_breaks = 55,
                                        width = 2000,
                                        height = 1500,
                                        res = 300,
                                        unit = "px",
                                        image_format = "png",
                                        save_data = TRUE, env, tumor, data_base,
                                        work_dir = "~/Desktop",
                                        name) {

    # local functions ####
    # plotMix adapted
    plot_local <- function(mc, data, nb_breaks, trace_density = TRUE,
                        title = "", xlim, ylim, ...) {
        if (length(data) <= 5000) {
            shapiro_pval <- shapiro.test(data)$p.value
        }
        p <- seq(min(data), to = max(data), length = 1000)
        d <- mclust::cdens(
            modelName = mc$modelName, data = p,
            parameters = mc$parameters
        )
        title <- paste0(mc$G, "-component Gaussian Mixture Model")

        par(mar = c(5, 5, 2, 5))
        h_freq <- hist(as.numeric(data),
            breaks = nb_breaks, freq = TRUE,
            axes = FALSE, xlab = "", ylab = "", main = ""
        )

        color <- vector("character", length = length(h_freq$mids))
        possible_colors <- RColorBrewer::brewer.pal(8, "Set2")
        last_delimiter <- 0
        classification_number <- length(unique(mc$classification))

        for (groups in seq(classification_number, 1)) {
            delimiter <- min(data[mc$classification == groups])
            if (groups == max(unique(mc$classification))) {
                color[h_freq$mids >= delimiter] <- possible_colors[groups]
            } else if (groups == 1) {
                color[color == ""] <- possible_colors[groups]
            } else {
                tmp1 <- h_freq$mids >= delimiter
                tmp2 <- h_freq$mids <= last_delimiter
                color[tmp1 & tmp2] <- possible_colors[groups]
            }
            last_delimiter <- delimiter
        }

        suppressWarnings(rug(as.numeric(data)))
        axis(4, las = 1)
        mtext("Frequency", side = 4, line = 2.5, cex.lab = 1, las = 3)
        par(new = TRUE)

        h_dens <- hist(as.numeric(data),
            breaks = nb_breaks, col = color, ylab = "",
            freq = FALSE, las = 1,
            xlab = expression("Log"[2] * "(expression+1)"),
            main = (paste(title, "\n", if (length(data) <= 5000) {
                paste("shapiro_pval = ",
                    signif(shapiro_pval, digits = 3),
                    sep = ""
                )
            }))
        )
        mtext("Density", side = 2, line = 2.5, cex.lab = 1, las = 3)

        if (trace_density) {
            par(new = TRUE)
            plot(mfr,
                what = "density", axes = FALSE, col = scales::alpha("red", .5),
                ann = FALSE, main = "", lwd = 3.5
            )
        }

        legend("topright",
            legend = paste0("G", seq(1, mc[["G"]])),
            title = "Groups",
            fill = unique(color), border = FALSE, bty = "n",
            y.intersp = 0.7, cex = 0.7
        )
    }

    # Code ####
    data_type <- gsub(" ", "_", data_type)
    name <- gsub("-", "_", name)

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    work_dir <- ifelse(
        missing("work_dir"), string_vars[["envir_link"]]$work_dir, work_dir
    )

    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    path <- file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses")
    assign("path", path, envir = get(envir_link))

    if (exists("name_e", envir = get(envir_link))) {
        dir.create(
            path = file.path(path, string_vars[["envir_link"]]$name_e),
            showWarnings = FALSE
        )
        dir.create(
            path = file.path(
                path, string_vars[["envir_link"]]$name_e,
                paste0(
                    "Mclust_Results_", tolower(data_type),
                    "_", toupper(name)
                )
            ),
            showWarnings = FALSE
        )
        dir <- file.path(
            path, string_vars[["envir_link"]]$name_e,
            paste0(
                "Mclust_Results_", tolower(data_type),
                "_", toupper(name)
            )
        )
    } else {
        dir.create(
            path = file.path(path, paste0(
                "Mclust_Results_",
                tolower(data_type), "_",
                toupper(name)
            )),
            showWarnings = FALSE
        )
        dir <- file.path(path, paste0(
            "Mclust_Results_", tolower(data_type),
            "_", toupper(name)
        ))
    }

    filtered_results <- eval(parse(text = paste0(
        "string_vars[['envir_link']]$",
        ifelse(
            tolower(data_type) == "methylation", "methylation_",
            strsplit(tolower(data_type), "_")[[1]][1]
        ),
        "_tumor_normalized_selected_",
        toupper(name)
    )))

    if (!rerun_plots) {
        unique_col_row <- unique(rownames(filtered_results))
        filtered_results <- filtered_results[unique_col_row, , drop = FALSE]
        # starting
        message("Running Model-Based Clustering...")

        # could not find function "mclustBIC" FIX
        mclustBIC <- mclust::mclustBIC

        if (is.numeric(group_number) && is.character(modelName)) {
            mfr <- mclust::Mclust(log2(as.numeric(filtered_results[, 1]) + 1),
                G = group_number, modelNames = modelName
            )
        } else if (is.character(modelName)) {
            mfr <- mclust::Mclust(log2(as.numeric(filtered_results[, 1]) + 1),
                modelNames = modelName
            )
        } else if (is.numeric(group_number)) {
            mfr <- mclust::Mclust(log2(as.numeric(filtered_results[, 1]) + 1),
                G = group_number
            )
        } else if (group_number == "" && is.null(modelName)) {
            mfr <- mclust::Mclust(log2(as.numeric(filtered_results[, 1]) + 1))
        }

        icl <- mclust::mclustICL(log2(as.numeric(filtered_results[, 1]) + 1))
        bic <- mclust::mclustBIC(log2(as.numeric(filtered_results[, 1]) + 1))

        if (mfr$G > 1) {
            lrt <- mclust::mclustBootstrapLRT(
                log2(as.numeric(filtered_results[, 1]) + 1),
                modelName = mfr[["modelName"]]
            )
        } else {
            stop(message(
                "Using ", name, " does not generate two or more groups.",
                " Please, try another 'name'!"
            ))
        }
        assign("mclust_result_final", mfr,
            envir = get(envir_link)
        )

        # could not find function "cdensE" FIX
        cdensE <- mclust::cdensE

        # PLots
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "mixture_log2_expression_%01d.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "mixture_log2_expression_%01d.svg"),
                width = width, height = height, onefile = FALSE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format!",
                " ('png' or 'svg')"
            ))
        }
        par(mar = c(5, 6, 4, 3.5))
        plot_local(
            mc = mfr,
            data = log2(as.numeric(filtered_results[, 1]) + 1),
            nb_breaks = n_breaks, las = 1,
            trace_density = TRUE,
            cex.main = 1.4,
            cex.lab = 1,
            cex.axis = 1.1
        )
        plot(bic, las = 1)
        plot(icl, las = 1)
        plot(mfr,
            what = "classification",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1
        )
        plot(mfr,
            what = "density",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1,
            ylab = "Density"
        )
        plot(mfr,
            what = "uncertain",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1
        )
        dev.off()

        # preparando output final
        groups <- as.data.frame(matrix(
            ncol = 2,
            nrow = length(mfr$classification)
        ))
        colnames(groups) <- c("Selected_classification", "Selected_uncertainty")
        rownames(groups) <- rownames(filtered_results)
        groups[, 1] <- mfr$classification
        groups[, 2] <- mfr$uncertainty
        write.csv(groups,
            file = file.path(dir, "groups_classifition.csv"),
            quote = TRUE, row.names = TRUE
        )
        assign("groups", groups, envir = get(envir_link))

        # Uncertainty possible cutoffs1
        uncertainty_range <- seq(0, 1, 0.0001)

        remaining_patients <- sapply(
            uncertainty_range,
            function(pnow) {
                length(which(groups[, "Selected_uncertainty"] < pnow))
            }
        )

        tmp <- uncertainty_range == 0.05
        remaning_example_cutoff <- remaining_patients[which(tmp)]

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "uncertainty_possible_cutoffs1.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "uncertainty_possible_cutoffs1.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format!",
                " ('png' or 'svg')"
            ))
        }
        plot(uncertainty_range, remaining_patients,
            main = "Uncertainty possible cutoffs",
            xlab = "Uncertainty cutoff", ylab = "Remaining Patients", lwd = 2,
            axes = FALSE, cex.lab = 1.1, col = "green4", type = "l"
        )

        axis(
            side = 2, seq(
                min(remaining_patients),
                (max(remaining_patients) + 30), 60
            ),
            lwd = 1, las = 1, cex.axis = 1.2
        )
        axis(
            side = 1, seq(0, 1, 0.05), lwd = 1, las = 1,
            cex.axis = 1.2
        )

        # example cutoff = 0.05
        abline(v = 0.05, lwd = 2)
        text(
            x = 0.05 + 0.01, y = 10, labels = "0.05",
            cex = 1.2, pos = 4, srt = 90
        )

        tmp <- remaining_patients == max(remaining_patients)
        abline(v = uncertainty_range[which(tmp)[1]], lwd = 2)
        text(
            x = 1.008 * uncertainty_range[which(tmp)[1]], y = 10,
            labels = paste0(
                "Plateau cut = ",
                uncertainty_range[which(tmp)[1]]
            ),
            cex = 1.2, pos = 4, srt = 90
        )

        abline(h = remaning_example_cutoff, lwd = 2, lty = 2)
        text(
            x = 0.9, y = remaning_example_cutoff,
            labels = paste0(
                remaning_example_cutoff, " (",
                scales::percent(
                    remaning_example_cutoff / max(remaining_patients)
                ), ")"
            ),
            cex = 1.2, pos = 1, srt = 0
        )

        dev.off()

        # Uncertainty possible cutoffs2
        uncertainty_range <- seq(0, 1, 0.01)

        remaining_patients <- sapply(
            uncertainty_range,
            function(pnow) {
                length(which(groups[, "Selected_uncertainty"] < pnow))
            }
        )

        to_plot <- data.frame(
            xaxis = uncertainty_range[seq(2, length(uncertainty_range))],
            yaxix = sapply(
                2:length(remaining_patients),
                function(n) {
                    remaining_patients[n] - remaining_patients[n - 1]
                }
            )
        )

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "uncertainty_possible_cutoffs2.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "uncertainty_possible_cutoffs2.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
            stop(message(tmp))
        }
        plot(
            x = to_plot$xaxis,
            y = to_plot$yaxix,
            ylim = c(0, 10), axes = FALSE, cex.lab = 1.1, type = "l",
            col = "blue2", lwd = 3, main = "Uncertainty possible cutoffs",
            xlab = "Uncertainty cutoff", ylab = "\U0394Patients"
        )

        legend("topright",
            bty = "n", legend = "Cutoff = 0.05",
            lwd = 1, lty = 3, cex = 0.9
        )

        axis(side = 2, lwd = 1, las = 1, cex.axis = 1.5)
        axis(side = 1, lwd = 1, las = 1, cex.axis = 1.5)

        # example cutoff  =  0.05
        abline(v = 0.05, lwd = 1, lty = 3)

        dev.off()

        assign("mclust_group_number", group_number, envir = get(envir_link))
        sink(paste0(dir, "/mclust_summary.txt"))
        print(summary(mfr))
        cat("\n\n\n")
        print(summary(bic))
        cat("\n\n\n")
        print(summary(icl))
        if (mfr$G > 1) {
            cat("\n\n\n")
            print(lrt)
        }
        sink()
    } else {
        before <- nrow(string_vars[["envir_link"]]$groups)
        groups <- subset(
            x = string_vars[["envir_link"]]$groups,
            string_vars[["envir_link"]]$groups[, 2] < uncertainty_cutoff
        )
        after <- nrow(groups)
        message(paste0("It was removed ", before - after, " patients!"))

        filtered_results_filtered <- filtered_results[rownames(groups), ]

        # could not find function "mclustBIC" FIX
        mclustBIC <- mclust::mclustBIC

        mclust_result <- string_vars[["envir_link"]]$mfr

        group_number <- mclust_result[["G"]]
        modelName <- mclust_result[["modelName"]]
        mfr <- mclust::Mclust(log2(filtered_results_filtered + 1),
            G = group_number, modelNames = modelName
        )

        # could not find function "cdensE" FIX
        cdensE <- mclust::cdensE

        if (tolower(image_format) == "png") {
            png(
                filename = paste0(
                    dir,
                    "/mixture_log2_expression_", uncertainty_cutoff,
                    "_filtered_%01d.png"
                ),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(
                    dir,
                    "/mixture_log2_expression_", uncertainty_cutoff,
                    "_filtered_%01d.svg"
                ),
                width = width, height = height, onefile = FALSE
            )
        } else {
            tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
            stop(message(tmp))
        }
        par(mar = c(5, 6, 4, 3.5))
        plot_local(
            mc = mfr,
            data = log2(filtered_results_filtered + 1),
            nb_breaks = n_breaks, las = 1,
            trace_density = TRUE,
            cex.main = 1.4,
            cex.lab = 1,
            cex.axis = 1.1
        )
        plot(mfr,
            what = "BIC",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1
        )
        plot(mfr,
            what = "classification",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1
        )
        plot(mfr,
            what = "density",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1,
            ylab = "Density"
        )
        plot(mfr,
            what = "uncertain",
            xlab = expression("Log"[2] * "(expression+1)"), las = 1
        )

        dev.off()

        groups[, 1] <- mfr$classification
        groups[, 2] <- mfr$uncertainty

        if (save_data) {
            # saving
            write.csv(groups,
                file = paste0(
                    dir, "/groups_classifition_",
                    uncertainty_cutoff,
                    "_filtered.csv"
                ),
                quote = TRUE, row.names = TRUE
            )

            assign("groups", groups, envir = get(envir_link))
        }
        assign("mclust_result_final_filtered", mfr,
            envir = get(envir_link)
        )
    }

    message("Done!")
}
