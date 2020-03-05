#' Search for clinical terms for all available tumors
#'
#' @param name
#' @param work_dir
#' @param tumor A character vector indicating which tumors should be used in
#'    the analysis. The default is "all".
#' @param data_type A character string indicating which data_type will be used
#'    to group the patients. Only needed when the environment was not created
#'    yet and the user does not have the \code{term} argument.
#' @param group_number Numerical value indicating how many groups should be
#'    generated. This argument is advisable if subcategory number is larger than
#'    5. The default is 2.
#' @param term A character string containing the clinical term to be used.
#' @param data_base
#' @param term_keyword A character string containing a possible fragment of the
#'    term. Used when the specific term is unknown.
#' @param tumor_with_term A numerical value indicating the minimum number of
#'    tumors containing a specific term with the keyword used.
#'    The default is 1.
#' @param p_cutoff
#' @param fdr_cutoff
#' @param confidence_level A numerical value containing the confidence interval
#'    to be used. The default is 0.95.
#' @param width,height,res,unit,image_format
#' @param env
#' @param cex_axix_x,cex_axix_y A numerical value giving the amount by which
#'    labels in X and Y axis respectively should be modified relative to the
#'    default size. The default value is 16 for both arguments.
#' @inheritParams concatenate_exon
#' @inheritParams download_gdc
#' @inheritParams dea_ebseq
#' @inheritParams groups_identification_mclust
#'
#' @return the groups generated after using the seleted clinical category and
#'    store it inside the determined environment name.
#' @export
#'
#' @importFrom data.table fread
#' @importFrom forcats fct_collapse
#'
#' @examples
#' library(DOAGDC)
#'
#' for (tumor in c("CHOL", "UCS", "OV")) {
#'     download_gdc(
#'         data_type = "clinical",
#'         data_base = "gdc",
#'         work_dir = "~/Desktop",
#'         tumor = tumor
#'     )
#' }
#'
#' # searching terms having the keyword inserted in common 'tumor_with_term = 3'
#' clinical_terms(
#'     tumor = c("CHOL", "UCS", "OV"),
#'     data_base = "gdc",
#'     term_keyword = "year",
#'     work_dir = "~/Desktop",
#'     tumor_with_term = 3
#' )
#'
#' # searching all terms available with at least three tumors in common
#' clinical_terms(
#'     tumor = c("CHOL", "UCS", "OV"),
#'     term_keyword = "",
#'     data_base = "gdc",
#'     work_dir = "~/Desktop",
#'     tumor_with_term = 3
#' )
#'
#' # using the analysis with a specified term
#'
#' download_gdc(
#'         data_type = "gene",
#'         data_base = "gdc",
#'         htseq = "FPKM",
#'         work_dir = "~/Desktop",
#'         tumor = "CHOL")
#'
#' concatenate_expression(
#'     name = "ENSG00000124440",
#'     data_type = "gene",
#'     data_base = "gdc",
#'     work_dir = "~/Desktop",
#'     htseq = "FPKM",
#'     tumor = "CHOL")
#'
#' clinical_terms(
#'     "ENSG00000124440",
#'     "~/Desktop",
#'     "CHOL",
#'     "gene",
#'     data_base = "gdc",
#'     env = CHOL_GDC_gene_tumor_data,
#'     term = "gender")
clinical_terms <- function(name,
                           work_dir,
                           tumor = "all",
                           data_type,
                           group_number = 2,
                           term,
                           data_base = "legacy",
                           term_keyword = NULL,
                           tumor_with_term = 1,
                           p_cutoff = 0.05,
                           fdr_cutoff = 0.05,
                           confidence_level = 0.95,
                           width = 7,
                           height = 7,
                           res = 300,
                           unit = "in",
                           image_format = "svg",
                           env,
                           cex_axix_x = 16,
                           cex_axix_y = 16) {

    # local functions ####
    extract_terms <- function(only_one_tumor) {
        folder <- ifelse(tolower(data_base) == "legacy", "clinical_data",
            "gdc_clinical_data"
        )

        dir <- file.path(
            work_dir, "DOAGDC", toupper(tumor),
            folder
        )

        count <- 0
        for (tumors in names(lista)) {
            count <- count + 1
            patient_file <- dir(
                path = dir[count],
                pattern = "clinical_patient"
            )
            clinical_df <- read.table(
                file.path(dir[count], patient_file), TRUE, "\t",
                stringsAsFactors = FALSE
            )

            clinical_df <- clinical_df[-c(1, 2), , drop = FALSE]
            rownames(clinical_df) <- as.character(unlist(clinical_df[, 2]))
            clinical_df <- clinical_df[, -2]

            lista[[count]] <- clinical_df
        }

        if (is.null(term_keyword)) {

            term_in_tumors <- apply(clinical_df, 2, table)

            # at least two groups
            terms_unique <- lapply(term_in_tumors, length) >= 2
            terms_with_keyword <- names(term_in_tumors)[terms_unique]
            term_in_tumors <- term_in_tumors[terms_unique]

            terms_polished <- cbind(
                terms_with_keyword,
                unlist(lapply(term_in_tumors, length))
            )
            colnames(terms_polished) <- c("Term", "Subcategories_number")
            terms_polished <- terms_polished[order(terms_polished[, 1]), ]

            message(paste0(
                "\nAll possible terms were saved in ",
                work_dir, " as '",
                "DOAGDC_possible_terms_for_",
                toupper(tumor), "_",
                tolower(data_base), ".txt"
            ))
            write.csv(terms_polished,
                file.path(
                    work_dir,
                    paste0(
                        "DOAGDC_possible",
                        "_terms_for_",
                        toupper(tumor), "_",
                        tolower(data_base),
                        ".csv"
                    )
                ),
                row.names = FALSE, quote = FALSE
            )

            assign("clinical_df", clinical_df, envir = get(envir_link))
        } else {

            # list all terms in the files
            terms_repeted <- unlist(lapply(lista, colnames))

            # list all possible terms
            terms_unique <- table(terms_repeted)[order(table(terms_repeted),
                decreasing = TRUE
            )]

            tmp <- names(terms_unique)[grep(
                pattern = term_keyword,
                x = names(terms_unique),
                ignore.case = TRUE,
                perl = TRUE
            )]
            if (!only_one_tumor) {
                terms_with_keyword <- terms_unique[tmp]
                message(paste("These are the possible terms:",
                    paste0(names(terms_with_keyword),
                        collapse = "\n-"
                    ),
                    sep = "\n-"
                ))
                message(paste0(
                    "\nAll of them were saved in ",
                    work_dir, " as '", "DOAGDC_possible_terms_for_",
                    tolower(term_keyword), "_",
                    tolower(data_base), ".txt"
                ))
                write.table(terms_with_keyword,
                    file.path(
                        work_dir,
                        paste0(
                            "DOAGDC_possible_terms_for_",
                            tolower(term_keyword), "_",
                            tolower(data_base), ".txt"
                        )
                    ),
                    row.names = FALSE
                )
            } else {
                terms_with_keyword <- tmp
                message(paste("These are the possible terms:",
                    paste0(terms_with_keyword,
                        collapse = "\n-"
                    ),
                    sep = "\n-"
                ))
                message(paste0(
                    "\nAll of them were saved in ",
                    work_dir, " as '",
                    "DOAGDC_possible_terms_for_",
                    tolower(term_keyword),
                    toupper(tumor),
                    tolower(data_base), ".txt"
                ))
                write.table(terms_with_keyword,
                    file.path(
                        work_dir,
                        paste0(
                            "DOAGDC_possible_terms_for_",
                            tolower(term_keyword),
                            toupper(tumor),
                            tolower(data_base), ".txt"
                        )
                    ),
                    row.names = FALSE,
                    col.names = "Term(s)"
                )
                assign("terms_with_keyword", terms_with_keyword,
                    envir = get(envir_link)
                )
                assign("clinical_df", lista[[tumor]], envir = get(envir_link))
            }
        }
    }

    # code ####
    if ("all" %in% tolower(tumor)) {
        lista <- vector("list", length(dir(path = file.path(
            work_dir,
            "DOAGDC"
        ))))
        names(lista) <- dir(path = file.path(work_dir, "DOAGDC"))
        extract_terms(only_one_tumor = FALSE)
    } else if (length(tumor) > 1) {
        lista <- vector("list", length(tumor))
        names(lista) <- tumor
        extract_terms(only_one_tumor = FALSE)
    } else {

        lista <- tumor

        data_type <- gsub(" ", "_", data_type)
        name <- gsub("-", "_", name)

        if (!(paste(toupper(tumor), toupper(data_base), gsub(
            " ", "_", data_type), "tumor_data", sep = "_") %in% ls(
            all.names = TRUE, envir = .GlobalEnv))) {
            assign(paste(toupper(tumor), toupper(data_base), gsub(
                " ", "_", data_type), "tumor_data", sep = "_"),
                new.env(parent = emptyenv()), envir = .GlobalEnv)

            envir_link <- paste(toupper(tumor), toupper(data_base),
                gsub(" ", "_", data_type), "tumor_data",
                sep = "_")
        } else {
            envir_link <- ifelse(missing(env), paste(toupper(tumor),
                toupper(data_base), gsub(" ", "_", data_type),
                "tumor_data", sep = "_"), deparse(substitute(env)))
        }

        string_vars <- list(envir_link = get(envir_link))

        dir.create(file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
            showWarnings = FALSE
        )

        path <- ifelse(exists("name_e", envir = get(envir_link)),
            file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses",
                                        string_vars[["envir_link"]]$name_e),
            file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses",
                                                            toupper(name))
        )

        assign("path", path, envir = get(envir_link))

        dir.create(path, showWarnings = FALSE)
        dir.create(paste0(
            path, "/Clinical_Results_", tolower(data_type), "_", toupper(name)
        ),
        showWarnings = FALSE
        )
        dir <- paste0(
            path, "/Clinical_Results_", tolower(data_type), "_", toupper(name)
        )

        lista <- vector("list", 1)
        names(lista) <- tumor

        extract_terms(only_one_tumor = TRUE)

        term_in_tumors <- string_vars[["envir_link"]]$clinical_df[, term,
            drop = FALSE
        ]
        message("Removing 'unusable subcategories...")
        term_in_tumors[, 1] <- gsub(
            "\\[Not Available\\]",
            NA, term_in_tumors[, 1]
        )
        term_in_tumors[, 1] <- gsub(
            "\\[Unknown\\]",
            NA, term_in_tumors[, 1]
        )
        term_in_tumors[, 1] <- gsub(
            "\\[Not Evaluated\\]",
            NA, term_in_tumors[, 1]
        )
        term_in_tumors <- na.exclude(term_in_tumors)
        term_category <- unique(term_in_tumors[, 1])
        category_number <- length(term_category)

        if (category_number < 2) {
            stop(message(
                "This term has less than 2 usable categories and",
                " cannot be used in this",
                " analysis!\n Please",
                " try another 'term'."
            ))
        }

        if (category_number > 5) {
            if (!is.na(as.numeric(term_in_tumors[1, 1]))) {

                # Frequency Tables/Grouped values
                term_in_tumors[, 1] <- as.numeric(term_in_tumors[, 1])

                range <- max(term_in_tumors[, 1]) - min(term_in_tumors[, 1])
                # > 2 and < 6
                groups_number <- group_number
                # You must round up, not off
                class_width <- ceiling(range / (groups_number + 1))
                x_breaks <- seq(
                    min(term_in_tumors[, 1]),
                    max(term_in_tumors[, 1]),
                    class_width
                )
                step_val <- cut(term_in_tumors[, 1],
                    breaks = x_breaks,
                    include.lowest = TRUE
                )

                # include.highest
                tmp <- length(levels(step_val))
                step_val[is.na(step_val)] <- levels(step_val)[tmp]
                levels(step_val)[length(levels(step_val))] <- gsub(
                    "\\,.*?\\]$",
                    paste0(",", max(term_in_tumors[, 1]), "]"),
                    levels(step_val)[length(levels(step_val))]
                )
                term_in_tumors[, 1] <- step_val

                term_category <- levels(step_val)
            } else if (term == "clinical_stage") {
                term_in_tumors[, 1] <- suppressWarnings(
                    forcats::fct_collapse(term_in_tumors[, 1],
                        Stage_I = c(
                            "Stage I", "Stage IA", "Stage IB",
                            "Stage IC"
                        ),
                        Stage_II = c(
                            "Stage II", "Stage IIA", "Stage IIB",
                            "Stage IIC"
                        ),
                        Stage_III = c(
                            "Stage III", "Stage IIIA", "Stage IIIB",
                            "Stage IIIC", "Stage IIIC1", "Stage IIIC2"
                        ),
                        Stage_IV = c(
                            "Stage IV", "Stage IVA", "Stage IVB",
                            "Stage IVC"
                        )
                    )
                )
            }
        }

        list_patients_per_category <- sapply(term_category, function(x) {
            rownames(term_in_tumors)[term_in_tumors[, 1] == x]
        })

        patients_per_category <- term_in_tumors
        patients_per_category <- table(patients_per_category)

        assign("list_patients_per_category", list_patients_per_category,
            envir = get(envir_link)
        )

        tmp <- ifelse(
            tolower(data_type) == "methylation",
            paste0("methylation_tumor_filtered_selected_", name), paste0(
                tolower(data_type), "_tumor_normalized_selected_",
                toupper(name)
            )
        )

        selected <- eval(parse(text = paste0(
            "string_vars[['envir_link']]$", tmp
        )))

        # split patient names
        patient_complete_name <- rownames(selected)
        # get the fisrt 3 elemetns inside all the keys
        patient_partial_name <- lapply(strsplit(
            x = patient_complete_name, split = "-", perl = TRUE), "[", 1:3)
        # collapse them in a vetor
        patient_partial_name <- unlist(lapply(
            patient_partial_name,
            function(x) {
                paste0(x,
                    collapse = "-"
                )
            }
        ))

        # box plot category X expresion
        tmp <- !duplicated(patient_partial_name)
        final_df <- as.data.frame(selected)[tmp, , drop = FALSE]

        tmp <- rownames(final_df) %in% rownames(term_in_tumors)
        final_df <- final_df[tmp, , drop = FALSE]
        final_df[term] <- "NA"

        # stacking the data
        final_df[, 2] <- as.character(term_in_tumors[rownames(final_df), 1])
        final_df[, 2] <- factor(final_df[, 2])

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, paste0(
                    term, "boxplot_Clinical",
                    "_category_by_expression_values.png"
                )),
                width = width, height = height,
                res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, paste0(
                    term, "_boxplot_Clinical",
                    "_category_by_expression_values.svg"
                )),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        final_df_plot <- cbind(final_df, rownames(final_df))
        colnames(final_df_plot) <- c("expression", "category", "patient")
        plot <- ggplot2::ggplot(
            final_df_plot,
            ggplot2::aes(
                x = category,
                y = log2(expression + 1),
                fill = category
            )
        ) +
            ggplot2::stat_boxplot(
                geom = "errorbar",
                position = ggplot2::position_dodge(1)
            ) +
            ggplot2::geom_boxplot(
                position = ggplot2::position_dodge(1),
                inherit.aes = TRUE,
                outlier.colour = NA
            ) +
            ggplot2::geom_jitter(
                position = ggplot2::position_jitter(0.2),
                alpha = 0.5, color = "black"
            ) +
            ggplot2::labs(
                x = "Categories", y = "log2(Expresion + 1)",
                title = term
            ) +
            ggplot2::theme(
                axis.line = ggplot2::element_line(colour = "black"),
                plot.title = ggplot2::element_text(
                    hjust = 0.5,
                    size = ggplot2::rel(1.5)
                ),
                axis.title.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(size = 18),
                axis.text.y = ggplot2::element_text(size = cex_axix_y),
                axis.text.x = ggplot2::element_text(size = cex_axix_x),
                legend.position = "none",
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank()
            )
        print(plot)
        dev.off()

        # anova
        values <- log2(final_df[, 1] + 1)
        categories <- final_df[, 2]
        aovobject <- aov(values ~ categories)
        anova <- anova(aovobject)

        sink(file.path(dir, paste0(
            term, "_ANOVA_Analysis_",
            "of_Variance_report.txt"
        )))
        cat("Used category")
        cat("\n")
        cat(patients_per_category)
        cat("\n\n\n")
        cat("Summarizing the Analysis of Variance Model")
        cat("\n")
        print(summary.aov(aovobject))
        cat("\n")
        print(summary.lm(aovobject))
        cat("\n\n")
        cat("Anova Tables")
        cat("\n")
        print(anova(aovobject))
        cat("\n\n\n")
        # To view the model coefficients, use the coef function:
        cat("Model Coefficients")
        cat("\n")
        print(coef(aovobject))
        cat("\n\n\n")
        # To view confidence intervals for the coefficients, use the
        # confint function:
        cat("Confidence Intervals")
        cat("\n")
        print(confint(aovobject))

        if (anova[["Pr(>F)"]][1] == "NaN") {
            sink()
            stop(message(
                "\nThis category '",
                term,
                "' does not have enough patients",
                " in each subcategory!"
            ))
        }

        if (anova[["Pr(>F)"]][1] <= p_cutoff) {
            cat("\n\n\n")
            cat(paste0(
                "The p-value of ", anova[["Pr(>F)"]][1],
                " is not greater than ", p_cutoff, "."
            ))
            cat("\n")
            cat("Hence pairwise t-tests will be performed.")
            cat("\n")
            cat("Tukey test")
            # tukey
            tukey <- TukeyHSD(aovobject, conf.level = confidence_level)
            print(tukey)
            cat("\n\n\n")
            sink()

            tmp <- which(tukey$categories[, "p adj"] < fdr_cutoff)
            groups_export <- strsplit(names(tmp), "-")

            clinical_groups <- vector("list", length(groups_export))
            names(clinical_groups) <- names(tmp)

            for (index in seq(1, length(groups_export))) {
                selector <- final_df$clinical_stage %in% groups_export[[index]]
                clinical_groups[[index]] <- final_df[selector, 2,
                    drop = FALSE
                ]
                tmp <- droplevels(clinical_groups[[index]][, 1])
                clinical_groups[[index]][, 1] <- tmp
            }

            assign("clinical_groups_clinical", clinical_groups,
                envir = get(envir_link)
            )
        } else {
            cat("\n\n\n")
            # to file
            cat(paste0(
                "\nThe p-value of ", anova[["Pr(>F)"]][1],
                " is greater than ", p_cutoff, "."
            ))
            cat("\nHence, pairwise t-tests will not be performed.")
            cat("\n\n\n")
            sink()
            # to CLI
            message(
                "The p-value of ", anova[["Pr(>F)"]][1],
                " is greater than ", p_cutoff, "."
            )
            message("\nHence pairwise t-tests will not be performed.")
        }
    }
    message("Done!\n")
}
