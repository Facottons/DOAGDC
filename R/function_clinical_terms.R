#' Search for clinical terms for all available tumors
#'
#' @param Name
#' @param workDir
#' @param tumor A character vector indicating which tumors should be used in
#'    the analysis. The default is "all".
#' @param dataType A character string indicating which dataType will be used to
#'    group the patients. Only needed when the environment was not created yet
#'    and the user does not have the \code{term} argument.
#' @param group.number Numerical value indicating how many groups should be
#'    generated. This argument is advisable if subcategory number is larger than
#'    5. The default is 2.
#' @param term A character string containing the clinical term to be used.
#' @param dataBase
#' @param term_keyword A character string containing a possible fragment of the
#'    term. Used when the specific term is unknown.
#' @param tumor_with_term A numerical value indicating the minimum number of
#'    tumors containing a specific term with the keyword used. The default is 1.
#' @param p_cutoff
#' @param FDR_cutoff
#' @param confidence.level A numerical value containing the confidence interval
#'    to be used. The default is 0.95.
#' @param Width,Height,Res,Unit,image_format
#' @param env
#' @param cexAxixX,cexAxixY A numerical value giving the amount by which labels
#'    in X and Y axis respectively should be modified relative to the default
#'    size. The default value is 16 for both arguments.
#' @inheritParams concatenate_files
#' @inheritParams dea_EBSeq
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
#' \dontrun{
#' #searching for terms having the keyword inserted in common 'tumor_with_term = 3'
#' clinical_terms(tumor = c("BRCA", "UCS", "OV"), term_keyword = "year", tumor_with_term = 3)
#'
#' #searching for all terms available with at least three tumors in common 'tumor_with_term = 3'
#' clinical_terms(tumor = c("BRCA", "UCS", "OV"), term_keyword = "", tumor_with_term = 3)
#'
#' #using the analysis with a specified term
#' clinical_terms("UCS", "isoform", term = "vital_status")
#' }
clinical_terms <- function(Name,
                        workDir,
                        tumor = "all",
                        dataType,
                        group.number = 2,
                        term,
                        dataBase = "legacy",
                        term_keyword = NULL,
                        tumor_with_term = 1,
                        p_cutoff = 0.05,
                        confidence.level = 0.95,
                        Width = 7,
                        Height = 7,
                        Res = 300,
                        Unit = "in",
                        image_format = "svg",
                        env,
                        cexAxixX = 16,
                        cexAxixY = 16) {

    # local functions ####
    for_couple_tumors <- function(only_one_tumor) {
        count <- 0

        if (tolower(dataBase) == "legacy") {
            folder <- "clinical_data"
        } else {
            folder <- "clinical"
        }

        DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                folder)

        for (tumors in names(lista)) {
            count <- count + 1
            patient_file <- dir(path = DIR, pattern = "clinical_patient")
            clinical_df <- data.table::fread(input = file.path(DIR, patient_file),
                                                    sep = "\t")

            clinical_df <- as.data.frame(clinical_df[-c(1, 2),])
            rownames(clinical_df) <- clinical_df[, 2]
            clinical_df <- clinical_df[, -2]

            lista[[count]] <- clinical_df
        }

        if (is.null(term_keyword)) {
            if (!exists("clinical_df", envir = get(envir_link))) {
                concatenate_files(dataType = dataType,
                                                    dataBase = dataBase,
                                                    workDir = workDir,
                                                    tumor = tumor)
            }

            term_in_tumors <- apply(string_vars[["envir_link"]]$clinical_df,
                                                2, table)

            # at least two groups
            terms_unique <- lapply(term_in_tumors, length) >= 2
            terms_with_keyword <- names(term_in_tumors)[terms_unique]
            term_in_tumors <- term_in_tumors[terms_unique]

            terms_polished <- cbind(terms_with_keyword, unlist(lapply(term_in_tumors, length)))
            colnames(terms_polished) <- c("Term", "Subcategories_number")
            terms_polished <- terms_polished[order(terms_polished[, 1]), ]

            message(paste("These are the possible terms:",
                                        paste0(terms_with_keyword[order(terms_with_keyword)],
                                            collapse = "\n-"), sep = "\n-"))
            message(paste0("\nAll of them were saved in ",
                                        workDir, " as '",
                                        "DOAGDC_possible_terms_for_",
                                        toupper(tumor), "_",
                                        tolower(dataBase), ".txt"))
            write.csv(terms_polished,
                                    file.path(workDir,
                                                paste0("DOAGDC_possible",
                                                    "_terms_for_",
                                                    toupper(tumor), "_",
                                                    tolower(dataBase),
                                                    ".csv")),
                                    row.names = FALSE, quote = FALSE)
        } else {

            #list all terms in the files
            terms_repeted <- unlist(lapply(lista, colnames))

            #list all possible terms
            terms_unique <- table(terms_repeted)[order(table(terms_repeted),
                                                            decreasing = TRUE)]

            tmp <- names(terms_unique)[grep(pattern = term_keyword,
                                                    x = names(terms_unique),
                                                    ignore.case = TRUE,
                                                    perl = TRUE)]
            if (!only_one_tumor) {
                terms_with_keyword <- terms_unique[tmp]
                message(paste("These are the possible terms:",
                                            paste0(names(terms_with_keyword),
                                                collapse = "\n-"),
                                                sep = "\n-"))
                message(paste0("\nAll of them were saved in ",
                                workDir, " as '", "DOAGDC_possible_terms_for_",
                                tolower(term_keyword), "_",
                                tolower(dataBase), ".txt"))
                write.table(terms_with_keyword,
                            file.path(workDir,
                                        paste0("DOAGDC_possible_terms_for_",
                                            tolower(term_keyword), "_",
                                            tolower(dataBase), ".txt")),
                                            row.names = FALSE)
            } else {
                terms_with_keyword <- tmp
                message(paste("These are the possible terms:",
                                                paste0(terms_with_keyword,
                                                collapse = "\n-"),
                                                sep = "\n-"))
                message(paste0("\nAll of them were saved in ",
                                                workDir, " as '",
                                                "DOAGDC_possible_terms_for_",
                                                tolower(term_keyword),
                                                toupper(tumor),
                                                tolower(dataBase), ".txt"))
                write.table(terms_with_keyword,
                            file.path(workDir,
                                        paste0("DOAGDC_possible_terms_for_",
                                            tolower(term_keyword),
                                            toupper(tumor),
                                            tolower(dataBase), ".txt")),
                                            row.names = FALSE,
                                            col.names = "Term(s)")
                assign("terms_with_keyword", terms_with_keyword,
                                            envir = get(envir_link))
                assign("clinical_df", lista[[tumor]], envir = get(envir_link))
            }
        }
    }

    #code ####
    dataType <- gsub(" ", "_", dataType)
    Name <- gsub("-", "_", Name)

    if ("all" %in% tolower(tumor)) {
        lista <- vector("list", length(dir(path = file.path(workDir,
                                                                "DOAGDC"))))
        names(lista) <- dir(path = file.path(workDir, "DOAGDC"))
        for_couple_tumors(only_one_tumor = FALSE)
    } else if (length(tumor) > 1) {
        lista <- vector("list", length(tumor))
        names(lista) <- tumor
        for_couple_tumors(only_one_tumor = FALSE)
    } else {

        if (!(paste(toupper(tumor), toupper(dataBase), gsub(" ", "_",
                                                                dataType),
                        "tumor_data", sep = "_") %in% ls(all.names = TRUE,
                                                    envir = .GlobalEnv))) {
            assign(paste(toupper(tumor), toupper(dataBase), gsub(" ", "_",
                                                                    dataType),
                            "tumor_data", sep = "_"),
                        new.env(parent = emptyenv()), envir = .GlobalEnv)
            envir_link <- paste(toupper(tumor), toupper(dataBase),
                                    gsub(" ", "_", dataType), "tumor_data",
                                    sep = "_")
        } else {
            if (missing(env)) {
                envir_link <- paste(toupper(tumor), toupper(dataBase),
                                            gsub(" ", "_", dataType),
                                            "tumor_data", sep = "_")
            } else {
                envir_link <- deparse(substitute(env))
            }
        }
        string_vars <- list(envir_link = get(envir_link))

        dir.create(file.path(workDir, "DOAGDC", toupper(tumor), "Analyses"),
                        showWarnings = FALSE)

        PATH <- file.path(workDir, "DOAGDC", toupper(tumor), "Analyses")

        if (exists("Name.e", envir = get(envir_link))) {
            PATH <- file.path(PATH, string_vars[["envir_link"]]$Name.e)
        }
        assign("PATH", PATH, envir = get(envir_link))

        dir.create(PATH, showWarnings = FALSE)
        dir.create(paste0(PATH, "/Clinical_Results.",
                            tolower(string_vars[["envir_link"]]$dataType), "_",
                            toupper(string_vars[["envir_link"]]$Name)),
                    showWarnings = FALSE)
        DIR <- paste0(PATH, "/Clinical_Results.",
                        tolower(string_vars[["envir_link"]]$dataType), "_",
                            toupper(string_vars[["envir_link"]]$Name))

        lista <- vector("list", 1)
        names(lista) <- tumor

        if (missing("term")) {
            for_couple_tumors(only_one_tumor = TRUE)
        } else {
            if (!exists("clinical_df", envir = get(envir_link))) {
                concatenate_files(dataType = dataType,
                                        dataBase = dataBase,
                                        workDir = workDir,
                                        tumor = tumor)
            }

            term_in_tumors <- string_vars[["envir_link"]]$clinical_df[, term,
                                                                drop = FALSE]
            message("Removing 'unusable subcategories...")
            term_in_tumors[, 1] <- gsub("\\[Not Available\\]",
                                                NA, term_in_tumors[, 1])
            term_in_tumors[, 1] <- gsub("\\[Unknown\\]",
                                                NA, term_in_tumors[, 1])
            term_in_tumors[, 1] <- gsub("\\[Not Evaluated\\]",
                                                NA, term_in_tumors[, 1])
            term_in_tumors <- na.exclude(term_in_tumors)
            term_category <- unique(term_in_tumors[, 1])
            category_number <- length(term_category)

            if (category_number < 2) {
                stop(message("This term has less than 2 usable categories and",
                                    " cannot be used in this",
                                    " analysis!\n Please",
                                    " try another 'term'."))
            } else if (category_number > 5) {
                if (!is.na(as.numeric(term_in_tumors[1, 1]))) {

                    # Frequency Tables/Grouped Values
                    term_in_tumors[, 1] <- as.numeric(term_in_tumors[, 1])

                    Range <- max(term_in_tumors[, 1]) - min(term_in_tumors[, 1])
                    # > 2 and < 6
                    Groups_number <- group.number
                    # You must round up, not off
                    class_width <- ceiling(Range / (Groups_number + 1))
                    x_breaks <- seq(min(term_in_tumors[, 1]),
                                                max(term_in_tumors[, 1]),
                                                class_width)
                    step_val <- cut(term_in_tumors[, 1], breaks = x_breaks,
                                                include.lowest = TRUE)

                    # include.highest
                    step_val[is.na(step_val)] <- levels(step_val)[length(levels(step_val))]
                    levels(step_val)[length(levels(step_val))] <- gsub(',.*?\\]$',
                                                                                paste0(",", max(term_in_tumors[, 1]), "]"),
                                                                                levels(step_val)[length(levels(step_val))])
                    term_in_tumors[, 1] <- step_val

                    term_category <- levels(step_val)
                } else if (term == "clinical_stage") {
                    term_in_tumors[, 1] <- suppressWarnings(forcats::fct_collapse(term_in_tumors[, 1],
                                                                                Stage_I = c("Stage I", "Stage IA", "Stage IB", "Stage IC"),
                                                                                Stage_II = c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"),
                                                                                Stage_III = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IIIC1", "Stage IIIC2"),
                                                                                Stage_IV = c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC")))
                }
            }

            list_patients_per_category <- sapply(term_category, function(x) {
                rownames(term_in_tumors)[term_in_tumors[, 1] == x]
            })

            patients_per_category <- term_in_tumors
            patients_per_category <- table(patients_per_category)
            print(c("there are", patients_per_category, "patients"))

            assign("list_patients_per_category", list_patients_per_category,
                        envir = get(envir_link))

            #get expression values
            if (tolower(dataType) == "methylation") {
                # all <- string_vars[["envir_link"]]$gene_tumor_normalized
                tmp <- paste0("methylation_tumor_filtered_selected_", Name)
            } else {
                tmp <- paste0(tolower(dataType), "_tumor_normalized_selected_",
                                    toupper(Name))
            }

            selected <- eval(parse(text = paste0('string_vars[[',
                                                    '"envir_link"]]$',
                                                    paste0(tmp))))

            #split patient names
            patient_complete_name <- rownames(selected)
            #get the fisrt 3 elemetns inside all the keys
            patient_partial_name <- lapply(strsplit(x = patient_complete_name,
                                        split = "-", perl = TRUE), "[", 1:3)
            #collapse them in a vetor
            patient_partial_name <- unlist(lapply(patient_partial_name,
                                                        function(x) {
                                                            paste0(x,
                                                                collapse = "-")
                                                        }))
            #box plot category X expresion
            tmp <- !duplicated(patient_partial_name)
            final_df <- as.data.frame(selected)[tmp, ,drop = FALSE]
            final_df <- final_df[rownames(term_in_tumors), , drop = FALSE]
            final_df[term] <- "NA"
            # rownames(final_df) <- unique(patient_partial_name)
            #stacking the data
            desired_levels <- levels(term_in_tumors[, 1])
            final_df[rownames(term_in_tumors), 2] <- as.character(term_in_tumors[, 1])
            final_df[, 2] <- factor(final_df[, 2], levels = desired_levels)

            if (tolower(image_format) == "png") {
                png(filename = file.path(DIR, paste0(term, "boxplot_Clinical",
                                        "_category_by_expression_values.png")),
                            width = Width, height = Height,
                            res = Res, units = Unit)
            } else if (tolower(image_format) == "svg") {
                svg(filename = file.path(DIR, paste0(term, "boxplot_Clinical",
                                        "_category_by_expression_values.svg")),
                            width = Width, height = Height, onefile = TRUE)
            } else {
                stop(message("Insert a valid image_format! ('png' or 'svg')"))
            }
            final_df_plot <- cbind(final_df, rownames(final_df))
            colnames(final_df_plot) <- c("expression", "category", "patient")
            dodge <- ggplot2::position_dodge(width = 0.9)
            plot <- ggplot2::ggplot(final_df_plot,
                                ggplot2::aes(x = category,
                                            y = log2(expression + 1),
                                            fill = category)) +
                        ggplot2::stat_boxplot(geom = 'errorbar',
                                            position = ggplot2::position_dodge(1)) +
                        ggplot2::geom_boxplot(position = ggplot2::position_dodge(1),
                                            inherit.aes = TRUE,
                                            outlier.colour = NA) +
                        ggplot2::geom_jitter(position = ggplot2::position_jitter(0.2),
                                            alpha = 0.5, color = "black") +
                        ggplot2::labs(x = "Categories", y = "log2(RSEM+1)",
                                    title = term) +
                        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                            plot.title = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1.5)),
                            axis.title.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_text(size = 18),
                            axis.text.y = ggplot2::element_text(size = cexAxixY),
                            axis.text.x = ggplot2::element_text(size = cexAxixX),
                            legend.position = "none",
                            panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank())
            print(plot)
            dev.off()

            #anova
            Values <- log2(final_df[, 1] + 1)
            Categories <- final_df[, 2]
            aovobject <- aov(Values ~ Categories)
            ANOVA <- anova(aovobject)

            sink(file.path(DIR, paste0(term, "_ANOVA_Analysis_",
                                            "of_Variance_report.txt")))
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

            if (ANOVA[["Pr(>F)"]][1] == "NaN") {
                sink()
                stop(message("\nThis category '",
                                    term,
                                    "' does not have enough patients",
                                    " in each subcategory!"))
            }

            if (ANOVA[["Pr(>F)"]][1] <= p_cutoff) {
                cat("\n\n\n")
                cat(paste0("The p-value of ", ANOVA[["Pr(>F)"]][1],
                                " is not greater than ", p_cutoff, "."))
                cat("\n")
                cat("Hence pairwise t-tests will be performed.")
                cat("\n")
                cat("Tukey test")
                #tukey
                tukey <- TukeyHSD(aovobject, conf.level = confidence.level)
                print(tukey)
                cat("\n\n\n")
                sink()

                groups_export <- strsplit(names(which(tukey$Categories[, "p adj"] < FDR_cutoff)), "-")

                clinical_groups <- vector('list', length(groups_export))
                names(clinical_groups) <- names(which(tukey$Categories[, "p adj"] < FDR_cutoff))

                for (index in seq(1, length(groups_export))) {
                    selector <- final_df$clinical_stage %in% groups_export[[index]]
                    clinical_groups[[index]] <- final_df[selector, 2,
                                                                drop = FALSE]
                    clinical_groups[[index]][, 1] <- droplevels(clinical_groups[[index]][, 1])
                }

                message("The selected groups after",
                                " statistical analysis are:\n",
                                paste(names(clinical_groups), "\n"),
                                "\n\nPlease insert these names to use them",
                                " in posterior analysis.")

                assign("clinical_groups_clinical", clinical_groups,
                            envir = get(envir_link))
            } else {
                cat("\n\n\n")
                #to file
                cat(paste0("\nThe p-value of ", ANOVA[["Pr(>F)"]][1],
                                    " is greater than ", p_cutoff, "."))
                cat("\nHence, pairwise t-tests will not be performed.")
                cat("\n\n\n")
                sink()
                #to CLI
                message("The p-value of ", ANOVA[["Pr(>F)"]][1],
                                " is greater than ", p_cutoff, ".")
                message("\nHence pairwise t-tests will not be performed.")
            }
        }
    }
    message("Done!\n")
}
