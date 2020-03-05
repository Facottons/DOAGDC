#' Separate patients in groups
#'
#' \code{groups_identification_cox} is a function designed to separate patients
#' in groups, based in clinical data, and coxHR.
#'
#' @param name
#' @param data_type
#' @param width,height,res,unit,image_format
#' @param save_data
#' @param env
#' @param tumor
#' @param data_base
#' @param work_dir
#' @inheritParams concatenate_exon
#' @inheritParams download_gdc
#' @inheritParams groups_identification_mclust
#'
#' @return the groups generated after using the coxHR analysis.
#' @export
#'
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom survminer ggsurvplot
#' @importFrom XML xmlParse
#' @importFrom XML xmlToList
#' @importFrom yarrr pirateplot
#' @importFrom pROC roc
#' @importFrom pROC plot.roc
#'
#' @examples
#' library(DOAGDC)
#'
#' concatenate_expression("gene",
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
#'
#' groups_identification_cox_hr("HIF3A", "gene",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop",
#'     data_base = "legacy",
#'     env = CHOL_LEGACY_gene_tumor_data
#' )
groups_identification_cox_hr <- function(name,
                                         data_type,
                                         width = 3000,
                                         height = 2000,
                                         res = 400,
                                         unit = "px",
                                         image_format = "png",
                                         save_data = TRUE, env,
                                         tumor, data_base,
                                         work_dir) {

    # Library load (xlsx tools httr jsonlite data.table XML ggplot2 survival
    # survminer ggthemes grid yarrr reshape extrafont scales cgdsr) (plyr png)

    # localFUN
    expander <- function(big_list) {
        tmp_list <- vector("list")

        for (j in seq_len(length(big_list))) {
            name <- as.character(names(big_list[j]))
            tmp_list <- append(tmp_list, name)
        }

        return(tmp_list)
    }

    ty <- function(a) {
        ifelse(is.null(a), "", a)
    }

    open_xml <- function() {
        message("Extracting clinical data from XML...")
        pb <- txtProgressBar(min = 0, max = length(clinical), style = 3)
        for (i in seq_len(length(clinical))) {
            if (i > 1) {
                list[[i]] <- XML::xmlToDataFrame(
                    clinical[i],
                    stringsAsFactors = FALSE
                )
            } else {
                list <- list(XML::xmlToDataFrame(
                    clinical[i],
                    stringsAsFactors = FALSE
                ))
            }
            setTxtProgressBar(pb, i)
        }
        close(pb)

        fds <- as.data.frame(matrix(nrow = length(clinical), ncol = 90))
        fd <- as.data.frame(matrix(nrow = 1, ncol = 74))

        message("Exporting to csv...")
        pb <- txtProgressBar(min = 0, max = length(clinical), style = 3)
        for (patient in seq_len(length(clinical))) {
            for (col in seq_len(74)) {
                fd[1, col] <- na.omit(list[[patient]][col])
            }

            # stage_event
            data <- XML::xmlParse(clinical[patient])
            xml_data <- XML::xmlToList(data)$patient
            stage_event <- xml_data$stage_event
            tnm <- stage_event$tnm_categories$pathologic_categories
            fd$system_version <- ty(stage_event$system_version[[1]])
            fd$pathologic_stage <- ty(stage_event$pathologic_stage[[1]])
            fd$pathologic_t <- ty(tnm$pathologic_t[[1]])
            fd$pathologic_n <- ty(tnm$pathologic_n[[1]])
            fd$pathologic_m <- ty(tnm$pathologic_m[[1]])

            # follow_ups
            follow_ups <- xml_data$follow_ups$follow_up
            fd$days_to_death_followup <- ty(follow_ups$days_to_death[[1]])
            fd$days_to_last_followup_followup <- ty(
                follow_ups$days_to_last_followup[[1]]
            )

            if (ty(
                follow_ups$days_to_new_tumor_event_after_initial_treatment[[1]]
            ) == "") {
                fd$days_to_new_tumor_event_after_initial_treatment_followup <- ty(
                    follow_ups[["new_tumor_events"]][["new_tumor_event"]][["days_to_new_tumor_event_after_initial_treatment"]][[1]]
                )
            } else {
                fd$days_to_new_tumor_event_after_initial_treatment_followup <- ty(
                    follow_ups$days_to_new_tumor_event_after_initial_treatment[[1]]
                )
            }

            fd$day_of_form_completion2 <- ty(
                follow_ups$day_of_form_completion[[1]]
            )
            fd$month_of_form_completion2 <- ty(
                follow_ups$month_of_form_completion[[1]]
            )
            fd$year_of_form_completion2 <- ty(
                follow_ups$year_of_form_completion[[1]]
            )

            if (ty(follow_ups$new_neoplasm_occurrence_anatomic_site_text[[1]]) == "") {
                fd$new_neoplasm_event_occurrence_anatomic_site_followup <- ty(
                    follow_ups[["new_tumor_events"]][["new_tumor_event"]][["new_neoplasm_event_occurrence_anatomic_site"]][[1]]
                )
            } else {
                fd$new_neoplasm_event_occurrence_anatomic_site_followup <- ty(
                    follow_ups$new_neoplasm_occurrence_anatomic_site_text[[1]]
                )
            }

            fd$new_neoplasm_event_type_followup <- ty(
                follow_ups[["new_tumor_events"]][["new_tumor_event"]][["new_neoplasm_event_type"]][[1]]
            )

            if (ty(follow_ups$new_tumor_event_after_initial_treatment[[1]]) == "") {
                fd$new_tumor_event_after_initial_treatment_followup <- ty(
                    follow_ups$new_tumor_events$new_tumor_event_after_initial_treatment[[1]]
                )
            } else {
                fd$new_tumor_event_after_initial_treatment_followup <- ty(
                    follow_ups$new_tumor_event_after_initial_treatment[[1]]
                )
            }
            fd$person_neoplasm_cancer_status_followup <- ty(
                follow_ups$person_neoplasm_cancer_status[[1]]
            )
            fd$vital_status_followup <- ty(follow_ups$vital_status[[1]])

            fds[patient, ] <- fd

            setTxtProgressBar(pb, patient)
        }
        close(pb)

        colnames(fds)[1:74] <- names(list[[1]])
        colnames(fds)[75:90] <- colnames(fd)[75:90]
        rownames(fds) <- fds$bcr_patient_barcode

        fds$form_completion_ymd <- paste(fds$year_of_form_completion,
            fds$month_of_form_completion, fds$day_of_form_completion,
            sep = "/"
        )

        fds$form_completion_ymd_followup <- paste(
            fds$day_of_form_completion2,
            fds$month_of_form_completion2,
            fds$day_of_form_completion2,
            sep = "/"
        )

        fds$new_tumor_event_after_initial_treatment_initial <- fds$new_tumor_events
        fds$days_to_new_tumor_event_after_initial_treatment_initial <- NA
        for (patient in seq_len(nrow(fds))) {
            if (!(fds$new_tumor_events[patient] %in% c("", "NO"))) {
                tmp <- stringr::str_extract_all(
                    fds$new_tumor_events[patient],
                    "^(YES)\\d+"
                )[[1]][1]
                fds[patient, 93] <- "YES"
                fds[patient, 94] <- as.numeric(gsub("YES", "", tmp))
            }
        }

        # final fix
        tmp <- fds$days_to_new_tumor_event_after_initial_treatment_followup == "day"
        fds$days_to_new_tumor_event_after_initial_treatment_followup[tmp] <- ""
        tmp <- fds$days_to_new_tumor_event_after_initial_treatment_followup == "false"
        fds$days_to_new_tumor_event_after_initial_treatment_followup[tmp] <- ""

        fds$final_vital_status <- fds$vital_status
        tmp <- !(fds$vital_status_followup %in% c("", "vital_status"))
        fds$final_vital_status[tmp] <- fds$vital_status_followup[tmp]

        fds$final_vital_status_times <- ifelse(
            fds$final_vital_status == "Dead",
            ifelse(
                fds$days_to_death == "",
                fds$days_to_death_followup,
                fds$days_to_death
            ),
            fds$days_to_last_followup
        )

        fds$final_vital_status_times <- as.numeric(
            fds$final_vital_status_times
        )
        fds$final_vital_status_times[is.na(
            fds$final_vital_status_times
        )] <- 0

        fds$final_rfs_status <- ifelse(
            fds$new_tumor_event_after_initial_treatment_initial == "YES",
            "relapse",
            ifelse(
                fds$new_tumor_event_after_initial_treatment_followup == "YES",
                "relapse",
                "censored"
            )
        )

        tmp <- is.na(
            fds$days_to_new_tumor_event_after_initial_treatment_initial
        )
        fds$days_to_new_tumor_event_after_initial_treatment_initial[tmp] <- ""

        fds$final_rfs_status_times <- ifelse(
            fds$days_to_new_tumor_event_after_initial_treatment_initial == "",
            fds$days_to_new_tumor_event_after_initial_treatment_followup, ""
        )

        fds$final_rfs_status_times <- as.numeric(
            fds$final_rfs_status_times
        )
        fds$final_rfs_status_times[is.na(fds$final_rfs_status_times)] <- 0

        # sort before save
        fds <- fds[, order(colnames(fds))]

        return(fds)
    }

    roc_auc <- function(dir, name, width, height, res, unit, image_format) {
        values <- framel[, "log2p1"]
        groups <- framel[, "group"]

        # area under the curve...
        save_plot("coxHR", "_AUC_LR")
        par(pty = "s")
        pROC::roc(groups, glm_fit$fitted.values,
            plot = TRUE, legacy.axes = TRUE, percent = TRUE, quiet = TRUE,
            xlab = "False Positive Percentage",
            ylab = "True Postive Percentage",
            col = "#377eb8", lwd = 4, print.auc = TRUE
        )
        dev.off()

        rf_model <- randomForest::randomForest(groups ~ values)

        # ROC for random forest
        save_plot("coxHR", "_AUC_RF")
        par(pty = "s")
        pROC::roc(groups, rf_model$votes[, 1],
            plot = TRUE, legacy.axes = TRUE, percent = TRUE, quiet = TRUE,
            xlab = "False Positive Percentage",
            ylab = "True Postive Percentage",
            col = "#4daf4a", lwd = 4, print.auc = TRUE
        )
        dev.off()

        save_plot("coxHR", "_LR_RF")
        par(pty = "s")
        pROC::roc(groups, glm_fit$fitted.values,
            plot = TRUE, legacy.axes = TRUE, percent = TRUE, quiet = TRUE,
            xlab = "False Positive Percentage",
            ylab = "True Postive Percentage",
            col = "#377eb8", lwd = 4, print.auc = TRUE
        )
        pROC::plot.roc(groups, rf_model$votes[, 1],
            percent = TRUE, col = "#4daf4a", lwd = 4,
            print.auc = TRUE, add = TRUE, print.auc.y = 40
        )
        par(
            fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
            new = TRUE
        )
        plot(0, 0,
            type = "n", bty = "n", xaxt = "n", yaxt = "n",
            ylab = "", xlab = ""
        )
        legend("top",
            legend = c("Logisitic Regression", "Random Forest"),
            col = c("#377eb8", "#4daf4a"), lwd = 4,
            xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n",
            cex = 1, lty = 1
        )
        dev.off()
    }

    save_plot <- function(dir_name, file_name) {
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(
                    dir, dir_name,
                    paste0(
                        name,
                        file_name, ".png"
                    )
                ),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(
                    dir, dir_name,
                    paste0(
                        name,
                        file_name, ".svg"
                    )
                ),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format! ('png' or 'svg')"
            ))
        }
    }

    cox_ph <- function(step) {
        path_capture <- file.path(
            dir, "kaplan_meier", paste(
                name, toupper(step), "coxHR_summary.txt",
                sep = "_"
            )
        )

        if (step == "overall") {
            coxmodel <- summary(survival::coxph(survival::Surv(
                final_vital_status_times,
                final_vital_status
            ) ~ group_factor,
            data = kaplan_now
            ))

            surv <- summary(survival::Surv(
                final_vital_status_times,
                final_vital_status
            ) ~ group_factor,
            data = kaplan_now
            )
        } else if (step == "rfs") {
            coxmodel <- summary(survival::coxph(survival::Surv(
                final_rfs_status_times,
                final_rfs_status
            ) ~ group_factor,
            data = kaplan_now
            ))

            surv <- summary(survival::Surv(
                final_rfs_status_times,
                final_rfs_status
            ) ~ group_factor,
            data = kaplan_now
            )
        } else {
            coxmodel <- summary(survival::coxph(survival::Surv(
                final_dmfs_status_times,
                final_dmfs_status
            ) ~ group_factor,
            data = kaplan_now
            ))

            surv <- summary(survival::Surv(
                final_dmfs_status_times,
                final_dmfs_status
            ) ~ group_factor,
            data = kaplan_now
            )
        }

        cat("CoxHR summary\n", file = path_capture)
        capture.output(print(coxmodel), file = path_capture, append = TRUE)

        cat("\n********************\n\n", file = path_capture, append = TRUE)
        cat("CoxHR value\n", file = path_capture, append = TRUE)
        capture.output(print(coxmodel$conf.int[1]),
            file = path_capture,
            append = TRUE
        )

        # logrank pvalue
        cat("\n********************\n\n", file = path_capture, append = TRUE)
        cat("logrank pvalue\n", file = path_capture, append = TRUE)
        capture.output(print(coxmodel$logtest[3]),
            file = path_capture,
            append = TRUE
        )

        cat("\n\nSurv only\n", file = path_capture, append = TRUE)
        capture.output(print(surv), file = path_capture, append = TRUE)
        cat("\n\nlogrank_p\n", file = path_capture, append = TRUE)
        capture.output(as.numeric(surv$sctest["pvalue"]))
    }

    surv_plot <- function(dir, name, step,
                                    width, height, res, unit, image_format) {
        amount <- table(kaplan_now$group_factor)
        percent <- amount / nrow(kaplan_now)

        if (step == "overall") {
            fit <- survival::survfit(survival::Surv(
                final_vital_status_times,
                final_vital_status
            ) ~ group_factor,
            data = kaplan_now
            )
        } else if (step == "rfs") {
            fit <- survival::survfit(survival::Surv(
                final_rfs_status_times,
                final_rfs_status
            ) ~ group_factor,
            data = kaplan_now
            )
        } else {
            fit <- survival::survfit(survival::Surv(
                final_dmfs_status_times,
                final_dmfs_status
            ) ~ group_factor,
            data = kaplan_now
            )
        }

        # Plot data
        color <- RColorBrewer::brewer.pal(
            length(levels(kaplan_now$group_factor)), "Set2"
        )

        save_plot("kaplan_meier", paste0("Kaplan_", step))
        par(mar = c(3, 3, 3, 6))
        print(survminer::ggsurvplot(fit,
            data = kaplan_now, risk.table = TRUE,
            linetype = rep(1, length(color)),
            # conf.int = TRUE,
            palette = color, xlab = "Days",
            ggtheme = ggplot2::theme_bw(),
            size = 2,
            font.main = 30,
            font.x = 26,
            font.y = 26,
            font.tickslab = 24
        ))
        dev.off()

        # Create 5yr status
        kaplan_now$yr5_status <- "unknown"
        tmp <- kaplan_now$final_vital_status_times >= (365 * 5)
        kaplan_now$yr5_status[tmp] <- "alive_5yr"
        tmp1 <- kaplan_now$final_vital_status_times < (365 * 5)
        tmp2 <- kaplan_now$final_vital_status == 1
        kaplan_now$yr5_status[tmp1 & tmp2] <- "dead_5yr"
        kaplan_now$yr5_status <- as.factor(kaplan_now$yr5_status)

        # define the formula now
        pirate_formula <- formula(expression ~ yr5_status)

        message("Pirate plotting \n")
        if (nrow(kaplan_now) > 0) {
            save_plot(
                "kaplan_meier", paste0("PiratePlot_log2_expression_5yr_", step)
            )
            par(mar = c(3, 4.8, 4, 2))
            yarrr::pirateplot(
                formula = pirate_formula,
                data = kaplan_now,
                ylab = expression("Log"[2] * "(Expression + 1)"),
                inf.method = "iqr", cex.lab = 1.2, cex.axis = 1.8,
                inf.disp = "rect", jitter.val = 0.08, theme = 2,
                cex.names = 1.5,
                pal = "pony", avg.line.fun = median, point.cex = 1.3,
                inf.f.col = "#969696", xlab = ""
            )
            dev.off()
        }

        # stats
        path_capture <- file.path(
            dir, "kaplan_meier",
            paste0(
                name, "PiratePlot_log2Expression_5yr_", step, "_stats.txt"
            )
        )

        fit <- aov(pirate_formula, data = kaplan_now)

        cat("********************\n", file = path_capture, append = TRUE)
        cat("A - ANOVA\n", file = path_capture, append = TRUE)
        capture.output(summary(fit), file = path_capture, append = TRUE)

        cat("\n********************\n", file = path_capture, append = TRUE)
        cat("B - Shapiro Test for Normality\n\n",
            file = path_capture,
            append = TRUE
        )
        cat("*B1 - all values test\n", file = path_capture, append = TRUE)
        tmp <- shapiro.test(kaplan_now$expression)
        capture.output(print(tmp), file = path_capture, append = TRUE)

        cat("\n*B2 - fit residuals test\n", file = path_capture, append = TRUE)
        tmp <- shapiro.test(residuals(fit))
        capture.output(print(tmp), file = path_capture, append = TRUE)

        cat("\n********************\n\n", file = path_capture, append = TRUE)
        cat("C - Tukey Honest Significant Differences post-test\n",
            file = path_capture, append = TRUE
        )
        tmp <- TukeyHSD(fit)
        capture.output(print(tmp), file = path_capture, append = TRUE)

        if (length(levels(kaplan_now$group_factor)) == 2) {
            cat("\n********************\n", file = path_capture, append = TRUE)
            cat("D - unpaired t-test\n", file = path_capture, append = TRUE)
            test <- var.test(expression ~ yr5_status, kaplan_now,
                alternative = "two.sided"
            )
            tmp <- t.test(pirate_formula,
                data = kaplan_now, alternative = "two.sided",
                mu = 0, paired = FALSE, var.equal = as.numeric(test[3]) <= 0.05,
                conf.level = 0.95
            )
            capture.output(print(tmp), file = path_capture, append = TRUE)
        }

        cat("\n********************\n", file = path_capture, append = TRUE)
        cat("non-parametric assumption:\n", file = path_capture, append = TRUE)
        cat("********************\n", file = path_capture, append = TRUE)
        cat("A - Kruskal-Wallis test\n", file = path_capture, append = TRUE)
        tmp <- kruskal.test(data = kaplan_now, pirate_formula)
        capture.output(print(tmp), file = path_capture, append = TRUE)

        cat("\n********************\n", file = path_capture, append = TRUE)
        cat("B - Non-adjusted Wilcoxon (Mann-Whitney U)\n",
            file = path_capture,
            append = TRUE
        )
        tmp <- pairwise.wilcox.test(
            x = kaplan_now$expression,
            g = kaplan_now$yr5_status,
            p.adjust.method = "none",
            paired = FALSE,
            exact = TRUE,
            correct = TRUE,
            alternative = "two.sided"
        )
        capture.output(print(tmp), file = path_capture, append = TRUE)

        cat("\n********************\n\n", file = path_capture, append = TRUE)
        cat("C - Benjamini FDR adjusted Wilcoxon (Mann-Whitney U)\n",
            file = path_capture, append = TRUE
        )
        tmp <- pairwise.wilcox.test(
            x = kaplan_now$expression,
            g = kaplan_now$yr5_status,
            p.adjust.method = "fdr",
            paired = FALSE,
            exact = TRUE,
            correct = TRUE,
            alternative = "two.sided"
        )
        capture.output(print(tmp), file = path_capture, append = TRUE)
    }

    # code
    data_type <- gsub(" ", "_", data_type)
    name <- gsub("-", "_", name)

    data_base_boo <- tolower(data_base) == "legacy"

    path <- file.path(
        work_dir, "DOAGDC", toupper(tumor), "Analyses",
        toupper(name)
    )

    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    dir.create(
        path = file.path(
            work_dir, "DOAGDC", toupper(tumor), "Analyses",
            toupper(name)
        ),
        showWarnings = FALSE
    )

    dir.create(file.path(
        path,
        paste0("/survival_Results_", toupper(name))
    ),
    showWarnings = FALSE
    )
    dir <- paste0(path, "/survival_Results_", toupper(name))
    dir.create(file.path(dir, "coxHR"), showWarnings = FALSE)
    dir.create(file.path(dir, "kaplan_meier"), showWarnings = FALSE)

    # Download clinical
    download_gdc(
        data_type = "clinical_supplement", tumor = tumor,
        data_base = ifelse(
            data_base_boo, "legacy", "gdc"
        ), work_dir = work_dir
    )

    folder_name <- ifelse(
        data_base_boo,
        "clinical_supplement_data",
        "gdc_clinical_supplement_data"
    )

    manifest <- data.frame(read.table(
        file = file.path(
            work_dir, "DOAGDC",
            toupper(tumor),
            folder_name,
            "manifest.sdrf"
        ),
        stringsAsFactors = FALSE, header = TRUE,
        sep = "\t"
    ))

    available_files <- manifest$file_name

    clinical <- file.path(
        workDir, "DOAGDC", tumor, folder_name, available_files
    )

    final_df  <- open_xml()

    # Save table
    write.csv(final_df, file.path(dir, paste0(tumor, "_clinical_xml", ".csv")))

    # Table File position
    if (missing(env)) {
        stop(message("Insert the Environment name as 'env' argument!"))
    } else {
        ev <- deparse(substitute(env))
    }

    assign("clinical_xml", final_df, envir = get(ev))

    # get expression values
    if (tolower(data_type) == "methylation") {
        tmp <- paste0("methylation_tumor_filtered_selected_", name)
    } else {
        tmp <- paste0(
            tolower(data_type), "_tumor_normalized_selected_",
            toupper(name)
        )
    }

    string_vars <- list(ev = get(ev))
    framel <- eval(parse(text = paste0(
        'string_vars[["ev"]]$',
        paste0(tmp)
    )))

    framel <- as.data.frame(framel)

    # remove undisered tissues
    tmp <- c("Hypopharynx", "Larynx", "Lip")
    boolean_tmp <- final_df$anatomic_neoplasm_subdivision %in% tmp
    final_df <- final_df[!boolean_tmp, ]

    # NOTE Start Cutoff Finder
    # https://www.ncbi.nlm.nih.gov/pubmed/23251644

    # clinical
    rownames(final_df) <- final_df$bcr_patient_barcode

    best_z_table <- NULL

    # add patient code
    framel$patient_code <- rownames(framel)

    # collect ordering
    matching_ordering <- match(
        framel[, "patient_code"],
        rownames(final_df)
    )

    # Add clinical values to the frames list
    framel <- cbind(framel, final_df[
        matching_ordering,
        c(
            "race", "gender",
            "pathologic_stage",
            "pathologic_t",
            "pathologic_n",
            "pathologic_m",
            "final_vital_status",
            "final_vital_status_times",
            "final_rfs_status",
            "final_rfs_status_times",
            "final_dmfs_status",
            "final_dmfs_status_times"
        )
    ])

    tissue_type_simple <- unname(sapply(rownames(framel), function(w) {
        paste(unlist(strsplit(w, "-"))[4], collapse = "-")
    }))
    tissue_type_simple <- unname(sapply(tissue_type_simple, function(w) {
        paste(unlist(strsplit(w, ""))[1:2], collapse = "")
    }))
    framel[, "tissue_type_simple"] <- tissue_type_simple

    # Calculate z-scores
    tumor_valueszscore <- as.numeric(scale(framel[, name]))

    # Save values to table
    framel[, "log2p1"] <- log2(framel[, name] + 1)
    framel[, "zscore"] <- NA
    framel[, "zscore"] <- tumor_valueszscore

    # Collect gene with
    table <- framel

    # Sort by expression level
    table <- table[order(table[, "zscore"]), ]

    ### Remove missing data
    table <- table[(!table$final_rfs_status == "unknown"), ]
    tmp <- table$final_rfs_status == "dead_other_cause_before_relapse"
    table <- table[!tmp, ]
    table <- table[(!table$final_rfs_status == "vital_status"), ]
    table <- table[!is.na(table$final_rfs_status), ]
    table <- table[!is.nan(table$final_rfs_status), ]
    table <- table[(!table$final_rfs_status_times == "unknown"), ]
    tmp <- table$final_rfs_status_times == "dead_other_cause_before_relapse"
    table <- table[!tmp, ]
    table <- table[(!table$final_rfs_status_times == "day"), ]
    table <- table[!is.na(table$final_rfs_status_times), ]
    table <- table[!is.nan(table$final_rfs_status_times), ]
    table$final_rfs_status_times <- as.numeric(
        table$final_rfs_status_times
    )

    # remove patients without relapse information
    table$final_rfs_status <- ifelse(
        table$final_rfs_status == "relapse", 1, 0
    )

    nmin <- round(0.04 * dim(table)[1])
    if (nmin < 10) {
        nmin <- 10
    }
    nmax <- dim(table)[1] - nmin
    patient_cuts_sequence <- nmin:nmax

    # create matrix to receive the results
    cutoff_opt <- matrix(
        nrow = length(patient_cuts_sequence),
        ncol = 6
    )
    colnames(cutoff_opt) <- c(
        "patient_no", "zscorecut", "HR",
        "HR_min95", "HR_max95", "logrank_p"
    )
    cutoff_opt <- as.data.frame(cutoff_opt)
    cutoff_opt$patient_no <- patient_cuts_sequence

    for (cut in patient_cuts_sequence) {
        if (table[, 1][cut] == 0) {
            # check if it isnt the last zero
            if (table[, 1][cut + 1] == 0) {
                # CoxHR
                cutoff_opt[cut, "HR"] <- Inf

                # CoxHR min .95
                cutoff_opt[cut, "HRmin95"] <- Inf

                # CoxHR max .95
                cutoff_opt[cut, "HRmax95"] <- Inf

                # logrank pvalue
                cutoff_opt[cut, "logrank_p"] <- 1

                # if is the last zero, do normally
            } else {
                # Create the classification vector for now
                classification_vector <- c(
                    rep("low", cut),
                    rep(
                        "high",
                        dim(table)[1] - cut
                    )
                )
                classification_vector <- factor(classification_vector,
                    levels = c("low", "high")
                )

                # create the model summary
                coxmodel <- suppressWarnings(summary(
                    survival::coxph(survival::Surv(
                        final_rfs_status_times, final_rfs_status
                    ) ~ classification_vector,
                    data = table
                    )
                ))

                # collect line now for table
                line_now <- match(cut, cutoff_opt[, "patient_no"])

                cutoff_opt[line_now, "zscorecut"] <- table[cut, "zscore"]

                # CoxHR
                cutoff_opt[line_now, "HR"] <- coxmodel$conf.int[1]

                # CoxHR min .95
                cutoff_opt[line_now, "HR_min95"] <- coxmodel$conf.int[3]

                # CoxHR max .95
                cutoff_opt[line_now, "HR_max95"] <- coxmodel$conf.int[4]

                # logrank pvalue
                cutoff_opt[line_now, "logrank_p"] <- coxmodel$logtest[3]
            }
        } else {

            # Create the classification vector for now
            classification_vector <- c(
                rep("low", cut),
                rep(
                    "high",
                    dim(table)[1] - cut
                )
            )
            classification_vector <- factor(classification_vector,
                levels = c("low", "high")
            )

            # create the model summary
            coxmodel <- suppressWarnings(
                summary(
                    survival::coxph(survival::Surv(
                        final_rfs_status_times, final_rfs_status
                    ) ~ classification_vector,
                    data = table
                    )
                )
            )

            # collect line now for table
            line_now <- match(cut, cutoff_opt[, "patient_no"])

            cutoff_opt[line_now, "zscorecut"] <- table[
                cut,
                "zscore"
            ]

            # CoxHR
            cutoff_opt[line_now, "HR"] <- coxmodel$conf.int[1]
            # CoxHR min .95
            cutoff_opt[line_now, "HR_min95"] <- coxmodel$conf.int[3]
            # CoxHR max .95
            cutoff_opt[line_now, "HR_max95"] <- coxmodel$conf.int[4]
            # logrank pvalue
            cutoff_opt[line_now, "logrank_p"] <- coxmodel$logtest[3]
        }
    }

    # remove infinite rows
    # pick infinite rows
    remove_rows <- -sort(unique(c(
        which(is.infinite(cutoff_opt$HR)),
        which(is.infinite(cutoff_opt$HR_min95)),
        which(is.infinite(cutoff_opt$HR_max95)),
        which(is.na(cutoff_opt$logrank_p))
    )))

    cutoff_opt <- cutoff_opt[!is.na(cutoff_opt$HR_min95), ]

    # save optim table
    write.csv(cutoff_opt,
        file = file.path(dir, "coxHR", paste0(
            "Optimization_",
            name, ".csv"
        ))
    )

    if (nrow(cutoff_opt) == 0) {
        stop(message(
            "\nLow expression values. Please, ",
            "try another expression data!\n"
        ))
    }

    # Collect best z score
    best_z_table <- suppressWarnings(
        cutoff_opt$zscorecut[which(min(cutoff_opt$logrank_p,
            na.rm = TRUE
        ) == cutoff_opt$logrank_p)][1]
    )

    if (nrow(cutoff_opt[remove_rows, ]) > 0) {
        save_plot("coxHR", "_with_confidence_interval")
        par(mar = c(5, 5, 3, 5))
        plot(
            x = cutoff_opt$zscorecut,
            y = cutoff_opt$HR, type = "l",
            ylim = c(
                min(cutoff_opt$HR_min95, na.rm = TRUE),
                max(cutoff_opt$HR_max95, na.rm = TRUE)
            ),
            xlim = c(
                min(cutoff_opt$zscorecut, na.rm = TRUE),
                max(cutoff_opt$zscorecut, na.rm = TRUE)
            ),
            lwd = 2, lty = 1, xlab = "Expression z-score",
            ylab = "Cox HR RFS", axes = FALSE,
            col = hsv(.58, .95, .62)
        )

        axis(
            side = 1, lwd = 2, cex.axis = 1.5, cex = 1.5,
            at = round(seq(min(cutoff_opt$zscorecut, na.rm = TRUE),
                max(cutoff_opt$zscorecut, na.rm = TRUE),
                by = 0.1
            ), 2)
        )
        axis(
            side = 2, lwd = 2, cex.axis = 1.5, cex = 1.5,
            col = hsv(.58, .95, .62), las = 1
        )

        lines(
            x = cutoff_opt$zscorecut,
            y = cutoff_opt$HR_min95,
            lwd = 1, lty = 3, col = hsv(.34, .88, .60, 0.78)
        )
        lines(
            x = cutoff_opt$zscorecut,
            y = cutoff_opt$HR_max95,
            lwd = 1, lty = 3, col = hsv(.34, .88, .60, 0.78)
        )

        text(best_z_table, max(cutoff_opt$HR, na.rm = TRUE) + 0.12,
            paste0("z-score = ", round(best_z_table, 2)),
            cex = 0.8
        )
        points(best_z_table, max(cutoff_opt$HR, na.rm = TRUE))

        points(table[, "zscore"],
            rep(
                min(cutoff_opt$HR_min95, na.rm = TRUE),
                length(table[, "zscore"])
            ),
            pch = "|", col = hsv(.0, .0, .15, .35)
        )

        par(new = TRUE)
        plot(
            x = cutoff_opt$zscorecut,
            y = cutoff_opt$logrank_p, type = "l",
            ylim = c(
                max(cutoff_opt$logrank_p, na.rm = TRUE),
                min(cutoff_opt$logrank_p, na.rm = TRUE) - 0.1
            ),
            xlim = c(
                min(cutoff_opt$zscorecut, na.rm = TRUE),
                max(cutoff_opt$zscorecut, na.rm = TRUE)
            ),
            lwd = 2, lty = 1, xlab = NA, ylab = NA,
            axes = FALSE, cex = 1.5, col = hsv(.028, .83, .99, .59)
        )
        axis(
            side = 4, lwd = 2, cex.axis = 1.5, cex = 1.5,
            col = hsv(.028, .83, .99), las = 1
        )
        mtext(side = 4, line = 3.8, "p logrank", cex = 1.2)

        min_p <- min(cutoff_opt$logrank_p, na.rm = TRUE)

        text(best_z_table, min_p - 0.01, paste0(
            "p logrank = ",
            round(min_p, 4)
        ), cex = 0.8)
        points(best_z_table, min_p)
        dev.off()

        save_plot("coxHR", "_log2_with_confidence_interval")
        par(mar = c(5, 5, 3, 5))
        plot(
            x = cutoff_opt$zscorecut,
            y = log2(cutoff_opt$HR), type = "l",
            ylim = c(
                min(log2(cutoff_opt$HR_min95), na.rm = TRUE),
                max(log2(cutoff_opt$HR_max95), na.rm = TRUE)
            ),
            xlim = c(
                min(cutoff_opt$zscorecut, na.rm = TRUE),
                max(cutoff_opt$zscorecut, na.rm = TRUE)
            ),
            lwd = 2, lty = 1, xlab = "Expression z-score",
            ylab = expression("Log"[2] * "(Cox HR RFS)"), axes = FALSE,
            col = hsv(.58, .95, .62)
        )

        axis(
            side = 1, lwd = 2, cex.axis = 1.5, cex = 1.5,
            at = round(seq(min(cutoff_opt$zscorecut, na.rm = TRUE),
                max(cutoff_opt$zscorecut, na.rm = TRUE),
                by = 0.1
            ), 2)
        )
        axis(
            side = 2, lwd = 2, cex.axis = 1.5, cex = 1.5,
            col = hsv(.58, .95, .62), las = 1
        )

        lines(
            x = cutoff_opt$zscorecut,
            y = log2(cutoff_opt$HR_min95),
            lwd = 1, lty = 3, col = hsv(.34, .88, .60, 0.78)
        )
        lines(
            x = cutoff_opt$zscorecut,
            y = log2(cutoff_opt$HR_max95),
            lwd = 1, lty = 3, col = hsv(.34, .88, .60, 0.78)
        )

        text(best_z_table, max(log2(cutoff_opt$HR),
            na.rm = TRUE
        ) + 0.105,
        paste0("z-score = ", round(best_z_table, 2)),
        cex = 0.8
        )
        points(best_z_table, max(log2(cutoff_opt$HR), na.rm = TRUE))

        points(table[, "zscore"],
            rep(
                min(log2(cutoff_opt$HR_min95), na.rm = TRUE),
                length(table[, "zscore"])
            ),
            pch = "|", col = hsv(.0, .0, .15, .35)
        )
        par(new = TRUE)

        plot(
            x = cutoff_opt$zscorecut,
            y = cutoff_opt$logrank_p, type = "l",
            ylim = c(
                max(cutoff_opt$logrank_p, na.rm = TRUE),
                min(cutoff_opt$logrank_p, na.rm = TRUE)
            ),
            xlim = c(
                min(cutoff_opt$zscorecut, na.rm = TRUE),
                max(cutoff_opt$zscorecut, na.rm = TRUE)
            ),
            lwd = 2, lty = 1, xlab = NA, ylab = NA,
            axes = FALSE, cex = 1.5, col = hsv(.028, .83, .99, .59)
        )
        axis(
            side = 4, lwd = 2, cex.axis = 1.5, cex = 1.5,
            col = hsv(.028, .83, .99), las = 1
        )
        mtext(side = 4, line = 3.8, "p logrank", cex = 1.2)

        min_p <- min(cutoff_opt$logrank_p, na.rm = TRUE)

        text(best_z_table, min_p - 0.015, paste0(
            "p logrank = ",
            round(min_p, 4)
        ), cex = 0.8)
        points(best_z_table, min_p)

        dev.off()
    }

    save_plot("coxHR", "_hide_confidence_interval")
    par(mar = c(5, 5, 3, 5))
    plot(
        x = cutoff_opt$zscorecut,
        y = cutoff_opt$HR, type = "l",
        lwd = 2, lty = 1, xlab = "Expression z-score", ylab = "Cox HR RFS",
        axes = FALSE,
        col = hsv(.58, .95, .62, .39)
    )

    axis(
        side = 1, lwd = 2, cex.axis = 1.5, cex = 1.5,
        at = round(seq(min(cutoff_opt$zscorecut, na.rm = TRUE),
            max(cutoff_opt$zscorecut, na.rm = TRUE),
            by = 0.1
        ), 2)
    )
    axis(
        side = 2, lwd = 2, cex.axis = 1.5, cex = 1.5,
        col = hsv(.58, .95, .62), las = 1
    )
    text(best_z_table, max(cutoff_opt$HR, na.rm = TRUE) + 0.12,
        paste0("z-score = ", round(best_z_table, 2)),
        cex = 0.8
    )
    points(best_z_table, max(cutoff_opt$HR, na.rm = TRUE))

    points(cutoff_opt$zscorecut,
        rep(
            min(cutoff_opt$HR, na.rm = TRUE),
            length(cutoff_opt$zscorecut)
        ),
        pch = "|", col = hsv(.0, .0, .15, .35)
    )

    par(new = TRUE)
    plot(
        x = cutoff_opt$zscorecut,
        y = cutoff_opt$logrank_p, type = "l",
        ylim = c(
            max(cutoff_opt$logrank_p, na.rm = TRUE),
            min(cutoff_opt$logrank_p, na.rm = TRUE)
        ),
        lwd = 3, lty = 1, xlab = NA, ylab = NA,
        axes = FALSE, cex = 1.5, col = hsv(.028, .83, .99, .39)
    )
    axis(
        side = 4, lwd = 2, cex.axis = 1.5, cex = 1.5,
        col = hsv(.028, .83, .99), las = 1
    )
    mtext(side = 4, line = 3.8, "p logrank", cex = 1.2)

    min_p <- min(cutoff_opt$logrank_p, na.rm = TRUE)
    text(best_z_table, min_p - 0.01, paste0("p logrank = ", round(min_p, 4)),
        cex = 0.8
    )
    points(best_z_table, min_p)

    dev.off()

    save_plot("coxHR", "_log2_hide_confidence_interval")
    par(mar = c(5, 5, 3, 5))
    plot(
        x = cutoff_opt$zscorecut,
        y = log2(cutoff_opt$HR), type = "l",
        lwd = 2, lty = 1, xlab = "Expression z-score",
        ylab = expression("Log"[2] * "(Cox HR RFS)"), axes = FALSE,
        col = hsv(.58, .95, .62, .39)
    )

    axis(
        side = 1, lwd = 2, cex.axis = 1.5, cex = 1.5,
        at = round(seq(min(cutoff_opt$zscorecut, na.rm = TRUE),
            max(cutoff_opt$zscorecut, na.rm = TRUE),
            by = 0.1
        ), 2)
    )
    axis(
        side = 2, lwd = 2, cex.axis = 1.5, cex = 1.5,
        col = hsv(.58, .95, .62), las = 1
    )

    text(best_z_table, max(log2(cutoff_opt$HR), na.rm = TRUE) + 0.105,
        paste0("z-score = ", round(best_z_table, 2)),
        cex = 0.8
    )
    points(best_z_table, max(log2(cutoff_opt$HR), na.rm = TRUE))

    points(cutoff_opt$zscorecut,
        rep(
            min(log2(cutoff_opt$HR), na.rm = TRUE),
            length(cutoff_opt$zscorecut)
        ),
        pch = "|", col = hsv(.0, .0, .15, .35)
    )

    par(new = TRUE)

    plot(
        x = cutoff_opt$zscorecut,
        y = cutoff_opt$logrank_p, type = "l",
        ylim = c(
            max(cutoff_opt$logrank_p, na.rm = TRUE),
            min(cutoff_opt$logrank_p, na.rm = TRUE)
        ),
        lwd = 3, lty = 1, xlab = NA, ylab = NA,
        axes = FALSE, cex = 1.5, col = hsv(.028, .83, .99, .39)
    )
    axis(
        side = 4, lwd = 2, cex.axis = 1.5, cex = 1.5,
        col = hsv(.028, .83, .99), las = 1
    )
    mtext(side = 4, line = 3.8, "p logrank", cex = 1.2)

    min_p <- min(cutoff_opt$logrank_p, na.rm = TRUE)
    text(best_z_table, min_p - 0.01, paste0(
        "p logrank = ",
        round(min_p, 4)
    ), cex = 0.8)
    points(best_z_table, min_p)

    dev.off()

    # Save best z table
    # save optim table
    write.csv(best_z_table, file = file.path(
        dir, "coxHR",
        paste0("Best_z_Table", ".csv")
    ))

    # NOTE A
    # add unknown to everyone
    framel[, "classification"] <- "unknown"
    # add low expressing tag
    framel[
        which(framel[, "zscore"] <= best_z_table),
        "classification"
    ] <- "low"
    # add high expressing tag
    framel[
        which(framel[, "zscore"] > best_z_table),
        "classification"
    ] <- "high"
    # create group tissue_type
    framel[, "group"] <- framel[, "classification"]

    move <- which(framel[, "classification"] == "unknown")
    framel[move, "group"] <- framel[move, "tissue_type_simple"]
    framel[, "group"] <- factor(framel[, "group"],
        levels = c("low", "high")
    )

    roc_auc()

    color <- RColorBrewer::brewer.pal(length(levels(framel$group)), "Set2")

    plot_pie <- data.frame(
        Group = names(table(framel$group)),
        n = as.numeric(table(framel$group)),
        prop = round(as.numeric(table(framel$group) / nrow(framel)) * 100, 2)
    )
    plot_pie <- plot_pie[order(plot_pie$Group, decreasing = TRUE), ]
    plot_pie$lab.ypos <- cumsum(plot_pie$prop) - (0.5 * plot_pie$prop)

    save_plot("coxHR", "_patients_all_conditions")
    a <- ggplot2::ggplot(plot_pie, ggplot2::aes(2, y = prop, fill = Group)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar("y", start = 0) +
        ggplot2::geom_text(ggplot2::aes(y = lab.ypos, label = n),
            color = "white"
        ) +
        ggplot2::scale_fill_manual(values = color) +
        ggplot2::theme_void() +
        ggplot2::xlim(0.5, 2.5)
    print(a)
    dev.off()

    # define the formula now
    formula_now <- formula(framel[, "log2p1"] ~ framel[, "group"])

    save_plot("coxHR", "_PiratePlot_log2_expression_all_conditions")
    par(mar = c(3, 4.8, 4, 2))
    yarrr::pirateplot(
        formula = formula_now,
        data = framel, gl.lwd = 0,
        ylab = expression("Log"[2] * "(Expression + 1)"),
        inf.method = "iqr", cex.lab = 1.2, cex.axis = 1.8,
        inf.disp = "rect", jitter.val = 0.08, theme = 2,
        cex.names = 1.5,
        pal = "pony", avg.line.fun = median, point.cex = 1.3,
        inf.f.col = hsv(0, 0, .78), xlab = "",
        inf.b.col = hsv(0, 0, .47)
    )
    dev.off()

    path_capture <- file.path(
        dir, "coxHR",
        paste0(
            name, "_PiratePlot_log2_expression",
            "_all_conditions_stats.txt"
        )
    )

    fit <- aov(formula_now, data = framel)
    cat("********************\n", file = path_capture, append = TRUE)
    cat("A - ANOVA\n", file = path_capture, append = TRUE)
    capture.output(summary(fit), file = path_capture, append = TRUE)

    cat("\n********************\n", file = path_capture, append = TRUE)
    cat("B - Shapiro Test for Normality\n\n",
        file = path_capture,
        append = TRUE
    )
    cat("*B1 - all values test\n", file = path_capture, append = TRUE)
    tmp <- shapiro.test(framel[, "log2p1"])
    capture.output(print(tmp), file = path_capture, append = TRUE)

    cat("\n*B2 - fit residuals test\n", file = path_capture, append = TRUE)
    tmp <- shapiro.test(residuals(fit))
    capture.output(print(tmp), file = path_capture, append = TRUE)

    cat("\n********************\n\n", file = path_capture, append = TRUE)
    cat("C - Tukey Honest Significant Differences post-test\n",
        file = path_capture, append = TRUE
    )
    tmp <- TukeyHSD(fit)
    capture.output(print(tmp), file = path_capture, append = TRUE)

    if (length(levels(framel[, "group"])) == 2) {
        cat("\n********************\n", file = path_capture, append = TRUE)
        cat("D - unpaired t-test\n", file = path_capture, append = TRUE)
        test <- var.test(expression ~ yr5_status, framel,
            alternative = "two.sided"
        )
        tmp <- t.test(pirate_formula,
            data = framel, alternative = "two.sided",
            mu = 0, paired = FALSE, var.equal = as.numeric(test[3]) <= 0.05,
            conf.level = 0.95
        )
        capture.output(print(tmp), file = path_capture, append = TRUE)
    }

    cat("\n********************\n", file = path_capture, append = TRUE)
    cat("non-parametric assumption:\n", file = path_capture, append = TRUE)
    cat("********************\n", file = path_capture, append = TRUE)
    cat("A - Kruskal-Wallis test\n", file = path_capture, append = TRUE)
    tmp <- kruskal.test(data = framel, pirate_formula)
    capture.output(print(tmp), file = path_capture, append = TRUE)

    cat("\n********************\n", file = path_capture, append = TRUE)
    cat("B - Non-adjusted Wilcoxon (Mann-Whitney U)\n",
        file = path_capture,
        append = TRUE
    )
    tmp <- pairwise.wilcox.test(
        x = framel[, "log2p1"],
        g = framel[, "group"],
        p.adjust.method = "none",
        paired = FALSE,
        exact = TRUE,
        correct = TRUE,
        alternative = "two.sided"
    )
    capture.output(print(tmp), file = path_capture, append = TRUE)

    cat("\n********************\n\n", file = path_capture, append = TRUE)
    cat("C - Benjamini FDR adjusted Wilcoxon (Mann-Whitney U)\n",
        file = path_capture, append = TRUE
    )
    tmp <- pairwise.wilcox.test(
        x = framel[, "log2p1"],
        g = framel[, "group"],
        p.adjust.method = "fdr",
        paired = FALSE,
        exact = TRUE,
        correct = TRUE,
        alternative = "two.sided"
    )
    capture.output(print(tmp), file = path_capture, append = TRUE)


    # NOTE Overall
    # Create table for patient numbers
    patient_amount_table <- matrix(ncol = 4, nrow = 8)
    colnames(patient_amount_table) <- c(
        "type", "classification",
        "amount", "percent"
    )
    patient_amount_table <- data.frame(patient_amount_table)

    # use previously defined markers

    # Copy patient_amount_table
    patient_n <- patient_amount_table

    # OVERALL SURVIVAL ANALYSIS
    # Kaplan meyer using best scenario
    kaplan_now <- framel[, c(
        name,
        "classification",
        "final_vital_status_times",
        "final_vital_status",
        "final_rfs_status",
        "final_rfs_status_times",
        "final_dmfs_status",
        "final_dmfs_status_times"
    )]

    # Kaplan
    colnames(kaplan_now)[2] <- "classification"

    # Collect just rows with classification
    kaplan_now <- kaplan_now[c(
        which(kaplan_now[, "classification"] == "low"),
        which(kaplan_now[, "classification"] == "high")
    ), ]


    # Factorize classification
    kaplan_now$group_factor <- factor(
        kaplan_now$classification,
        levels = c("low", "high")
    )

    ############# COLLECT DATA

    # Fill table with patient amount - ALL
    rows_now <- c(1, 2)
    patient_n[rows_now, "type"] <- "All"
    patient_n[rows_now, "classification"] <- c("low", "high")
    patient_n[rows_now[1], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["low"]
    )
    patient_n[rows_now[2], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["high"]
    )

    tmp <- patient_n[rows_now[1], "amount"] + patient_n[rows_now[2], "amount"]
    patient_n[rows_now[1], "percent"] <- patient_n[rows_now[1], "amount"] / tmp
    patient_n[rows_now[2], "percent"] <- patient_n[rows_now[2], "amount"] / tmp


    # Retrieve non-NA entries for final vital status
    kaplan_now <- kaplan_now[(!kaplan_now$final_vital_status == "unknown"), ]
    tmp <- kaplan_now$final_vital_status == "vital_status"
    kaplan_now <- kaplan_now[!tmp, ]
    kaplan_now <- kaplan_now[!is.na(kaplan_now$final_vital_status), ]
    kaplan_now <- kaplan_now[!is.nan(kaplan_now$final_vital_status), ]
    tmp <- kaplan_now$final_vital_status_times == "unknown"
    kaplan_now <- kaplan_now[!tmp, ]
    kaplan_now <- kaplan_now[(!kaplan_now$final_vital_status_times == "day"), ]
    kaplan_now <- kaplan_now[!is.na(kaplan_now$final_vital_status_times), ]
    kaplan_now <- kaplan_now[!is.nan(kaplan_now$final_vital_status_times), ]
    kaplan_now$final_vital_status_times <- as.numeric(
        kaplan_now$final_vital_status_times
    )

    kaplan_now$final_vital_status <- ifelse(
        kaplan_now$final_vital_status == "Dead", 1, 0
    )

    kaplan_now$expression <- log2(kaplan_now[, name] + 1)

    # Fill table with patient amount - OSurv
    rows_now <- c(3, 4)
    patient_n[rows_now, "type"] <- "OSurv"
    patient_n[rows_now, "classification"] <- c("low", "high")
    patient_n[rows_now[1], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["low"]
    )
    patient_n[rows_now[2], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["high"]
    )

    patient_n[rows_now[1], "percent"] <- patient_n[rows_now[1], "amount"] / tmp
    patient_n[rows_now[2], "percent"] <- patient_n[rows_now[2], "amount"] / tmp

    cox_ph(step = "overall")
    surv_plot(step = "overall")


    # NOTE RFS
    kaplan_now <- framel[, c(
        name,
        "classification",
        "final_vital_status_times",
        "final_vital_status",
        "final_rfs_status",
        "final_rfs_status_times",
        "final_dmfs_status",
        "final_dmfs_status_times"
    )]

    # Kaplan
    colnames(kaplan_now)[2] <- "classification"

    # Collect just rows with classification
    kaplan_now <- kaplan_now[c(
        which(kaplan_now[, "classification"] == "low"),
        which(kaplan_now[, "classification"] == "high")
    ), ]

    # Factorize classification
    kaplan_now$group_factor <- factor(
        kaplan_now$classification,
        levels = c("low", "high")
    )

    kaplan_now <- kaplan_now[(!kaplan_now$final_rfs_status == "unknown"), ]
    kaplan_now <- kaplan_now[!is.na(kaplan_now$final_rfs_status), ]
    kaplan_now <- kaplan_now[!is.nan(kaplan_now$final_rfs_status), ]

    kaplan_now <- kaplan_now[!kaplan_now$final_rfs_status_times == "unknown", ]
    tmp <- kaplan_now$final_rfs_status_times == "dead_other_cause_before_relapse"
    kaplan_now <- kaplan_now[!tmp, ]
    kaplan_now <- kaplan_now[(!kaplan_now$final_rfs_status_times == "day"), ]
    kaplan_now <- kaplan_now[!is.na(kaplan_now$final_rfs_status_times), ]
    kaplan_now <- kaplan_now[!is.nan(kaplan_now$final_rfs_status_times), ]
    kaplan_now <- kaplan_now[!is.infinite(kaplan_now$final_rfs_status_times), ]
    kaplan_now$final_rfs_status_times <- as.numeric(
        kaplan_now$final_rfs_status_times
    )

    kaplan_now$final_rfs_status <- ifelse(
        kaplan_now$final_rfs_status == "relapse", 1, 0
    )

    kaplan_now$expression <- log2(kaplan_now[, name] + 1)

    # Fill table with patient amount - Relapse
    rows_now <- c(5, 6)
    patient_n[rows_now, "type"] <- "Relapse"
    patient_n[rows_now, "classification"] <- c("low", "high")
    patient_n[rows_now[1], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["low"]
    )
    patient_n[rows_now[2], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["high"]
    )

    tmp <- patient_n[rows_now[1], "amount"] + patient_n[rows_now[2], "amount"]
    patient_n[rows_now[1], "percent"] <- patient_n[rows_now[1], "amount"] / tmp
    patient_n[rows_now[2], "percent"] <- patient_n[rows_now[2], "amount"] / tmp
    cox_ph(step = "rfs")
    surv_plot(step = "rfs")


    # NOTE dmfs

    # dmfs
    # Kaplan meyer using best scenario
    kaplan_now <- framel[, c(
        name,
        "classification",
        "final_vital_status_times",
        "final_vital_status",
        "final_rfs_status",
        "final_rfs_status_times",
        "final_dmfs_status",
        "final_dmfs_status_times"
    )]

    # Kaplan
    colnames(kaplan_now)[2] <- "classification"

    # Collect just rows with classification
    kaplan_now <- kaplan_now[c(
        which(kaplan_now[, "classification"] == "low"),
        which(kaplan_now[, "classification"] == "high")
    ), ]

    # Factorize classification
    kaplan_now$group_factor <- factor(
        kaplan_now$classification,
        levels = c("low", "high")
    )

    # Retrieve non-NA entries for final vital status
    kaplan_now <- kaplan_now[(!kaplan_now$final_dmfs_status == "unknown"), ]
    kaplan_now <- kaplan_now[!is.na(kaplan_now$final_dmfs_status), ]
    kaplan_now <- kaplan_now[!is.nan(kaplan_now$final_dmfs_status), ]
    kaplan_now <- kaplan_now[!kaplan_now$final_dmfs_status_times == "unknown", ]
    tmp <- kaplan_now$final_dmfs_status_times == "dead_other_cause_before_relapse"
    kaplan_now <- kaplan_now[!tmp, ]
    kaplan_now <- kaplan_now[(!kaplan_now$final_dmfs_status_times == "day"), ]
    kaplan_now <- kaplan_now[!is.na(kaplan_now$final_dmfs_status_times), ]
    kaplan_now <- kaplan_now[!is.nan(kaplan_now$final_dmfs_status_times), ]
    kaplan_now <- kaplan_now[!is.infinite(kaplan_now$final_dmfs_status_times), ]
    kaplan_now$final_dmfs_status_times <- as.numeric(
        kaplan_now$final_dmfs_status_times
    )

    # Overall survival
    kaplan_now$final_dmfs_status <- as.integer(gsub(
        "metastasis", 1, gsub("censored", 0, kaplan_now$final_dmfs_status)
    ))

    kaplan_now$expression <- log2(kaplan_now[, name] + 1)


    # Fill table with patient amount - DMFS
    rows_now <- c(7, 8)
    patient_n[rows_now, "type"] <- "DMFS"
    patient_n[rows_now, "classification"] <- c("low", "high")
    patient_n[rows_now[1], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["low"]
    )
    patient_n[rows_now[2], "amount"] <- as.numeric(
        table(kaplan_now[, "classification"])["high"]
    )

    tmp <- patient_n[rows_now[1], "amount"] + patient_n[rows_now[2], "amount"]
    patient_n[rows_now[1], "percent"] <- patient_n[rows_now[1], "amount"] / tmp
    patient_n[rows_now[2], "percent"] <- patient_n[rows_now[2], "amount"] / tmp

    if (sum(kaplan_now$final_dmfs_status) > 0) {
        cox_ph(step = "dmfs")
        surv_plot(step = "dmfs")
    }


    # NOTE PatientAmount
    patient_n$type <- factor(
        patient_n$type,
        levels = c("All", "OSurv", "Relapse", "DMFS")
    )
    patient_n$classification <- factor(
        patient_n$classification,
        levels = c("high", "low")
    )

    # Add position for plotting
    patient_n$pos <- 0
    for (type in levels(patient_n$type)) {
        lines <- patient_n$type == type
        numbers <- patient_n$amount[lines]
        patient_n$pos[lines] <- cumsum(numbers) - (0.5 * numbers)
    }

    # Add percents
    patient_n$percent <- (patient_n$percent * 100)
    patient_n$percent <- round(patient_n$percent, 0)

    save_plot("kaplan_meier", "_PatientAmount")
    p4 <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_bar(ggplot2::aes(
            y = amount, x = type,
            fill = classification
        ),
        data = patient_n,
        stat = "identity"
        ) +
        ggplot2::geom_text(
            data = patient_n,
            ggplot2::aes(
                x = type, y = pos,
                label = paste0(percent, "%")
            ),
            size = 3.6
        ) +
        ggplot2::scale_fill_manual(values = c(
            hsv(.008, .82, .84, .86),
            hsv(.61, .66, .78, .70)
        )) +
        ggplot2::theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(
                colour = "black",
                size = 14
            ),
            axis.text.y = ggplot2::element_text(
                colour = "black",
                size = 13
            )
        )
    print(p4)
    dev.off()


    # Save table with patient amounts
    write.csv(patient_n, file.path(
        dir, "kaplan_meier",
        paste0(
            "PatientAmount_",
            name, ".csv"
        )
    ))

    clinical_groups <- framel[, c("group", "classification")]
    rownames(clinical_groups) <- framel[, "patient_code"]

    assign("clinical_groups", clinical_groups, envir = get(ev))

    write.csv(clinical_groups, file = file.path(dir, "clinical_groups"))

    gc()
    message("Done!\n")
}
