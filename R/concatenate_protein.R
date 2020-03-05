#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_protein} is a function designed to concatenate GDC files
#' into a single matrix, where the columns stand for patients code and rows
#' stand for data names.
#'
#' @param name A character string indicating the desired values to be used in
#'    next analysis. For instance, "HIF3A" in the legacy gene expression matrix,
#'    "mir-1307" in the miRNA quantification matrix, or "HER2" in the protein
#'    quantification matrix.
#' @param data_base
#' @param work_dir
#' @param tumor
#' @param tumor_data Logical value where \code{TRUE} specifies the desire to
#'    work
#'    with tumor tissue files only. When set to FALSE, it creates two matrices,
#'    one containing tumor data and other containing data from not-tumor tissue.
#'    The default is \code{TRUE}.
#' @param only_filter Logical value where \code{TRUE} indicates that the matrix
#'    is already concatenate and the function should choose a different
#'    \code{name}, without concatenate all the files again. The default is
#'    FALSE.
#' @param tumor_type Numerical value(s) correspondent to barcode data types:
#'    \itemize{Tumor codes: \item{1: Primary Solid Tumor}\item{2: Recurrent
#'    Solid Tumor}\item{3: Primary Blood Derived Cancer - Peripheral
#'    Blood}\item{4: Recurrent Blood Derived Cancer - Bone Marrow}\item{5:
#'    Additional - New Primary}\item{6: Metastatic}\item{7: Additional
#'    Metastatic}\item{8: Human Tumor Original Cells}\item{9: Primary Blood
#'    Derived Cancer - Bone Marrow}} The default is 1.
#' @param normal_type Numerical value(s) correspondent to barcode data types:
#'    \itemize{Normal codes: \item{10: Blood Derived Normal}\item{11: Solid
#'    Tissue Normal}\item{12: Buccal Cell Normal}\item{13: EBV Immortalized
#'    Normal}\item{14: Bone Marrow Normal}\item{15: sample type 15}\item{16-19:
#'    sample type 16}}{ or}\itemize{Control codes: \item{use '20:29' without
#'    quotes}} The default is 11.
#' @param env A character string containing the environment name that should be
#'    used. If none has been set yet, the function will create one in global
#'    environment following the standard criteria:
#'    \itemize{\item{'<tumor>_<data_base>_protein_tumor_data'}{ or}
#'    \item{'<tumor>_<data_base>_protein_both_data'}{ (for tumor and not tumor
#'    data
#'    in separated matrices).}}
#' @param save_data Logical value where \code{TRUE} indicates that the
#'    concatenate and filtered matrix should be saved in local storage. The
#'    default is FALSE.
#' @inheritParams download_gdc
#'
#' @return A matrix with data names in row and patients code in column.
#' @export
#'
#' @importFrom data.table fread
#'
#' @examples
#' library(DOAGDC)
#'
#' # Concatenating gene expression data into a single matrix
#' # data already downloaded using the 'download_gdc' function
#' concatenate_protein(
#'     name = "Caspase-8-M-E",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
concatenate_protein <- function(name, data_base,
                                work_dir, tumor,
                                tumor_data = TRUE,
                                only_filter = FALSE,
                                tumor_type = 1,
                                normal_type = 11,
                                env,
                                save_data = FALSE) {

    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    if (!only_filter) {
        if (tolower(data_base) == "gdc") {
            stop(message(
                "\nThrere is no protein expression ",
                "data in GDC data base!!",
                "\nPlease use 'legacy' data base"
            ))
        }

        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), "protein_data")

        if (!dir.exists(file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "mage_data"
        ))) {
            message("Downloading magetab data...")
            suppressWarnings(download_gdc(
                data_type = "mage",
                data_base = data_base,
                tumor = tumor,
                work_dir = work_dir
            ))
        }

        desing_array_dir <- file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "mage_data"
        )

        to_gunzip <- file.path(
            desing_array_dir,
            paste0(
                "mdanderson.org_", toupper(tumor),
                ".MDA_RPPA_Core.mage-tab.1.1.0.tar.gz"
            )
        )
        ## or, if you just want to extract the target file:
        untar(to_gunzip, exdir = desing_array_dir)

        design_array <- data.table::fread(file.path(
            desing_array_dir,
            paste0(
                "mdanderson.org_", toupper(tumor),
                ".MDA_RPPA_Core.mage-tab.1.1.0"
            ),
            paste0(
                "mdanderson.org_", toupper(tumor),
                ".MDA_RPPA_Core.array_design.txt"
            )
        ),
        select = c(6, 7)
        )

        design_array <- unique(design_array)

        rownames(design_array) <- design_array[[2]]

        patient_files <- dir(path = dir, pattern = "protein_expression")

        pb <- txtProgressBar(
            min = 0, max = length(patient_files),
            style = 3
        )
        count <- 0
        for (i in file.path(dir, patient_files)) {
            count <- count + 1
            if (count == 1) {
                patient_id <- character(length = length(patient_files))
                pr <- data.table::fread(i)
                uuid <- colnames(pr)[2]
                protein_table <- matrix(
                    ncol = length(patient_files),
                    nrow = nrow(pr) - 1
                )
                colnames(protein_table) <- seq_len(ncol(protein_table))
                patient_id[1] <- as.character(design_array[[1]][grep(
                    uuid,
                    design_array[[2]]
                )])
                pr <- pr[-1, ]
                rownames(protein_table) <- as.character(pr[[1]])
                protein_table[, 1] <- as.numeric(pr[[2]])
                colnames(protein_table)[1] <- patient_id[1]
            } else {
                pr <- data.table::fread(i, select = 2)
                uuid <- colnames(pr)[1]
                patient_id[count] <- as.character(
                    design_array[[1]][grep(uuid, design_array[[2]])]
                )
                pr <- pr[-1, ]
                protein_table[, count] <- as.numeric(pr[[1]])
                colnames(protein_table)[count] <- patient_id[count]
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)

        # separing files
        seletor <- unname(sapply(
            colnames(protein_table),
            function(w) {
                paste(unlist(strsplit(w, "-"))[4], collapse = "-")
            }
        ))

        regex <- ifelse(length(tumor_type) > 1,
            paste0("(", paste(formatC(tumor_type, width = 2, flag = "0"),
                collapse = "|"
            ), ")", "[A-Z]"),
            paste0(formatC(tumor_type, width = 2, flag = "0"), "[A-Z]")
        )

        only_tumor_tissue <- grepl(
            pattern = regex,
            seletor
        )

        # not related the control samples
        patients <- as.character(colnames(protein_table)[only_tumor_tissue])
        protein_table <- protein_table[, patients]

        patients_short <- unname(sapply(patients, function(w) {
            paste(unlist(strsplit(w, "-"))[1:3], collapse = "-")
        }))
        # duplicate fix
        if (length(patients) != length(unique(patients_short))) {
            coluns_2rename <- numeric()
            names_2rename <- character()
            for (Patients in unique(patients_short)) {
                selector <- Patients == patients_short
                if (sum(selector) > 1) {
                    coluns_2rename <- c(
                        coluns_2rename,
                        grep(
                            Patients,
                            patients_short
                        )[-1]
                    )
                    names_2rename <- c(
                        names_2rename,
                        paste0(
                            Patients,
                            seq((sum(selector) - 1))
                        )
                    )
                }
            }
            colnames(protein_table)[coluns_2rename] <- names_2rename
            colnames(protein_table)[-coluns_2rename] <- unique(
                patients_short
            )
        } else {
            colnames(protein_table) <- patients_short
        }

        protein_table <- gsub("-.-.", "", rownames(protein_table))

        if (save_data) {
            message("\nSaving your data...\n")
            # saving
            write.table(
                x = protein_table,
                file = file.path(
                    dir,
                    "protein_tumor_normalized.tsv"
                ),
                row.names = TRUE, quote = FALSE, sep = "\t"
            )
        }


        assign(paste(toupper(tumor), toupper(data_base),
            "methylation",
            ifelse(tumor_data, "tumor_data", "both_data"),
            sep = "_"
        ),
        new.env(parent = emptyenv()),
        envir = .GlobalEnv
        )

        ev <- paste(toupper(tumor), toupper(data_base),
            "methylation",
            ifelse(tumor_data, "tumor_data", "both_data"),
            sep = "_"
        )

        sv <- list(ev = get(ev))
        attr(ev, "name") <- paste0(
            "Environment created by DOAGDC",
            " package, use its name in ",
            "'env' argument"
        )

        assign("protein_patients", patients, envir = get(ev))
        assign("protein_tumor_normalized", protein_table,
            envir = get(ev)
        )

        if (!missing(name)) {
            protein_selector <- rownames(
                sv[["ev"]]$protein_tumor_normalized
            )
            protein_selector_final <- grepl(paste0("^", name, "$"),
                protein_selector,
                perl = TRUE
            )

            protein_exp_selected <- t(
                sv[["ev"]]$protein_tumor_normalized[protein_selector_final, ]
            )
            name <- gsub("-", "_", name)
            colnames(protein_exp_selected) <- name
            assign(paste0(
                "protein_tumor_normalized_selected_",
                toupper(name)
            ),
            protein_exp_selected,
            envir = get(ev)
            )
        }

        if (!tumor_data) {
            seletor <- unname(sapply(
                colnames(protein_table),
                function(w) {
                    paste(unlist(strsplit(w, "-"))[4], collapse = "-")
                }
            ))

            regex <- ifelse(length(normal_type) > 1,
                paste0("(", paste(formatC(normal_type, width = 2, flag = "0"),
                    collapse = "|"
                ), ")", "[A-Z]"),
                paste0(formatC(normal_type, width = 2, flag = "0"), "[A-Z]")
            )

            not_tumor_tissue <- grepl(
                pattern = regex,
                seletor
            )

            if (sum(not_tumor_tissue) == 0) {
                stop(message(
                    "There is no 'Normal' data in ", tumor,
                    " tumor folder! Please, rerun after",
                    " set 'tumor_data = TRUE'.",
                    "\n\n"
                ))
            }
            not_tumor_patients <- as.character(
                colnames(protein_table)[not_tumor_tissue]
            )
            protein_not_tumor <- protein_table[, not_tumor_patients]
            protein_not_tumor <- gsub(
                "-.-.", "",
                rownames(protein_not_tumor)
            )

            assign("not_tumor_patients", not_tumor_patients,
                envir = get(ev)
            )
            assign("protein_not_tumor_normalized", protein_not_tumor,
                envir = get(ev)
            )

            if (save_data) {
                message("\nSaving your data...\n")
                # saving
                write.table(
                    x = tumor_protein_table,
                    file = file.path(
                        dir,
                        "protein_tumor_normalized.tsv"
                    ),
                    row.names = TRUE, quote = FALSE, sep = "\t"
                )
                write.table(
                    x = protein_not_tumor,
                    file = file.path(
                        dir,
                        "protein_not_tumor_normalized.tsv"
                    ),
                    row.names = TRUE, quote = FALSE, sep = "\t"
                )
            }

            if (!missing(name)) {
                protein_selector <- rownames(
                    sv[["ev"]]$protein_not_tumor_normalized
                )
                protein_selector_final <- grepl(paste0("^", name, "$"),
                    protein_selector,
                    perl = TRUE
                )

                protein_exp_selected <- t(sv[["ev"]]$protein_not_tumor_normalized[protein_selector_final, ])
                name <- gsub("-", "_", name)
                colnames(protein_exp_selected) <- name
                assign(paste0(
                    "protein_not_tumor_normalized_selected_",
                    toupper(name)
                ),
                protein_exp_selected,
                envir = get(ev)
                )
            }
        }
    } else {

        if (missing(env)) {
            stop(message(
                "Please, before using 'only_filter'",
                " argument, insert the Environment name."
            ))
        }

        ev <- deparse(substitute(env))

        sv <- list(ev = get(ev))
        attr(ev, "name") <- paste0(
            "Environment created by DOAGDC",
            " package, use its name in ",
            "'env' argument"
        )

        if (!missing(name)) {
            protein_selector <- rownames(
                sv[["ev"]]$protein_tumor_normalized
            )
            protein_selector_final <- grepl(paste0("^", name, "$"),
                protein_selector,
                perl = TRUE
            )

            protein_exp_selected <- t(
                sv[["ev"]]$protein_tumor_normalized[protein_selector_final, , drop = FALSE]
            )
            name <- gsub("-", "_", name)
            colnames(protein_exp_selected) <- name
            assign(paste0(
                "protein_tumor_normalized_selected_",
                toupper(name)
            ),
            protein_exp_selected,
            envir = get(ev)
            )
            if (!tumor_data) {
                protein_selector <- rownames(
                    sv[["ev"]]$protein_not_tumor_normalized
                )
                protein_selector_final <- grepl(paste0("^", name, "$"),
                    protein_selector,
                    perl = TRUE
                )

                protein_exp_selected <- t(
                    sv[["ev"]]$protein_not_tumor_normalized[protein_selector_final, , drop = FALSE]
                )
                name <- gsub("-", "_", name)
                colnames(protein_exp_selected) <- name
                assign(paste0(
                    "protein_not_tumor_normalized_selected_",
                    toupper(name)
                ),
                protein_exp_selected,
                envir = get(ev)
                )
            }
        }
    }
}