#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_expression} is a function designed to concatenate GDC
#' files into a single matrix, where the columns stand for patients code and
#' rows stand for data names.
#'
#' @param data_type
#' @param normalization Logical value where \code{TRUE} specify the desire to
#'    work with normalized files only. When FALSE, in the second run, do not
#'    forget to set env argument. This argument is only applyable to gene and
#'    isoform expression data from GDC Legacy Archive. The default is
#'    \code{TRUE}.
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
#'    \itemize{\item{'<tumor>_<data_base>_<data_type>_tumor_data'}{ or}
#'    \item{'<tumor>_<data_base>_<data_type>__both_data'}{ (for tumor and not
#'    tumor data in separated matrices).}}
#' @param save_data Logical value where \code{TRUE} indicates that the
#'    concatenate and filtered matrix should be saved in local storage. The
#'    default is FALSE.
#' @inheritParams download_gdc
#'
#' @return A matrix with data names in row and patients code in column.
#' @export
#'
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom reshape2 melt
#' @importFrom XML xmlToDataFrame
#'
#' @examples
#' library(DOAGDC)
#'
#' # Concatenating gene expression data into a single matrix
#' # data already downloaded using the 'download_gdc' function
#' concatenate_expression("gene",
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
concatenate_expression <- function(data_type,
                                    normalization = TRUE,
                                    name, data_base,
                                    work_dir, tumor,
                                    tumor_data = TRUE,
                                    only_filter = FALSE,
                                    tumor_type = 1,
                                    normal_type = 11,
                                    env,
                                    save_data = FALSE) {

    # local functions ####
    gene_isoform <- function(df) {
        # selecting files by type (normalized or not)
        regex <- ifelse(
            normalization,
            ifelse(
                data_type_boo, "^.+(genes.normalized).+$",
                "^.+(isoforms.normalized).+$"
            ),
            ifelse(
                data_type_boo, "^.+(.rsem.genes.results)$",
                "^.+(isoforms.results)$",
            )
        )

        return(df[grepl(regex, df$file_name), ])
    }

    open_genexpress <- function(files) {
        message("Reading data...")
        # reading files in one file
        pb <- txtProgressBar(min = 0, max = length(files$file_name), style = 3)
        count <- 0
        for (i in files$file_name) {
            count <- count + 1
            setTxtProgressBar(pb, count)
            actual_file <- data.table::fread(file.path(direc, i))
            if (count == 1) {
                actual_file <- data.table::fread(file.path(direc, i))
                completed_table1 <- matrix(
                    nrow = nrow(actual_file),
                    ncol = length(files$file_name)
                )
                rownames(completed_table1) <- actual_file[[1]]
                completed_table1[, 1] <- as.numeric(actual_file[[2]])
            } else {
                actual_file <- data.table::fread(file.path(direc, i), select = 2)
                completed_table1[, count] <- as.numeric(actual_file[[1]])
            }
        }
        close(pb)

        return(completed_table1)
    }

    name_finder <- function(completed_table) {
        to_google <- ifelse(suppressWarnings(is.na(as.numeric(name))),
            paste0("^", toupper(name)),
            paste0(name, "$")
        )
        gene_selector_final <- grep(to_google,
            rownames(completed_table),
            perl = TRUE
        )

        selecionado <- completed_table[ifelse(
            data_base_boo && data_type_boo, gene_selector_final, name
        ), , drop = FALSE]

        stop_caller(
            nrow(selecionado),
            paste0(message("'", name, "' not found!!!\n"))
        )

        express_transpo <- t(selecionado)
        return(express_transpo)
    }

    stop_caller <- function(size, texto) {
        if (size == 0) {
            stop(message(texto))
        }
    }

    fpkm2_tpm <- function(fpkm) {
        exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    }

    # code ####
    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    data_base_boo <- tolower(data_base) == "legacy"
    data_type_boo <- tolower(data_type) == "gene"

    if (!only_filter) {
        assign(paste(toupper(tumor), toupper(data_base),
            gsub(" ", "_", tolower(data_type)),
            ifelse(tumor_data, "tumor_data", "both_data"),
            sep = "_"
        ),
        new.env(parent = emptyenv()),
        envir = .GlobalEnv
        )

        ev <- paste(toupper(tumor), toupper(data_base),
            gsub(" ", "_", tolower(data_type)),
            ifelse(tumor_data, "tumor_data", "both_data"),
            sep = "_"
        )

        attr(ev, "name") <- paste0(
            "Environment created by DOAGDC",
            " package, use its name in ",
            "'env' argument"
        )

        direc <- file.path(
            work_dir, "DOAGDC", toupper(tumor),
            ifelse(
                data_base_boo, ifelse(
                    data_type_boo, "gene_data", "isoform_data"
                ), "gdc_gene_data"
            )
        )

        codes <- dir(direc, ".sdrf$")

        codigos <- read.table(file.path(direc, codes[1]),
            stringsAsFactors = FALSE,
            header = TRUE, sep = "\t"
        )

        codigos <- codigos[, c("file_name", ifelse(
            data_base_boo, "submitter_id", "cases"
        ))]

        regex <- ifelse(length(tumor_type) > 1,
            paste0("(", paste(formatC(tumor_type, width = 2, flag = "0"),
                collapse = "|"
            ), ")", "[A-Z]"),
            paste0(formatC(tumor_type, width = 2, flag = "0"), "[A-Z]")
        )

        if (data_base_boo) {

            seletor <- unname(sapply(
                codigos$submitter_id,
                function(w) {
                    paste0(unlist(strsplit(w, "-"))[4])
                }
            ))

            codigos_t <- gene_isoform(codigos)
            codigos_t <- codigos_t[grepl(regex, x = seletor), ]
            patients <- codigos_t$submitter_id
        } else {

            seletor <- unname(sapply(
                codigos$cases,
                function(w) {
                    paste0(unlist(strsplit(w, "-"))[4])
                }
            ))

            codigos_t <- codigos[grepl(regex, x = seletor), ]

            file_end <- ifelse(normalization, "FPKM", "count")
            codigos_t <- codigos_t[grep(
                file_end, codigos_t$file_name,
                perl = TRUE
            ), ]
            patients <- codigos_t$cases
        }

        remove_duplicated <- duplicated(codigos_t$file_name)
        codigos_t <- codigos_t[!remove_duplicated, ]
        rownames(codigos_t) <- codigos_t$file_name

        completed_table <- open_genexpress(files = codigos_t)
        colnames(completed_table) <- patients
        rownames(completed_table) <- gsub(
            "\\.{1}\\d+$", "", rownames(completed_table)
        )

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
            colnames(completed_table)[coluns_2rename] <- names_2rename
            colnames(completed_table)[-coluns_2rename] <- unique(patients_short)
        } else {
            colnames(completed_table) <- patients_short
        }

        if (!tumor_data) {

            regex <- ifelse(length(normal_type) > 1,
                paste0("(", paste(formatC(normal_type, width = 2, flag = "0"),
                    collapse = "|"
                ), ")", "[A-Z]"),
                paste0(formatC(normal_type, width = 2, flag = "0"), "[A-Z]")
            )

            if (data_base_boo) {
                seletor <- unname(sapply(
                    codigos$submitter_id,
                    function(w) {
                        paste0(unlist(strsplit(w, "-"))[4])
                    }
                ))

                codigos_nt <- codigos[grepl(regex, x = seletor), ]
                patients_nt <- codigos_nt$submitter_id
            } else {
                seletor <- unname(sapply(
                    codigos$cases,
                    function(w) {
                        paste0(unlist(strsplit(w, "-"))[4])
                    }
                ))

                codigos_nt <- codigos[grepl(regex, x = seletor), ]
                patients_nt <- codigos_nt$cases
            }

            remove_duplicated_nt <- duplicated(
                codigos_nt$file_name
            )
            codigos_nt <- codigos_nt[!remove_duplicated_nt, ]
            rownames(codigos_nt) <- codigos_nt$file_name

            stop_caller(
                nrow(codigos_nt),
                paste0(
                    "There is no 'Normal' data in ",
                    tumor,
                    " tumor folder! Please, rerun after set",
                    " 'tumor_data = TRUE'.\n\n"
                )
            )

            completed_nt_table <- open_genexpress(codigos_nt)
            colnames(completed_nt_table) <- patients_nt
            rownames(completed_nt_table) <- gsub(
                "\\.{1}\\d+$", "", rownames(completed_nt_table)
            )

            patients_short <- unname(sapply(
                patients_nt,
                function(w) {
                    paste(unlist(strsplit(w, "-"))[1:3],
                        collapse = "-"
                    )
                }
            ))

            # duplicate fix
            if (length(patients_nt) != length(unique(
                patients_short
            ))) {
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
                colnames(
                    completed_nt_table
                )[coluns_2rename] <- names_2rename
                colnames(
                    completed_nt_table
                )[-coluns_2rename] <- unique(patients_short)
            } else {
                colnames(completed_nt_table) <- patients_short
            }

            if (save_data) {
                message("\nSaving your data...\n")
                write.table(table_tumor,
                    paste0(
                        direc, "/", tumor,
                        "_tumor_data.tsv"
                    ),
                    sep = "\t"
                )
                write.table(
                    completed_nt_table,
                    paste0(
                        direc, "/",
                        tumor, "_nt_data.tsv"
                    ),
                    sep = "\t"
                )
            }
            assign("patients_nt", patients_nt, envir = get(ev))

            # export TPM values
            if (normalization) {
                completed_nt_table <- apply(completed_nt_table, 2, fpkm2_tpm)
            }

            assign(paste0(
                tolower(data_type), "_nt_",
                ifelse(
                    normalization,
                    "normalized",
                    "nn"
                )
            ),
            completed_nt_table,
            envir = get(ev)
            )

            if (!missing(name)) {
                nt_selected <- name_finder(completed_nt_table)

                assign(paste0(
                    tolower(data_type), "_nt_",
                    ifelse(
                        normalization,
                        "normalized_selected_",
                        "nn_selected_"
                    ),
                    toupper(name)
                ),
                nt_selected,
                envir = get(ev)
                )
            }
        }

        assign("patients", patients, envir = get(ev))

        # export TPM values
        if (normalization) {
            completed_table <- apply(completed_table, 2, fpkm2_tpm)
        }

        assign(paste0(
            tolower(data_type), "_tumor_",
            ifelse(
                normalization,
                "normalized",
                "nn"
            )
        ),
        completed_table,
        envir = get(ev)
        )

        assign("work_dir", work_dir, envir = get(ev))

        if (!missing(name)) {

            express_transpo <- name_finder(completed_table)

            assign(paste0(
                tolower(data_type), "_tumor_",
                ifelse(
                    normalization,
                    "normalized_selected_",
                    "nn_selected_"
                ),
                toupper(name)
            ),
            express_transpo,
            envir = get(ev)
            )

            assign("name", name, envir = get(ev))
            assign("name_e", name, envir = get(ev))
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

        if (data_type_boo) {
            if (normalization) {
                table_tumor <- name_finder(sv[["ev"]]$gene_tumor_normalized)
            } else {
                table_tumor <- name_finder(sv[["ev"]]$gene_tumor_nn)
            }

            assign(paste0(
                ifelse(
                    normalization,
                    "gene_tumor_normalized_selected_",
                    "gene_tumor_nn_selected_"
                ),
                toupper(name)
            ),
            table_tumor,
            envir = get(ev)
            )

            if (!tumor_data) {
                # not tumor
                if (normalization) {
                    table_nt <- name_finder(sv[["ev"]]$gene_nt_normalized)
                } else {
                    table_nt <- name_finder(sv[["ev"]]$gene_nt_nn)
                }

                assign(paste0(
                    ifelse(
                        normalization,
                        "gene_nt_normalized_selected_",
                        "gene_nt_nn_selected_"
                    ),
                    toupper(name)
                ),
                table_nt,
                envir = get(ev)
                )

            }
        } else if (tolower(data_type) == "isoform") {
            if (normalization) {
                table_tumor <- name_finder(
                    sv[["ev"]]$isoform_tumor_normalized
                )
            } else {
                table_tumor <- name_finder(sv[["ev"]]$isoform_tumor_nn)
            }

            assign(paste0(
                ifelse(
                    normalization,
                    "isoform_tumor_normalized_selected_",
                    "isoform_tumor_nn_selected_"
                ),
                toupper(name)
            ),
            table_tumor,
            envir = get(ev)
            )

            if (!tumor_data) {
                # not tumor
                if (normalization) {
                    table_nt <- name_finder(sv[["ev"]]$isoform_nt_normalized)
                } else {
                    table_nt <- name_finder(sv[["ev"]]$isoform_nt_nn)
                }

                assign(paste0(
                    ifelse(
                        normalization,
                        "isoform_nt_normalized_selected_",
                        "isoform_nt_nn_selected_"
                    ),
                    toupper(name)
                ),
                table_nt,
                envir = get(ev)
                )
            }
        }

        assign("name", name, envir = get(ev))
        assign("name_e", name, envir = get(ev))
    }
}
