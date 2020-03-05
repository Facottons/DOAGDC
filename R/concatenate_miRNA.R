#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_mirna} is a function designed to concatenate GDC files into
#' a single matrix, where the columns stand for patients code and rows stand for
#' data names.
#'
#' @param data_type A character string indicating the desired miRNA file:
#     \itemize{\item{"mirna gene quantification"}, \item{"mirna expression
#     quantification"}}
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
#' @param platform
#' @param use_hg19_mirbase20 Logical value where \code{TRUE} indicates that only
#'    hg19.mirbase20 should be used. This parameter is needed when using
#'    \code{data_base = "legacy"} and one of the available miRNA
#'    \code{data_type}
#'    in "legacy" ("miRNA gene quantification" and "miRNA isoform
#'    quantification"). The default is FALSE.
#' @param env A character string containing the environment name that should be
#'    used. If none has been set yet, the function will create one in global
#'    environment following the standard criteria:
#'    \itemize{\item{'tumor_data_base_data_type_tumor_data'}{ or}
#'    \item{'tumor_data_base_data_type_both_data'}{ (for tumor and not tumor
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
#' concatenate_mirna(
#'     name = "hsa-mir-21",
#'     use_hg19_mirbase20 = TRUE,
#'     data_base = "gdc",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
concatenate_mirna <- function(data_type,
                            normalization = TRUE,
                            name, data_base,
                            work_dir, tumor,
                            tumor_data = TRUE,
                            only_filter = FALSE,
                            tumor_type = 1,
                            normal_type = 11,
                            platform = "",
                            use_hg19_mirbase20 = FALSE,
                            env,
                            save_data = FALSE) {

    # local functions ####
    open_mirna <- function(mirna_files_arg) {
        if (normalization) {
            select <- c(3, 4)
        } else {
            select <- c(2, 4)
        }
        # reading files in one file
        pb <- txtProgressBar(min = 0, max = length(mirna_files_arg), style = 3)
        count <- 0
        # try to open with lapply(list, function)
        for (file in file.path(dir, mirna_files_arg)) {
            count <- count + 1
            if (count > 1) {
                actual_file <- data.table::fread(file,
                    sep = "\t",
                    showProgress = FALSE,
                    select = select
                )
                data_mirna_table[, count] <- as.numeric(actual_file[[1]])
                tmp <- as.character(actual_file[[2]])
                cross_mapped_mirna_table[, count] <- tmp
            } else if (count == 1) {
                actual_file <- data.table::fread(file,
                    sep = "\t",
                    showProgress = FALSE
                )
                tmp <- as.character(actual_file[["miRNA_ID"]])
                data_mirna_table <- matrix(
                    nrow = nrow(actual_file),
                    dimnames = list(
                        tmp,
                        manifest$cases
                    ),
                    ncol = length(mirna_files_arg)
                )
                cross_mapped_mirna_table <- matrix(
                    nrow = nrow(actual_file),
                    dimnames = list(
                        tmp,
                        manifest$cases
                    ),
                    ncol = length(mirna_files_arg)
                )
                data_mirna_table[, 1] <- as.numeric(actual_file[["read_count"]])
                tmp <- as.character(actual_file[["cross-mapped"]])
                cross_mapped_mirna_table[, 1] <- tmp
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)

        assign("data_mirna_table", data_mirna_table, envir = get(ev))
        assign("cross_mapped_mirna_table", cross_mapped_mirna_table,
            envir = get(ev)
        )
    }

    # code ####
    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    tmp1 <- "mirna" %in% strsplit(tolower(data_type), split = " ")[[1]][1]
    tmp2 <- "isoform expression quantification" %in% tolower(data_type)
    if (tmp1 || tmp2) {
        if (!only_filter) {
            if (tolower(data_base) == "legacy") {
                dir <- file.path(
                    work_dir, "DOAGDC", toupper(tumor),
                    paste0(tolower(data_type), "_data")
                )

                manifest <- data.table::fread(file.path(dir, "manifest.sdrf"),
                    select = c(
                        "file_name", "cases",
                        "platform"
                    )
                )

                tmp <- !grepl(
                    pattern = paste(".sdrf", "Data_access_time.txt",
                        sep = "|"
                    ),
                    x = dir(path = dir)
                )

                mirna_files_downloaded <- dir(path = dir)[tmp]
                if (use_hg19_mirbase20) {
                    mirna_files_downloaded <- mirna_files_downloaded[grepl(
                        "hg19.mirbase20", mirna_files_downloaded,
                        fixed = TRUE
                    )]
                } else {
                    mirna_files_downloaded <- mirna_files_downloaded[!grepl(
                        "hg19.mirbase20", mirna_files_downloaded,
                        fixed = TRUE
                    )]
                }
                tmp <- manifest$file_name %in% mirna_files_downloaded
                manifest <- manifest[tmp, ]

                # only mirna data!!
                if (sum(tmp) > 0) {
                    if (tolower(platform) %in% "illumina hiseq") {
                        tmp <- manifest$platform == "Illumina HiSeq"
                        mirna_files <- manifest[tmp, ]
                    } else if (tolower(platform) %in% "illumina ga") {
                        tmp <- manifest$platform == "Illumina GA"
                        mirna_files <- manifest[tmp, ]
                    } else if (tolower(platform) %in% "h-mirna_8x15kv2") {
                        tmp <- manifest$platform == "H-miRNA_8x15Kv2"
                        mirna_files <- manifest[tmp, ]
                    } else if (tolower(platform) %in% "h-mirna_8x15kv") {
                        tmp <- manifest$platform == "H-miRNA_8x15Kv"
                        mirna_files <- manifest[tmp, ]
                    }
                } else {
                    stop(message("No miRNA data was downloaded!"))
                }
            } else if (tolower(data_base) == "gdc") {
                dir <- file.path(
                    work_dir, "DOAGDC", toupper(tumor),
                    paste0("gdc_", tolower(data_type), "_data")
                )

                manifest <- data.table::fread(file.path(dir, "manifest.sdrf"),
                    select = c("file_name", "cases")
                )
                tmp <- !grepl(
                    pattern = paste(".sdrf", "Data_access_time.txt",
                        sep = "|"
                    ),
                    x = dir(path = dir)
                )

                mirna_files_downloaded <- dir(path = dir)[tmp]
                tmp <- manifest$file_name %in% mirna_files_downloaded
                manifest <- manifest[tmp, ]
                mirna_files <- manifest
            }

            # checking for unmatched data
            if (nrow(mirna_files) > 0) {
                message("Reading data...")
            } else {
                message(paste0(
                    "Wrong Plataform '", tolower(platform),
                    "' chosen!"
                ))
                stop()
            }

            # separing files
            seletor <- unname(sapply(
                mirna_files$cases,
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
                pattern = regex, seletor
            )
            # not related the control samples
            mirna_tumor_files <- as.character(
                mirna_files[[1]][only_tumor_tissue]
            )
            patients <- as.character(as.matrix(
                mirna_files[only_tumor_tissue, "cases"]
            ))
            assign("patients", patients, envir = get(ev))
            open_mirna(mirna_tumor_files)
            # saving
            # write.table(x = data_mirna_table,
            # file = "./mirnaylation_table/mirna_table_completed.csv",
            #           row.names = TRUE, quote = FALSE)

            tumor_data_mirna_table <- sv[["ev"]]$data_mirna_table
            tumor_cross_mapped_mirna_table <- sv[["ev"]]$cross_mapped_mirna_table


            assign(
                "mirna_tumor_cross_mapped",
                sv[["ev"]]$cross_mapped_mirna_table,
                envir = get(ev)
            )

            assign(
                paste0(
                    "mirna_tumor_",
                    ifelse(normalization, "normalized", "nn")
                ),
                sv[["ev"]]$data_mirna_table,
                envir = get(ev)
            )
            if (save_data) {
                message("\nSaving your data...\n")
                # saving
                write.table(
                    x = sv[["ev"]]$data_mirna_table,
                    file = file.path(
                        dir,
                        paste0(
                            "mirna_tumor_",
                            ifelse(normalization, "normalized", "nn"),
                            ".tsv"
                        )
                    ),
                    row.names = TRUE, quote = FALSE,
                    sep = "\t"
                )
                write.table(
                    x = sv[["ev"]]$cross_mapped_mirna_table,
                    file = file.path(
                        dir,
                        "mirna_tumor_cross_mapped.tsv"
                    ),
                    row.names = TRUE, quote = FALSE,
                    sep = "\t"
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

            assign(
                paste0(
                    "mirna_tumor_",
                    ifelse(normalization, "normalized", "nn")
                ),
                sv[["ev"]]$data_mirna_table,
                envir = get(ev)
            )

            if (!missing(name)) {
                mirna_selector <- rownames(
                    sv[["ev"]]$mirna_tumor_normalized
                )
                mirna_selector_final <- grepl(paste0(
                    "^",
                    paste0(
                        "hsa-",
                        tolower(name)
                    ),
                    "$"
                ),
                mirna_selector,
                perl = TRUE
                )

                mirna_exp_selected <- t(sv[["ev"]]$mirna_tumor_normalized[mirna_selector_final, , drop = FALSE])
                cross_mapped_selected <- t(sv[["ev"]]$mirna_tumor_normalized[mirna_selector_final, , drop = FALSE])

                if (sum(mirna_selector_final) == 0) {
                    stop(message("\nmiRNA not found!!!\n\n"))
                }

                colnames(mirna_exp_selected) <- paste0(
                    "hsa-",
                    tolower(name)
                )
                colnames(cross_mapped_selected) <- paste0(
                    "hsa-",
                    tolower(name)
                )
                name <- gsub("-", "_", paste0("hsa-", tolower(name)))
                assign(paste0(
                    "mirna_tumor_normalized_selected_",
                    toupper(name)
                ), mirna_exp_selected,
                envir = get(ev)
                )
                assign(paste0(
                    "mirna_tumor_cross_mapped_selected_",
                    toupper(name)
                ), cross_mapped_selected,
                envir = get(ev)
                )
            }

            if (!tumor_data) {
                # separing files
                seletor <- unname(sapply(
                    mirna_files$cases,
                    function(w) {
                        paste(unlist(strsplit(w, "-"))[4], collapse = "-")
                    }
                ))

                regex <- ifelse(length(normal_type) > 1,
                    paste0(
                        "(", paste(formatC(normal_type, width = 2, flag = "0"),
                            collapse = "|"
                        ), ")", "[A-Z]"
                    ),
                    paste0(formatC(normal_type, width = 2, flag = "0"), "[A-Z]")
                )

                not_tumor_tissue <- grepl(
                    pattern = regex, seletor
                )

                # not tumor
                mirna_not_tumor_files <- as.character(
                    mirna_files[[1]][not_tumor_tissue]
                )
                if (length(mirna_not_tumor_files) == 0) {
                    stop(message(
                        "There is no 'Normal' data in ", tumor,
                        " tumor folder! Please, rerun after",
                        " set 'tumor_data = TRUE'.",
                        "\n\n"
                    ))
                }
                patients <- as.character(as.matrix(
                    mirna_files[not_tumor_tissue, "cases"]
                ))
                open_mirna(mirna_not_tumor_files)
                nt_data_mirna_table <- sv[["ev"]]$data_mirna_table
                nt_cross_mapped_mirna_table <- sv[["ev"]]$cross_mapped_mirna_table

                assign(
                    paste0(
                        "mirna_nt_",
                        ifelse(normalization, "normalized", "nn")
                    ),
                    nt_data_mirna_table,
                    envir = get(ev)
                )

                if (save_data) {
                    message("\nSaving your data...\n")
                    write.table(
                        x = nt_data_mirna_table,
                        file = file.path(
                            dir,
                            "mirna_nt_",
                            ifelse(normalization, "normalized", "nn"),
                            ".tsv"
                        ),
                        row.names = TRUE, quote = FALSE,
                        sep = "\t"
                    )
                    write.table(
                        x = nt_cross_mapped_mirna_table,
                        file = file.path(
                            dir,
                            "mirna_nt_cross_mapped.tsv"
                        ),
                        row.names = TRUE, quote = FALSE,
                        sep = "\t"
                    )
                }

                if (!missing(name)) {
                    # not tumor
                    mirna_selector <- rownames(sv[["ev"]]$mirna_nt_normalized)
                    mirna_selector_final <- grepl(paste0(
                        "^",
                        paste0(
                            "hsa-",
                            tolower(name)
                        ),
                        "$"
                    ),
                    mirna_selector,
                    perl = TRUE
                    )

                    mirna_exp_selected <- t(sv[["ev"]]$mirna_nt_normalized[mirna_selector_final, , drop = FALSE])
                    cross_mapped_selected <- t(sv[["ev"]]$mirna_nt_normalized[mirna_selector_final, , drop = FALSE])

                    colnames(mirna_exp_selected) <- paste0(
                        "hsa-",
                        tolower(name)
                    )
                    colnames(cross_mapped_selected) <- paste0(
                        "hsa-",
                        tolower(name)
                    )
                    name <- gsub("-", "_", paste0("hsa-", tolower(name)))
                    assign(paste0(
                        "mirna_nt_normalized_selected_",
                        toupper(name)
                    ), mirna_exp_selected,
                    envir = get(ev)
                    )
                    assign(paste0(
                        "mirna_nt_cross_mapped_selected_",
                        toupper(name)
                    ), cross_mapped_selected,
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
                mirna_selector <- rownames(
                    sv[["ev"]]$mirna_tumor_normalized
                )
                mirna_selector_final <- grepl(paste0(
                    "^",
                    paste0(
                        "hsa-",
                        tolower(name)
                    ),
                    "$"
                ),
                mirna_selector,
                perl = TRUE
                )

                mirna_exp_selected <- t(sv[["ev"]]$mirna_tumor_normalized[mirna_selector_final, , drop = FALSE])
                cross_mapped_selected <- t(sv[["ev"]]$mirna_tumor_normalized[mirna_selector_final, , drop = FALSE])

                if (sum(mirna_selector_final) == 0) {
                    stop(message("\nmiRNA not found!!!\n\n"))
                }

                colnames(mirna_exp_selected) <- paste0(
                    "hsa-",
                    tolower(name)
                )
                colnames(cross_mapped_selected) <- paste0(
                    "hsa-",
                    tolower(name)
                )
                name <- gsub("-", "_", paste0(
                    "hsa-",
                    tolower(name)
                ))
                assign(paste0(
                    "mirna_tumor_normalized_selected_",
                    toupper(name)
                ), mirna_exp_selected,
                envir = get(ev)
                )
                assign(paste0(
                    "mirna_tumor_cross_mapped_selected_",
                    toupper(name)
                ), cross_mapped_selected,
                envir = get(ev)
                )
                if (!tumor_data) {
                    # not tumor
                    mirna_selector <- rownames(
                        sv[["ev"]]$mirna_nt_normalized
                    )
                    mirna_selector_final <- grepl(paste0(
                        "^",
                        paste0(
                            "hsa-",
                            tolower(name)
                        ),
                        "$"
                    ),
                    mirna_selector,
                    perl = TRUE
                    )

                    mirna_exp_selected <- t(sv[["ev"]]$mirna_nt_normalized[mirna_selector_final, , drop = FALSE])
                    cross_mapped_selected <- t(sv[["ev"]]$mirna_nt_normalized[mirna_selector_final, , drop = FALSE])

                    colnames(mirna_exp_selected) <- paste0(
                        "hsa-",
                        tolower(name)
                    )
                    colnames(cross_mapped_selected) <- paste0(
                        "hsa-",
                        tolower(name)
                    )
                    name <- gsub("-", "_", paste0("hsa-", tolower(name)))
                    assign(paste0(
                        "mirna_nt_normalized_selected_",
                        toupper(name)
                    ), mirna_exp_selected,
                    envir = get(ev)
                    )
                    assign(paste0(
                        "mirna_nt_cross_mapped_selected_",
                        toupper(name)
                    ), cross_mapped_selected,
                    envir = get(ev)
                    )
                }
            }
        }
        suppressWarnings(remove(data_mirna_table,
            cross_mapped_mirna_table,
            envir = get(ev)
        ))
    }
}