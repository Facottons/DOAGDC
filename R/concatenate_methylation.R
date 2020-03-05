#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_methylation} is a function designed to concatenate GDC
#' files into a single matrix, where the columns stand for patients code and
#' rows stand for data names.
#'
#' @param name A character string indicating the desired values to be used in
#'    next analysis. For instance, "HIF3A" in the legacy gene expression matrix,
#'    "mir-1307" in the miRNA quantification matrix, or "HER2" in the protein
#'    quantification matrix.
#' @param data_base
#' @param work_dir
#' @param tumor
#' @param tumor_data Logical value where \code{TRUE} specifies the desire to
#'    work with tumor tissue files only. When set to FALSE, it creates two
#'    matrices, one containing tumor data and other containing data from
#'    not-tumor tissue. The default is \code{TRUE}.
#' @param cutoff_beta_na Numerical value indicating the maximum threshold
#'    percentage (in decimal form) to tolerate and to remove rows containing NA
#'    for beta values (methylation data). The default is 0.25.
#' @param cutoff_betasd Numerical value indicating the standard deviation
#'    threshold of beta values (methylation data). It keeps only rows that have
#'    standard deviation of beta values higher than the threshold. The default
#'    is \code{0.005}.
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
#' @param env A character string containing the environment name that should be
#'    used. If none has been set yet, the function will create one in global
#'    environment following the standard criteria:
#'    \itemize{\item{'<tumor>_<data_base>_methylation_tumor_data'}{ or}
#'    \item{'<tumor>_<data_base>_methylation_both_data'}{ (for tumor and
#'    not tumor data in separated matrices).}}
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
#' concatenate_methylation(
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     platform = "Illumina Human Methylation 450",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
concatenate_methylation <- function(name, data_base,
                                    work_dir, tumor,
                                    tumor_data = TRUE,
                                    cutoff_beta_na = 0.25,
                                    cutoff_betasd = 0.005,
                                    only_filter = FALSE,
                                    tumor_type = 1,
                                    normal_type = 11,
                                    platform = "",
                                    env,
                                    save_data = FALSE) {

    # local functions ####
    open_meth <- function(meth_files_arg) {
        # reading files in one file
        pb <- txtProgressBar(min = 0, max = length(meth_files_arg), style = 3)
        count <- 0
        # try to open with lapply(list, function)
        for (file in file.path(dir, meth_files_arg)) {
            count <- count + 1
            if (count > 1) {
                actual_file <- data.table::fread(file,
                    sep = "\t",
                    showProgress = FALSE,
                    select = 2
                )
                meth_table[, count] <- suppressWarnings(
                    as.numeric(actual_file[[1]])
                )
                colnames(meth_table)[count] <- names(actual_file)
            } else if (count == 1) {
                actual_file <- data.table::fread(file,
                    sep = "\t",
                    showProgress = FALSE
                )
                refererence_meth_table <- as.data.frame(
                    actual_file[, select_ref, with = FALSE]
                )
                rownames(refererence_meth_table) <- as.character(
                    as.matrix(actual_file[, 1])
                )
                meth_table <- matrix(
                    nrow = nrow(actual_file),
                    dimnames = list(rownames(refererence_meth_table), c()),
                    ncol = length(meth_files_arg)
                )
                meth_table[, 1] <- suppressWarnings(
                    as.numeric(actual_file[[2]])
                )
                colnames(meth_table) <- seq_len(length(meth_files_arg))
                colnames(meth_table)[1] <- names(actual_file)[2]
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)
        if (tolower(data_base) == "legacy") {
            colnames(refererence_meth_table) <- refererence_meth_table[1, ]
            refererence_meth_table <- refererence_meth_table[-1, , drop = FALSE]
            meth_table <- meth_table[-1, , drop = FALSE]
        }
        return(cbind(meth_table, refererence_meth_table))
    }

    filter_meth_fun <- function(meth_table, files) {

        ### preparing data          FASTER WITH MATRIX
        message("Please wait... Filtering data!")
        pb <- txtProgressBar(min = 0, max = 3, style = 3)
        count <- 0
        setTxtProgressBar(pb, count)
        initial_row <- nrow(meth_table)
        # The beta value (β) estimates the methylation level of the CpG locus.
        # beta value (β) = methylated/unmethylated alleles
        # Remove sites containing NA for beta values based in cutoff
        beta_na_filter <- rowMeans(is.na(
            meth_table[, seq_len(length(files))]
        )) < cutoff_beta_na
        meth_table <- meth_table[beta_na_filter, ]

        count <- 1
        setTxtProgressBar(pb, count)

        # Keep sites for which the beta values have standard deviation value
        # higher than cutoff
        beta_na_filter <- apply(
            meth_table[, seq_len(length(files))],
            1, sd,
            na.rm = TRUE
        ) > cutoff_betasd
        meth_table <- meth_table[beta_na_filter, ]

        count <- 2
        setTxtProgressBar(pb, count)

        # Remove CpGs on sex chromosomes
        tmp1 <- meth_table$Chromosome != "X"
        tmp2 <- meth_table$Chromosome != "Y"
        sex_filter <- tmp1 & tmp2
        sex_filter[is.na(sex_filter)] <- as.logical("FALSE")
        meth_table <- meth_table[sex_filter, ]

        # only with gene symbol info (can be changed after released)
        tmp1 <- meth_table$Gene_Symbol != ""
        tmp2 <- meth_table$Gene_Symbol != "."
        gene_symbol_filter <- tmp1 & tmp2
        meth_table <- meth_table[gene_symbol_filter, ]

        count <- 3
        setTxtProgressBar(pb, count)

        # Number of removed lines after filtering
        message(paste("It was filtered", initial_row - nrow(meth_table),
            "lines from data!",
            sep = " "
        ))
        close(pb)
        # normalization (3 options) Batch effect: systematic differences across
        # groups of samp The background was subtracted using the methylumi
        # package (method "noob") [1]. The signal intensity values were
        # normalized using the SWAN normalization method, as implemented in the
        # minfi package. or lumi R package or minfi Bioconductor packag(SWAN)
        # library(methylumi) ?? wateRmellow??
        if (nrow(meth_table) == 0) {
            stop(message(
                "All data were filtered out! Please, choose a different",
                " cutoff_beta_na and cutoff_betasd parameters"
            ))
        }

        return(meth_table)
    }

    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    if (!only_filter) {
        dir <- ifelse(tolower(data_base) == "legacy", file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "methylation_data"
        ), file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "gdc_methylation_data"
        ))

        manifest <- data.table::fread(file.path(dir, "manifest.sdrf"),
            select = c(
                "file_name",
                "cases",
                "platform"
            )
        )

        manifest <- as.data.frame(as.matrix(manifest),
            stringsAsFactors = FALSE
        )

        tmp <- !grepl(
            pattern = paste(".sdrf", "Data_access_time.txt",
                sep = "|"
            ),
            x = dir(path = dir)
        )

        meth_files_downloaded <- dir(path = dir)[tmp]
        tmp <- manifest$file_name %in% meth_files_downloaded
        manifest <- manifest[tmp, ]

        dir.create(file.path(
            work_dir, "DOAGDC",
            toupper(tumor), "Analyses"
        ),
        showWarnings = FALSE
        )

        dir.create(file.path(
            work_dir, "DOAGDC",
            toupper(tumor), "Analyses", name
        ),
        showWarnings = FALSE
        )

        dir.create(file.path(
            work_dir, "DOAGDC",
            toupper(tumor), "Analyses", name,
            "Methylation"
        ),
        showWarnings = FALSE
        )

        # only meth data!!
        if (sum(tmp) > 0) {
            tmp <- tolower(manifest$platform) == tolower(platform)
            meth_files <- manifest[tmp, ]

            # checking for unmatched data
            if (nrow(meth_files) > 0) {
                message("Reading data...")
            } else {
                stop(message(
                    "There is no data from ", tolower(platform),
                    " in your local storage!\n"
                ))
            }
        } else {
            stop(message("No Methylation data was downloaded!\n"))
        }

        # separing files
        if (tolower(data_base) == "legacy") {
            tmp <- meth_files[, 1]
            select_ref <- c(3:5)
        } else {
            tmp <- meth_files[, 2]
            select_ref <- c(3:11)
        }

        seletor <- unname(sapply(
            tmp,
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

        only_tumor_tissue <- grepl(pattern = regex, seletor, perl = TRUE)
        meth_tumor_files <- as.character(meth_files[, 1][only_tumor_tissue])
        meth_table <- open_meth(meth_tumor_files)
        tmp <- (length(meth_tumor_files) + 1):(ncol(meth_table))

        # filtering
        meth_table <- filter_meth_fun(meth_table, meth_tumor_files)

        reference_table_filtered <- meth_table[, tmp]
        tmp <- seq_len(length(meth_tumor_files))
        tumor_meth_table_filtered <- meth_table[, tmp]

        if (tolower(data_base) == "gdc") {
            reference_table_filtered[, 4] <- unlist(lapply(strsplit(
                x = reference_table_filtered[, 4],
                split = ";", perl = TRUE
            ), "[", 1))
        }

        if (save_data) {
            message("\nSaving data, this could take a while...\n")
            # saving
            tmp <- "methylation_tumor_filtered.csv"
            write.csv(
                x = cbind(
                    reference_table_filtered,
                    tumor_meth_table_filtered
                ),
                file = file.path(
                    work_dir, "DOAGDC",
                    toupper(tumor), "Analyses", toupper(name),
                    "Methylation", tmp
                ),
                row.names = TRUE, quote = FALSE
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

        assign("methylation_tumor_filtered", tumor_meth_table_filtered,
            envir = get(ev)
        )
        assign("reference_table_filtered", reference_table_filtered,
            envir = get(ev)
        )

        if (!missing(name)) {
            select_gene <- grepl(
                pattern = name,
                x = reference_table_filtered$Gene_Symbol
            )
            if (sum(select_gene) == 0) {
                message(
                    name, " not found...\n\n",
                    "Please insert a valid gene symbol!!"
                )
            }
            assign(paste0("reference_table_filtered_selected_", name),
                reference_table_filtered[select_gene, ],
                envir = get(ev)
            )

            message("\nSaving required data...\n")
            # saving
            tmp <- "methylation_tumor_filtered_selected_"
            write.csv(
                x = cbind(
                    reference_table_filtered[select_gene, ],
                    tumor_meth_table_filtered[select_gene, ]
                ),
                file = file.path(
                    work_dir, "DOAGDC",
                    toupper(tumor), "Analyses", toupper(name),
                    "Methylation", paste0(tmp, name, ".csv")
                ),
                row.names = TRUE, quote = FALSE
            )

            assign(paste0("methylation_tumor_filtered_selected_", name),
                tumor_meth_table_filtered[select_gene, ],
                envir = get(ev)
            )
            assign("name", name, envir = get(ev))
        }

        if (!tumor_data) {
            # separing files

            seletor <- unname(sapply(
                meth_files,
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

            not_tumor_tissue <- grepl(pattern = regex, seletor)

            if (sum(not_tumor_tissue) == 0) {
                stop(message(
                    "\nThere is no 'normal' data available to ",
                    tumor, " tumor!!\n\n"
                ))
            }
            meth_nt_files <- as.character(meth_files[, 1][not_tumor_tissue])
            meth_table <- open_meth(meth_nt_files)
            tmp <- (length(meth_nt_files) + 1):(ncol(meth_table))

            ## filtering
            meth_table <- filter_meth_fun(meth_table, meth_nt_files)

            reference_table_filtered <- meth_table[, tmp]
            tmp <- seq_len(length(meth_nt_files))
            nt_meth_table_filtered <- meth_table[, tmp]

            assign("methylation_nt_filtered",
                nt_meth_table_filtered,
                envir = get(ev)
            )
            assign("reference_nt_table_filtered", reference_table_filtered,
                envir = get(ev)
            )

            # saving
            if (save_data) {
                message("\nSaving data...\n")
                tmp <- "not_tumor_meth_table_completed_filtered.tsv"
                write.table(
                    x = cbind(
                        reference_table_filtered,
                        nt_meth_table_filtered
                    ),
                    file = file.path(dir, tmp),
                    row.names = TRUE, quote = FALSE, sep = "\t"
                )
            }

            if (!missing(name)) {
                assign(paste0(
                    "methylation_nt_filtered_selected_",
                    name
                ),
                nt_meth_table_filtered[select_gene, ],
                envir = get(ev)
                )

                message("\nSaving required data...\n")
                # saving
                tmp <- "methylation_nt_filtered_selected_"
                write.table(
                    x = cbind(
                        reference_table_filtered[select_gene, ],
                        nt_meth_table_filtered[select_gene, ]
                    ),
                    file = file.path(
                        work_dir, "DOAGDC",
                        toupper(tumor), "Analyses", toupper(name),
                        "Methylation", paste0(tmp, name, ".csv")
                    ),
                    row.names = TRUE, quote = FALSE
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

        if (!missing(name)) {
            select_gene <- grepl(
                pattern = name,
                x = sv[["ev"]]$reference_table_filtered$Gene_Symbol
            )
            if (sum(select_gene) == 0) {
                stop(message(
                    name, " not found...\n\n",
                    "Please insert a valid gene symbol!!"
                ))
            }
            assign(paste0("reference_table_filtered_selected_", name),
                sv[["ev"]]$reference_table_filtered[select_gene, ],
                envir = get(ev)
            )

            assign(paste0("methylation_tumor_filtered_selected_", name),
                sv[["ev"]]$methylation_tumor_filtered[select_gene, ],
                envir = get(ev)
            )

            message("\nSaving required data...\n")
            # saving
            write.table(
                x = cbind(
                    sv[["ev"]]$reference_table_filtered[select_gene, ],
                    sv[["ev"]]$methylation_tumor_filtered[select_gene, ]
                ),
                file = file.path(
                    work_dir, "DOAGDC", toupper(tumor),
                    "Analyses", name,
                    "Methylation",
                    paste0(
                        "methylation_tumor",
                        "_filtered_selected_",
                        name,
                        ".tsv"
                    )
                ),
                row.names = TRUE, quote = FALSE, sep = "\t"
            )

            # not tumor
            if (!tumor_data) {
                assign(paste0(
                    "methylation_nt_filtered_selected_",
                    name
                ),
                sv[["ev"]]$methylation_nt_filtered[select_gene, ],
                envir = get(ev)
                )

                message("\nSaving required data...\n")
                # saving
                tmp <- "methylation_nt_filtered_selected_"
                write.table(
                    x = cbind(
                        sv[["ev"]]$reference_nt_table_filtered[select_gene, ],
                        sv[["ev"]]$methylation_nt_filtered[select_gene, ]
                    ),
                    file = file.path(
                        dir,
                        paste0(tmp, name, ".tsv")
                    ),
                    row.names = TRUE, quote = FALSE, sep = "\t"
                )
            }
        }
    }
}