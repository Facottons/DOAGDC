#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_exon} is a function designed to concatenate GDC files into
#' a single matrix, where the columns stand for patients code and rows stand for
#' data names.
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
#' @param htseq A character string indicating which htseq workflow data should
#'    be downloaded (only applied to "GDC" gene expression): "Counts", "FPKM"
#'    or "FPKM-UQ".
#' @param work_dir
#' @param tumor
#' @param workflow_type A character string specifying the workflow type for
#'    mutation data in "gdc". Where: \itemize{\item{"varscan" stands for}{
#'    VarScan2 Variant Aggregation and Masking} \item{"mutect" stands for}{
#'    MuTect2 Variant Aggregation and Masking}\item{"muse" stands for}{ MuSE
#'    Variant Aggregation and Masking}\item{"somaticsniper" stands for}{
#'    SomaticSniper Variant Aggregation and Masking}\item{"all" means to}{
#'    concatenate all workflows into a single matrix.}}
#' @param tumor_data Logical value where \code{TRUE} specifies the desire to
#'    work
#'    with tumor tissue files only. When set to FALSE, it creates two matrices,
#'    one containing tumor data and other containing data from not-tumor tissue.
#'    The default is \code{TRUE}.
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
#' @importFrom R.utils gunzip
#' @importFrom reshape2 melt
#' @importFrom XML xmlToDataFrame
#'
#' @examples
#' library(DOAGDC)
#'
#' # Concatenating gene expression data into a single matrix
#' # data already downloaded using the 'download_gdc' function
#' concatenate_exon("gene",
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
concatenate_exon <- function(data_type,
                            normalization = TRUE,
                            name, data_base,
                            htseq = NULL,
                            work_dir, tumor,
                            workflow_type,
                            tumor_data = TRUE,
                            only_filter = FALSE,
                            tumor_type = 1,
                            normal_type = 11,
                            platform = "",
                            env,
                            save_data = FALSE) {

    # create env if not created yet
    if (missing(env) && tumor_data && !only_filter) {
        assign(paste(toupper(tumor), toupper(data_base),
            "exon", "tumor_data", sep = "_"
        ),
        new.env(parent = emptyenv()),
        envir = .GlobalEnv
        )
        envir_link <- paste(toupper(tumor), toupper(data_base),
            "exon", "tumor_data", sep = "_"
        )
    } else if (missing(env) && !tumor_data && !only_filter) {
        assign(paste(toupper(tumor), toupper(data_base),
            "exon", "both_data", sep = "_"
        ),
        new.env(parent = emptyenv()),
        envir = .GlobalEnv
        )
        envir_link <- paste(toupper(tumor), toupper(data_base),
            "exon", "both_data", sep = "_"
        )
    } else if (missing(env) && only_filter) {
        message(
            "Please, before using 'only_filter'",
            " argument, insert the Environment name."
        )
    } else {
        envir_link <- deparse(substitute(env))
    }

    string_vars <- list(envir_link = get(envir_link))
    attr(envir_link, "name") <- paste0(
        "Environment created by DOAGDC",
        " package, use its name in ",
        "'env' argument"
    )

    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    if (!only_filter) {
        if (tolower(data_base) == "gdc") {
            stop(message(
                "\nThrere is no Exon quantification",
                " data in GDC data base!!",
                "\nPlease use 'Legacy' data base\n"
            ))
        }
        dir <- file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "exon quantification_data"
        )


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

        if (tolower(platform) == "") {
            message(
                "Please, insert 'illumina hiseq' or 'illumina ga' ",
                "for Exon quantification data."
            )
            stop()
        } else if (tolower(platform) == "illumina ga") {
            to_gunzip <- file.path(
                desing_array_dir,
                paste0(
                    "unc.edu_", toupper(tumor),
                    ".IlluminaHiSeq_RNASeqV2.",
                    "mage-tab.1.1.0.tar.gz"
                )
            )
            ## or, if you just want to extract the target file:
            untar(to_gunzip, exdir = desing_array_dir)

            design_array <- data.table::fread(file.path(
                desing_array_dir,
                paste0(
                    "unc.edu_",
                    toupper(tumor),
                    ".IlluminaHiSeq_",
                    "RNASeqV2.mage-",
                    "tab.1.1.0"
                ),
                paste0(
                    "unc.edu_",
                    toupper(tumor),
                    ".IlluminaHiSeq_",
                    "RNASeqV2.1.1.",
                    "0.sdrf.txt"
                )
            ),
            select = c(2, 22)
            )
        } else if (tolower(platform) == "illumina hiseq") {
            to_gunzip <- file.path(
                desing_array_dir,
                paste0(
                    "unc.edu_", toupper(tumor),
                    ".IlluminaHiSeq_RNASeqV2",
                    ".mage-tab.1.1.0.tar.gz"
                )
            )
            ## or, if you just want to extract the target file:
            untar(to_gunzip, exdir = desing_array_dir)

            design_array <- data.table::fread(file.path(
                desing_array_dir,
                paste0(
                    "unc.edu_",
                    toupper(tumor),
                    ".IlluminaHiSeq_",
                    "RNASeqV2.mage-",
                    "tab.1.1.0"
                ),
                paste0(
                    "unc.edu_",
                    toupper(tumor),
                    ".IlluminaHiSeq_",
                    "RNASeqV2.1.",
                    "1.0.sdrf.txt"
                )
            ),
            select = c(2, 22)
            )
        }

        exon_files <- dir(dir, pattern = "unc.edu")
        tmp <- unlist(design_array[design_array[[2]] %in% exon_files, 1])
        patient_code <- as.character(tmp)

        if (normalization) {
            select_this_column <- 4 # RPKM
        } else {
            select_this_column <- 2 # raw counts
        }

        pb <- txtProgressBar(min = 0, max = length(exon_files), style = 3)
        count <- 0
        for (i in file.path(dir, exon_files)) {
            count <- count + 1
            if (count == 1) {
                exon <- data.table::fread(i)
                exon_table <- matrix(
                    ncol = length(exon_files),
                    nrow = nrow(exon)
                )
                colnames(exon_table) <- seq(length(exon_files))
                rownames(exon_table) <- as.character(exon[[1]])
                exon_table[, 1] <- as.numeric(exon[[select_this_column]])
            } else {
                exon <- data.table::fread(i, select = select_this_column)
                exon_table[, count] <- as.numeric(exon[[1]])
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)

        colnames(exon_table) <- patient_code

        # separing files
        seletor <- unname(sapply(
            patient_code,
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
        patients_tumor <- patient_code[only_tumor_tissue]
        assign("patients", patients_tumor, envir = get(envir_link))
        assign("exon_tumor_", exon_table, envir = get(envir_link))

        if (save_data) {
            message("\nSaving your data...\n")
            # saving
            write.table(
                x = exon_table,
                file = file.path(dir, "tumor_exon_table.tsv"),
                row.names = TRUE, quote = FALSE, sep = "\t"
            )
        }
        if (!tumor_data) {
            # separing files
            seletor <- unname(sapply(
                patient_code,
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
                    " set 'tumor_data = TRUE'.\n\n"
                ))
            }

            not_tumor_patients <- as.character(
                colnames(exon_table)[not_tumor_tissue]
            )
            exon_not_tumor <- exon_table[, not_tumor_patients]
            assign("not_tumor_patients", not_tumor_patients,
                envir = get(envir_link)
            )
            assign("exon_not_tumor_", exon_not_tumor,
                envir = get(envir_link)
            )

            if (save_data) {
                message("\nSaving your data...\n")
                # saving
                write.table(
                    x = tumor_protein_table,
                    file = file.path(
                        dir,
                        "tumor_protein_table.tsv"
                    ),
                    row.names = TRUE, quote = FALSE, sep = "\t"
                )
                write.table(
                    x = exon_not_tumor,
                    file = file.path(dir, "exon_not_tumor.tsv"),
                    row.names = TRUE, quote = FALSE, sep = "\t"
                )
            }
        }
    } else {
        if (normalization) {
            protein_selector <- rownames(
                string_vars[["envir_link"]]$tumor_protein_table
            )
        } else {
            stop(message(
                "Only normalized data can ",
                "be used to separate patients!"
            ))
        }
        protein_selector_final <- grepl(name,
            protein_selector,
            perl = TRUE
        )

        protein_exp_selected <- t(string_vars[["envir_link"]]$tumor_protein_table[protein_selector_final, ])
        colnames(protein_exp_selected) <- name
        colnames(cross_mapped_selected) <- name
        assign(paste0("exon_exp_selected_", name), protein_exp_selected,
            envir = get(envir_link)
        )
    }

    # NOTE Mage
    if (FALSE) {
        if (tolower(data_base) == "legacy") {
            dir <- file.path(work_dir, "DOAGDC", toupper(tumor), "mage_data")
        } else if (tolower(data_base) == "gdc") {
            dir <- file.path(
                work_dir, "DOAGDC", toupper(tumor),
                "gdc_protein_data"
            )
        }
        dir(path = dir, pattern = "mage-tab")
    }

    # end ####
    assign("tumor", tumor, envir = get(envir_link))
    assign("data_base", data_base, envir = get(envir_link))
    assign("data_type", data_type, envir = get(envir_link))
    assign("work_dir", work_dir, envir = get(envir_link))
    if (!missing(name)) {
        assign("name", name, envir = get(envir_link))
        assign("name_e", name, envir = get(envir_link))
    }
    message("\nDone!")
}
