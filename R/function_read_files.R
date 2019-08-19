#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_files} is a function designed to concatenate GDC files into
#' a single matrix, where the columns stand for patients code and rows stand for
#' data names.
#'
#' @param dataType
#' @param normalization Logical value where \code{TRUE} specify the desire to
#'    work with normalized files only. When FALSE, in the second run, do not
#'    forget to set env argument. This argument is only applyable to gene and
#'    isoform expression data from GDC Legacy Archive. The default is
#'    \code{TRUE}.
#' @param Name A character string indicating the desired values to be used in
#'    next analysis. For instance, "HIF3A" in the legacy gene expression matrix,
#'    "mir-1307" in the miRNA quantification matrix, or "HER2" in the protein
#'    quantification matrix.
#' @param dataBase
#' @param HTSeq A character string indicating which HTSeq workflow data should
#'    be downloaded (only applied to "GDC" gene expression): "Counts", "FPKM",
#'    or "FPKM-UQ".
#' @param workDir
#' @param tumor
#' @param workflowType A character string specifying the workflow type for
#'    mutation data in "gdc". Where: \itemize{\item{"varscan" stands for}{
#'    VarScan2 Variant Aggregation and Masking} \item{"mutect" stands for}{
#'    MuTect2 Variant Aggregation and Masking}\item{"muse" stands for}{ MuSE
#'    Variant Aggregation and Masking}\item{"somaticsniper" stands for}{
#'    SomaticSniper Variant Aggregation and Masking}\item{"all" means to}{
#'    concatenate all workflows into a single matrix.}}
#' @param tumorData Logical value where \code{TRUE} specifies the desire to work
#'    with tumor tissue files only. When set to FALSE, it creates two matrices,
#'    one containing tumor data and other containing data from not-tumor tissue.
#'    The default is \code{TRUE}.
#' @param cutoffBetaNA Numerical value indicating the maximum threshold
#'    percentage (in decimal form) to tolerate and to remove rows containing NA
#'    for beta values (methylation data). The default is 0.25.
#' @param cutoffBetasd Numerical value indicating the standard deviation
#'    threshold of beta values (methylation data). It keeps only rows that have
#'    standard deviation of beta values higher than the threshold. The default
#'    is \code{0.005}.
#' @param onlyFilter Logical value where \code{TRUE} indicates that the matrix
#'    is already concatenate and the function should choose a different
#'    \code{Name}, without concatenate all the files again. The default is
#'    FALSE.
#' @param Platform
#' @param use_hg19_mirbase20 Logical value where \code{TRUE} indicates that only
#'    hg19.mirbase20 should be used. This parameter is needed when using
#'    \code{dataBase = "legacy"} and one of the available miRNA \code{dataType}
#'    in "legacy" ("miRNA gene quantification" and "miRNA isoform
#'    quantification"). The default is FALSE.
#' @param env A character string containing the environment name that should be
#'    used. If none has been set yet, the function will create one in global
#'    environment following the standard criteria:
#'    \itemize{\item{'tumor_dataBase_dataType_tumor_data'}{ or}
#'    \item{'tumor_dataBase_dataType_both_data'}{ (for tumor and not tumor data
#'    in separated matrices).}}
#' @param saveData Logical value where \code{TRUE} indicates that the
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
#' # Concatenating gene expression data into a single matrix
#' # data already downloaded using the 'download_gdc' function
#' concatenate_files("gene",
#'                    Name = "HIF3A",
#'                    dataBase = "legacy",
#'                    tumor = "CHOL",
#'                    workDir = "~/Desktop")
#'
concatenate_files <- function(dataType,
                    normalization = TRUE,
                    Name, dataBase,
                    HTSeq = NULL,
                    workDir, tumor,
                    workflowType,
                    tumorData = TRUE,
                    cutoffBetaNA = 0.25,
                    cutoffBetasd = 0.005,
                    onlyFilter = FALSE,
                    Platform = "",
                    use_hg19_mirbase20 = FALSE,
                    env,
                    saveData = FALSE) {

    #  c("data.table", "R.utils")

    #local functions ####
    SW <- function(x) {suppressWarnings(x)}

    gene_isoform <- function() {
        codes <- dir(path = DIR, pattern = ".sdrf$")
        codigos <- read.table(file.path(DIR, codes[1]),
                            stringsAsFactors = FALSE,
                            header = TRUE, sep = "\t")
        codigos <- codigos[, c("file_name", "submitter_id")]
        tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                            collapse = "[aAbBcCdD]-|-")
        not_tumor_regex <- paste(10:19, collapse = "[aAbBcCdD]-|-")

        if (tumorData) {
            codigos <- codigos[grepl(paste0("-", tumor_regex, "-"),
                                    x = codigos$submitter_id), ]
            remove_duplicated <- duplicated(codigos$file_name)
            codigos <- codigos[!remove_duplicated, ]
            rownames(codigos) <- codigos$file_name

            #selecting files by type (normalized or not)
            if (normalization) {
                genes_normalized <- grepl("^.+(genes.normalized).+$",
                                        codigos$file_name)
                isoforms_normalized <- grepl("^.+(isoforms.normalized).+$",
                                            codigos$file_name)
                gene_files <- codigos[genes_normalized, ]
                isoform_files <- codigos[isoforms_normalized, ]
                assign("genes_normalized", genes_normalized,
                    envir = get(envir_link))
                assign("isoforms_normalized", isoforms_normalized,
                    envir = get(envir_link))
                assign("gene_files", gene_files, envir = get(envir_link))
                assign("isoform_files", isoform_files, envir = get(envir_link))
            } else {
                genes_not_normalized <- grepl("^.+(.rsem.genes.results)$",
                                            codigos$file_name)
                isoforms_not_normalized <- grepl("^.+(isoforms.results)$",
                                                codigos$file_name)
                gene_files <- codigos[genes_not_normalized, ]
                isoform_files <- codigos[isoforms_not_normalized, ]
                assign("genes_not_normalized", genes_not_normalized,
                    envir = get(envir_link))
                assign("isoforms_not_normalized", isoforms_not_normalized,
                    envir = get(envir_link))
                assign("gene_files", gene_files, envir = get(envir_link))
                assign("isoform_files", isoform_files, envir = get(envir_link))
            }
            assign("codigos", codigos, envir = get(envir_link))
        } else {
            codigos_tumor <- codigos[grepl(paste0("-", tumor_regex, "-"),
                                    x = codigos$submitter_id), ]
            codigos_not_tumor <- codigos[grepl(paste0("-", not_tumor_regex,
                                                    "[aAbBcCdD]-"),
                                    x = codigos$submitter_id), ]
            remove_duplicated_tumor <- duplicated(codigos_tumor$file_name)
            remove_duplicated_not_tumor <- duplicated(codigos_not_tumor$file_name)
            codigos_tumor <- codigos_tumor[!remove_duplicated_tumor, ]
            codigos_not_tumor <- codigos_not_tumor[!remove_duplicated_not_tumor, ]

            rownames(codigos_tumor) <- codigos_tumor$file_name
            rownames(codigos_not_tumor) <- codigos_not_tumor$file_name

            #selecting files by type (normalized or not)
            if (normalization) {
                #tumor
                genes_normalized <- grepl("^.+(genes.normalized).+$",
                                        codigos_tumor$file_name)
                isoforms_normalized <- grepl("^.+(isoforms.normalized).+$",
                                            codigos_tumor$file_name)
                gene_files <- codigos_tumor[genes_normalized, ]
                isoform_files <- codigos_tumor[isoforms_normalized, ]
                assign("genes_normalized", genes_normalized,
                    envir = get(envir_link))
                assign("isoforms_normalized", isoforms_normalized,
                    envir = get(envir_link))
                assign("gene_files", gene_files, envir = get(envir_link))
                assign("isoform_files", isoform_files, envir = get(envir_link))

                #not tumor
                genes_normalized_not_tumor <- grepl("^.+(genes.normalized).+$",
                                                    codigos_not_tumor$file_name)
                isoforms_normalized_not_tumor <- grepl("^.+(isoforms.normalized).+$", codigos_not_tumor$file_name)
                tmp <- codigos_not_tumor[genes_normalized_not_tumor, ]
                gene_files_not_tumor <- tmp
                tmp <- codigos_not_tumor[isoforms_normalized_not_tumor, ]
                isoform_files_not_tumor <- tmp
                assign("genes_normalized_not_tumor", genes_normalized_not_tumor,
                    envir = get(envir_link))
                assign("isoforms_normalized_not_tumor",
                    isoforms_normalized_not_tumor, envir = get(envir_link))
                assign("gene_files_not_tumor", gene_files_not_tumor,
                    envir = get(envir_link))
                assign("isoform_files_not_tumor", isoform_files_not_tumor,
                    envir = get(envir_link))
            } else {
                #tumor
                genes_not_normalized <- grepl("^.+(.rsem.genes.results)$",
                                            codigos_tumor$file_name)
                isoforms_not_normalized <- grepl("^.+(isoforms.results)$",
                                                codigos_tumor$file_name)
                gene_files <- codigos_tumor[genes_not_normalized, ]
                isoform_files <- codigos_tumor[isoforms_not_normalized, ]
                assign("genes_not_normalized", genes_not_normalized,
                    envir = get(envir_link))
                assign("isoforms_not_normalized", isoforms_not_normalized,
                    envir = get(envir_link))
                assign("gene_files", gene_files, envir = get(envir_link))
                assign("isoform_files", isoform_files, envir = get(envir_link))

                #not tumor
                genes_not_normalized_not_tumor <- grepl("^.+(.rsem.genes.results)$", codigos_not_tumor$file_name)
                isoforms_not_normalized_not_tumor <- grepl("^.+(isoforms.results)$", codigos_not_tumor$file_name)
                gene_files_not_tumor <- codigos_not_tumor[genes_not_normalized_not_tumor, ]
                isoform_files_not_tumor <- codigos_not_tumor[isoforms_not_normalized_not_tumor, ]
                assign("genes_not_normalized_not_tumor",
                        genes_not_normalized_not_tumor,
                        envir = get(envir_link))
                assign("isoforms_not_normalized_not_tumor",
                        isoforms_not_normalized_not_tumor,
                        envir = get(envir_link))
                assign("gene_files_not_tumor", gene_files_not_tumor,
                    envir = get(envir_link))
                assign("isoform_files_not_tumor", isoform_files_not_tumor,
                    envir = get(envir_link))
            }
            assign("codigos_tumor", codigos_tumor, envir = get(envir_link))
            assign("codigos_not_tumor", codigos_not_tumor,
                envir = get(envir_link))
        }
    }

    OPEN_genexpress <- function(files, Var) {
        message("Reading data...")
        #reading files in one file
        pb <- txtProgressBar(min = 0, max = length(files$file_name), style = 3)
        count <- 0
        for (i in files$file_name) {
            count <- count + 1
            setTxtProgressBar(pb, count)
            actual_file <- data.table::fread(file.path(DIR, i))
            if (count == 1) {
                actual_file <- data.table::fread(file.path(DIR, i))
                completed_table1 <- matrix(nrow = nrow(actual_file),
                                        ncol = length(files$file_name))
                rownames(completed_table1) <- actual_file[[1]]
                completed_table1[, 1] <- as.numeric(actual_file[[2]])
            } else {
                actual_file <- data.table::fread(file.path(DIR, i), select = 2)
                completed_table1[, count] <- as.numeric(actual_file[[1]])
            }
        }
        close(pb)
        tmp1 <- Var == "gene_files_not_tumor"
        tmp2 <- Var == "isoform_files_not_tumor"
        if (Var == "gene_files" || Var == "isoform_files") {
            assign("name_row", rownames(completed_table1),
                envir = get(envir_link))
            assign("completed_table1", completed_table1,
                envir = get(envir_link))
        } else if (tmp1 || tmp2) {
            assign("name_row_not_tumor", rownames(completed_table1),
                envir = get(envir_link))
            assign("completed_table1_not_tumor", completed_table1,
                envir = get(envir_link))
        }
    }

    export <- function() {
        tmp <- tolower(HTSeq) != "counts" && !is.null(HTSeq)
        if (tumorData) {
            if (normalization || tmp) {
                assign("gene_tumor_normalized",
                    string_vars[["envir_link"]]$completed_table1,
                    envir = get(envir_link))
            } else {
                assign("gene_tumor_not_normalized",
                    string_vars[["envir_link"]]$completed_table1,
                    envir = get(envir_link))
            }
            if (saveData) {
                message("\nSaving your data...\n")
                write.table(string_vars[["envir_link"]]$completed_table1,
                            paste0(DIR, "/", tumor, "_tumor_data.tsv"),
                            sep = "\t")
            }
            assign("patients", patients, envir = get(envir_link))
        } else {
            #tumor
            if (normalization || tmp) {
                assign("gene_tumor_normalized",
                        string_vars[["envir_link"]]$completed_table1,
                        envir = get(envir_link))
            } else {
                assign("gene_tumor_not_normalized",
                        string_vars[["envir_link"]]$completed_table1,
                        envir = get(envir_link))

            }
            #not tumor
            row.names(string_vars[["envir_link"]]$completed_table1_not_tumor) <- string_vars[["envir_link"]]$name_row_not_tumor
            if (normalization || tolower(HTSeq) != "counts") {
                colnames(string_vars[["envir_link"]]$completed_table1_not_tumor) <- patients_short
                assign("gene_not_tumor_normalized",
                    string_vars[["envir_link"]]$completed_table1_not_tumor,
                    envir = get(envir_link))
            } else {
                assign("gene_not_tumor_not_normalized",
                    string_vars[["envir_link"]]$completed_table1_not_tumor,
                    envir = get(envir_link))
            }
            if (saveData) {
                message("\nSaving your data...\n")
                write.table(string_vars[["envir_link"]]$completed_table1,
                            paste0(DIR, "/", tumor, "_tumor_data.tsv"),
                            sep = "\t")
                write.table(string_vars[["envir_link"]]$completed_table1_not_tumor,
                            paste0(DIR, "/", tumor, "_not_tumor_data.tsv"),
                            sep = "\t")
            }
            assign("patients", patients, envir = get(envir_link))
            assign("patients_not_tumor", patients_not_tumor,
                envir = get(envir_link))
        }
    }

    gene_finder <- function(completed_table, variable_name) {
        to_google <- ifelse(suppressWarnings(is.na(as.numeric(Name))),
                            paste0("(^", toupper(Name), "\\|", ")"),
                            paste0("(", "\\|", Name, ")$"))
        gene_selector_final <- grepl(to_google,
                                    rownames(completed_table), perl = TRUE)

        genes_1 <- completed_table[gene_selector_final, ]
        if (length(genes_1) == 0) {
            message("Gene '", Name, "' not found!!!\n")
            stop()
        }
        genes_transpo1 <- as.matrix(genes_1)
        colnames(genes_transpo1) <- rownames(completed_table)[gene_selector_final]
        assign(paste0(variable_name, toupper(Name)), genes_transpo1,
            envir = get(envir_link))
    }

    isoform_finder <- function(completed_table, variable_name) {
        Isoforms_1 <- completed_table[Name, , drop = FALSE]
        if (length(Isoforms_1) == 0) {
            message("'", Name, "' not found!!!\n")
            stop()
        }
        Isoforms_transpo1 <- t(Isoforms_1)
        Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)

        assign(paste0(variable_name, toupper(Name)), Isoforms_transpo1,
            envir = get(envir_link))
    }

    OPEN_meth <- function(meth_files_arg) {
        #reading files in one file
        pb <- txtProgressBar(min = 0, max = length(meth_files_arg), style = 3)
        count <- 0
        #try to open with lapply(list, function)
        for (file in file.path(DIR, meth_files_arg)) {
            count <- count + 1
            if (count > 1) {
                actual_file <- data.table::fread(file, sep = "\t",
                                                showProgress = FALSE,
                                                select = 2)
                meth_table[, count] <- as.numeric(actual_file[[1]])
            } else if (count == 1) {
                actual_file <- data.table::fread(file, sep = "\t",
                                                showProgress = FALSE)
                refererence_meth_table <- as.data.frame(actual_file[, Select_ref, with = FALSE])
                rownames(refererence_meth_table) <- as.character(as.matrix(actual_file[, 1]))
                meth_table <- matrix(nrow = nrow(actual_file),
                                    dimnames = list(rownames(refererence_meth_table), patients),
                                                ncol = length(meth_files_arg))
                meth_table[, 1] <- as.numeric(actual_file[[2]])
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)
        if (tolower(dataBase) == "legacy") {
            colnames(refererence_meth_table) <- refererence_meth_table[1, ]
            refererence_meth_table <- refererence_meth_table[-1, , drop = FALSE]
            meth_table <- meth_table[-1, , drop = FALSE]
        }
        assign("meth_table", meth_table, envir = get(envir_link))
        assign("refererence_meth_table", refererence_meth_table,
            envir = get(envir_link))
    }

    filter_meth_FUN <- function(meth_table, refererence_meth_table) {

        ### preparing data          FASTER WITH MATRIX
        message("Please wait... Filtering data!")
        pb <- txtProgressBar(min = 0, max = 3, style = 3)
        count <- 0
        setTxtProgressBar(pb, count)
        initial_row <- nrow(meth_table)
        # The beta value (β) estimates the methylation level of the CpG locus.
        # beta value (β) = methylated/unmethylated alleles
        #Remove sites containing NA for beta values based in cutoff
        beta_na_filter <- rowMeans(is.na(meth_table)) < cutoffBetaNA
        meth_table <- meth_table[beta_na_filter, ]
        refererence_meth_table <- refererence_meth_table[beta_na_filter, ]

        count <- 1
        setTxtProgressBar(pb, count)

        #Keep sites for which the beta values have standard deviation value
        #higher than cutoff
        beta_na_filter <- apply(meth_table, 1, sd, na.rm = TRUE) > cutoffBetasd
        meth_table <- meth_table[beta_na_filter, ]
        refererence_meth_table <- refererence_meth_table[beta_na_filter, ]

        count <- 2
        setTxtProgressBar(pb, count)

        #Remove CpGs on sex chromosomes
        tmp1 <- refererence_meth_table$Chromosome != "X"
        tmp2 <- refererence_meth_table$Chromosome != "Y"
        sex_filter <- tmp1 & tmp2
        sex_filter[is.na(sex_filter)] <- as.logical("FALSE")
        meth_table <- meth_table[sex_filter, ]
        refererence_meth_table <- refererence_meth_table[sex_filter, ]

        #only with gene symbol info (can be changed after released)
        tmp1 <- refererence_meth_table$Gene_Symbol != ""
        tmp2 <- refererence_meth_table$Gene_Symbol != "."
        gene_symbol_filter <- tmp1 & tmp2
        meth_table <- meth_table[gene_symbol_filter, ]
        refererence_meth_table <- refererence_meth_table[gene_symbol_filter, ]

        count <- 3
        setTxtProgressBar(pb, count)

        #Number of removed lines after filtering
        message(paste("It was filtered", initial_row - nrow(meth_table),
                        "lines from data!", sep = " "))
        close(pb)
        #normalization (3 options) Batch effect: systematic differences across
        #groups of samp The background was subtracted using the methylumi
        #package (method "noob") [1]. The signal intensity values were
        #normalized using the SWAN normalization method, as implemented in the
        #minfi package. or lumi R package or minfi Bioconductor packag(SWAN)
        # library(methylumi) ?? wateRmellow??
        if (nrow(meth_table) == 0) {
            message("All data were filtered out! Please, choose a different",
                    " cutoffBetaNA and cutoffBetasd parameters")
            stop()
        }

        #export to global enviroment
        assign("meth_table", meth_table, envir = get(envir_link))
        assign("refererence_meth_table", refererence_meth_table,
            envir = get(envir_link))
    }

    OPEN_miRNA <- function(mirna_files_arg) {
        tmp1 <- tolower(dataType) == "mirna gene quantification"
        tmp2 <- tolower(dataType) == "mirna expression quantification"
        tmp3 <- tolower(dataType) == "mirna isoform quantification"
        tmp4 <- tolower(dataType) == "isoform expression quantification"
        if (normalization) {
            if (tmp1 || tmp2) {
                Select <- c(3,4)
                first_col <- "reads_per_million_miRNA_mapped"
            } else if (tmp3 || tmp4) {
                message("Under development for this dataType!\n")
                stop()
            }

        } else {
            if (tmp1 || tmp2) {
                Select <- c(2,4)
                first_col <- "read_count"
            } else if (tmp3 || tmp4) {
                message("Under development for this dataType!\n")
                stop()
            }
        }
        #reading files in one file
        pb <- txtProgressBar(min = 0, max = length(mirna_files_arg), style = 3)
        count <- 0
        #try to open with lapply(list, function)
        for (file in file.path(DIR, mirna_files_arg)) {
            count <- count + 1
            if (count > 1) {
                actual_file <- data.table::fread(file, sep = "\t",
                                                showProgress = FALSE,
                                                select = Select)
                data_mirna_table[, count] <- as.numeric(actual_file[[1]])
                tmp <- as.character(actual_file[[2]])
                cross_mapped_mirna_table[, count] <- tmp
            } else if (count == 1) {
                actual_file <- data.table::fread(file, sep = "\t",
                                                showProgress = FALSE)
                tmp <- as.character(actual_file[["miRNA_ID"]])
                data_mirna_table <- matrix(nrow = nrow(actual_file),
                                        dimnames = list(tmp,
                                                        manifest$cases),
                                        ncol = length(mirna_files_arg))
                cross_mapped_mirna_table <- matrix(nrow = nrow(actual_file),
                                                dimnames = list(tmp,
                                                                manifest$cases),
                                                ncol = length(mirna_files_arg))
                data_mirna_table[, 1] <- as.numeric(actual_file[["read_count"]])
                tmp <- as.character(actual_file[["cross-mapped"]])
                cross_mapped_mirna_table[, 1] <- tmp
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)

        assign("data_mirna_table", data_mirna_table, envir = get(envir_link))
        assign("cross_mapped_mirna_table", cross_mapped_mirna_table,
            envir = get(envir_link))
    }
    # code ####
    #create env if not created yet
    if (missing(env) && tumorData && !onlyFilter) {
        assign(paste(toupper(tumor), toupper(dataBase),
                    gsub(" ", "_", tolower(dataType)),
                    "tumor_data", sep = "_"),
            new.env(parent = emptyenv()), envir = .GlobalEnv)
        envir_link <- paste(toupper(tumor), toupper(dataBase),
                            gsub(" ", "_", tolower(dataType)),
                            "tumor_data", sep = "_")
    } else if (missing(env) && !tumorData && !onlyFilter) {
        assign(paste(toupper(tumor), toupper(dataBase),
                    gsub(" ", "_", tolower(dataType)),
                    "both_data", sep = "_"),
            new.env(parent = emptyenv()), envir = .GlobalEnv)
        envir_link <- paste(toupper(tumor), toupper(dataBase),
                            gsub(" ", "_", tolower(dataType)),
                            "both_data", sep = "_")
    } else if (missing(env) && onlyFilter) {
        message("Please, before using 'onlyFilter'",
                " argument, insert the Environment name.")
    } else {
        envir_link <- deparse(substitute(env))
    }

    string_vars <- list(envir_link = get(envir_link))
    attr(envir_link, "name" ) <- paste0("Environment created by DOAGDC",
                                        " package, use its name in ",
                                        "'env' argument")

    dir.create(path = file.path(workDir, "DOAGDC", toupper(tumor), "Analyses"),
            showWarnings = FALSE)

    # mutation ####
    if (tolower(dataType) == "mutation") {
        if (tolower(dataBase) == "legacy") {
            #selecting platform file
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor), "mutation_data")
            manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"),
                                        select = c("file_name", "platform"))
            if (tolower(Platform) %in% "illumina ga") {
                tmp <- manifest[[1]][manifest$platform == "Illumina GA"]
                maf_files <- as.character(tmp)
            } else if (tolower(Platform) %in% "illumina hiseq") {
                tmp <- manifest[[1]][manifest$platform == "Illumina HiSeq"]
                maf_files <- as.character(tmp)
            }

            #checking for unmatched data
            if (length(maf_files) > 0) {
                message("Loading mutation file...")
            } else {
                stop(message(paste0("There is no data from ",
                                    tolower(Platform),
                                    " in your local storage!")))
            }

            message("\nIt will be open ", length(maf_files), " file(s)...\n")
            pb <- txtProgressBar(min = 0, max = length(maf_files), style = 3)
            count <- 0
            for (file in maf_files) {
                count <- count + 1
                assign(file, read.delim(file = file.path(DIR, file),
                                        comment.char = "#", sep = "\t",
                                        header = TRUE,
                                        stringsAsFactors = FALSE),
                    envir = get(envir_link))
                message("\nThe file is load as '", file, "' in ",
                        envir_link, " Environment.")

                if (saveData) {
                    write.table(x = get(file, envir = get(envir_link)),
                                file = file.path(DIR, paste0(file, ".tsv")),
                                row.names = FALSE, quote = FALSE, sep = "\t")
                }
                setTxtProgressBar(pb, count)
            }
            close(pb)

        } else if (tolower(dataBase) == "gdc") {
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                            "gdc_mutation_data")
            manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"),
                                        select = c(11,12))
            manifest$submitter_id <- unlist(lapply(strsplit(x = manifest$submitter_id,
                                                            split = "-",
                                                            perl = TRUE), "[[", 3))
            if (tolower(workflowType) %in% "all") {
                to_gunzip <- dir(path = DIR, pattern = "maf.gz$")
            } else if (tolower(workflowType) %in% "varscan") {
                tmp <- manifest[[1]][manifest$submitter_id == "varscan"]
                to_gunzip <- as.character(tmp)
            } else if (tolower(workflowType) %in% "mutect") {
                tmp <- manifest[[1]][manifest$submitter_id == "mutect"]
                to_gunzip <- as.character(tmp)
            } else if (tolower(workflowType) %in% "muse") {
                tmp <- manifest[[1]][manifest$submitter_id == "muse"]
                to_gunzip <- as.character(tmp)
            } else if (tolower(workflowType) %in% "somaticsniper") {
                tmp <- manifest[[1]][manifest$submitter_id == "somaticsniper"]
                to_gunzip <- as.character(tmp)
            }

            tryCatch(R.utils::gunzip(file.path(DIR, to_gunzip),
                    remove = FALSE),
                    error = function(e) {})
            import_maf <- gsub(".gz", "", to_gunzip)
            #openning files
            message("\nLoading mutation file...\n")
            maf <- read.delim(file = paste(DIR, import_maf, sep = "/"),
                            comment.char = "#", sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE)#, quote = FALSE)

            assign(import_maf, maf, envir = get(envir_link))
            message("\nThe file is open as '", import_maf, "' in ",
                    envir_link, " environment.")
            if (saveData) {
                write.table(x = get(import_maf, envir = get(envir_link)),
                            file = file.path(DIR, paste0(import_maf, ".tsv")),
                            row.names = FALSE, quote = FALSE, sep = "\t")
            }
        }
    }

    # methylation ####
    if (tolower(dataType) == "methylation") {
        if (!onlyFilter) {
            if (tolower(dataBase) == "legacy") {
                DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                "methylation_data")
                Select_ref <- c(3:5)
                patient_selector <- "submitter_id"
            } else if (tolower(dataBase) == "gdc") {
                DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                "gdc_methylation_data")
                Select_ref <- c(3:11)
                patient_selector <- "cases"
            }

            manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"),
                                        select = c("file_name",
                                                    patient_selector,
                                                    "platform"))

            tmp <- !grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                        sep = "|"),
                        x = dir(path = DIR))

            meth_files_downloaded <- dir(path = DIR)[tmp]
            manifest <- manifest[manifest$file_name %in% meth_files_downloaded, ]

            #only meth data!!
            if (sum(tmp) > 0) {
                tmp <- tolower(Platform) %in% "illumina human methylation 27"
                if (tolower(Platform) %in% "illumina human methylation 450") {
                    tmp <- manifest$platform == "Illumina Human Methylation 450"
                    meth_files <- manifest[tmp, ]
                } else if (tmp) {
                    tmp <- manifest$platform == "Illumina Human Methylation 27"
                    meth_files <- manifest[tmp, ]
                }

                #checking for unmatched data
                if (nrow(meth_files) > 0) {
                    message("Reading data...")
                } else {
                    message(paste0("There is no data from ", tolower(Platform),
                                " in your local storage!\n"))
                    stop()
                }

            } else {
                message("No Methylation data was downloaded!\n")
                stop()
            }

            if (tumorData) {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-",
                                                            tumor_regex,
                                                            "[aAbB]-"),
                                        meth_files[[2]])
                #not related the control samples
                meth_tumor_files <- as.character(meth_files[[1]][only_tumor_tissue])
                patients <- as.character(as.matrix(meth_files[only_tumor_tissue, 2]))
                assign("patients", patients, envir = get(envir_link))
                suppressWarnings(OPEN_meth(meth_tumor_files))
                #saving
                # write.table(x = meth_table,
                        # file = "./methylation_table/meth_table_completed.csv",
                #           row.names = TRUE, quote = FALSE)

                ##filtering
                meth_table <- string_vars[["envir_link"]]$meth_table
                refererence_meth_table <- string_vars[["envir_link"]]$refererence_meth_table
                filter_meth_FUN(meth_table, refererence_meth_table)
                tumor_meth_table_filtered <- string_vars[["envir_link"]]$meth_table
                reference_table_filtered <- string_vars[["envir_link"]]$refererence_meth_table

                if (tolower(dataBase) == "gdc") {
                    reference_table_filtered[, 4] <- unlist(lapply(strsplit(x = reference_table_filtered[, 4],
                                                                    split = ";", perl = TRUE), "[", 1))
                }

                assign("methylation_tumor_filtered", tumor_meth_table_filtered,
                    envir = get(envir_link))
                assign("reference_table_filtered", reference_table_filtered,
                    envir = get(envir_link))

                if (saveData) {
                    message("\nSaving data, this could take a while...\n")
                    #saving
                    tmp <- "methylation_tumor_filtered.tsv"
                    write.table(x = cbind(refererence_meth_table,
                                        tumor_meth_table_filtered),
                                file = file.path(DIR, tmp),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }

                if (!missing(Name)) {
                    select_gene <- grepl(pattern = Name,
                                        x = reference_table_filtered$Gene_Symbol)
                    if (sum(select_gene) == 0) {
                        message(cat(paste0(Name, " not found...\n\n"),
                                    "Please insert a valid gene symbol!!"))
                    }
                    assign(paste0("reference_table_filtered_selected_", Name),
                        reference_table_filtered[select_gene, ],
                        envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    tmp <- "methylation_tumor_filtered_selected_"
                    write.table(x = cbind(reference_table_filtered[select_gene, ],
                                        tumor_meth_table_filtered[select_gene, ]),
                                file = file.path(DIR,
                                                paste0(tmp ,Name, ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")

                    assign(paste0("methylation_tumor_filtered_selected_", Name),
                        tumor_meth_table_filtered[select_gene, ],
                        envir = get(envir_link))
                    assign("Name", Name, envir = get(envir_link))
                }

            } else {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                not_tumor_regex <- paste(10:19, collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-", tumor_regex,
                                                            "[aAbB]-"),
                                        meth_files)
                not_tumor_tissue <- grepl(pattern = paste0("-", not_tumor_regex,
                                                        "[aAbB]-"),
                                        meth_files)
                #not related the control samples

                if (sum(not_tumor_tissue) == 0) {
                    message("\nThere is no 'normal' data available to ",
                            tumor, " tumor!!\n\n")
                    stop()
                }

                ##filtering
                #tumor
                meth_tumor_files <- meth_files[only_tumor_tissue]
                patients <- as.character(as.matrix(meth_files[only_tumor_tissue, 2]))
                patients_tumor <- patients
                suppressWarnings(OPEN_meth(meth_tumor_files))
                meth_table <- string_vars[["envir_link"]]$meth_table
                refererence_meth_table <- string_vars[["envir_link"]]$refererence_meth_table
                filter_meth_FUN(meth_table, refererence_meth_table)
                tumor_meth_table_filtered <- string_vars[["envir_link"]]$meth_table
                reference_table_filtered <- string_vars[["envir_link"]]$refererence_meth_table

                ##filtering
                #not tumor
                meth_not_tumor_files <- string_vars[["envir_link"]]$meth_files[not_tumor_tissue]
                patients <- as.character(as.matrix(meth_files[not_tumor_tissue, 2]))
                patients_not_tumor <- patients
                suppressWarnings(OPEN_meth(meth_tumor_files))
                meth_table <- string_vars[["envir_link"]]$meth_table
                refererence_meth_table <- string_vars[["envir_link"]]$refererence_meth_table
                filter_meth_FUN(meth_table, refererence_meth_table)
                not_tumor_meth_table_filtered <- string_vars[["envir_link"]]$meth_table

                assign("methylation_tumor_filtered", tumor_meth_table_filtered,
                    envir = get(envir_link))
                assign("methylation_not_tumor_filtered",
                    not_tumor_meth_table_filtered,
                    envir = get(envir_link))
                assign("reference_table_filtered", reference_table_filtered,
                    envir = get(envir_link))

                #saving
                if (saveData) {
                    message("\nSaving data...\n")
                    tmp <- "tumor_meth_table_completed_filtered.tsv"
                    write.table(x = cbind(refererence_meth_table,
                                        tumor_meth_table_filtered),
                                file = file.path(DIR, tmp),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                    tmp <- "not_tumor_meth_table_completed_filtered.tsv"
                    write.table(x = cbind(refererence_meth_table,
                                        not_tumor_meth_table_filtered),
                                file = file.path(DIR, tmp),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }

                if (!missing(Name)) {
                    select_gene <- grepl(pattern = Name,
                                        x = reference_table_filtered$Gene_Symbol)
                    if (sum(select_gene) == 0) {
                        message(cat(paste0(Name, " not found...\n\n"),
                                    "Please insert a valid gene symbol!!"))
                        suppressWarnings(remove(meth_table,
                                                refererence_meth_table,
                                                envir = get(envir_link)))
                        stop()
                    }
                    assign(paste0("reference_table_filtered_selected_", Name),
                        reference_table_filtered[select_gene, ],
                        envir = get(envir_link))

                    assign(paste0("methylation_tumor_filtered_selected_", Name),
                        tumor_meth_table_filtered[select_gene, ],
                        envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    write.table(x = cbind(reference_table_filtered[select_gene, ],
                                        tumor_meth_table_filtered[select_gene, ]),
                                file = file.path(DIR,
                                                paste0("methylation_tumor_",
                                                "filtered_selected_",
                                                Name,".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")

                    # not tumor
                    assign(paste0("methylation_not_tumor_filtered_selected_",
                                Name),
                        not_tumor_meth_table_filtered[select_gene, ],
                        envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    tmp <- "methylation_not_tumor_filtered_selected_"
                    write.table(x = cbind(reference_table_filtered[select_gene, ],
                                        not_tumor_meth_table_filtered[select_gene, ]),
                                file = file.path(DIR,paste0(tmp, Name, ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")

                    assign("Name", Name, envir = get(envir_link))
                }
            }
            #removing useless variables
            suppressWarnings(remove(meth_table, refererence_meth_table,
                                    envir = get(envir_link)))
            invisible(gc())
        } else {
            if (!missing(Name)) {
                select_gene <- grepl(pattern = Name,
                                    x = string_vars[["envir_link"]]$reference_table_filtered$Gene_Symbol)
                if (sum(select_gene) == 0) {
                    message(cat(paste0(Name, " not found...\n\n"),
                                "Please insert a valid gene symbol!!"))
                    stop()
                }
                assign(paste0("reference_table_filtered_selected_", Name),
                    string_vars[["envir_link"]]$reference_table_filtered[select_gene, ],
                    envir = get(envir_link))

                assign(paste0("methylation_tumor_filtered_selected_", Name),
                    string_vars[["envir_link"]]$methylation_tumor_filtered[select_gene, ],
                    envir = get(envir_link))

                message("\nSaving required data...\n")
                #saving
                write.table(x = cbind(string_vars[["envir_link"]]$reference_table_filtered[select_gene, ],
                                    string_vars[["envir_link"]]$methylation_tumor_filtered[select_gene, ]),
                            file = file.path(DIR,
                                            paste0("methylation_tumor",
                                                    "_filtered_selected_",
                                                    Name,
                                                    ".tsv")),
                            row.names = TRUE, quote = FALSE, sep = "\t")

                # not tumor
                if (!tumorData) {
                    assign(paste0("methylation_not_tumor_filtered_selected_",
                                Name),
                        string_vars[["envir_link"]]$methylation_not_tumor_filtered[select_gene, ],
                        envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    tmp <- "methylation_tumor_filtered_selected_"
                    write.table(x = cbind(string_vars[["envir_link"]]$reference_table_filtered[select_gene, ],
                                        string_vars[["envir_link"]]$methylation_not_tumor_filtered[select_gene, ]),
                                file = file.path(DIR,
                                                paste0(tmp, Name, ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            }
        }
    }

    # gene and isoform ####
    if ("gene" == tolower(dataType) || "isoform" == tolower(dataType)) {
        if (!onlyFilter) {
            if (tolower(dataType) == "gene") {
                if (tolower(dataBase) == "legacy") {
                    DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                    "gene_data")
                    #prepare data selectors
                    gene_isoform()
                    message(tumor)
                    OPEN_genexpress(files = get(envir_link)$gene_files,
                                    Var = "gene_files")
                    patients <- string_vars[["envir_link"]]$gene_files$submitter_id
                    colnames(string_vars[["envir_link"]]$completed_table1) <- patients

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
                                coluns_2rename <- c(coluns_2rename,
                                                    grep(Patients,
                                                        patients_short)[-1])
                                names_2rename <- c(names_2rename,
                                                paste0(Patients,
                                                        seq((sum(selector) - 1))))
                            }
                        }
                        colnames(string_vars[["envir_link"]]$completed_table1)[coluns_2rename] <- names_2rename
                        colnames(string_vars[["envir_link"]]$completed_table1)[-coluns_2rename] <- unique(patients_short)
                    } else {
                        colnames(string_vars[["envir_link"]]$completed_table1) <- patients_short
                    }

                    # Search gene
                    if (tumorData) {

                        #Getting and export the final results
                        export()

                        if (!missing(Name)) {
                            if (normalization) {
                                to_google <- ifelse(SW(is.na(as.numeric(Name))),
                                                    paste0("(^", toupper(Name),
                                                        "\\|", ")"),
                                                    paste0("(", "\\|",
                                                        Name, ")$"))
                                gene_selector_final <- grepl(to_google,
                                                            rownames(string_vars[["envir_link"]]$completed_table1),
                                                            perl = TRUE)
                                genes_1 <- string_vars[["envir_link"]]$completed_table1[gene_selector_final, , drop = FALSE]

                                if (nrow(genes_1) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1 <- t(genes_1)
                                colnames(genes_transpo1) <- Name

                                assign(paste0("gene_tumor_normalized_selected_",
                                            toupper(Name)),
                                    genes_transpo1, envir = get(envir_link))

                            } else {
                                to_google <- ifelse(SW(is.na(as.numeric(Name))),
                                                    paste0("(^",
                                                                toupper(Name),
                                                                "\\|", ")"),
                                                    paste0("(", "\\|",
                                                                Name, ")$"))
                                gene_selector_final <- grepl(to_google,
                                                            rownames(string_vars[["envir_link"]]$completed_table1),
                                                            perl = TRUE)
                                genes_1 <- string_vars[["envir_link"]]$completed_table1[gene_selector_final, , drop = FALSE]

                                if (nrow(genes_1) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1 <- t(genes_1)
                                colnames(genes_transpo1) <- Name

                                assign(paste0("gene_tumor_not_",
                                    "normalized_selected_", toupper(Name)),
                                    genes_transpo1, envir = get(envir_link))
                            }
                        }
                    } else {
                        # #tumor and not tumor data
                        if (!missing(Name)) {
                            if (normalization) {
                                to_google <- ifelse(suppressWarnings(is.na(as.numeric(Name))),
                                                    paste0("(^", toupper(Name),
                                                                "\\|", ")"),
                                                    paste0("(", "\\|", Name,
                                                                        ")$"))
                                gene_selector_final <- grepl(to_google,
                                                            rownames(string_vars[["envir_link"]]$completed_table1), perl = TRUE)
                                genes_1 <- string_vars[["envir_link"]]$completed_table1[gene_selector_final, , drop = FALSE]
                                if (nrow(genes_1) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1 <- t(genes_1)
                                colnames(genes_transpo1) <- Name

                                assign(paste0("gene_tumor_",
                                    "normalized_selected_", toupper(Name)),
                                    genes_transpo1, envir = get(envir_link))
                            } else {
                                to_google <- ifelse(SW(is.na(as.numeric(Name))),
                                                    paste0("(^",
                                                                toupper(Name),
                                                                "\\|",
                                                                ")"),
                                                    paste0("(",
                                                                "\\|",
                                                                Name, ")$"))
                                gene_selector_final <- grepl(to_google,
                                                            rownames(string_vars[["envir_link"]]$completed_table1),
                                                            perl = TRUE)
                                genes_1 <- string_vars[["envir_link"]]$completed_table1[gene_selector_final, , drop = FALSE]
                                if (nrow(genes_1) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1 <- t(genes_1)
                                colnames(genes_transpo1) <- Name

                                assign(paste0("gene_tumor_not_",
                                            "normalized_selected_",
                                            toupper(Name)), genes_transpo1,
                                    envir = get(envir_link))

                            }
                        }

                        #not tumor
                        if (length(string_vars[["envir_link"]]$gene_files_not_tumor$file_name) == 0) {
                            message("There is no 'Normal' data in ", tumor,
                                    " tumor folder! Please, rerun after set",
                                    " 'tumorData = TRUE'.\n\n")
                            stop()
                        }
                        OPEN_genexpress(files = string_vars[["envir_link"]]$gene_files_not_tumor,
                                        Var = "gene_files_not_tumor")
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$completed_table1_not_tumor
                        patients_not_tumor <- string_vars[["envir_link"]]$gene_files_not_tumor$submitter_id
                        colnames(completed_table1_not_tumor) <- patients_not_tumor



                        patients_short <- unname(sapply(patients_not_tumor,
                                                                function(w) {
                            paste(unlist(strsplit(w, "-"))[1:3],
                                                                collapse = "-")
                        }))

                        #Getting and export the final results
                        export()

                        # duplicate fix
                        if (length(patients_not_tumor) != length(unique(patients_short))) {
                            coluns_2rename <- numeric()
                            names_2rename <- character()
                            for (Patients in unique(patients_short)) {
                                selector <- Patients == patients_short
                                if (sum(selector) > 1) {
                                    coluns_2rename <- c(coluns_2rename,
                                            grep(Patients, patients_short)[-1])
                                    names_2rename <- c(names_2rename,
                                                    paste0(Patients,
                                                    seq((sum(selector) - 1))))
                                }
                            }
                            colnames(completed_table1_not_tumor)[coluns_2rename] <- names_2rename
                            colnames(completed_table1_not_tumor)[-coluns_2rename] <- unique(patients_short)
                        } else {
                            colnames(completed_table1_not_tumor) <- patients_short
                        }

                        rownames(completed_table1_not_tumor) <- string_vars[["envir_link"]]$name_row_not_tumor

                        if (!missing(Name)) {
                            if (normalization) {
                                to_google_not_tumor <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                            paste0("(",
                                                                        Name,
                                                                        "\\|",
                                                                        ")"))
                                gene_selector_final_not_tumor <- grepl(to_google_not_tumor,
                                                                    rownames(string_vars[["envir_link"]]$completed_table1),
                                                                    perl = TRUE)
                                genes_1_not_tumor <- completed_table1_not_tumor[gene_selector_final_not_tumor, , drop = FALSE]

                                if (nrow(genes_1_not_tumor) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1_not_tumor <- t(genes_1_not_tumor)
                                colnames(genes_transpo1_not_tumor) <- Name

                                FILTERED_RESULTS_not_tumor <- genes_transpo1_not_tumor
                                assign(paste0("gene_not_tumor_",
                                                        "normalized_selected_",
                                                        toupper(Name)),
                                    genes_transpo1_not_tumor,
                                    envir = get(envir_link))
                            } else {
                                to_google_not_tumor <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                            paste0("(", Name,
                                                                "\\|", ")"))
                                gene_selector_final_not_tumor <- grepl(to_google_not_tumor,
                                                                    rownames(string_vars[["envir_link"]]$completed_table1),
                                                                    perl = TRUE)
                                genes_1_not_tumor <- completed_table1_not_tumor[gene_selector_final_not_tumor, , drop = FALSE]

                                if (nrow(genes_1_not_tumor) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1_not_tumor <- t(genes_1_not_tumor)
                                colnames(genes_transpo1_not_tumor) <- Name

                                FILTERED_RESULTS_not_tumor <- genes_transpo1_not_tumor
                                tmp <- "gene_not_tumor_not_normalized_selected_"
                                assign(paste0(tmp, Name),
                                        genes_transpo1_not_tumor,
                                        envir = get(envir_link))
                            }
                        }
                    }
                } else if (tolower(dataBase) == "gdc") {
                    DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                    "gdc_gene_data")

                    #open manifest
                    codes <- dir(DIR, paste0(toupper(HTSeq), "_manifest.sdrf$"))
                    codigos <- read.table(file.path(DIR, codes[1]),
                                        stringsAsFactors = FALSE,
                                        header = TRUE, sep = "\t")
                    codigos <- codigos[, c("file_name", "cases")]

                    tumor_regex <- paste(formatC(seq(1:9), width = 2,
                                                flag = "0"),
                                        collapse = "[aAbB]-|-")

                    # Search gene
                    if (tumorData) {
                        codigos <- codigos[grepl(paste0("-", tumor_regex, "-"),
                                                x = codigos$cases), ]
                        remove_duplicated <- duplicated(codigos$file_name)
                        codigos <- codigos[!remove_duplicated, ]
                        rownames(codigos) <- codigos$file_name
                        patients <- codigos$cases

                        #concatenate and finishing table
                        OPEN_genexpress(files = codigos, Var = "gene_files")
                        rownames(string_vars[["envir_link"]]$completed_table1) <- string_vars[["envir_link"]]$name_row
                        colnames(string_vars[["envir_link"]]$completed_table1) <- patients

                        patients_short <- unname(sapply(patients, function(w) {
                            paste(unlist(strsplit(w, "-"))[1:3],
                                                                collapse = "-")
                        }))

                        #Getting and export the final results
                        export()

                        # duplicate fix
                        if (length(patients) != length(unique(patients_short))) {
                            coluns_2rename <- numeric()
                            names_2rename <- character()
                            for (Patients in unique(patients_short)) {
                                selector <- Patients == patients_short
                                if (sum(selector) > 1) {
                                    coluns_2rename <- c(coluns_2rename,
                                            grep(Patients, patients_short)[-1])
                                    names_2rename <- c(names_2rename,
                                                    paste0(Patients,
                                                    seq((sum(selector) - 1))))
                                }
                            }
                            colnames(string_vars[["envir_link"]]$completed_table1)[coluns_2rename] <- names_2rename
                            colnames(string_vars[["envir_link"]]$completed_table1)[-coluns_2rename] <- unique(patients_short)
                        } else {
                            colnames(string_vars[["envir_link"]]$completed_table1) <- patients_short
                        }

                        #saving completed_table
                        # write.table(string_vars[["envir_link"]]$completed_table1,
                        #             paste0(DIR, "/concatenate_",
                        #                                         HTSeq, ".tsv"),
                        #             sep = "\t")


                        if (!missing(Name)) {
                            if (normalization) {
                                genes_1 <- string_vars[["envir_link"]]$completed_table1[Name, , drop = FALSE]

                                if (nrow(genes_1) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1 <- t(genes_1)
                                colnames(genes_transpo1) <- Name

                                FILTERED_RESULTS <- genes_transpo1
                                assign(paste0("gene_tumor_normalized_selected_",
                                            toupper(Name)), genes_transpo1,
                                    envir = get(envir_link))
                            } else {
                                genes_1 <- string_vars[["envir_link"]]$completed_table1[Name, , drop = FALSE]

                                if (nrow(genes_1) == 0) {
                                    stop("Gene not found")
                                }
                                genes_transpo1 <- t(genes_1)
                                colnames(genes_transpo1) <- Name

                                FILTERED_RESULTS <- genes_transpo1
                                assign(paste0("gene_tumor_not_",
                                        "normalized_selected_", toupper(Name)),
                                    genes_transpo1, envir = get(envir_link))
                            }
                        }
                    } else {
                        not_tumor_regex <- paste(10:19, collapse = "[aAbB]-|-")

                        codigos_not_tumor <- codigos[grepl(paste0("-",
                                                        not_tumor_regex, "-"),
                                                        x = codigos$cases), ]
                        remove_duplicated_not_tumor <- duplicated(codigos_not_tumor$file_name)
                        codigos_not_tumor <- codigos_not_tumor[!remove_duplicated_not_tumor, ]

                        rownames(codigos_not_tumor) <- codigos_not_tumor$file_name

                        #concatenate
                        if (length(codigos_not_tumor$file_name) == 0) {
                            message("There is no 'Normal' data in ", tumor,
                                    " tumor folder! Please, rerun after set",
                                    " 'tumorData = TRUE'.",
                                    "\n\n")
                            stop()
                        }
                        OPEN_genexpress(files = codigos_not_tumor,
                                                Var = "gene_files_not_tumor")
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$completed_table1_not_tumor
                        patients_not_tumor <- string_vars[["envir_link"]]$gene_files_not_tumor$submitter_id
                        colnames(completed_table1_not_tumor) <- patients_not_tumor

                        patients_short <- unname(sapply(patients_not_tumor,
                                        function(w) {
                                            paste(unlist(strsplit(w, "-"))[1:3],
                                                                collapse = "-")
                        }))
                        # duplicate fix
                        if (length(patients_not_tumor) != length(unique(patients_short))) {
                            coluns_2rename <- numeric()
                            names_2rename <- character()
                            for (Patients in unique(patients_short)) {
                                selector <- Patients == patients_short
                                if (sum(selector) > 1) {
                                    coluns_2rename <- c(coluns_2rename,
                                                        grep(Patients,
                                                        patients_short)[-1])
                                    names_2rename <- c(names_2rename,
                                                    paste0(Patients,
                                                    seq((sum(selector) - 1))))
                                }
                            }
                            colnames(completed_table1_not_tumor)[coluns_2rename] <- names_2rename
                            colnames(completed_table1_not_tumor)[-coluns_2rename] <- unique(patients_short)
                        } else {
                            colnames(completed_table1_not_tumor) <- patients_short
                        }
                    }
                    if (tolower(HTSeq) != "counts") {
                        assign("HTSeq_normalized", toupper(HTSeq),
                                                    envir = get(envir_link))
                    }
                }
            }

            if (tolower(dataType) == "isoform") {
                DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                                                "isoform_data")
                #prepare data selectors
                gene_isoform()
                message(tumor)
                OPEN_genexpress(files = string_vars[["envir_link"]]$isoform_files, Var = "isoform_files")
                completed_table1_tumor <- string_vars[["envir_link"]]$completed_table1
                patients <- string_vars[["envir_link"]]$isoform_files$submitter_id
                colnames(completed_table1_tumor) <- patients

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
                            coluns_2rename <- c(coluns_2rename, grep(Patients,
                                                        patients_short)[-1])
                            names_2rename <- c(names_2rename, paste0(Patients,
                                                    seq((sum(selector) - 1))))
                        }
                    }
                    colnames(completed_table1_tumor)[coluns_2rename] <- names_2rename
                    colnames(completed_table1_tumor)[-coluns_2rename] <- unique(patients_short)
                } else {
                    colnames(completed_table1_tumor) <- patients_short
                }

                if (!tumorData) {
                    if (nrow(string_vars[["envir_link"]]$isoform_files_not_tumor) == 0) {
                        message("There is no 'Normal' data in ", tumor,
                                " tumor folder! Please, rerun after set ",
                                "'tumorData = TRUE'.",
                                "\n\n")
                        stop()
                    }

                    OPEN_genexpress(files = string_vars[["envir_link"]]$isoform_files_not_tumor,
                                    Var = "isoform_files_not_tumor")
                    patients_not_tumor <- string_vars[["envir_link"]]$isoform_files_not_tumor$submitter_id
                    colnames(string_vars[["envir_link"]]$completed_table1_not_tumor) <- patients_not_tumor

                    patients_short <- unname(sapply(patients_not_tumor,
                                    function(w) {
                                        paste(unlist(strsplit(w, "-"))[1:3],
                                                                collapse = "-")
                    }))
                    # duplicate fix
                    if (length(patients_not_tumor) != length(unique(patients_short))) {
                        coluns_2rename <- numeric()
                        names_2rename <- character()
                        for (Patients in unique(patients_short)) {
                            selector <- Patients == patients_short
                            if (sum(selector) > 1) {
                                coluns_2rename <- c(coluns_2rename,
                                                    grep(Patients,
                                                        patients_short)[-1])
                                names_2rename <- c(names_2rename,
                                                    paste0(Patients,
                                                    seq((sum(selector) - 1))))
                            }
                        }
                        colnames(string_vars[["envir_link"]]$completed_table1_not_tumor)[coluns_2rename] <- names_2rename
                        colnames(string_vars[["envir_link"]]$completed_table1_not_tumor)[-coluns_2rename] <- unique(patients_short)
                    } else {
                        colnames(string_vars[["envir_link"]]$completed_table1_not_tumor) <- patients_not_tumor
                    }

                    #Getting and export the final results
                    if (tumorData) {
                        if (normalization) {
                            assign("isoform_tumor_normalized",
                                completed_table1_tumor,
                                envir = get(envir_link))
                        } else {
                            assign("isoform_tumor_not_normalized",
                                completed_table1_tumor,
                                envir = get(envir_link))
                        }
                        if (saveData) {
                            message("\nSaving your data...\n")
                            write.table(completed_table1_tumor,
                                        paste0(DIR, "/", tumor,
                                            "_tumor_data.tsv"), sep = "\t")
                        }
                        assign("patients", patients, envir = get(envir_link))

                    } else {
                        #tumor
                        if (normalization) {
                            assign("isoform_tumor_normalized",
                                completed_table1_tumor,
                                envir = get(envir_link))
                        } else {
                            assign("isoform_tumor_not_normalized",
                                completed_table1_tumor,
                                envir = get(envir_link))

                        }
                        #not tumor
                        row.names(string_vars[["envir_link"]]$completed_table1_not_tumor) <- string_vars[["envir_link"]]$name_row_not_tumor
                        if (normalization) {
                            assign("isoform_not_tumor_normalized", string_vars[["envir_link"]]$completed_table1_not_tumor,
                                envir = get(envir_link))
                        } else {
                            assign("isoform_not_tumor_not_normalized", string_vars[["envir_link"]]$completed_table1_not_tumor,
                                envir = get(envir_link))
                        }
                        if (saveData) {
                            message("\nSaving your data...\n")
                            write.table(completed_table1_tumor,
                                        paste0(DIR, "/", tumor,
                                            "_tumor_data.tsv"), sep = "\t")
                            write.table(string_vars[["envir_link"]]$completed_table1_not_tumor,
                                        paste0(DIR, "/",
                                            tumor, "_not_tumor_data.tsv"),
                                        sep = "\t")
                        }
                        assign("patients", patients, envir = get(envir_link))
                        assign("patients_not_tumor", patients_not_tumor,
                            envir = get(envir_link))
                    }

                    # Search gene
                    if (!missing(Name)) {
                        if (normalization) {

                            #tumor
                            Isoforms_1 <- completed_table1_tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms_1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms_transpo1 <- t(Isoforms_1)
                            Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)
                            colnames(Isoforms_transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms_transpo1
                            assign(paste0("isoform_tumor_normalized_selected_",
                                        toupper(Name)), Isoforms_transpo1,
                                envir = get(envir_link))


                            #not tumor
                            Isoforms_1 <- string_vars[["envir_link"]]$completed_table1_not_tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms_1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms_transpo1 <- t(Isoforms_1)
                            Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)
                            colnames(Isoforms_transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms_transpo1
                            assign(paste0("isoform_not_tumor_",
                                                    "normalized_selected_",
                                        toupper(Name)), Isoforms_transpo1,
                                envir = get(envir_link))

                        } else {
                            #tumor
                            Isoforms_1 <- completed_table1_tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms_1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms_transpo1 <- t(Isoforms_1)
                            Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)
                            colnames(Isoforms_transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms_transpo1
                            assign(paste0("isoform_tumor_not_",
                                                    "normalized_selected_",
                                                    toupper(Name)),
                                                    Isoforms_transpo1,
                                                    envir = get(envir_link))


                            #not tumor
                            Isoforms_1 <- string_vars[["envir_link"]]$completed_table1_not_tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms_1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms_transpo1 <- t(Isoforms_1)
                            Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)
                            colnames(Isoforms_transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms_transpo1
                            assign(paste0("isoform_not_tumor_not_",
                                                    "normalized_selected_",
                                                    toupper(Name)),
                                                    Isoforms_transpo1,
                                                    envir = get(envir_link))
                        }
                    }
                } else {
                    if (!missing(Name)) {
                        if (normalization) {
                            Isoforms_1 <- completed_table1_tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms_1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms_transpo1 <- t(Isoforms_1)
                            Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)
                            colnames(Isoforms_transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms_transpo1
                            assign(paste0("isoform_tumor_normalized_selected_",
                                        toupper(Name)), Isoforms_transpo1,
                                envir = get(envir_link))
                        } else {
                            Isoforms_1 <- completed_table1_tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms_1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms_transpo1 <- t(Isoforms_1)
                            Isoforms_transpo1 <- as.data.frame(Isoforms_transpo1)
                            colnames(Isoforms_transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms_transpo1
                            assign(paste0("isoform_tumor_not_",
                                                    "normalized_selected_",
                                                    toupper(Name)),
                                                    Isoforms_transpo1,
                                                    envir = get(envir_link))
                        }
                    }
                }
            }

            suppressWarnings(remove(genes_normalized, isoforms_normalized,
                                    genes_not_normalized,
                                    isoforms_not_normalized,
                                    genes_normalized_not_tumor,
                                    isoforms_normalized_not_tumor,
                                    genes_not_normalized_not_tumor,
                                    isoforms_not_normalized_not_tumor,
                                    envir = get(envir_link)))
            suppressWarnings(remove(completed_table1,
                                    completed_table1_not_tumor,
                                    codigos_not_tumor, isoform_files,
                                    isoform_files_not_tumor, gene_files,
                                    codigos_tumor, gene_files_not_tumor,
                                    envir = get(envir_link)))
            invisible(gc(verbose = FALSE))
        } else {
            if (tolower(dataType) == "gene" && tolower(dataBase) == "legacy") {

                if (tumorData) {
                    if (normalization) {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        gene_finder(completed_table1_tumor,
                            variable_name = "gene_tumor_normalized_selected_")
                    } else {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        gene_finder(completed_table1_tumor,
                                    variable_name = paste0("gene_tumor_not_",
                                                    "normalized_selected_"))
                    }
                } else {
                    #tumor
                    if (normalization) {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        gene_finder(completed_table1_tumor,
                                        variable_name = paste0("gene_tumor_",
                                        "normalized_selected_"))
                    } else {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        gene_finder(completed_table1_tumor,
                                    variable_name = paste0("gene_tumor_not_",
                                    "normalized_selected_"))
                    }
                    #not tumor
                    if (normalization) {
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$gene_not_tumor_normalized
                        gene_finder(completed_table1_not_tumor,
                                    variable_name = paste0("gene_not_tumor",
                                    "_normalized_selected_"))
                    } else {
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$gene_not_tumor_not_normalized
                        gene_finder(completed_table1_not_tumor,
                                    variable_name = paste0("gene_not_tumor_",
                                    "not_normalized_selected_"))
                    }
                }
            } else if (tolower(dataType) == "gene" && tolower(dataBase) == "gdc") {
                if (tumorData) {
                    if (normalization) {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        isoform_finder(completed_table1_tumor,
                                        variable_name = paste0("gene_tumor_",
                                        "normalized_selected_"))
                    } else {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        isoform_finder(completed_table1_tumor,
                                    variable_name = paste0("gene_tumor_not_",
                                    "normalized_selected_"))
                    }
                } else {
                    #tumor
                    if (normalization) {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        isoform_finder(completed_table1_tumor,
                                        variable_name = paste0("gene_tumor_",
                                        "normalized_selected_"))
                    } else {
                        completed_table1_tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        isoform_finder(completed_table1_tumor,
                                        variable_name = paste0("gene_tumor_",
                                        "not_normalized_selected_"))
                    }
                    #not tumor
                    if (normalization) {
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$gene_not_tumor_normalized
                        isoform_finder(completed_table1_not_tumor,
                                            variable_name = paste0("gene_not_",
                                            "tumor_normalized_selected_"))
                    } else {
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$gene_not_tumor_not_normalized
                        isoform_finder(completed_table1_not_tumor,
                                            variable_name = paste0("gene_not_",
                                            "tumor_not_normalized_selected_"))
                    }
                }
            } else if (tolower(dataType) == "isoform") {
                if (tumorData) {
                    if (normalization) {
                        completed_table1_tumor <- string_vars[["envir_link"]]$isoform_tumor_normalized
                        isoform_finder(completed_table1_tumor,
                                            variable_name = paste0("isoform_",
                                            "tumor_normalized_selected_"))
                    } else {
                        completed_table1_tumor <- string_vars[["envir_link"]]$isoform_tumor_not_normalized
                        isoform_finder(completed_table1_tumor,
                                            variable_name = paste0("isoform",
                                            "_tumor_not_normalized_selected_"))
                    }
                } else {
                    #tumor
                    if (normalization) {
                        completed_table1_tumor <- string_vars[["envir_link"]]$isoform_tumor_normalized
                        isoform_finder(completed_table1_tumor,
                                            variable_name = paste0("isoform_",
                                            "tumor_normalized_selected_"))
                    } else {
                        completed_table1_tumor <- string_vars[["envir_link"]]$isoform_tumor_not_normalized
                        isoform_finder(completed_table1_tumor,
                                            variable_name = paste0("isoform",
                                            "_tumor_not_normalized_selected_"))
                    }
                    #not tumor
                    if (normalization) {
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$isoform_not_tumor_normalized
                        isoform_finder(completed_table1_not_tumor,
                                            variable_name = paste0("isoform",
                                            "_not_tumor_normalized_selected_"))
                    } else {
                        completed_table1_not_tumor <- string_vars[["envir_link"]]$isoform_not_tumor_not_normalized
                        isoform_finder(completed_table1_not_tumor,
                                        variable_name = paste0("isoform",
                                        "_not_tumor_not_normalized_selected_"))
                    }
                }
            }
        }
    }

    # clinical ####
    if ("clinical" == tolower(dataType)) {
        if (tolower(dataBase) == "legacy") {
            #selecting platform file
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                                            "clinical_data")
            patient_file <- dir(path = DIR, pattern = "clinical_patient")

            clinical_df <- data.table::fread(input = file.path(DIR,
                                                                patient_file))

            clinical_df <- as.data.frame(clinical_df[-c(1,2),])
            rownames(clinical_df) <- clinical_df[, 2]
            clinical_df <- clinical_df[, -2]

            assign("clinical_df", clinical_df, envir = get(envir_link))

        } else if (tolower(dataBase) == "gdc") {
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                            "gdc_biospecimen_data")
            patient_files <- dir(path = DIR, pattern = "xml$")

            pb <- txtProgressBar(min = 0, max = length(patient_files),
                                style = 3)
            count <- 0
            #concatenating all the XML files at once
            for (i in file.path(DIR, patient_files)) {
                data_df <- XML::xmlToDataFrame(file.path(DIR, i))
                data_df <- na.omit(reshape2::melt(t(data_df)))
                rownames(data_df) <- data_df[, 1, drop = TRUE]
                data_df <- t(data_df[, 3, drop = FALSE])
                if (count == 0) {
                    clinical_df <- as.data.frame(matrix(ncol = ncol(data_df),
                                                        nrow = length(patient_files)))
                    RowNames <- character(length = length(patient_files))
                    colnames(clinical_df) <- colnames(data_df)
                    clinical_df[1, ] <- data_df
                    RowNames[1] <- data_df[, "bcr_patient_barcode"]
                } else if (count <= length(file.path(DIR, patient_files))) {
                    clinical_df[count, ] <- data_df
                    RowNames[count] <- data_df[, "bcr_patient_barcode"]
                } else{
                    clinical_df[count, ] <- data_df
                    RowNames[count] <- data_df[, "bcr_patient_barcode"]
                    row.names(clinical_df) <- RowNames
                    assign("clinical_df", clinical_df, envir = get(envir_link))
                }
                count <- count + 1
                setTxtProgressBar(pb, count)
            }
            close(pb)
        }
    }

    # protein ####
    if ("protein" == tolower(dataType)) {
        if (!onlyFilter) {
            if (tolower(dataBase) == "gdc") {
                stop(message("\nThrere is no protein expression ",
                            "data in GDC data base!!",
                            "\nPlease use 'legacy' data base"))
            }

            DIR <- file.path(workDir, "DOAGDC", toupper(tumor), "protein_data")

            if (!dir.exists(file.path(workDir, "DOAGDC", toupper(tumor),
                                    "mage_data"))) {
                message("Downloading magetab data...")
                suppressWarnings(download_gdc(dataType = "mage",
                                            dataBase = dataBase,
                                            tumor = tumor,
                                            workDir = workDir))
            }

            desing_array_DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                        "mage_data")

            to_gunzip <- file.path(desing_array_DIR,
                                paste0("mdanderson.org_", toupper(tumor),
                                    ".MDA_RPPA_Core.mage-tab.1.1.0.tar.gz"))
            ## or, if you just want to extract the target file:
            untar(to_gunzip, exdir = desing_array_DIR)

            design_array <- data.table::fread(file.path(desing_array_DIR,
                                                        paste0("mdanderson.org_", toupper(tumor), ".MDA_RPPA_Core.mage-tab.1.1.0"),
                                                        paste0("mdanderson.org_", toupper(tumor), ".MDA_RPPA_Core.array_design.txt")),
                                            select = c(6,7))

            design_array <- unique(design_array)

            rownames(design_array) <- design_array[[2]]

            patient_files <- dir(path = DIR, pattern = "protein_expression")

            pb <- txtProgressBar(min = 0, max = length(patient_files),
                                style = 3)
            count <- 0
            for (i in file.path(DIR, patient_files)) {
                count <- count + 1
                if (count == 1) {
                    RowNames <- character(length = length(patient_files))
                    patient_id <- character(length = length(patient_files))
                    pr <- data.table::fread(i)
                    uuid <- colnames(pr)[2]
                    protein_table <- matrix(ncol = length(patient_files),
                                            nrow = nrow(pr) - 1)
                    colnames(protein_table) <- 1:ncol(protein_table)
                    patient_id[1] <- as.character(design_array[[1]][grep(uuid,
                                                        design_array[[2]])])
                    pr <- pr[-1, ]
                    rownames(protein_table) <- as.character(pr[[1]])
                    protein_table[, 1] <- as.numeric(pr[[2]])
                    colnames(protein_table)[1] <- patient_id[1]
                } else {
                    pr <- data.table::fread(i, select = 2)
                    uuid <- colnames(pr)[1]
                    patient_id[count] <- as.character(design_array[[1]][grep(uuid, design_array[[2]])])
                    pr <- pr[-1, ]
                    protein_table[, count] <- as.numeric(pr[[1]])
                    colnames(protein_table)[count] <- patient_id[count]
                }
                setTxtProgressBar(pb, count)
            }
            close(pb)

            if (tumorData) {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-", tumor_regex,
                                                            "[aAbB]-"),
                                        colnames(protein_table))
                #not related the control samples
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
                            coluns_2rename <- c(coluns_2rename,
                                                grep(Patients,
                                                    patients_short)[-1])
                            names_2rename <- c(names_2rename,
                                            paste0(Patients,
                                                    seq((sum(selector) - 1))))
                        }
                    }
                    colnames(protein_table)[coluns_2rename] <- names_2rename
                    colnames(protein_table)[-coluns_2rename] <- unique(patients_short)
                } else {
                    colnames(protein_table) <- patients_short
                }

                protein_table <- gsub("-.-.", "", rownames(protein_table))

                assign("protein_patients", patients, envir = get(envir_link))
                assign("protein_tumor_normalized", protein_table,
                    envir = get(envir_link))

                if (saveData) {
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = protein_table,
                                file = file.path(DIR,
                                            "protein_tumor_normalized.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            } else {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                not_tumor_regex <- paste(10:19, collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-", tumor_regex,
                                                            "[aAbB]-"),
                                        colnames(protein_table))
                not_tumor_tissue <- grepl(pattern = paste0("-",
                                                            not_tumor_regex,
                                                        "[aAbB]-"),
                                        colnames(protein_table))

                tumor_patients <- as.character(colnames(protein_table)[only_tumor_tissue])
                tumor_protein_table <- protein_table[, tumor_patients]
                tumor_protein_table <- gsub("-.-.", "",
                                            rownames(tumor_protein_table))

                assign("tumor_patients", tumor_patients,
                    envir = get(envir_link))
                assign("protein_tumor_normalized", tumor_protein_table,
                    envir = get(envir_link))

                if (sum(not_tumor_tissue) == 0) {
                    message("There is no 'Normal' data in ", tumor,
                            " tumor folder! Please, rerun after",
                            " set 'tumorData = TRUE'.",
                            "\n\n")
                    stop()
                }
                not_tumor_patients <- as.character(colnames(protein_table)[not_tumor_tissue])
                protein_not_tumor <- protein_table[, not_tumor_patients]
                protein_not_tumor <- gsub("-.-.", "",
                                                rownames(protein_not_tumor))

                assign("not_tumor_patients", not_tumor_patients,
                    envir = get(envir_link))
                assign("protein_not_tumor_normalized", protein_not_tumor,
                    envir = get(envir_link))

                if (saveData) {
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = tumor_protein_table,
                                file = file.path(DIR,
                                            "protein_tumor_normalized.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                    write.table(x = protein_not_tumor,
                                file = file.path(DIR,
                                        "protein_not_tumor_normalized.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            }

            if (!missing(Name)) {
                if (tumorData) {
                    if (normalization) {
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is only normalized protein data!")
                        stop()
                    }
                    protein_selector_final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector,
                                                    perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein_selector_final, ])
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_",
                                toupper(Name)),
                        protein_exp_selected, envir = get(envir_link))
                } else {
                    # tumor
                    if (normalization) {
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein_selector_final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector,
                                                    perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein_selector_final, ])
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_",
                                toupper(Name)),
                        protein_exp_selected, envir = get(envir_link))

                    # not tumor
                    protein_selector <- rownames(string_vars[["envir_link"]]$protein_not_tumor_normalized)
                    protein_selector_final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector,
                                                    perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_not_tumor_normalized[protein_selector_final, ])
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_not_tumor_normalized_selected_",
                                toupper(Name)),
                        protein_exp_selected, envir = get(envir_link))
                }
            }
        } else {
            if (!missing(Name)) {
                if (tumorData) {
                    if (normalization) {
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein_selector_final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector,
                                                    perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein_selector_final, , drop = FALSE])
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_",
                                toupper(Name)),
                        protein_exp_selected, envir = get(envir_link))
                } else {
                    # tumor
                    if (normalization) {
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein_selector_final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector,
                                                    perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein_selector_final, , drop = FALSE])
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_",
                                toupper(Name)),
                        protein_exp_selected, envir = get(envir_link))

                    # not tumor
                    protein_selector <- rownames(string_vars[["envir_link"]]$protein_not_tumor_normalized)
                    protein_selector_final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector,
                                                    perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_not_tumor_normalized[protein_selector_final, , drop = FALSE])
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_not_tumor_normalized_selected_",
                                toupper(Name)),
                        protein_exp_selected, envir = get(envir_link))
                }
            }
        }
    }

    # mirna ####
    tmp1 <- "mirna" %in% strsplit(tolower(dataType), split = " ")[[1]][1]
    tmp2 <- "isoform expression quantification" %in% tolower(dataType)
    if (tmp1 || tmp2) {
        if (!onlyFilter) {
            if (tolower(dataBase) == "legacy") {
                DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                paste0(tolower(dataType), "_data"))

                manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"),
                                            select = c("file_name", "cases",
                                                        "platform"))

                tmp <- !grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                            sep = "|"),
                            x = dir(path = DIR))

                mirna_files_downloaded <- dir(path = DIR)[tmp]
                if (use_hg19_mirbase20) {
                    mirna_files_downloaded <- mirna_files_downloaded[grepl("hg19.mirbase20", mirna_files_downloaded, fixed = TRUE)]
                } else {
                    mirna_files_downloaded <- mirna_files_downloaded[!grepl("hg19.mirbase20", mirna_files_downloaded, fixed = TRUE)]
                }
                manifest <- manifest[manifest$file_name %in% mirna_files_downloaded, ]

                #only mirna data!!
                if (sum(tmp) > 0) {
                    if (tolower(Platform) %in% "illumina hiseq") {
                        tmp <- manifest$platform == "Illumina HiSeq"
                        mirna_files <- manifest[tmp, ]
                    } else if (tolower(Platform) %in% "illumina ga") {
                        tmp <- manifest$platform == "Illumina GA"
                        mirna_files <- manifest[tmp, ]
                    } else if (tolower(Platform) %in% "h-mirna_8x15kv2") {
                        tmp <- manifest$platform == "H-miRNA_8x15Kv2"
                        mirna_files <- manifest[tmp, ]
                    } else if (tolower(Platform) %in% "h-mirna_8x15kv") {
                        tmp <- manifest$platform == "H-miRNA_8x15Kv"
                        mirna_files <- manifest[tmp, ]
                    }
                } else {
                    stop(message("No miRNA data was downloaded!"))
                }

            } else if (tolower(dataBase) == "gdc") {
                DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                paste0("gdc_", tolower(dataType), "_data"))

                manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"),
                                            select = c("file_name", "cases"))
                tmp <- !grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                            sep = "|"),
                            x = dir(path = DIR))

                mirna_files_downloaded <- dir(path = DIR)[tmp]
                tmp <- manifest$file_name %in% mirna_files_downloaded
                manifest <- manifest[tmp, ]
                mirna_files <- manifest
            }

            #checking for unmatched data
            if (nrow(mirna_files) > 0) {
                message("Reading data...")
            } else {
                message(paste0("Wrong Plataform '", tolower(Platform),
                            "' chosen!"))
                stop()
            }

            if (tumorData) {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-", tumor_regex,
                                                            "[aAbB]-"),
                                        mirna_files$cases)
                #not related the control samples
                mirna_tumor_files <- as.character(mirna_files[[1]][only_tumor_tissue])
                patients <- as.character(as.matrix(mirna_files[only_tumor_tissue, "cases"]))
                assign("patients", patients, envir = get(envir_link))
                OPEN_miRNA(mirna_tumor_files)
                #saving
                # write.table(x = data_mirna_table,
                    # file = "./mirnaylation_table/mirna_table_completed.csv",
                #           row.names = TRUE, quote = FALSE)

                assign("mirna_tumor_cross_mapped", string_vars[["envir_link"]]$cross_mapped_mirna_table, envir = get(envir_link))

                if (normalization) {
                    assign("mirna_tumor_normalized", string_vars[["envir_link"]]$data_mirna_table, envir = get(envir_link))
                    if (saveData) {
                        message("\nSaving your data...\n")
                        #saving
                        write.table(x = string_vars[["envir_link"]]$data_mirna_table,
                                    file = file.path(DIR,
                                                "mirna_tumor_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = string_vars[["envir_link"]]$cross_mapped_mirna_table,
                                    file = file.path(DIR,
                                                "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                    }
                } else {
                    assign("mirna_tumor_not_normalized", string_vars[["envir_link"]]$data_mirna_table, envir = get(envir_link))
                    if (saveData) {
                        message("\nSaving your data...\n")
                        #saving
                        write.table(x = string_vars[["envir_link"]]$data_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_tumor_not_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = string_vars[["envir_link"]]$cross_mapped_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                    }
                }
            } else {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                not_tumor_regex <- paste(10:19, collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-",
                                                            tumor_regex, "-"),
                                        mirna_files$cases)
                not_tumor_tissue <- grepl(pattern = paste0("-",
                                                        not_tumor_regex, "-"),
                                        mirna_files$cases)

                # tumor
                mirna_tumor_files <- as.character(mirna_files[[1]][only_tumor_tissue])
                patients <- as.character(as.matrix(mirna_files[only_tumor_tissue, "cases"]))
                OPEN_miRNA(mirna_tumor_files)
                tumor_data_mirna_table <- string_vars[["envir_link"]]$data_mirna_table
                tumor_cross_mapped_mirna_table <- string_vars[["envir_link"]]$cross_mapped_mirna_table
                patients_tumor <- patients

                if (normalization) {
                    assign("mirna_tumor_normalized", tumor_data_mirna_table,
                        envir = get(envir_link))
                } else {
                    assign("mirna_tumor_not_normalized", tumor_data_mirna_table,
                        envir = get(envir_link))
                }

                suppressWarnings(remove(data_mirna_table,
                                        cross_mapped_mirna_table,
                                        envir = get(envir_link)))

                # not tumor
                mirna_not_tumor_files <- as.character(mirna_files[[1]][not_tumor_tissue])
                if (length(mirna_not_tumor_files) == 0) {
                    message("There is no 'Normal' data in ", tumor,
                            " tumor folder! Please, rerun after",
                            " set 'tumorData = TRUE'.",
                            "\n\n")
                    stop()
                }
                patients <- as.character(as.matrix(mirna_files[not_tumor_tissue, "cases"]))
                OPEN_miRNA(mirna_not_tumor_files)
                not_tumor_data_mirna_table <- string_vars[["envir_link"]]$data_mirna_table
                not_tumor_cross_mapped_mirna_table <- string_vars[["envir_link"]]$cross_mapped_mirna_table
                patients_not_tumor <- patients

                if (normalization) {
                    assign("mirna_not_tumor_normalized",
                        not_tumor_data_mirna_table, envir = get(envir_link))
                } else {
                    assign("mirna_not_tumor_not_normalized",
                        not_tumor_data_mirna_table, envir = get(envir_link))
                }

                if (normalization) {
                    if (saveData) {
                        message("\nSaving your data...\n")
                        write.table(x = tumor_data_mirna_table,
                                    file = file.path(DIR,
                                                "mirna_tumor_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = tumor_cross_mapped_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = not_tumor_data_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_not_tumor_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = not_tumor_cross_mapped_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_not_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                    }
                } else {
                    if (saveData) {
                        message("\nSaving your data...\n")
                        write.table(x = tumor_data_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_tumor_not_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = tumor_cross_mapped_mirna_table,
                                    file = file.path(DIR,
                                                "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = not_tumor_data_mirna_table,
                                    file = file.path(DIR,
                                        "mirna_not_tumor_not_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                        write.table(x = not_tumor_cross_mapped_mirna_table,
                                    file = file.path(DIR,
                                            "mirna_not_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE,
                                                                    sep = "\t")
                    }
                }
            }

            if (!missing(Name)) {
                if (tumorData) {
                    if (normalization) {
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA_selector_final <- grepl(paste0("^",
                                                            paste0("hsa-",
                                                                tolower(Name)),
                                                                "$"),
                                                    miRNA_selector,
                                                    perl = TRUE)


                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])

                        if (sum(miRNA_selector_final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        colnames(cross_mapped_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        Name <- gsub("-", "_", paste0("hsa-", tolower(Name)))
                        assign(paste0("mirna_tumor_normalized_selected_",
                                        toupper(Name)), miRNA_exp_selected,
                                envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_",
                                        toupper(Name)), cross_mapped_selected,
                                envir = get(envir_link))
                    } else {
                        message("Only normalized data can be",
                                                " used to separate patients!")
                        stop()
                    }
                } else {
                    # tumor
                    if (normalization) {
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA_selector_final <- grepl(paste0("^",
                                                        paste0("hsa-",
                                                                tolower(Name)),
                                                                "$"),
                                                    miRNA_selector,
                                                    perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])

                        if (sum(miRNA_selector_final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        colnames(cross_mapped_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        Name <- gsub("-", "_", paste0("hsa-", tolower(Name)))
                        assign(paste0("mirna_tumor_normalized_selected_",
                                                    toupper(Name)),
                                                    miRNA_exp_selected,
                                                    envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_",
                                                    toupper(Name)),
                                                    cross_mapped_selected,
                                                    envir = get(envir_link))
                    } else {
                        message("Only normalized data can be",
                                                " used to separate patients!")
                        stop()
                    }

                    # not tumor
                    if (normalization) {
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_not_tumor_normalized)
                        miRNA_selector_final <- grepl(paste0("^",
                                                            paste0("hsa-",
                                                                tolower(Name)),
                                                                "$"),
                                                    miRNA_selector,
                                                    perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA_selector_final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA_selector_final, , drop = FALSE])

                        colnames(miRNA_exp_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        colnames(cross_mapped_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        Name <- gsub("-", "_", paste0("hsa-", tolower(Name)))
                        assign(paste0("mirna_not_tumor_normalized_selected_",
                                    toupper(Name)), miRNA_exp_selected,
                            envir = get(envir_link))
                        assign(paste0("mirna_not_tumor_cross_mapped_selected_",
                                    toupper(Name)), cross_mapped_selected,
                            envir = get(envir_link))
                    } else {
                        message("Only normalized data can be",
                                                " used to separate patients!")
                        stop()
                    }
                }
            }
        } else {
            if (!missing(Name)) {
                if (tumorData) {
                    if (normalization) {
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA_selector_final <- grepl(paste0("^",
                                                            paste0("hsa-",
                                                                tolower(Name)),
                                                                "$"),
                                                    miRNA_selector,
                                                    perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])

                        if (sum(miRNA_selector_final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- paste0("hsa-",
                                                            tolower(Name))
                        colnames(cross_mapped_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        Name <- gsub("-", "_", paste0("hsa-",
                                                                tolower(Name)))
                        assign(paste0("mirna_tumor_normalized_selected_",
                                    toupper(Name)), miRNA_exp_selected,
                            envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_",
                                    toupper(Name)), cross_mapped_selected,
                            envir = get(envir_link))
                    } else {
                        message("Only normalized data can be",
                                            " used to separate patients!\n")
                        stop()
                    }
                } else {
                    # tumor
                    if (normalization) {
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA_selector_final <- grepl(paste0("^",
                                                            paste0("hsa-",
                                                                tolower(Name)),
                                                                "$"),
                                                    miRNA_selector,
                                                    perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA_selector_final, , drop = FALSE])

                        if (sum(miRNA_selector_final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        colnames(cross_mapped_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        Name <- gsub("-", "_", paste0("hsa-", tolower(Name)))
                        assign(paste0("mirna_tumor_normalized_selected_",
                                    toupper(Name)), miRNA_exp_selected,
                            envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_",
                                    toupper(Name)), cross_mapped_selected,
                            envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to",
                                " separate patients!\n")
                        stop()
                    }

                    # not tumor
                    if (normalization) {
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_not_tumor_normalized)
                        miRNA_selector_final <- grepl(paste0("^",
                                                            paste0("hsa-",
                                                                tolower(Name)),
                                                                "$"),
                                                    miRNA_selector,
                                                    perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA_selector_final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA_selector_final, , drop = FALSE])

                        colnames(miRNA_exp_selected) <- paste0("hsa-",
                                                            tolower(Name))
                        colnames(cross_mapped_selected) <- paste0("hsa-",
                                                                tolower(Name))
                        Name <- gsub("-", "_", paste0("hsa-", tolower(Name)))
                        assign(paste0("mirna_not_tumor_normalized_selected_",
                                    toupper(Name)), miRNA_exp_selected,
                            envir = get(envir_link))
                        assign(paste0("mirna_not_tumor_cross_mapped_selected_",
                                    toupper(Name)), cross_mapped_selected,
                            envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to ",
                                "separate patients!\n")
                        stop()
                    }
                }
            }
        }
        suppressWarnings(remove(data_mirna_table,
                                cross_mapped_mirna_table,
                                envir = get(envir_link)))
    }

    # Exon quantification ####
    if ("exon" == strsplit(tolower(dataType), split = " ")[[1]][1]) {
        if (!onlyFilter) {
            if (tolower(dataBase) == "gdc") {
                stop(message("\nThrere is no Exon quantification",
                            " data in GDC data base!!",
                            "\nPlease use 'Legacy' data base\n"))
            }
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                            "exon quantification_data")


            if (!dir.exists(file.path(workDir, "DOAGDC", toupper(tumor),
                                    "mage_data"))) {
                message("Downloading magetab data...")
                suppressWarnings(download_gdc(dataType = "mage",
                                            dataBase = dataBase,
                                            tumor = tumor,
                                            workDir = workDir))
            }

            desing_array_DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                                        "mage_data")

            if (tolower(Platform) == "") {
                message("Please, insert 'illumina hiseq' or 'illumina ga' ",
                        "for Exon quantification data.")
                stop()
            } else if (tolower(Platform) == "illumina ga") {
                to_gunzip <- file.path(desing_array_DIR,
                                    paste0("unc.edu_", toupper(tumor),
                                            ".IlluminaHiSeq_RNASeqV2.",
                                            "mage-tab.1.1.0.tar.gz"))
                ## or, if you just want to extract the target file:
                untar(to_gunzip, exdir = desing_array_DIR)

                design_array <- data.table::fread(file.path(desing_array_DIR,
                                                        paste0("unc.edu_",
                                                            toupper(tumor),
                                                            ".IlluminaHiSeq_",
                                                            "RNASeqV2.mage-",
                                                            "tab.1.1.0"),
                                                        paste0("unc.edu_",
                                                            toupper(tumor),
                                                            ".IlluminaHiSeq_",
                                                            "RNASeqV2.1.1.",
                                                            "0.sdrf.txt")),
                                                select = c(2, 22))
            } else if (tolower(Platform) == "illumina hiseq") {
                to_gunzip <- file.path(desing_array_DIR,
                                    paste0("unc.edu_", toupper(tumor),
                                            ".IlluminaHiSeq_RNASeqV2",
                                            ".mage-tab.1.1.0.tar.gz"))
                ## or, if you just want to extract the target file:
                untar(to_gunzip, exdir = desing_array_DIR)

                design_array <- data.table::fread(file.path(desing_array_DIR,
                                                        paste0("unc.edu_",
                                                            toupper(tumor),
                                                            ".IlluminaHiSeq_",
                                                            "RNASeqV2.mage-",
                                                            "tab.1.1.0"),
                                                        paste0("unc.edu_",
                                                            toupper(tumor),
                                                            ".IlluminaHiSeq_",
                                                            "RNASeqV2.1.",
                                                            "1.0.sdrf.txt")),
                                                select = c(2, 22))
            }

            exon_files <- dir(DIR, pattern = "unc.edu")
            tmp <- unlist(design_array[design_array[[2]] %in% exon_files, 1])
            patient_code <- as.character(tmp)

            if (normalization) {
                select_this_column <- 4 #RPKM
            } else {
                select_this_column <- 2 #raw counts
            }

            pb <- txtProgressBar(min = 0, max = length(exon_files), style = 3)
            count <- 0
            for (i in file.path(DIR, exon_files)) {
                count <- count + 1
                if (count == 1) {
                    exon <- data.table::fread(i)
                    exon_table <- matrix(ncol = length(exon_files),
                                            nrow = nrow(exon))
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

            if (tumorData) {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-", tumor_regex,
                                                            "[aAbB]-"),
                                        patient_code)
                #not related the control samples
                patients_tumor <- patient_code[only_tumor_tissue]
                assign("patients", patients_tumor, envir = get(envir_link))
                assign("exon_tumor_", exon_table, envir = get(envir_link))

                if (saveData) {
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = exon_table,
                                file = file.path(DIR, "tumor_exon_table.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            } else {
                #separing files
                tumor_regex <- paste(formatC(seq(1:9), width = 2, flag = "0"),
                                    collapse = "[aAbB]-|-")
                not_tumor_regex <- paste(10:19, collapse = "[aAbB]-|-")
                only_tumor_tissue <- grepl(pattern = paste0("-", tumor_regex,
                                                            "[aAbB]-"),
                                        patient_code)
                not_tumor_tissue <- grepl(pattern = paste0("-", not_tumor_regex,
                                                        "[aAbB]-"),
                                        patient_code)

                tumor_patients <- as.character(colnames(exon_table)[only_tumor_tissue])
                tumor_protein_table <- exon_table[, tumor_patients]
                assign("tumor_patients", tumor_patients,
                    envir = get(envir_link))
                assign("exon_tumor_", tumor_protein_table,
                    envir = get(envir_link))

                if (sum(not_tumor_tissue) == 0) {
                    message("There is no 'Normal' data in ", tumor,
                            " tumor folder! Please, rerun after",
                            " set 'tumorData = TRUE'.\n\n")
                    stop()
                }
                not_tumor_patients <- as.character(colnames(exon_table)[not_tumor_tissue])
                exon_not_tumor <- exon_table[, not_tumor_patients]
                assign("not_tumor_patients", not_tumor_patients,
                    envir = get(envir_link))
                assign("exon_not_tumor_", exon_not_tumor,
                    envir = get(envir_link))

                if (saveData) {
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = tumor_protein_table,
                                file = file.path(DIR,
                                                "tumor_protein_table.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                    write.table(x = exon_not_tumor,
                                file = file.path(DIR, "exon_not_tumor.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            }
        } else {
            if (normalization) {
                protein_selector <- rownames(string_vars[["envir_link"]]$tumor_protein_table)
            } else {
                stop(message("Only normalized data can ",
                            "be used to separate patients!"))
            }
            protein_selector_final <- grepl(Name,
                                            protein_selector, perl = TRUE)

            protein_exp_selected <- t(string_vars[["envir_link"]]$tumor_protein_table[protein_selector_final, ])
            colnames(protein_exp_selected) <- Name
            colnames(cross_mapped_selected) <- Name
            assign(paste0("exon_exp_selected_", Name), protein_exp_selected,
                envir = get(envir_link))
        }
    }

    # mage ####
    if ("mage" == tolower(dataType)) {
        if (tolower(dataBase) == "legacy") {
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor), "mage_data")
        } else if (tolower(dataBase) == "gdc") {
            DIR <- file.path(workDir, "DOAGDC", toupper(tumor),
                            "gdc_protein_data")
        }
        patient_file <- dir(path = DIR, pattern = "mage-tab")
        count <- 0

    }

    # end ####
    assign("tumor", tumor, envir = get(envir_link))
    assign("dataBase", dataBase, envir = get(envir_link))
    assign("dataType", dataType, envir = get(envir_link))
    assign("workDir", workDir, envir = get(envir_link))
    if (!missing("Name")) {
        assign("Name", Name, envir = get(envir_link))
        assign("Name.e", Name, envir = get(envir_link))
    }
    message("\nDone!")
}
