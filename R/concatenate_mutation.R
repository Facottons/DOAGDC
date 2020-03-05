#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_mutation} is a function designed to concatenate GDC files into
#' a single matrix, where the columns stand for patients code and rows stand for
#' data names.
#'
#' @param name A character string indicating the desired values to be used in
#'    next analysis. For instance, "HIF3A" in the legacy gene expression matrix,
#'    "mir-1307" in the miRNA quantification matrix, or "HER2" in the protein
#'    quantification matrix.
#' @param data_base
#' @param work_dir
#' @param tumor
#' @param workflow_type A character string specifying the workflow type for
#'    mutation data in "gdc". Where: \itemize{\item{"varscan" -}{
#'    VarScan2 Variant Aggregation and Masking} \item{"mutect" -}{
#'    MuTect2 Variant Aggregation and Masking}\item{"muse" -}{ MuSE
#'    Variant Aggregation and Masking}\item{"somaticsniper" -}{
#'    SomaticSniper Variant Aggregation and Masking}\item{"all" means to}{
#'    concatenate all workflows into a single matrix.}}
#' @param tumor_data Logical value where \code{TRUE} specifies the desire to
#'    work
#'    with tumor tissue files only. When set to FALSE, it creates two matrices,
#'    one containing tumor data and other containing data from not-tumor tissue.
#'    The default is \code{TRUE}.
#' @param platform
#' @param env A character string containing the environment name that should be
#'    used. If none has been set yet, the function will create one in global
#'    environment following the standard criteria:
#'    \itemize{\item{'<tumor>_<data_base>_mutation_tumor_data'}{ or}
#'    \item{'<tumor>_<data_base>_mutation_both_data'}{ (for tumor and not tumor
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
#'
#' @examples
#' library(DOAGDC)
#'
#' # Concatenating gene expression data into a single matrix
#' # data already downloaded using the 'download_gdc' function
#' concatenate_mutation(
#'     name = "HIF3A",
#'     workflow_type = "varscan",
#'     platform = "Illumina GA",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
concatenate_mutation <- function(name, data_base,
                                 work_dir, tumor,
                                 workflow_type,
                                 tumor_data = TRUE,
                                 platform = "",
                                 env,
                                 save_data = FALSE) {


    # TODO: create matrix for mutation data from normal tissue

    # create env if not created yet
    if (missing(env) && tumor_data) {
        assign(paste(toupper(tumor), toupper(data_base),
            "mutation",
            "tumor_data",
            sep = "_"
        ),
        new.env(parent = emptyenv()),
        envir = .GlobalEnv
        )
        ev <- paste(toupper(tumor), toupper(data_base),
            "mutation",
            "tumor_data",
            sep = "_"
        )
    } else if (missing(env) && !tumor_data) {
        assign(paste(toupper(tumor), toupper(data_base),
            "mutation",
            "both_data",
            sep = "_"
        ),
        new.env(parent = emptyenv()),
        envir = .GlobalEnv
        )
        ev <- paste(toupper(tumor), toupper(data_base),
            "mutation",
            "both_data",
            sep = "_"
        )
    } else {
        ev <- deparse(substitute(env))
    }

    attr(ev, "name") <- paste0(
        "Environment created! Use this name in 'env' argument"
    )

    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor), "Analyses"),
        showWarnings = FALSE
    )

    if (tolower(data_base) == "legacy") {
        # selecting platform file
        dir <- file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "mutation_data"
        )
        manifest <- data.table::fread(file.path(dir, "manifest.sdrf"),
            select = c("file_name", "platform")
        )

        test <- tolower(platform)
        tmp <- manifest[[1]][tolower(manifest$platform) == test]
        maf_files <- as.character(tmp)

        # checking for unmatched data
        if (length(maf_files) > 0) {
            message("Loading mutation file...")
        } else {
            stop(message(paste0(
                "There is no data from ",
                tolower(platform),
                " in your local storage!"
            )))
        }

        message("\nIt will be open ", length(maf_files), " file(s)...\n")
        pb <- txtProgressBar(min = 0, max = length(maf_files), style = 3)
        count <- 0
        for (file in maf_files) {
            count <- count + 1
            assign(file, read.delim(
                file = file.path(dir, file),
                comment.char = "#", sep = "\t",
                header = TRUE,
                stringsAsFactors = FALSE
            ),
            envir = get(ev)
            )
            message(
                "\nThe file is load as '", file, "' in ",
                ev, " Environment."
            )

            if (save_data) {
                write.table(
                    x = get(file, envir = get(ev)),
                    file = file.path(dir, paste0(file, ".tsv")),
                    row.names = FALSE, quote = FALSE, sep = "\t"
                )
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)
    } else if (tolower(data_base) == "gdc") {
        dir <- file.path(
            work_dir, "DOAGDC", toupper(tumor),
            "gdc_mutation_data"
        )
        manifest <- data.table::fread(file.path(dir, "manifest.sdrf"),
            select = c(11, 12)
        )
        manifest$submitter_id <- unlist(lapply(strsplit(
            x = manifest$submitter_id,
            split = "-",
            perl = TRUE
        ), "[[", 3))
        if (tolower(workflow_type) %in% "all") {
            to_gunzip <- dir(path = dir, pattern = "maf.gz$")
        } else {
            test <- tolower(workflow_type)
            tmp <- manifest[[1]][manifest$submitter_id == test]
            to_gunzip <- as.character(tmp)
        }

        tryCatch(R.utils::gunzip(file.path(dir, to_gunzip),
            remove = FALSE
        ),
        error = function(e) {

        }
        )
        import_maf <- gsub(".gz", "", to_gunzip)
        # openning files
        message("\nLoading mutation file...\n")
        maf <- read.delim(
            file = paste(dir, import_maf, sep = "/"),
            comment.char = "#", sep = "\t",
            header = TRUE,
            stringsAsFactors = FALSE
        ) # , quote = FALSE)

        assign(import_maf, maf, envir = get(ev))
        message(
            "\nThe file is open as '", import_maf, "' in ",
            ev, " environment."
        )
        if (save_data) {
            write.table(
                x = get(import_maf, envir = get(ev)),
                file = file.path(dir, paste0(import_maf, ".tsv")),
                row.names = FALSE, quote = FALSE, sep = "\t"
            )
        }
    }
}