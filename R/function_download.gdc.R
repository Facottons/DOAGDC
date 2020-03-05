#' Download data from GDC Data Portal and GDC Legacy Archive
#'
#' \code{download_gdc} is a function designed to download methylation, mutation,
#' clinical data, protein expression, MAGETAB, gene expression, isoform
#' expression, miRNA expression and clinical images data from GDC Data Portal
#' and GDC Legacy Archive.
#'
#' @param data_type Type of data. It could be \code{"methylation", "mutation",
#'    "clinical_supplement", "biospecimen", "gene", or "clinical"(biotab)}.
#'    \itemize{ \item{Only present
#'    in "Legacy" database:}{\code{"protein", "Exon quantification", "miRNA gene
#'    quantification", "miRNA isoform quantification", "isoform", and "image"}.}
#'    \item{Only present in "GDC" database:}{\code{"miRNA
#'    Expression Quantification", and "Isoform Expression Quantification"
#'    (miRNA)}.}}
#' @param tumor A character string contaning one of the 33 tumors available in
#'    the TCGA project. For instance, the \code{"BRCA"} stands for breast
#'    cancer.
#' @param data_base A character string specifying \code{"GDC"} for GDC Data
#'    Portal or \code{"legacy"} for GDC Legacy Archive.
#' @param htseq A character string indicating which HTSeq workflow data should
#'    be downloaded: \code{"Counts", "FPKM", or "all"}. The default is
#'    \code{"all"}.
#' @param work_dir A character string specifying the path to work directory.
#' @param all_files A logical value. Set \code{FALSE} to avoid the download of
#'    not used data to reduce download size, e.g. quantification files. The
#'    default is \code{FALSE}.
#' @param platform A character string indicating the platform name for
#'    methylation, exon quantificaton, miRNA, and mutation data. \itemize{
#'    \item{For mutation and exon quantificaton data:}{\code{"Illumina GA",
#'    "Illumina HiSeq" or "all"}.} \item{For methylation data}{\code{"Illumina
#'    Human Methylation 450", "Illumina Human Methylation 27" or "all"}.}
#'    \item{For miRNA data:}{\code{"Illumina GA", "Illumina HiSeq",
#'    "H-miRNA_8x15K" (for GBM tumor), "H-miRNA_8x15Kv2" (for OV tumor), or
#'    "all"}.} }The default for all data_type cited is \code{"all"} (when
#'    downloading data).
#'
#' @return the files download are stored inside the determined folders in the
#'    user machine.
#'
#' @import AnnotationDbi clusterProfiler devtools DOSE ggbiplot ggplot2 methods
#'    stringi survminer yarrr
#' @export
#'
#' @importFrom curl curl
#' @importFrom httr content
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' @importFrom tools md5sum
#' @importFrom grDevices dev.off hsv png svg
#' @importFrom graphics abline axis hist image layout legend lines matplot
#'    mtext par plot.new points rect rug text title
#' @importFrom stats TukeyHSD anova aov as.dendrogram as.dist chisq.test coef
#'    confint cor density dist formula hclust kruskal.test median model.matrix
#'    na.exclude na.omit order.dendrogram p.adjust pairwise.wilcox.test prcomp
#'    reorder residuals sd shapiro.test summary.aov summary.lm
#' @importFrom utils combn read.csv read.delim read.table setTxtProgressBar
#'   txtProgressBar untar write.csv write.table
#'
#' @examples
#' library(DOAGDC)
#'
#' # Downloading gene expression data from GDC Legacy Archive
#' download_gdc("gene", "CHOL", "legacy", work_dir = "~/Desktop")
download_gdc <- function(data_type = "gene",
                        tumor,
                        data_base = "legacy",
                        htseq = "",
                        work_dir,
                        all_files = FALSE,
                        platform = "all") {

    # local functions ####
    download_httr <- function(url, destfile) {
        first <- httr::GET(url = url)
        second <- httr::content(x = first, as = "raw")
        writeBin(object = second, con = destfile)
    }

    size_par <- function(tumor, type_of_data, db) {
        if (db == "legacy") {
            first_part <- "https://api.gdc.cancer.gov/legacy/projects/TCGA-"
        } else {
            first_part <- "https://api.gdc.cancer.gov/projects/TCGA-"
        }
        url <- paste0(
            first_part, toupper(tumor),
            "?expand=summary,summary.data_categories&pretty=true"
        )
        jason <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
        jason <- jason$data
        jason <- jason$summary
        jason <- jason$data_categories
        jason$data_category <- tolower(jason$data_category)
        size <- as.numeric(subset(
            x = jason,
            subset = data_category == tolower(type_of_data),
            "file_count"
        ))
        return(size)
    }

    # old api url https://gdc-api.nci.nih.gov/
    # selecting the right API

    # code ####
    if (tolower(data_base) == "legacy") {
        inicio <- "https://api.gdc.cancer.gov/legacy/data/"
        url_inicio <- "https://api.gdc.cancer.gov/legacy/files/"
        status <- "https://api.gdc.cancer.gov/legacy/status"
        folder_name <- paste0(tolower(data_type), "_data")
    } else if (tolower(data_base) == "gdc") {
        inicio <- "https://api.gdc.cancer.gov/data/"
        url_inicio <- "https://api.gdc.cancer.gov/files/"
        status <- "https://api.gdc.cancer.gov/status"
        folder_name <- paste(tolower(data_base), tolower(data_type),
            "data",
            sep = "_"
        )
    } else {
        stop("Please insert a data base name!")
    }

    message("Please wait, accessing GDC server...")
    tryCatch(tmp <- read.csv(status),
        error = function(e) {
            stop(message(cat(
                "GDC server or your internet conection is",
                " off. \n Please try again later!"
            )))
        }
    )

    dir.create(path = file.path(work_dir, "DOAGDC"), showWarnings = FALSE)
    dir.create(
        path = file.path(work_dir, "DOAGDC", toupper(tumor)),
        showWarnings = FALSE
    )

    # legacy ####
    if (tolower(data_base) == "legacy") {
        # NOTE gene and isoform ####
        if ("gene" %in% tolower(data_type) || "isoform" %in% tolower(data_type)) {
            size_par_rsem <- function(tumor) {
                url <- paste0(
                    "https://api.gdc.cancer.gov/legacy/",
                    "projects/TCGA-", toupper(tumor),
                    "?expand=summary,summary.data_categories&",
                    "pretty=true"
                )
                jason <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
                jason <- jason$data
                jason <- jason$summary
                jason <- jason$data_categories
                size <- subset(
                    x = jason,
                    subset = data_category == "Gene expression",
                    "file_count"
                )
                return(size)
            }

            size <- size_par_rsem(tumor = tumor)

            if ("gene" %in% tolower(data_type)) {
                dir.create(
                    path = file.path(
                        work_dir, "DOAGDC",
                        toupper(tumor), "gene_data"
                    ),
                    showWarnings = FALSE
                )
                dir <- file.path(work_dir, "DOAGDC", toupper(tumor),
                                                                "gene_data")

                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22:%22and%22,%22",
                    "content%22:%5B%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:%22",
                    "files.access%22,%22value%22:%5B",
                    "%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:",
                    "%22cases.project.project_id%22,%22value%22:",
                    "%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22files",
                    ".data_category%22,%22value%22:%5B%22Gene%20",
                    "expression%22%5D%7D%7D,%7B%",
                    "22op%22:%22in%22,%22content%22:%7B%22field",
                    "%22:%22files.data_type%22,%",
                    "22value%22:%5B%22Gene%20expression%20",
                    "quantification%22%5D%7D%7D%5D%7D&",
                    "pagination=%7B%22files%22:%7B%22from%22:0,%22",
                    "size%22:100,%22sort%22:",
                    "%22cases.project.project_id:asc%22%7D%7D&",
                    "format=JSON"
                )
            } else if ("isoform" %in% tolower(data_type)) {
                dir.create(
                    path = file.path(
                        work_dir, "DOAGDC", toupper(tumor),
                        "isoform_data"
                    ),
                    showWarnings = FALSE
                )
                dir <- file.path(
                    work_dir, "DOAGDC", toupper(tumor),
                    "isoform_data"
                )

                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22:%22and%22,%22",
                    "content%22:%5B%7B%22op%22:%22in",
                    "%22,%22content%22:%7B%22field%22:%22files.",
                    "access%22,%22value%22:%5B%22",
                    "open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:%22",
                    "cases.project.project_id%22,%22value%22:",
                    "%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22",
                    "files.data_category%22,%22value%22:%5B%22",
                    "Gene%20expression%22%5D%7D%",
                    "7D,%7B%22op%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22files.data_",
                    "type%22,%22value%22:%5B%22Isoform%20expression",
                    "%20quantification%22%5D",
                    "%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22",
                    "from%22:0,%22size%22:100,%",
                    "22sort%22:%22cases.project.project_id:asc",
                    "%22%7D%7D&format=JSON"
                )
            }

            message("\n\nDownloading manifest...\n")
            json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

            manifest_df <- json$data$hits

            manipular <- manifest_df[, "cases"]
            # 16(submitter_id) will be cases
            manifest_df[, c(
                "center", "acl", "state_comment", "cases",
                "tags"
            )] <- NULL

            cases <- matrix(nrow = nrow(manifest_df), ncol = 1)
            for (index in seq_len(length(manipular))) {
                patient_code <- manipular[[index]][2]
                tmp <- as.character(unlist(patient_code))
                cases[index, 1] <- tmp[grep("TCGA", tmp)]
            }

            manifest_df$submitter_id <- cases

            if (!all_files) {
                # no need to download these files
                manifest_df <- manifest_df[!grepl(
                    "^.+(.pergene.txt)$",
                    manifest_df$file_name
                ), ]
                manifest_df <- manifest_df[!grepl(
                    "^.+(.quantification.txt)$",
                    manifest_df$file_name
                ), ]
                manifest_df <- manifest_df[!grepl(
                    "^.+(level3.data.txt)$",
                    manifest_df$file_name
                ), ]
                manifest_df <- manifest_df[!grepl(
                    "^.+(tags.txt)$",
                    manifest_df$file_name
                ), ]
                manifest_df <- manifest_df[!grepl(
                    "^.+(genes.txt)$",
                    manifest_df$file_name
                ), ]
            }

            write.table(
                x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
                quote = FALSE, row.names = FALSE, sep = "\t"
            )

            manifest_df <- manifest_df[, c("file_name", "md5sum", "file_id")]

            colnames(manifest_df) <- c("filename", "md5", "id")

            id_matrix <- manifest_df[, "id"]
        }
    } else if (tolower(data_base) == "gdc") {
        # NOTE gene GDC ####
        if ("gene" %in% tolower(data_type)) {
            dir.create(
                path = file.path(
                    work_dir, "DOAGDC", toupper(tumor),
                    "gdc_gene_data"
                ),
                showWarnings = FALSE
            )
            dir <- file.path(work_dir, "DOAGDC", toupper(tumor),
                                                            "gdc_gene_data")

            size_par <- function(tumor) {
                url <- paste0(
                    "https://api.gdc.cancer.gov/projects/TCGA-",
                    toupper(tumor),
                    "?expand=summary,summary.",
                    "data_categories&pretty=true"
                )
                jason <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
                jason <- jason$data
                jason <- jason$summary
                jason <- jason$data_categories
                tmp <- jason$data_category == "Transcriptome Profiling"
                size <- subset(x = jason, subset = tmp, "file_count")
                return(as.numeric(size))
            }

            size <- size_par(tumor = tumor)
            # selecting which HTSeq data to download
            if (tolower(htseq) == "counts") {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22:%22and%22,%22",
                    "content%22:%5B%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:%22",
                    "files.data_type%22,%22value%22:%5B%22",
                    "Gene%20Expression%20Quantification%22%5D%7D%7D",
                    ",%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22cases.project.project_id",
                    "%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                    "%22%",
                    "5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22files.access%22,",
                    "%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22",
                    ":%22in%22,%22content%22:%7B%22field",
                    "%22:%22files.data_format%22,%22value%22:%5B",
                    "%22TXT%22%5D%7D%7D,%7B%22op%22:",
                    "%22in%22,%22content%22:%7B%22field%22:%22files",
                    ".experimental_strategy%22,%22value%22",
                    ":%5B%22RNA-Seq%22%5D%7D%7D,%7B%22op%22:%22in",
                    "%22,%22content%22:%7B%22field%22:",
                    "%22files.data_category%22,%22value%22:%5B%22",
                    "Transcriptome%20Profiling%22%5D%7D%7D,",
                    "%7B%22op%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22files.analysis.workflow_type%22",
                    ",%22value%22:%5B%22HTSeq%20-%20Counts%22%5D%",
                    "7D%7D%5D%7D&format=JSON"
                )
            } else if (toupper(htseq) == "FPKM") {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22:%22and%22,%22",
                    "content%22:%5B%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:%22files.",
                    "data_type%22,%22value%22:%5B%22",
                    "Gene%20Expression%20Quantification%22%5D%7D",
                    "%7D,%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22cases.project.project_",
                    "id%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                    "%22%",
                    "5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:",
                    "%7B%22field%22:%22files.access%22,",
                    "%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%",
                    "22:%22in%22,%22content%22:%7B%22field",
                    "%22:%22files.data_format%22,%22value%22:%5B%",
                    "22TXT%22%5D%7D%7D,%7B%22op%22:",
                    "%22in%22,%22content%22:%7B%22field%22:%22",
                    "files.experimental_strategy%22,%22value%22",
                    ":%5B%22RNA-Seq%22%5D%7D%7D,%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:",
                    "%22files.data_category%22,%22value%22:%5B%22",
                    "Transcriptome%20Profiling%22%5D%7D%7D,",
                    "%7B%22op%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22files.analysis.workflow_type%22",
                    ",%22value%22:%5B%22HTSeq%20-%20", toupper(htseq),
                    "%22%5D%7D%7D%5D%7D&format=JSON"
                )
            } else if (toupper(htseq) == "ALL") {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22:%22and%22,%22",
                    "content%22:%5B%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:%22files.",
                    "data_type%22,%22value%22:",
                    "%5B%22Gene%20Expression%20Quantification%22",
                    "%5D%7D%7D,%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:%22files.",
                    "access%22,%22value%22:%5B%22",
                    "open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:%22",
                    "files.data_format%22,%22value%22:%5B%22TXT%22%",
                    "5D%7D%7D,%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22:%22files.",
                    "analysis.workflow_type%22,%22",
                    "value%22:%5B%22HTSeq%20-%20Counts%22,",
                    "%22HTSeq%20-%20FPKM%22%5D%7D%7D,%7B%22op%22:%22in%22,",
                    "%22content%22:%7B%22field%22:%22",
                    "files.experimental_strategy%22,%22value%22:",
                    "%5B%22RNA-Seq%22%5D%7D%7D,%7B%22",
                    "op%22:%22in%22,%22content%22:%7B%22field%22:",
                    "%22cases.project.project_id%22,",
                    "%22value%22:%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,",
                    "%22content%22:%7B%22field%22:%22files.data_",
                    "category%22,%22value%22:%5B%22",
                    "Transcriptome%20Profiling%22%5D%7D%7D%5D%7D&",
                    "format=JSON"
                )
            } else {
                stop(message('\nPlease insert a HTSeq value! ("Counts", ",
                            "FPKM", or "all")\n'))
            }

            message("\n\nDownloading manifest...\n")
            sw <- function(x) {
                suppressWarnings(x)
            }
            sw(json <- jsonlite::fromJSON(readLines(curl::curl(url))))

            manifest_df <- json$data$hits

            manipular <- manifest_df[, "cases"]
            # 16(submitter_id) will be cases
            manifest_df[, c("analysis", "acl")] <- NULL

            cases <- matrix(nrow = nrow(manifest_df), ncol = 1)
            for (index in seq_len(length(manipular))) {
                patient_code <- manipular[[index]][2]
                tmp <- as.character(unlist(patient_code))
                cases[index, 1] <- tmp[grep("TCGA", tmp)]
            }

            manifest_df$cases <- cases

            manifest_df[, "cases"] <- as.character(manifest_df[, "cases"])

            write.table(
                x = manifest_df, file = paste0(
                    dir, "/", toupper(htseq),
                    "_manifest.sdrf"
                ),
                quote = FALSE, row.names = FALSE, sep = "\t"
            )

            manifest_df <- manifest_df[, c("file_name", "md5sum", "file_id")]

            colnames(manifest_df) <- c("filename", "md5", "id")

            id_matrix <- manifest_df[, "id"]
        } else if ("isoform" %in% tolower(data_type)) {
            stop(message(
                "\nThere is no isoform data in GDC data",
                " base! (unfortunately...)\n"
            ))
        }
    }

    # NOTE mutation ####
    if ("mutation" %in% tolower(data_type)) {
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        if (tolower(data_base) == "legacy") {
            platform <- paste0(
                "%22value%22:%5B%22Illumina%20GA%22,",
                "%22Illumina%20HiSeq%22%5D%7D%7D"
            )
            # legacy exclusive
            if (tolower(platform) %in% "all") {
                tmp <- "Illumina GA|Illumina HiSeq"
            } else if (tolower(platform) %in% "illumina ga") {
                tmp <- "Illumina GA"
            } else if (tolower(platform) %in% "illumina hiseq") {
                tmp <- "Illumina HiSeq"
            }

            url <- paste0(
                url_inicio, "?pretty=true&expand=cases.samples.",
                "portions.analytes.aliquots,cases.project,center,",
                "analysis&filters=%7B%22op%22",
                ":%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,",
                "%22content%22:%7B%22field%22",
                ":%22files.data_format%22,%22value%22:%5B%22MAF%22%",
                "5D%7D%7D,%7B%22op%22:%22in%",
                "22,%22content%22:%7B%22field%22:%22files.access%22",
                ",%22value%22:%5B%22open%22%",
                "5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%",
                "22field%22:%22files.experimental",
                "_strategy%22,%22value%22:%5B%22DNA-Seq%22%5D%7D%7D",
                ",%7B%22op%22:%22in%22,%22content",
                "%22:%7B%22field%22:%22files.platform%22,", platform,
                ",%7B%22op%22:%22in%22,%22",
                "content%22:%7B%22field%22:%22cases.project.",
                "project_id%22,%22value%22:%5B%22TCGA-",
                toupper(tumor), "%22%5D%7D%7D%5D%7D&format=JSON"
            )
        } else if (tolower(data_base) == "gdc") {
            url <- paste0(
                url_inicio, "?pretty=true&expand=cases.samples.",
                "portions.analytes.aliquots,cases.project,center,",
                "analysis&filters=%7B%",
                "22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22",
                "in%22,%22content%22:",
                "%7B%22field%22:%22files.data_format%22,%22value%22",
                ":%5B%22MAF%22%5D%7D%",
                "7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field",
                "%22:%22files.access%",
                "22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:",
                "%22in%22,%22content%",
                "22:%7B%22field%22:%22files.data_category%22,%22",
                "value%22:%5B%22Simple%20",
                "Nucleotide%20Variation%22%5D%7D%7D,%7B%22op%22:%22",
                "in%22,%22content%22:",
                "%7B%22field%22:%22files.data_type%22,%22value%22:",
                "%5B%22Masked%20Somatic",
                "%20Mutation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                "content%22:%7B%22field%22",
                ":%22files.experimental_strategy%22,%22value%22:",
                "%5B%22WXS%22%5D%7D%7D,%7B",
                "%22op%22:%22in%22,%22content%22:%7B%22field%22:",
                "%22cases.project.project_id",
                "%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                "%22%5D%7D%7D%5D%7D&facetTab=cases&format=JSON"
            )
        }


        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
            error = function(e) stop(e),
            finally = message(
                "The tumor ", toupper(tumor),
                " does not have ", tolower(data_type),
                " data!"
            )
        )

        manifest_df <- json$data$hits

        # checkin' if there are available data to download (e.g. LAML)
        if (length(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        # center only available at legacy db
        if (!is.null(manifest_df$center)) {
            center <- manifest_df$center
            manifest_df[, c("tags", "acl", "center")] <- NULL
            manifest_df <- cbind(center, manifest_df)
        } else {
            analysis <- manifest_df$analysis
            manifest_df[, c("analysis", "acl")] <- NULL
            manifest_df <- cbind(analysis, manifest_df)
        }

        manipular <- manifest_df[, "cases"]

        pre_cases <- manipular[[1]][2]

        patient_code <- matrix(
            data = as.character(unlist(pre_cases)),
            nrow = 10,
            ncol = length(as.character(unlist(pre_cases))) / 10
        )

        manifest_df[, "cases"] <- paste(patient_code[6, ], collapse = ", ")

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[
            grep(tmp, head(manifest_df)$platform),
            c("file_name", "md5sum", "file_id")
        ]

        colnames(manifest_df) <- c("filename", "md5", "id")

        id_matrix <- manifest_df[, "id"]
    }

    # NOTE methylation ####
    if ("methylation" %in% tolower(data_type)) {
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        platform <- paste0(
            "%22value%22:%5B%22Illumina%20Human%20Methylation",
            "%20450%22,%22Illumina%20Human%20Methylation%2027%22%5D%7D%7D"
        )

        if (tolower(platform) %in% "all") {
            tmp <- paste0(
                "Illumina Human Methylation 450|Illumina ",
                "Human Methylation 27"
            )
        } else if (tolower(platform) %in% "illumina human methylation 450") {
            # platform <- paste0("%22value%22:%5B%22Illumina%20Human%20",
            # "Methylation%20450%22%5D%7D%7D")
            tmp <- "Illumina Human Methylation 450"
        } else if (tolower(platform) %in% "illumina human methylation 27") {
            # platform <- paste0("%22value%22:%5B%22Illumina%20Human%20",
            # "Methylation%2027%22%5D%7D%7D")
            tmp <- "Illumina Human Methylation 27"
        }

        size <- size_par(
            type_of_data = "DNA Methylation", tumor = tumor,
            db = tolower(data_base)
        )

        if (tolower(data_base) == "legacy") {
            url <- paste0(
                url_inicio, "?pretty=true&expand=cases.samples.",
                "portions.analytes.aliquots,cases.project,center,",
                "analysis&size=", size,
                "&filters=%7B%22op%22:%22and%22,%22content%22:",
                "%5B%7B%22op%22:%22in%22,",
                "%22content%22:%7B%22field%22:%22files.data_format",
                "%22,%22value%22:%5B",
                "%22TXT%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                "content%22:%7B%22field%22",
                ":%22files.access%22,%22value%22:%5B%22open%22",
                "%5D%7D%7D,%7B%22op%22:%",
                "22in%22,%22content%22:%7B%22field%22:%22files",
                ".platform%22,", platform,
                ",%7B%22op%22:%22in%22,%22content%22:%7B%22",
                "field%22:%22files.data_category%22,%22value%22:",
                "%5B%22DNA%20methylation",
                "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:",
                "%7B%22field%22:%22cases",
                ".project.project_id%22,%22value%22:%5B%22TCGA-",
                toupper(tumor),
                "%22%5D%7D%7D%5D%7D&format=JSON"
            )
        } else if (tolower(data_base) == "gdc") {
            url <- paste0(
                url_inicio, "?pretty=true&expand=cases.samples.",
                "portions.analytes.aliquots,cases.project,center,",
                "analysis&size=",
                size, "&filters=%7B%22op%22:%22and%22,%22content%22",
                ":%5B%7B%22op%22:%22in",
                "%22,%22content%22:%7B%22field%22:%22files.data_",
                "category%22,%22value%22:",
                "%5B%22DNA%20Methylation%22%5D%7D%7D,%7B%22op%22:",
                "%22in%22,%22content%22:",
                "%7B%22field%22:%22files.platform%22,", platform,
                ",%7B%22op%22:%22in%22,%22content%22:%7B%",
                "22field%22:%22files.access%22,%22value%22:%5B%22",
                "open%22%5D%7D%7D,%7B%",
                "22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
                "files.data_format%22,%",
                "22value%22:%5B%22TXT%22%5D%7D%7D,%7B%22op%22:%22",
                "in%22,%22content%22:%7B",
                "%22field%22:%22cases.project.project_id%22,%22",
                "value%22:%5B%22TCGA-", toupper(tumor),
                "%22%5D%7D%7D%5D%7D&format=JSON"
            )
        }
        message("\n\nDownloading manifest...\n")

        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)$data$hits

        manifest_df <- json

        # checkin' if there are available data to download (e.g. LAML)
        if (length(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        # center only available at legacy db
        if (!is.null(manifest_df$center)) {
            center <- manifest_df$center
            manifest_df[, c("tags", "acl", "center")] <- NULL
            manifest_df <- cbind(center, manifest_df)
        } else {
            analysis <- manifest_df$analysis
            manifest_df[, c("analysis", "acl")] <- NULL
            manifest_df <- cbind(analysis, manifest_df)
        }

        manipular <- manifest_df[, "cases"]

        cases <- matrix(nrow = nrow(manifest_df), ncol = 1)

        for (index in seq_len(length(manipular))) {
            patient_code <- manipular[[index]][2]
            tmp <- as.character(unlist(patient_code))
            cases[index, 1] <- tmp[grep("TCGA", tmp)]
        }

        manifest_df$cases <- cases

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[
            grep(tmp, manifest_df$platform),
            c("file_name", "md5sum", "file_id")
        ]

        colnames(manifest_df) <- c("filename", "md5", "id")

        id_matrix <- manifest_df[, "id"]
    }

    # NOTE clinical and image ####
    tmp <- c("clinical", "biospecimen", "clinical_supplement", "image")
    if (tolower(data_type) %in% tmp) {
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        if (tolower(data_base) == "legacy") {
            size <- size_par(
                tumor = tumor, type_of_data = "Clinical",
                db = "legacy"
            )
            if ("image" == tolower(data_type)) {
                message(
                    "A lot of data are going to be downloaded,",
                    " please wait..."
                )
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=10000&filters=%7B%22op%22",
                    ":%22and%22,%22content%22:%5B%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22",
                    ":%22files.data_category%22,%22value%22:%5B%22",
                    "Clinical%22%5D%7D%7D,%7B%22op%22",
                    ":%22in%22,%22content%22:%7B%22field%22:%22",
                    "files.access%22,%22value%22:%5B%22",
                    "open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:%22files.",
                    "platform%22,%22value%22:%5B%22Clinical%22%",
                    "5D%7D%7D,%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:%22cases.project.",
                    "project_id%22,%22value%22:%5B%22TCGA-",
                    toupper(tumor), "%22%5D%7D%7D,%7B%22op%22:%22",
                    "in%22,%22content%22:%7B%22field%22",
                    ":%22files.data_format%22,%22value%22:%5B%22",
                    "SVS%22%5D%7D%7D%5D%7D&format=JSON"
                )
            } else if ("clinical_supplement" == tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size, "&filters=",
                    "%7B%22op%22:%22and%22,%22content%22:%5B%7B",
                    "%22op%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22files.access%22,%22value%22:%5B",
                    "%22open%22%5D%7D%7D,%7B%22op%22:%22in%",
                    "22,%22content%22:%7B%22field%22:%22cases.",
                    "project.program.name%22,%22value%22:%5B%",
                    "22TCGA%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:%22cases.pro",
                    "ject.project_id%22,%22value%22:%5B%22TCGA-",
                    toupper(tumor),
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,",
                    "%22content%22:%7B%22field%22:%22files.data_",
                    "format%22,%22value%22:%5B%22BCR%20XML",
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22files.data_type",
                    "%22,%22value%22:%5B%22Clinical%20Supplement",
                    "%22%5D%7D%7D%5D%7D&format=JSON"
                )
            } else if ("clinical" == tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size, "&filters=%7B%",
                    "22op%22:%22and%22,%22content%22:%5B%7B%22op",
                    "%22:%22in%22,%22content%22:%7B%22field%22",
                    ":%22files.access%22,%22value%22:%5B%22open",
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                    "%22:%7B%22field%22:%22cases.project.program",
                    ".name%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,",
                    "%7B%22op%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22cases.project.project_id%22,%22",
                    "value%22:%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%",
                    "22:%7B%22field%22:%22files.data_category%22,",
                    "%22value%22:%5B%22Clinical%22%5D%7D%7D,",
                    "%7B%22op%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22files.data_format%22,%22value%",
                    "22:%5B%22Biotab%22%5D%7D%7D%5D%7D&format=JSON"
                )
            } else if ("biospecimen" == tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size, "&filters=%7B%22",
                    "op%22:%22and%22,%22content%22:%5B%7B%22op%22:",
                    "%22in%22,%22content%22:%7B%22field%22:%22",
                    "files.access%22,%22value%22:%5B%22open%22%5D",
                    "%7D%7D,%7B%22op%22:%22in%22,%22content%22:",
                    "%7B%22field%22:%22cases.project.program.name",
                    "%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,%7B%",
                    "22op%22:%22in%22,%22content%22:%7B%22field%22",
                    ":%22cases.project.project_id%22,%22value%",
                    "22:%5B%22TCGA-", toupper(tumor), "%22%5D%7D%",
                    "7D,%7B%22op%22:%22in%22,%22content%22:%7B%",
                    "22field%22:%22files.data_format%22,%22value",
                    "%22:%5B%22BCR%20XML%22%5D%7D%7D,%7B%22op%22",
                    ":%22in%22,%22content%22:%7B%22field%22:%22",
                    "files.data_type%22,%22value%22:%5B%22Biospeci",
                    "men%20Supplement%22%5D%7D%7D%5D%7D&format=JSON"
                )
            }
        } else if (tolower(data_base) == "gdc") {
            size <- size_par(
                tumor = tumor, type_of_data = "Clinical",
                db = "gdc"
            )

            if ("clinical_supplement" == tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size,
                    "&filters=%7B%22op%",
                    "22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%",
                    "22%3A%22in%22%2C%22content%22%3A%7B%22field%22",
                    "%3A%22cases.project.project_id%22%2C%22value",
                    "%22%3A%5B%22TCGA-", toupper(tumor), "%22%5D%",
                    "7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content",
                    "%22%3A%7B%22field%22%3A%22files.access%22%2C%22",
                    "value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22",
                    "op%22%3A%22in%22%2C%22content%22%3A%7B%22field",
                    "%22%3A%22files.data_format%22%2C%22value%",
                    "22%3A%5B%22BCR%20XML%22%5D%7D%7D%2C%7B%22op%22%",
                    "3A%22in%22%2C%22content%22%3A%7B%22field%22",
                    "%3A%22files.data_type%22%2C%22value%22%3A%5B%22",
                    "Clinical%20Supplement%22%5D%7D%7D%5D%7D&",
                    "format=JSON"
                )
            } else if ("clinical" == tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size, "&filters=",
                    "%7B%22op%22%3A%22and%22%2C%22content%22",
                    "%3A%5B%7B%22op%22%3A%22in%22%2C%22content",
                    "%22%3A%7B%22field%22%3A%22cases.project.",
                    "program.name%22%2C%22value%22%3A%5B%22",
                    "TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22",
                    "%2C%22content%22%3A%7B%22field%22%3A%22",
                    "cases.project.project_id%22%2C%22value%22%3A%5B",
                    "%22TCGA-", toupper(tumor), "%22%5D%7D%7D%2C",
                    "%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B",
                    "%22field%22%3A%22files.access%22%2C%22value",
                    "%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22",
                    "%3A%22in%22%2C%22content%22%3A%7B%22field",
                    "%22%3A%22files.data_category%22%2C%22value",
                    "%22%3A%5B%22Clinical%22%5D%7D%7D%2C%7B%22",
                    "op%22%3A%22in%22%2C%22content%22%3A%7B%22",
                    "field%22%3A%22files.data_format%22%2C%22",
                    "value%22%3A%5B%22BCR%20Biotab%22%5D%7D",
                    "%7D%5D%7D&format=JSON"
                )
            } else if ("biospecimen" == tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.proje",
                    "ct,center,analysis&size=", size,
                    "&filters=%7B%22op%",
                    "22%3A%22and%22%2C%22content%22%3A%5B%7B%22op",
                    "%22%3A%22in%22%2C%22content%22%3A%7B%22field%",
                    "22%3A%22cases.project.project_id%22%2C%22value%",
                    "22%3A%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D%2C%7B%22op",
                    "%22%3A%22in%22%2C%22content%22%3A%7B%22field",
                    "%22%3A%22files.access%22%2C%22value%22%3A%5B%22",
                    "open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%",
                    "22content%22%3A%7B%22field%22%3A%22files.",
                    "data_format%22%2C%22value%22%3A%5B%22BCR%20XML",
                    "%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22",
                    "content%22%3A%7B%22field%22%3A%22files.data_",
                    "type%22%2C%22value%22%3A%5B%22Biospecimen%20",
                    "Supplement%22%5D%7D%7D%5D%7D&format=JSON"
                )
            }
        }

        message("\n\nDownloading manifest...\n")
        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        manifest_df <- json$data$hits

        # checkin' if there are available data to download (e.g. LAML)
        if (length(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        # center only available at legacy db
        if (!is.null(manifest_df$center)) {
            center <- manifest_df$center
            manifest_df[, c("tags", "acl", "center", "cases")] <- NULL
            manifest_df <- cbind(center, manifest_df)
        } else {
            manifest_df[, c("tags", "acl", "cases")] <- NULL
        }

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[, c("file_name", "md5sum", "file_id")]
        colnames(manifest_df) <- c("filename", "md5", "id")
        id_matrix <- manifest_df[, "id"]
    }

    # NOTE protein ####
    if ("protein" == tolower(data_type)) {
        if (tolower(data_base) == "gdc") {
            stop(message(
                "\nThrere is no protein expression data in",
                " GDC data base!!",
                "\nPlease use 'legacy' data base"
            ))
        }
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        size <- size_par(
            tumor = tumor, type_of_data = "Protein expression",
            db = "legacy"
        )

        url <- paste0(
            url_inicio, "?pretty=true&expand=cases.samples.",
            "portions.analytes.aliquots,cases.project,center,",
            "analysis&size=", size, "&filters=%7B%22op%22",
            ":%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22",
            "content%22:%7B%22field%22",
            ":%22files.data_category%22,%22value%22:%5B%22",
            "Protein%20expression%22%5D%7D%",
            "7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field",
            "%22:%22files.access%22,%22",
            "value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22",
            ",%22content%22:%7B%22",
            "field%22:%22cases.project.project_id%22,%22value%22",
            ":%5B%22TCGA-",
            toupper(tumor), "%22%5D%7D%7D%5D%7D&format=JSON"
        )

        message("\n\nDownloading manifest...\n")
        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        manifest_df <- json$data$hits

        # checkin' if there are available data to download (e.g. LAML)
        if (length(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        center <- manifest_df$center
        manifest_df[, c("center", "acl", "cases")] <- NULL

        manifest_df <- cbind(center, manifest_df)

        # for some reason sometimes- 16(submitter_id) will be cases
        patient_code <- manifest_df$submitter_id

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[, c("file_name", "md5sum", "file_id")]

        colnames(manifest_df) <- c("filename", "md5", "id")

        id_matrix <- manifest_df[, "id"]
    }

    # NOTE mage ####
    if ("mage" == tolower(data_type)) {
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        url <- paste0(
            url_inicio, "?pretty=true&expand=cases.samples.",
            "portions.analytes.aliquots,cases.project,center,",
            "analysis&size=1000000&filters=%7B%22op%22:",
            "%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22",
            "content%22:%7B%22field%22:",
            "%22files.data_format%22,%22value%22:%5B%22MAGETAB%",
            "22%5D%7D%7D,%7B%22op%22:%22",
            "in%22,%22content%22:%7B%22field%22:%22files.access",
            "%22,%22value%22:%5B%22open%",
            "22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B",
            "%22field%22:%22files.data_type",
            "%22,%22value%22:%5B%22TCGA%20DCC%20Archive%22%5D%",
            "7D%7D,%7B%22op%22:%22in%22,%22",
            "content%22:%7B%22field%22:%22files.data_category",
            "%22,%22value%22:%5B%22Archive%22",
            "%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from",
            "%22:0,%22size%22:100,%22sort",
            "%22:%22cases.project.project_id:asc%22%7D",
            "%7D&format=JSON"
        )

        message("\n\nDownloading manifest...\n")
        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        manifest_df <- json$data$hits

        # checkin' if there are available data to download (e.g. LAML)
        if (nrow(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        # 16(submitter_id) will be cases
        manifest_df[, c("acl", "cases")] <- NULL

        manifest_df <- manifest_df[grepl(toupper(tumor), manifest_df$file_name,
            ignore.case = FALSE, perl = TRUE
        ), ]

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[, c("file_name", "md5sum", "file_id")]

        colnames(manifest_df) <- c("filename", "md5", "id")

        id_matrix <- manifest_df[, "id"]
    }

    # NOTE mirna ####
    is_mirna <- "mirna" %in% strsplit(tolower(data_type), split = " ")[[1]][1]
    is_isoform <- "isoform expression quantification" %in% tolower(data_type)
    if (is_mirna || is_isoform) {
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        platform <- ""

        if (tolower(platform) %in% "all") {
            tmp <- "Illumina HiSeq|Illumina GA|H-miRNA_8x15Kv2|H-miRNA_8x15Kv"
        } else if (tolower(platform) %in% "illumina hiseq") {
            platform <- paste0(
                ",%7B%22op%22:%22in%22,%22content%22:%7B%22",
                "field%22:%22files.platform%22,%22",
                "value%22:%5B%22Illumina%20HiSeq%22%5D%7D%7D"
            )
            tmp <- "Illumina HiSeq"
        } else if (tolower(platform) %in% "illumina ga") {
            platform <- paste0(
                ",%7B%22op%22:%22in%22,%22content%22:%7B%22",
                "field%22:%22files.platform%22,%22value%22",
                ":%5B%22Illumina%20GA%22%5D%7D%7D"
            )
            tmp <- "Illumina GA"
        } else if (tolower(platform) %in% "h-mirna_8x15kv2") {
            platform <- paste0(
                ",%7B%22op%22:%22in%22,%22content%22:%7B%22",
                "field%22:%22files.platform%22,%22",
                "value%22:%5B%22H-miRNA_8x15Kv2%22%5D%7D%7D"
            )
            tmp <- "H-miRNA_8x15Kv2"
        } else if (tolower(platform) %in% "h-mirna_8x15kv") {
            platform <- paste0(
                ",%7B%22op%22:%22in%22,%22content%22:",
                "%7B%22field%22:%22files.platform%22,",
                "%22value%22:%5B%22H-miRNA_8x15Kv%22%5D%7D%7D"
            )
            tmp <- "H-miRNA_8x15Kv"
        }

        if (tolower(data_base) == "legacy") {
            size <- size_par(
                type_of_data = "Gene expression", tumor = tumor,
                db = tolower(data_base)
            )
            if ("mirna gene quantification" %in% tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size,
                    "&filters=%7B%22op%22:%22and%22,%22content%22:",
                    "%5B%7B%22op%22:%22in%22,%",
                    "22content%22:%7B%22field%22:%22files.data_type",
                    "%22,%22value%22:%5B%22miRNA%",
                    "20gene%20quantification%22%5D%7D%7D,%7B%22op",
                    "%22:%22in%22,%22content%22:%7B%22",
                    "field%22:%22files.access%22,%22value%22:%5B%",
                    "22open%22%5D%7D%7D,%7B%22op%22:",
                    "%22in%22,%22content%22:%7B%22field%22:%22",
                    "cases.project.project_id%22,%22",
                    "value%22:%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D", platform, "%5D%7D&format=JSON"
                )
                print(url)
            } else if ("mirna isoform quantification" %in% tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=", size,
                    "&filters=%7B%22op%22:%22and%22,%22content%22",
                    ":%5B%7B%22op%22:%22in%22,%22",
                    "content%22:%7B%22field%22:%22files.data_type",
                    "%22,%22value%22:%5B%22miRNA%",
                    "20isoform%20quantification%22%5D%7D%7D,%7B%22",
                    "op%22:%22in%22,%22content%22",
                    ":%7B%22field%22:%22files.access%22,%22value",
                    "%22:%5B%22open%22%5D%7D%7D,%7B",
                    "%22op%22:%22in%22,%22content%22:%7B%22field",
                    "%22:%22cases.project.project_id",
                    "%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D", platform, "%5D%7D&format=JSON"
                )
            }
        } else if (tolower(data_base) == "gdc") {
            size <- size_par(
                type_of_data = "Transcriptome Profiling",
                tumor = tumor, db = tolower(data_base)
            )
            if ("mirna expression quantification" %in% tolower(data_type)) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22%3A%22and%22%2C%",
                    "22content%22%3A%5B%7B",
                    "%22op%22%3A%22in%22%2C%22content%22%3A%7B%",
                    "22field%22%3A%22cases.",
                    "project.project_id%22%2C%22value%22%3A%",
                    "5B%22TCGA-", toupper(tumor),
                    "%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%",
                    "2C%22content%22%3A%7B%22",
                    "field%22%3A%22files.access%22%2C%22value",
                    "%22%3A%5B%22open%22%5D%7D%",
                    "7D%2C%7B%22op%22%3A%22in%22%2C%22content%",
                    "22%3A%7B%22field%22%3A%22",
                    "files.data_type%22%2C%22value%22%3A%5B%22",
                    "miRNA%20Expression%20Quan",
                    "tification%22%5D%7D%7D%5D%7D&format=JSON"
                )
            } else if (is_isoform) {
                url <- paste0(
                    url_inicio, "?pretty=true&expand=cases.samples.",
                    "portions.analytes.aliquots,cases.project,",
                    "center,analysis&size=",
                    size, "&filters=%7B%22op%22%3A%22and%22%2C%22",
                    "content%22%3A%5B%7B%22",
                    "op%22%3A%22in%22%2C%22content%22%3A%7B%22",
                    "field%22%3A%22cases.project",
                    ".project_id%22%2C%22value%22%3A%5B%22TCGA-",
                    toupper(tumor),
                    "%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22",
                    "content%22%3A%7B%22field",
                    "%22%3A%22files.access%22%2C%22value%22%3A%5B",
                    "%22open%22%5D%7D%7D%2C%",
                    "7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%",
                    "22field%22%3A%22files.",
                    "data_type%22%2C%22value%22%3A%5B%22Isoform%",
                    "20Expression%20Quantifi",
                    "cation%22%5D%7D%7D%5D%7D&format=JSON"
                )
            }
        }
        message("\n\nDownloading manifest...\n")

        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
            error = function(e) stop(e),
            finally = message(
                "The tumor ", toupper(tumor),
                " does not have ",
                tolower(data_type), "data!"
            )
        )

        manifest_df <- json$data$hits

        # checkin' if there are available data to download (e.g. LAML)
        if (length(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        # center only available at legacy db
        if (!is.null(manifest_df$center)) {
            center <- manifest_df$center
            manifest_df[, c("tags", "acl", "center")] <- NULL
            manifest_df <- cbind(center, manifest_df)
        } else {
            analysis <- manifest_df$analysis
            manifest_df[, c("analysis", "acl")] <- NULL
            manifest_df <- cbind(analysis, manifest_df)
        }

        manipular <- manifest_df[, "cases"]

        cases <- matrix(nrow = nrow(manifest_df), ncol = 1)

        for (index in seq_len(length(manipular))) {
            patient_code <- manipular[[index]][2]
            cases[index, 1] <- as.character(unlist(patient_code))[6]
        }

        manifest_df$cases <- cases

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[
            grep(tmp, head(manifest_df)$platform),
            c("file_name", "md5sum", "file_id")
        ]

        colnames(manifest_df) <- c("filename", "md5", "id")

        id_matrix <- manifest_df[, "id"]
    }

    # NOTE Exon quantification ####
    if ("exon" == strsplit(tolower(data_type), split = " ")[[1]][1]) {
        if (tolower(data_base) == "gdc") {
            stop(message(
                "\nThrere is no Exon quantification ",
                "data in GDC data base!!",
                "\nPlease use 'legacy' data base"
            ))
        }
        dir.create(
            path = file.path(
                work_dir, "DOAGDC", toupper(tumor),
                folder_name
            ),
            showWarnings = FALSE
        )
        dir <- file.path(work_dir, "DOAGDC", toupper(tumor), folder_name)

        platform <- paste0(
            "%22value%22:%5B%22Illumina%20GA%22,",
            "%22Illumina%20HiSeq%22%5D%7D%7D"
        )

        if (tolower(platform) %in% "all") {
            tmp <- "Illumina GA|Illumina HiSeq"
        } else if (tolower(platform) %in% "illumina ga") {
            tmp <- "Illumina GA"
        } else if (tolower(platform) %in% "illumina hiseq") {
            tmp <- "Illumina HiSeq"
        }

        size <- size_par(
            tumor = tumor, type_of_data = "Gene expression",
            db = "legacy"
        )

        url <- paste0(
            url_inicio, "?pretty=true&expand=cases.samples.",
            "portions.analytes.aliquots,cases.project,center,",
            "analysis&size=", size, "&filters=",
            "%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:",
            "%22in%22,%22content%22:%7B%22",
            "field%22:%22cases.project.program.name%22,%22value",
            "%22:%5B%22TCGA%22%5D%7D%7D,%7B",
            "%22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
            "files.data_type%22,%22value%22:",
            "%5B%22Exon%20quantification%22%5D%7D%7D,%7B%22op%22:",
            "%22in%22,%22content%22:%7B%22",
            "field%22:%22files.access%22,%22value%22:%5B%22open%22%",
            "5D%7D%7D,%7B%22op%22:%22in%",
            "22,%22content%22:%7B%22field%22:%22files.platform%22,",
            platform, ",%7B%22op%22:",
            "%22in%22,%22content%22:%7B%22field%22:%22cases.project",
            ".project_id%22,%22value%22",
            ":%5B%22TCGA-", toupper(tumor),
            "%22%5D%7D%7D%5D%7D&format=JSON"
        )

        message("\n\nDownloading manifest...\n")

        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
            error = function(e) stop(e),
            finally = message(
                "The tumor ", toupper(tumor),
                " does not have ", tolower(data_type),
                "data!"
            )
        )

        manifest_df <- json$data$hits

        # checkin' if there are available data to download (e.g. LAML)
        if (length(manifest_df) == 0) {
            stop("There're not data to be downloaded for this cancer type!")
        }

        center <- manifest_df$center
        manifest_df[, c("center", "acl", "cases", "tags")] <- NULL

        manifest_df <- cbind(center, manifest_df)

        # for some reason sometimes- 16(submitter_id) will be cases
        patient_code <- manifest_df$submitter_id # same protein problem

        write.table(
            x = manifest_df, file = paste0(dir, "/manifest.sdrf"),
            quote = FALSE, row.names = FALSE, sep = "\t"
        )

        manifest_df <- manifest_df[
            grep(tmp, head(manifest_df)$platform),
            c("file_name", "md5sum", "file_id")
        ]

        colnames(manifest_df) <- c("filename", "md5", "id")

        id_matrix <- manifest_df[, "id"]
    }

    # NOTE Download PPD ####
    if (length(dir(dir)) > 1) {
        # verifying if the data is already downloaded
        pattern <- paste(".sdrf", "Data_access_time.txt", sep = "|")
        already_downloaded <- dir(
            path = dir, include.dirs = FALSE,
            recursive = FALSE,
            full.names = TRUE
        )[!grepl(
            pattern,
            dir(path = dir)
        )]
        message("Checking md5 from downloaded files\n")
        already_downloaded_md5 <- as.vector(tools::md5sum(already_downloaded))
        selector <- manifest_df[, "md5"] %in% already_downloaded_md5
        id_matrix <- manifest_df[!selector, "id"]
    }
    if (length(id_matrix) != 0) {
        message("OK!\n")
        ### download data
        url <- paste0(inicio, id_matrix)
        pb <- txtProgressBar(min = 0, max = length(url), style = 3)
        cont <- 0
        for (id in url) {
            cont <- cont + 1
            message(paste("\nDownloading", tumor, data_type, cont, "of",
                length(url),
                sep = " "
            ))
            setTxtProgressBar(pb, cont)
            tmp <- manifest_df[manifest_df$id == id_matrix[cont], "filename"]
            download_httr(url = id, destfile = paste0(dir, "/", tmp))
            md5 <- tools::md5sum(dir(
                path = dir, pattern = tmp,
                full.names = TRUE
            ))
            tmp_md5 <- manifest_df[manifest_df$id == id_matrix[cont], "md5"]
            while (md5[[1]] != tmp_md5) {
                message(paste0(
                    "The md5 of file '", tmp,
                    "' is wrong. Downloading again...\n"
                ))
                download_httr(url = id, destfile = paste0(dir, "/", tmp))
                md5 <- tools::md5sum(dir(dir, tmp, full.names = TRUE))
            }
        }
        close(pb)
        # from python
        # file_endpt = 'https://api.gdc.cancer.gov/files/'
        # file_uuid = 'd853e541-f16a-4345-9f00-88e03c2dc0bc'
        # response = requests.get(file_endpt + file_uuid)
        # message(sprintf(
        # "On %s I realized %s was...\n%s by the street", Sys.Date(), person,
        # action))
        message("\n\nDownload is done!\n\n")
    } else {
        message(paste0(
            "There is nothing to download for ", tumor,
            ". You already have all data available",
            " in the selected data base!"
        ))
    }

    # saving accession data
    write.table(Sys.time(), paste0(dir, "/Data_access_time.txt"),
        quote = FALSE,
        row.names = FALSE, col.names = FALSE
    )
}
