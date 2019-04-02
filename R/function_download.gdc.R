#' Download data from GDC Data Portal and GDC Legacy Archive
#'
#' \code{download_gdc} is a function designed to download methylation, mutation,
#' clinical data, protein expression, MAGETAB, gene expression, isoform
#' expression, miRNA expression and clinical images data from GDC Data Portal
#' and GDC Legacy Archive.
#'
#' @param dataType Type of data. It could be \code{"methylation", "mutation",
#'   "clinical_supplement", "biospecimen", "gene"}. \itemize{ \item{Only present
#'   in "Legacy" database:}{\code{"protein", "Exon quantification", "miRNA gene
#'   quantification", "miRNA isoform quantification", "isoform", "image",
#'   "clinical"}.} \item{Only present in "GDC" database:}{\code{"miRNA
#'   Expression Quantification", "Isoform Expression Quantification"(miRNA)}.} }
#' @param tumor A character string contaning one of the 33 tumors available in
#'   the TCGA project. For instance, the \code{"BRCA"} stands for breast cancer.
#' @param dataBase A character string specifying \code{"GDC"} for GDC Data
#'   Portal or \code{"legacy"} for GDC Legacy Archive.
#' @param HTSeq A character string indicating which HTSeq workflow data should
#'   be downloaded: \code{"Counts", "FPKM", "FPKM-UQ" or "all"}. The default is
#'   \code{"all"}.
#' @param workDir A character string specifying the path to work directory.
#' @param all.files A logical value. Set \code{FALSE} to avoid the download of
#'   not used data to reduce download size, e.g. quantification files. The
#'   default is \code{FALSE}.
#' @param Platform A character string indicating the platform name for
#'   methylation, exon quantificaton, miRNA, and mutation data. \itemize{
#'   \item{For mutation and exon quantificaton data:}{\code{"Illumina GA",
#'   "Illumina HiSeq" or "all"}.} \item{For methylation data}{\code{"Illumina
#'   Human Methylation 450", "Illumina Human Methylation 27" or "all"}.}
#'   \item{For miRNA data:}{\code{"Illumina GA", "Illumina HiSeq",
#'   "H-miRNA_8x15K" (for GBM tumor), "H-miRNA_8x15Kv2" (for OV tumor), or
#'   "all"}.} }The default for all dataType cited is \code{"all"} (when
#'   downloading data).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Downloading gene expression data from GDC Data Portal
#' download_gdc(dataType = "gene", "BRCA", "GDC", workDir = "~/Desktop", HTSeq = "counts")
#'
#' # Downloading mutation data from GDC Legacy Archive
#' download_gdc(dataType = "mutation", "BRCA", "GDC", Platform = "Illumina HiSeq")
#' }
download_gdc <- function(dataType = "gene",
                         tumor,
                         dataBase = "legacy",
                         HTSeq = "",
                         workDir = "~/Desktop",
                         all.files = FALSE,
                         Platform = "all"){

    #verifying if the package is already installed - c("jsonlite")
    # to.load <- c("RCurl", "R.utils", "readr", "jsonlite", "stringr", "tools", "httr")

    # local functions ####
    download.httr <- function(URL, destfile){
        first <- httr::GET(url = URL)
        second <- httr::content(x = first, as = "raw")
        writeBin(object = second, con = destfile)
    }

    size.par <- function(tumor, typeOfData, DB){
        if (DB == "legacy"){
            first.part <- "https://api.gdc.cancer.gov/legacy/projects/TCGA-"
        } else {
            first.part <- "https://api.gdc.cancer.gov/projects/TCGA-"
        }
        url <- paste0(first.part, toupper(tumor),
                      "?expand=summary,summary.data_categories&pretty=true")
        jason <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
        jason <- jason$data
        jason <- jason$summary
        jason <- jason$data_categories
        size <- as.numeric(subset(x = jason, subset = data_category == typeOfData, "file_count"))
        return(size)
    }

    # old api url https://gdc-api.nci.nih.gov/
    # selecting the right API

    # code ####
    if(tolower(dataBase) == "legacy"){
        inicio <- "https://api.gdc.cancer.gov/legacy/data/"
        url.inicio <- "https://api.gdc.cancer.gov/legacy/files/"
        status <- "https://api.gdc.cancer.gov/legacy/status"
        folder.name <- paste0(tolower(dataType), "_data")
    } else if(tolower(dataBase) == "gdc"){
        inicio <- "https://api.gdc.cancer.gov/data/"
        url.inicio <- "https://api.gdc.cancer.gov/files/"
        status <- "https://api.gdc.cancer.gov/status"
        folder.name <- paste(tolower(dataBase), tolower(dataType), "data", sep = "_")
    } else {
        stop("Please insert a data base name!")
    }

    message("Please wait, accessing GDC server...")
    tryCatch(status.gdc <- read.csv(status),
             error = function(e){
                 stop(message(cat("GDC server or your internet conection is off. \n Please try again later!")))
             })

    #check for colunm classes (verification only in manifest)
    # sapply(ncol(manifest.df), function(x){class(manifest.df[, x])}); names() <- colnames(manifest.df)


    dir.create(path = file.path(workDir, "GDCtools"), showWarnings = FALSE)
    dir.create(path = file.path(workDir, "GDCtools", toupper(tumor)), showWarnings = FALSE)

    # legacy ####
    if(tolower(dataBase) == "legacy"){
        # gene and isoform ####
        if ("gene" %in% tolower(dataType) || "isoform" %in% tolower(dataType)){

            size.par.rsem <- function(tumor){
                url <- paste0("https://api.gdc.cancer.gov/legacy/projects/TCGA-", toupper(tumor),
                             "?expand=summary,summary.data_categories&pretty=true")
                jason <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
                jason <- jason$data
                jason <- jason$summary
                jason <- jason$data_categories
                size <- subset(x = jason, subset = data_category == "Gene expression", "file_count")
                return(size)
            }

            size <- size.par.rsem(tumor = tumor)

            if ("gene" %in% tolower(dataType)) {

                dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), "gene_data"),
                           showWarnings = FALSE)
                DIR <- file.path(workDir, "GDCtools", toupper(tumor), "gene_data")

                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22",
                              "in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B",
                              "%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:",
                              "%22cases.project.project_id%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                              "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files",
                              ".data_category%22,%22value%22:%5B%22Gene%20expression%22%5D%7D%7D,%7B%",
                              "22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_type%22,%",
                              "22value%22:%5B%22Gene%20expression%20quantification%22%5D%7D%7D%5D%7D&",
                              "pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:100,%22sort%22:",
                              "%22cases.project.project_id:asc%22%7D%7D&format=JSON")


            } else if ("isoform" %in% tolower(dataType)) {

                dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), "isoform_data"),
                           showWarnings = FALSE)
                DIR <- file.path(workDir, "GDCtools", toupper(tumor), "isoform_data")

                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in",
                              "%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22",
                              "open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
                              "cases.project.project_id%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                              "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
                              "files.data_category%22,%22value%22:%5B%22Gene%20expression%22%5D%7D%",
                              "7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_",
                              "type%22,%22value%22:%5B%22Isoform%20expression%20quantification%22%5D",
                              "%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:100,%",
                              "22sort%22:%22cases.project.project_id:asc%22%7D%7D&format=JSON")
            }

            message("\n\nDownloading manifest...\n")
            json  <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

            manifest.df <- json$data$hits

            manipular <- manifest.df[, "cases"]
            #16(submitter_id) will be cases
            manifest.df[, c("center", "acl", "state_comment", "cases", "tags")] <- NULL

            cases <- matrix(nrow = nrow(manifest.df), ncol = 1)
            for (index in 1:length(manipular)){
                patient.code <- manipular[[index]][2]
                cases[index, 1] <- as.character(unlist(patient.code))[6]
            }

            manifest.df$submitter_id <- cases

            if (!all.files){
                #no need to download these files
                manifest.df <- manifest.df[!grepl("^.+(.pergene.txt)$",
                                                  manifest.df$file_name), ]
                manifest.df <- manifest.df[!grepl("^.+(.quantification.txt)$",
                                                  manifest.df$file_name), ]
                manifest.df <- manifest.df[!grepl("^.+(level3.data.txt)$",
                                                  manifest.df$file_name), ]
                manifest.df <- manifest.df[!grepl("^.+(tags.txt)$",
                                                  manifest.df$file_name), ]
                manifest.df <- manifest.df[!grepl("^.+(genes.txt)$",
                                                  manifest.df$file_name), ]
            }

            write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                        quote = FALSE, row.names = FALSE, sep = "\t")

            manifest.df <- manifest.df[, c("file_name", "md5sum", "file_id")]

            colnames(manifest.df) <- c("filename", "md5", "id" )

            id.matrix <- manifest.df[, "id"]
        }

    } else if(tolower(dataBase) == "gdc"){
        # gene GDC ####
        if ("gene" %in% tolower(dataType)){

            dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), "gdc_gene_data"),
                       showWarnings = FALSE)
            DIR <- file.path(workDir, "GDCtools", toupper(tumor), "gdc_gene_data")

            size.par <- function(tumor){
                url <- paste0("https://api.gdc.cancer.gov/projects/TCGA-", toupper(tumor),
                             "?expand=summary,summary.data_categories&pretty=true")
                jason <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
                jason <- jason$data
                jason <- jason$summary
                jason <- jason$data_categories
                size <- subset(x = jason, subset = data_category == "Transcriptome Profiling", "file_count")
                return(as.numeric(size))
            }

            size <- size.par(tumor = tumor)
            # selecting which HTSeq data to download
            if (tolower(HTSeq) == "counts"){
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22",
                              "in%22,%22content%22:%7B%22field%22:%22files.data_type%22,%22value%22:%5B%22",
                              "Gene%20Expression%20Quantification%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                              "%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-", toupper(tumor), "%22%",
                              "5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,",
                              "%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field",
                              "%22:%22files.data_format%22,%22value%22:%5B%22TXT%22%5D%7D%7D,%7B%22op%22:",
                              "%22in%22,%22content%22:%7B%22field%22:%22files.experimental_strategy%22,%22value%22",
                              ":%5B%22RNA-Seq%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:",
                              "%22files.data_category%22,%22value%22:%5B%22Transcriptome%20Profiling%22%5D%7D%7D,",
                              "%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.analysis.workflow_type%22",
                              ",%22value%22:%5B%22HTSeq%20-%20Counts%22%5D%7D%7D%5D%7D&format=JSON")


            } else if (toupper(HTSeq) == "FPKM" || toupper(HTSeq) == "FPKM-UQ"){
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22",
                              "in%22,%22content%22:%7B%22field%22:%22files.data_type%22,%22value%22:%5B%22",
                              "Gene%20Expression%20Quantification%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                              "%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-", toupper(tumor), "%22%",
                              "5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,",
                              "%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field",
                              "%22:%22files.data_format%22,%22value%22:%5B%22TXT%22%5D%7D%7D,%7B%22op%22:",
                              "%22in%22,%22content%22:%7B%22field%22:%22files.experimental_strategy%22,%22value%22",
                              ":%5B%22RNA-Seq%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:",
                              "%22files.data_category%22,%22value%22:%5B%22Transcriptome%20Profiling%22%5D%7D%7D,",
                              "%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.analysis.workflow_type%22",
                              ",%22value%22:%5B%22HTSeq%20-%20", toupper(HTSeq), "%22%5D%7D%7D%5D%7D&format=JSON")

            }  else if (toupper(HTSeq) == "ALL") {
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22",
                              "in%22,%22content%22:%7B%22field%22:%22files.data_type%22,%22value%22:",
                              "%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D,%7B%22op%22:%22",
                              "in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22",
                              "open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
                              "files.data_format%22,%22value%22:%5B%22TXT%22%5D%7D%7D,%7B%22op%22:%22",
                              "in%22,%22content%22:%7B%22field%22:%22files.analysis.workflow_type%22,%22",
                              "value%22:%5B%22HTSeq%20-%20FPKM-UQ%22,%22HTSeq%20-%20Counts%22,%22HTSeq%20",
                              "-%20FPKM%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
                              "files.experimental_strategy%22,%22value%22:%5B%22RNA-Seq%22%5D%7D%7D,%7B%22",
                              "op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,",
                              "%22value%22:%5B%22TCGA-", toupper(tumor), "%22%5D%7D%7D,%7B%22op%22:%22in%22,",
                              "%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22",
                              "Transcriptome%20Profiling%22%5D%7D%7D%5D%7D&format=JSON")
            } else {
                stop(message('\nPlease insert a HTSeq value! ("Counts", "FPKM", "FPKM-UQ" or "all")\n'))
            }

            # json  <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
            message("\n\nDownloading manifest...\n")
            suppressWarnings(json <- jsonlite::fromJSON(readLines(curl::curl(url))))

            manifest.df <- json$data$hits

            manipular <- manifest.df[, "cases"]
            #16(submitter_id) will be cases
            manifest.df[, c("analysis", "acl")] <- NULL

            cases <- matrix(nrow = nrow(manifest.df), ncol = 1)
            for (index in 1:length(manipular)){
                patient.code <- manipular[[index]][2]
                cases[index, 1] <- as.character(unlist(patient.code))[6]
            }

            manifest.df$cases <- cases

            manifest.df[, "cases"] <- as.character(manifest.df[, "cases"])

            write.table(x = manifest.df, file = paste0(DIR, "/", toupper(HTSeq),
                                                       "_manifest.sdrf"),
                        quote = FALSE, row.names = FALSE, sep = "\t")

            manifest.df <- manifest.df[, c("file_name", "md5sum", "file_id")]

            colnames(manifest.df) <- c("filename", "md5", "id" )

            id.matrix <- manifest.df[, "id"]
        } else if ("isoform" %in% tolower(dataType)){
            stop(message("\nThere is no isoform data in GDC data base! (unfortunately...)\n"))
        }
    }

    # mutation ####
    if ("mutation" %in% tolower(dataType)){

        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        if(tolower(dataBase) == "legacy"){

            platform <- "%22value%22:%5B%22Illumina%20GA%22,%22Illumina%20HiSeq%22%5D%7D%7D"
            #legacy exclusive
            if (tolower(Platform) %in% "all"){
                tmp <- "Illumina GA|Illumina HiSeq"
            } else if (tolower(Platform) %in% "illumina ga"){
                # platform <- "%22value%22:%5B%22Illumina%20GA%22%5D%7D%7D"
                tmp <- "Illumina GA"
            } else if (tolower(Platform) %in% "illumina hiseq"){
                # platform <- "%22value%22:%5B%22Illumina%20HiSeq%22%5D%7D%7D"
                tmp <- "Illumina HiSeq"
            }

            url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                          "portions.analytes.aliquots,cases.project,center,analysis&filters=%7B%22op%22",
                          ":%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                          ":%22files.data_format%22,%22value%22:%5B%22MAF%22%5D%7D%7D,%7B%22op%22:%22in%",
                          "22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%",
                          "5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.experimental",
                          "_strategy%22,%22value%22:%5B%22DNA-Seq%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                          "%22:%7B%22field%22:%22files.platform%22,", platform,",%7B%22op%22:%22in%22,%22",
                          "content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-",
                          toupper(tumor), "%22%5D%7D%7D%5D%7D&format=JSON")

        } else if(tolower(dataBase) == "gdc"){
            url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                          "portions.analytes.aliquots,cases.project,center,analysis&filters=%7B%",
                          "22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:",
                          "%7B%22field%22:%22files.data_format%22,%22value%22:%5B%22MAF%22%5D%7D%",
                          "7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%",
                          "22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%",
                          "22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Simple%20",
                          "Nucleotide%20Variation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:",
                          "%7B%22field%22:%22files.data_type%22,%22value%22:%5B%22Masked%20Somatic",
                          "%20Mutation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                          ":%22files.experimental_strategy%22,%22value%22:%5B%22WXS%22%5D%7D%7D,%7B",
                          "%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id",
                          "%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                          "%22%5D%7D%7D%5D%7D&facetTab=cases&format=JSON")
        }


        json  <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
                 error = function(e) stop(e),
                 finally = message("The tumor ", toupper(tumor),
                                   " does not have ", tolower(dataType), "data!"))

        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (length(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        #center only available at legacy DB
        if (!is.null(manifest.df$center)){
            center <- manifest.df$center
            manifest.df[, c("tags", "acl", "center")] <- NULL
            manifest.df <- cbind(center, manifest.df)
        } else {
            analysis <- manifest.df$analysis
            manifest.df[, c("analysis", "acl")] <- NULL
            manifest.df <- cbind(analysis, manifest.df)
        }

        manipular <- manifest.df[, "cases"]

        pre.cases <- manipular[[1]][2]

        patient.code <- matrix(data = as.character(unlist(pre.cases)), nrow = 10,
                               ncol = length(as.character(unlist(pre.cases)))/10)

        manifest.df[, "cases"] <- paste(patient.code[6, ], collapse = ", ")

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[grep(tmp, head(manifest.df)$platform),
                                   c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # methylation ####
    if ("methylation" %in% tolower(dataType)){
        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        platform <- "%22value%22:%5B%22Illumina%20Human%20Methylation%20450%22,%22Illumina%20Human%20Methylation%2027%22%5D%7D%7D"

        if (tolower(Platform) %in% "all"){
            tmp <- "Illumina Human Methylation 450|Illumina Human Methylation 27"
        } else if (tolower(Platform) %in% "illumina human methylation 450"){
            # platform <- "%22value%22:%5B%22Illumina%20Human%20Methylation%20450%22%5D%7D%7D"
            tmp <- "Illumina Human Methylation 450"
        } else if (tolower(Platform) %in% "illumina human methylation 27"){
            # platform <- "%22value%22:%5B%22Illumina%20Human%20Methylation%2027%22%5D%7D%7D"
            tmp <- "Illumina Human Methylation 27"
        }

        if(tolower(dataBase) == "legacy"){
            size <- size.par(typeOfData = "DNA methylation", tumor = tumor, DB = tolower(dataBase))
            url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                          "portions.analytes.aliquots,cases.project,center,analysis&size=", size,
                          "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,",
                          "%22content%22:%7B%22field%22:%22files.data_format%22,%22value%22:%5B",
                          "%22TXT%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                          ":%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%",
                          "22in%22,%22content%22:%7B%22field%22:%22files.platform%22,", platform,
                          ",%7B%22op%22:%22in%22,%22content%22:%7B%22",
                          "field%22:%22files.data_category%22,%22value%22:%5B%22DNA%20methylation",
                          "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases",
                          ".project.project_id%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                          "%22%5D%7D%7D%5D%7D&format=JSON")
        } else if(tolower(dataBase) == "gdc"){
            size <- size.par(typeOfData = "DNA Methylation", tumor = tumor, DB = tolower(dataBase))
            url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                          "portions.analytes.aliquots,cases.project,center,analysis&size=",
                          size, "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in",
                          "%22,%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:",
                          "%5B%22DNA%20Methylation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:",
                          "%7B%22field%22:%22files.platform%22,", platform,
                          ",%7B%22op%22:%22in%22,%22content%22:%7B%",
                          "22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%",
                          "22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_format%22,%",
                          "22value%22:%5B%22TXT%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B",
                          "%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-", toupper(tumor),
                          "%22%5D%7D%7D%5D%7D&format=JSON")
        }
        message("\n\nDownloading manifest...\n")
        # json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)


        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
                 error = function(e) stop(e),
                 finally = message("The tumor ", toupper(tumor),
                                   " does not have ", tolower(dataType), "data!"))

        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (length(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        #center only available at legacy DB
        if (!is.null(manifest.df$center)){
            center <- manifest.df$center
            manifest.df[, c("tags", "acl", "center")] <- NULL
            manifest.df <- cbind(center, manifest.df)
        } else {
            analysis <- manifest.df$analysis
            manifest.df[, c("analysis", "acl")] <- NULL
            manifest.df <- cbind(analysis, manifest.df)
        }

        manipular <- manifest.df[, "cases"]

        cases <- matrix(nrow = nrow(manifest.df), ncol = 1)

        for (index in 1:length(manipular)){
            patient.code <- manipular[[index]][2]
            cases[index, 1] <- as.character(unlist(patient.code))[6]
        }

        manifest.df$cases <- cases

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[grep(tmp, head(manifest.df)$platform),
                                   c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # clinical and image ####
    if ("clinical" == tolower(dataType) || "biospecimen" == tolower(dataType) || "clinical_supplement" == tolower(dataType) || "image" == tolower(dataType)){

        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        if(tolower(dataBase) == "legacy"){
            if ("image" == tolower(dataType)){
                message("A lot of data are going to be downloaded, please wait...")
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=10000&filters=%7B%22op%22",
                              ":%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                              ":%22files.data_category%22,%22value%22:%5B%22Clinical%22%5D%7D%7D,%7B%22op%22",
                              ":%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22",
                              "open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.",
                              "platform%22,%22value%22:%5B%22Clinical%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                              "content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-",
                              toupper(tumor), "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                              ":%22files.data_format%22,%22value%22:%5B%22SVS%22%5D%7D%7D%5D%7D&format=JSON")
            } else if ("clinical_supplement" == tolower(dataType)) {
                size <- size.par(tumor = tumor, typeOfData = "Clinical", DB = "legacy")
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=",
                              "%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22",
                              "field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%",
                              "22,%22content%22:%7B%22field%22:%22cases.project.program.name%22,%22value%22:%5B%",
                              "22TCGA%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.pro",
                              "ject.project_id%22,%22value%22:%5B%22TCGA-",toupper(tumor), "%22%5D%7D%7D,%7B%22op%22:%22in%22,",
                              "%22content%22:%7B%22field%22:%22files.data_format%22,%22value%22:%5B%22BCR%20XML",
                              "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_type",
                              "%22,%22value%22:%5B%22Clinical%20Supplement%22%5D%7D%7D%5D%7D&format=JSON")
            } else if ("clinical" == tolower(dataType)) {
                size <- size.par(tumor = tumor, typeOfData = "Clinical", DB = "legacy")
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=%7B%",
                              "22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                              ":%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content",
                              "%22:%7B%22field%22:%22cases.project.program.name%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,",
                              "%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22",
                              "value%22:%5B%22TCGA-", toupper(tumor), "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%",
                              "22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Clinical%22%5D%7D%7D,",
                              "%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_format%22,%22value%",
                              "22:%5B%22Biotab%22%5D%7D%7D%5D%7D&format=JSON")
            } else if ("biospecimen" == tolower(dataType)) {
                size <- size.par(tumor = tumor, typeOfData = "Clinical", DB = "legacy")
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=%7B%22",
                              "op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22",
                              "files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:",
                              "%7B%22field%22:%22cases.project.program.name%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,%7B%",
                              "22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%",
                              "22:%5B%22TCGA-", toupper(tumor), "%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%",
                              "22field%22:%22files.data_format%22,%22value%22:%5B%22BCR%20XML%22%5D%7D%7D,%7B%22op%22",
                              ":%22in%22,%22content%22:%7B%22field%22:%22files.data_type%22,%22value%22:%5B%22Biospeci",
                              "men%20Supplement%22%5D%7D%7D%5D%7D&format=JSON")
            }
        } else if(tolower(dataBase) == "gdc"){
            if ("clinical_supplement" == tolower(dataType)) {
                size <- size.par(tumor = tumor, typeOfData = "Clinical", DB = "gdc")
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=%7B%22op%",
                              "22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22",
                              "%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-", toupper(tumor), "%22%5D%",
                              "7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22",
                              "value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field",
                              "%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22BCR%20XML%22%5D%7D%7D%2C%7B%22op%22%",
                              "3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22",
                              "Clinical%20Supplement%22%5D%7D%7D%5D%7D&format=JSON")
            } else if ("biospecimen" == tolower(dataType)) {
                size <- size.par(tumor = tumor, typeOfData = "Clinical", DB = "gdc")
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=%7B%22op%",
                              "22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%",
                              "22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-", toupper(tumor),
                              "%22%5D%7D%7D%2C%7B%22op",
                              "%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22",
                              "open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.",
                              "data_format%22%2C%22value%22%3A%5B%22BCR%20XML%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22",
                              "content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Biospecimen%20",
                              "Supplement%22%5D%7D%7D%5D%7D&format=JSON")
            }

        }

        message("\n\nDownloading manifest...\n")
        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (length(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        #center only available at legacy DB
        if (!is.null(manifest.df$center)){
            center <- manifest.df$center
            manifest.df[, c("tags", "acl", "center", "cases")] <- NULL
            manifest.df <- cbind(center, manifest.df)
        } else {
            manifest.df[, c("tags", "acl", "cases")] <- NULL
        }

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[, c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # protein ####
    if ("protein" == tolower(dataType)){
        if (tolower(dataBase) == "gdc") {
            stop(message("\nThrere is no protein expression data in GDC data base!!",
                         "\nPlease use 'legacy' data base"))
        }
        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        size <- size.par(tumor = tumor, typeOfData = "Protein expression", DB = "legacy")

        url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                      "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=%7B%22op%22",
                      ":%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22",
                      ":%22files.data_category%22,%22value%22:%5B%22Protein%20expression%22%5D%7D%",
                      "7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22",
                      "value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22",
                      "field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-",
                      toupper(tumor), "%22%5D%7D%7D%5D%7D&format=JSON")

        message("\n\nDownloading manifest...\n")
        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (length(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        center <- manifest.df$center
        manifest.df[, c("center", "acl", "cases")] <- NULL

        manifest.df <- cbind(center, manifest.df)

        #for some reason sometimes- 16(submitter_id) will be cases
        patient.code <- manifest.df$submitter_id

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[, c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # mage ####
    if ("mage" == tolower(dataType)){
        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                      "portions.analytes.aliquots,cases.project,center,analysis&size=1000000&filters=%7B%22op%22:",
                      "%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:",
                      "%22files.data_format%22,%22value%22:%5B%22MAGETAB%22%5D%7D%7D,%7B%22op%22:%22",
                      "in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%",
                      "22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_type",
                      "%22,%22value%22:%5B%22TCGA%20DCC%20Archive%22%5D%7D%7D,%7B%22op%22:%22in%22,%22",
                      "content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Archive%22",
                      "%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:100,%22sort",
                      "%22:%22cases.project.project_id:asc%22%7D%7D&format=JSON")

        message("\n\nDownloading manifest...\n")
        json <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)

        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (nrow(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        #16(submitter_id) will be cases
        manifest.df[, c("acl", "cases")] <- NULL

        manifest.df <- manifest.df[grepl(toupper(tumor), manifest.df$file_name,
                                         ignore.case = FALSE, perl = TRUE), ]

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[, c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # mirna ####
    if ("mirna" %in% strsplit(tolower(dataType), split = " ")[[1]][1] || "isoform expression quantification" %in% tolower(dataType)){
        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        platform <- ""

        if (tolower(Platform) %in% "all"){
            tmp <- "Illumina HiSeq|Illumina GA|H-miRNA_8x15Kv2|H-miRNA_8x15Kv"
        } else if (tolower(Platform) %in% "illumina hiseq"){
            platform <- ",%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.platform%22,%22value%22:%5B%22Illumina%20HiSeq%22%5D%7D%7D"
            tmp <- "Illumina HiSeq"
        } else if (tolower(Platform) %in% "illumina ga"){
            platform <- ",%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.platform%22,%22value%22:%5B%22Illumina%20GA%22%5D%7D%7D"
            tmp <- "Illumina GA"
        } else if (tolower(Platform) %in% "h-mirna_8x15kv2") {
            platform <- ",%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.platform%22,%22value%22:%5B%22H-miRNA_8x15Kv2%22%5D%7D%7D"
            tmp <- "H-miRNA_8x15Kv2"
        } else if (tolower(Platform) %in% "h-mirna_8x15kv") {
            platform <- ",%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.platform%22,%22value%22:%5B%22H-miRNA_8x15Kv%22%5D%7D%7D"
            tmp <- "H-miRNA_8x15Kv"
        }

        if(tolower(dataBase) == "legacy"){
            size <- size.par(typeOfData = "Gene expression", tumor = tumor, DB = tolower(dataBase))
            if ("mirna gene quantification" %in% tolower(dataType)) {
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,
                              "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%",
                              "22content%22:%7B%22field%22:%22files.data_type%22,%22value%22:%5B%22miRNA%",
                              "20gene%20quantification%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22",
                              "field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:",
                              "%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22",
                              "value%22:%5B%22TCGA-", toupper(tumor),"%22%5D%7D%7D",platform,
                              "%5D%7D&format=JSON")
                print(url)
            } else if ("mirna isoform quantification" %in% tolower(dataType)) {
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=", size,
                              "&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22",
                              "content%22:%7B%22field%22:%22files.data_type%22,%22value%22:%5B%22miRNA%",
                              "20isoform%20quantification%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22",
                              ":%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B",
                              "%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id",
                              "%22,%22value%22:%5B%22TCGA-", toupper(tumor),"%22%5D%7D%7D", platform,
                              "%5D%7D&format=JSON")
            }

        } else if(tolower(dataBase) == "gdc"){
            size <- size.par(typeOfData = "Transcriptome Profiling", tumor = tumor, DB = tolower(dataBase))
            if ("mirna expression quantification" %in% tolower(dataType)) {
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B",
                              "%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.",
                              "project.project_id%22%2C%22value%22%3A%5B%22TCGA-", toupper(tumor),
                              "%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22",
                              "field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%",
                              "7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22",
                              "files.data_type%22%2C%22value%22%3A%5B%22miRNA%20Expression%20Quan",
                              "tification%22%5D%7D%7D%5D%7D&format=JSON")
            } else if ("isoform expression quantification" %in% tolower(dataType)) {
                url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                              "portions.analytes.aliquots,cases.project,center,analysis&size=",
                              size, "&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22",
                              "op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project",
                              ".project_id%22%2C%22value%22%3A%5B%22TCGA-", toupper(tumor),
                              "%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field",
                              "%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%",
                              "7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.",
                              "data_type%22%2C%22value%22%3A%5B%22Isoform%20Expression%20Quantifi",
                              "cation%22%5D%7D%7D%5D%7D&format=JSON")
            }
        }
        message("\n\nDownloading manifest...\n")

        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
                 error = function(e) stop(e),
                 finally = message("The tumor ", toupper(tumor),
                                   " does not have ", tolower(dataType), "data!"))

        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (length(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        #center only available at legacy DB
        if (!is.null(manifest.df$center)){
            center <- manifest.df$center
            manifest.df[, c("tags", "acl", "center")] <- NULL
            manifest.df <- cbind(center, manifest.df)
        } else {
            analysis <- manifest.df$analysis
            manifest.df[, c("analysis", "acl")] <- NULL
            manifest.df <- cbind(analysis, manifest.df)
        }

        manipular <- manifest.df[, "cases"]

        cases <- matrix(nrow = nrow(manifest.df), ncol = 1)

        for (index in 1:length(manipular)){
            patient.code <- manipular[[index]][2]
            cases[index, 1] <- as.character(unlist(patient.code))[6]
        }

        manifest.df$cases <- cases

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[grep(tmp, head(manifest.df)$platform),
                                   c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # Exon quantification #####
    if ("exon" == strsplit(tolower(dataType), split = " ")[[1]][1]){
        if (tolower(dataBase) == "gdc") {
            stop(message("\nThrere is no Exon quantification data in GDC data base!!",
                         "\nPlease use 'legacy' data base"))
        }
        dir.create(path = file.path(workDir, "GDCtools", toupper(tumor), folder.name),
                   showWarnings = FALSE)
        DIR <- file.path(workDir, "GDCtools", toupper(tumor), folder.name)

        platform <- "%22value%22:%5B%22Illumina%20GA%22,%22Illumina%20HiSeq%22%5D%7D%7D"

        if (tolower(Platform) %in% "all"){
            tmp <- "Illumina GA|Illumina HiSeq"
        } else if (tolower(Platform) %in% "illumina ga"){
            # platform <- "%22value%22:%5B%22Illumina%20GA%22%5D%7D%7D"
            tmp <- "Illumina GA"
        } else if (tolower(Platform) %in% "illumina hiseq"){
            # platform <- "%22value%22:%5B%22Illumina%20HiSeq%22%5D%7D%7D"
            tmp <- "Illumina HiSeq"
        }

        size <- size.par(tumor = tumor, typeOfData = "Gene expression", DB = "legacy")

        url <- paste0(url.inicio, "?pretty=true&expand=cases.samples.",
                      "portions.analytes.aliquots,cases.project,center,analysis&size=", size,"&filters=",
                      "%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22",
                      "field%22:%22cases.project.program.name%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,%7B",
                      "%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_type%22,%22value%22:",
                      "%5B%22Exon%20quantification%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22",
                      "field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D,%7B%22op%22:%22in%",
                      "22,%22content%22:%7B%22field%22:%22files.platform%22,", platform, ",%7B%22op%22:",
                      "%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22",
                      ":%5B%22TCGA-", toupper(tumor), "%22%5D%7D%7D%5D%7D&format=JSON")


        message("\n\nDownloading manifest...\n")

        json <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
                 error = function(e) stop(e),
                 finally = message("The tumor ", toupper(tumor),
                                   " does not have ", tolower(dataType), "data!"))


        manifest.df <- json$data$hits

        #checkin' if there are available data to download (e.g. LAML)
        if (length(manifest.df) == 0){
            stop("There're not data to be downloaded for this cancer type!")
        }

        center <- manifest.df$center
        manifest.df[, c("center", "acl", "cases", "tags")] <- NULL

        manifest.df <- cbind(center, manifest.df)

        #for some reason sometimes- 16(submitter_id) will be cases
        patient.code <- manifest.df$submitter_id # same protein problem

        write.table(x = manifest.df, file = paste0(DIR, "/manifest.sdrf"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        manifest.df <- manifest.df[grep(tmp, head(manifest.df)$platform),
                                   c("file_name", "md5sum", "file_id")]

        colnames(manifest.df) <- c("filename", "md5", "id" )

        id.matrix <- manifest.df[, "id"]
    }

    # Download PPD ####
    if (length(dir(DIR)) > 1){
        #verifying if the data is already downloaded
        already.downloaded <- dir(path = DIR, include.dirs = FALSE, recursive = FALSE,
                                  full.names = TRUE)[!grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                                                                       sep = "|"),
                                                                       x = dir(path = DIR))]
        message("Checking md5 from downloaded files\n")
        already.downloaded.md5 <- as.vector(tools::md5sum(already.downloaded))
        selector <- manifest.df[, "md5"] %in% already.downloaded.md5
        id.matrix <- manifest.df[!selector, "id"]
    } else {
        # ("id.matrix", id.matrix, envir = get(paste(tumor, dataBase, dataType, sep = "_")))
    }

    if (length(id.matrix) != 0){
        message("OK!\n")
        ###download data
        url <- paste0(inicio, id.matrix)
        pb <- txtProgressBar(min = 0, max = length(url), style = 3)
        contador <- 0
        for (id in url){
            contador <- contador + 1
            message(paste("\nDownloading", tumor, dataType, contador, "of", length(url), sep = " "))
            setTxtProgressBar(pb, contador)
            download.httr(URL = id,
                                   destfile = paste0(DIR, "/",
                                                    manifest.df[manifest.df$id == id.matrix[contador], "filename"]))
            md5 <- tools::md5sum(dir(path = DIR, pattern = manifest.df[manifest.df$id == id.matrix[contador],
                                                                "filename"], full.names = TRUE))
            while(md5[[1]] != manifest.df[manifest.df$id == id.matrix[contador], "md5"]){
                message(paste0("The md5 of file '", manifest.df[manifest.df$id == id.matrix[contador],
                                                              "filename"], "' differs from the original file. Downloading again...\n"))
                download.httr(URL = id,
                              destfile = paste0(DIR, "/",
                                               manifest.df[manifest.df$id == id.matrix[contador],
                                                           "filename"]))
                md5 <- tools::md5sum(dir(path = DIR, pattern = manifest.df[manifest.df$id == id.matrix[contador],
                                                                     "filename"], full.names = TRUE))
            }
        }
        close(pb)
        # from python
        # file_endpt = 'https://api.gdc.cancer.gov/files/'
        # file_uuid = 'd853e541-f16a-4345-9f00-88e03c2dc0bc'
        # response = requests.get(file_endpt + file_uuid)
        # message(sprintf("On %s I realized %s was...\n%s by the street", Sys.Date(), person, action))
        message("\n\nDownload is done!\n\n")
    } else {
        message(paste0("There is nothing to download for ", tumor,
                       ". You already have all data available in the selected data base!"))
    }

    #saving accession data
    write.table(Sys.time(), paste0(DIR, "/Data_access_time.txt"), quote = FALSE,
                row.names = FALSE, col.names = FALSE)

}
