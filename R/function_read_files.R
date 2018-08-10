#' Concatenate GDC files into a single matrix and prepar the data
#'
#' \code{concatenate_files} is a function designed to concatenate GDC files into
#' a single matrix, where the columns stand for patients code and rows stand for
#' data names.
#'
#' @param dataType
#' @param normalization Logical value where \code{TRUE} specify the desire to
#'   work with normalized files only. When FALSE, in the second run, do not
#'   forget to set env argument. This argument is only applyable to gene and
#'   isoform expression data from GDC Legacy Archive. The default is
#'   \code{TRUE}.
#' @param Name A character string indicating which row name will be used in
#'   mclust separation, e.g. "HIF3A" in the gene expression concatenated matrix.
#' @param dataBase
#' @param HTSeq A character string indicating which HTSeq workflow data should
#'   be downloaded: "Counts", "FPKM", or "FPKM-UQ".
#' @param workDir
#' @param tumor
#' @param workflowType A character string specifying the workflow type for
#'   mutation data in "gdc". Where "varscan" stands for 'VarScan2 Variant
#'   Aggregation and Masking', "mutect" stands for 'MuTect2 Variant Aggregation
#'   and Masking', "muse" stands for 'MuSE Variant Aggregation and Masking',
#'   "somaticsniper" stands for 'SomaticSniper Variant Aggregation and Masking'
#'   and "all" stands for concatenate all workflows in a single matrix.
#' @param tumorData Logical value where \code{TRUE} specifies the desire to work
#'   with tumor tissue files only. When set to FALSE, it creates two matrices,
#'   one containing tumor data and other containing data from not-tumor tissue.
#'   The default is \code{TRUE}.
#' @param cutoffBetaNA Numerical value indicating the maximum threshold
#'   percentage (in decimal form) to tolerate and to remove rows containing NA
#'   for beta values (methylation data). The default is 0.25.
#' @param cutoffBetasd Numerical value indicating the standard deviation
#'   threshold of beta values (methylation data). It keeps only rows that have
#'   standard deviation of beta values higher than the threshold. The default is
#'   \code{0.005}.
#' @param onlyFilter Logical value where \code{TRUE} indicates that the matrix
#'   is already concatenate and the function should choose a different
#'   \code{Name}, without concatenate all the files again. The default is FALSE.
#' @param Platform
#' @param use.hg19.mirbase20 Logical value where \code{TRUE} indicates that only
#'   hg19.mirbase20 should be used. This parameter is needed when using
#'   \code{dataBase = "legacy"} and one of the available miRNA \code{dataType}
#'   in "legacy" ("miRNA gene quantification" and "miRNA isoform
#'   quantification"). The default is FALSE.
#' @param env A character string containing the environment name that should be
#'   used. If none has been set yet, the function will create one in global
#'   environment following the standard criteria:
#'   \itemize{\item{'tumor_dataBase_dataType_tumor_data'}{ or}
#'   \item{'tumor_dataBase_dataType_both_data'}{ (for tumor and not tumor data in
#'   separated matrices).}}
#' @param saveData Logical value where \code{TRUE} indicates that the
#'   concatenate and filtered matrix should be saved in local storage. The
#'   default is FALSE.
#' @inheritParams download_gdc
#'
#' @return A matrix with data names in row and patients code in column.
#' @export
#'
#' @examples
#' \dontrun{
#' # Concatenating isoform expression data into a single matrix
#' concatenate_files(dataType = "isoform", Name = "uc002peh.2", dataBase = "legacy", tumor = "BRCA")
#' }
concatenate_files <- function(dataType,
                      normalization = TRUE,
                      Name , dataBase,
                      HTSeq = NULL,
                      workDir = "~/Desktop", tumor,
                      workflowType,
                      tumorData = TRUE,
                      cutoffBetaNA = 0.25,
                      cutoffBetasd = 0.005,
                      onlyFilter = FALSE,
                      Platform = "",
                      use.hg19.mirbase20 = FALSE,
                      env,
                      saveData = FALSE) {

    # #verifying if the package is already installed - c("data.table" , "R.utils")

    #local functions ####
    gene.isoform <- function(){
        codes <- dir(path = DIR, pattern = ".sdrf$")
        codigos <- read.table(file.path(DIR, codes[1]), stringsAsFactors = FALSE,
                              header = TRUE, sep = "\t")
        codigos <- codigos[, c("file_name", "submitter_id")]
        tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbBcCdD]-|-")
        not.tumor.regex <- paste(10:19, collapse = "[aAbBcCdD]-|-")

        if (tumorData){
            codigos <- codigos[grepl(paste0("-", tumor.regex, "-"),
                                     x = codigos$submitter_id), ]
            remove.duplicated <- duplicated(codigos$file_name)
            codigos <- codigos[!remove.duplicated, ]
            rownames(codigos) <- codigos$file_name

            #selecting files by type (normalized or not)
            if (normalization) {
                genes.normalized <- grepl("^.+(genes.normalized).+$", codigos$file_name)
                isoforms.normalized <- grepl("^.+(isoforms.normalized).+$", codigos$file_name)
                gene.files <- codigos[genes.normalized, ]
                isoform.files <- codigos[isoforms.normalized, ]
                assign("genes.normalized", genes.normalized, envir = get(envir_link))
                assign("isoforms.normalized", isoforms.normalized, envir = get(envir_link))
                assign("gene.files", gene.files, envir = get(envir_link))
                assign("isoform.files", isoform.files, envir = get(envir_link))
            } else {
                genes.not.normalized <- grepl("^.+(.rsem.genes.results)$", codigos$file_name)
                isoforms.not.normalized <- grepl("^.+(isoforms.results)$", codigos$file_name)
                gene.files <- codigos[genes.not.normalized, ]
                isoform.files <- codigos[isoforms.not.normalized, ]
                assign("genes.not.normalized", genes.not.normalized, envir = get(envir_link))
                assign("isoforms.not.normalized", isoforms.not.normalized, envir = get(envir_link))
                assign("gene.files", gene.files, envir = get(envir_link))
                assign("isoform.files", isoform.files, envir = get(envir_link))
            }
            assign("codigos", codigos, envir = get(envir_link))
        } else {
            codigos.tumor <- codigos[grepl(paste0("-", tumor.regex, "-"),
                                     x = codigos$submitter_id), ]
            codigos.not.tumor <- codigos[grepl(paste0("-", not.tumor.regex, "[aAbBcCdD]-"),
                                     x = codigos$submitter_id), ]
            remove.duplicated.tumor <- duplicated(codigos.tumor$file_name)
            remove.duplicated.not.tumor <- duplicated(codigos.not.tumor$file_name)
            codigos.tumor <- codigos.tumor[!remove.duplicated.tumor, ]
            codigos.not.tumor <- codigos.not.tumor[!remove.duplicated.not.tumor, ]

            rownames(codigos.tumor) <- codigos.tumor$file_name
            rownames(codigos.not.tumor) <- codigos.not.tumor$file_name

            #selecting files by type (normalized or not)
            if (normalization) {
                #tumor
                genes.normalized <- grepl("^.+(genes.normalized).+$", codigos.tumor$file_name)
                isoforms.normalized <- grepl("^.+(isoforms.normalized).+$", codigos.tumor$file_name)
                gene.files <- codigos.tumor[genes.normalized, ]
                isoform.files <- codigos.tumor[isoforms.normalized, ]
                assign("genes.normalized", genes.normalized, envir = get(envir_link))
                assign("isoforms.normalized", isoforms.normalized, envir = get(envir_link))
                assign("gene.files", gene.files, envir = get(envir_link))
                assign("isoform.files", isoform.files, envir = get(envir_link))

                #not tumor
                genes.normalized.not.tumor <- grepl("^.+(genes.normalized).+$", codigos.not.tumor$file_name)
                isoforms.normalized.not.tumor <- grepl("^.+(isoforms.normalized).+$", codigos.not.tumor$file_name)
                gene.files.not.tumor <- codigos.not.tumor[genes.normalized.not.tumor, ]
                isoform.files.not.tumor <- codigos.not.tumor[isoforms.normalized.not.tumor, ]
                assign("genes.normalized.not.tumor", genes.normalized.not.tumor, envir = get(envir_link))
                assign("isoforms.normalized.not.tumor", isoforms.normalized.not.tumor, envir = get(envir_link))
                assign("gene.files.not.tumor", gene.files.not.tumor, envir = get(envir_link))
                assign("isoform.files.not.tumor", isoform.files.not.tumor, envir = get(envir_link))
            } else {
                #tumor
                genes.not.normalized <- grepl("^.+(.rsem.genes.results)$", codigos.tumor$file_name)
                isoforms.not.normalized <- grepl("^.+(isoforms.results)$", codigos.tumor$file_name)
                gene.files <- codigos.tumor[genes.not.normalized, ]
                isoform.files <- codigos.tumor[isoforms.not.normalized, ]
                assign("genes.not.normalized", genes.not.normalized, envir = get(envir_link))
                assign("isoforms.not.normalized", isoforms.not.normalized, envir = get(envir_link))
                assign("gene.files", gene.files, envir = get(envir_link))
                assign("isoform.files", isoform.files, envir = get(envir_link))

                #not tumor
                genes.not.normalized.not.tumor <- grepl("^.+(.rsem.genes.results)$", codigos.not.tumor$file_name)
                isoforms.not.normalized.not.tumor <- grepl("^.+(isoforms.results)$", codigos.not.tumor$file_name)
                gene.files.not.tumor <- codigos.not.tumor[genes.not.normalized.not.tumor, ]
                isoform.files.not.tumor <- codigos.not.tumor[isoforms.not.normalized.not.tumor, ]
                assign("genes.not.normalized.not.tumor", genes.not.normalized.not.tumor, envir = get(envir_link))
                assign("isoforms.not.normalized.not.tumor", isoforms.not.normalized.not.tumor, envir = get(envir_link))
                assign("gene.files.not.tumor", gene.files.not.tumor, envir = get(envir_link))
                assign("isoform.files.not.tumor", isoform.files.not.tumor, envir = get(envir_link))
            }
            assign("codigos.tumor", codigos.tumor, envir = get(envir_link))
            assign("codigos.not.tumor", codigos.not.tumor, envir = get(envir_link))
        }
    }

    OPEN.genexpress <- function(files, Var){
        message("Reading data...")
        #reading files in one file
        pb <- txtProgressBar(min = 0, max = length(files$file_name), style = 3)
        count <- 0
        for (i in files$file_name) {
            count <- count + 1
            setTxtProgressBar(pb, count)
            actual.file <- data.table::fread(file.path(DIR, i))
            if (count == 1) {
                actual.file <- data.table::fread(file.path(DIR, i))
                completed.table1 <- matrix(nrow = nrow(actual.file),
                                           ncol = length(files$file_name))
                rownames(completed.table1) <- actual.file[[1]]
                completed.table1[, 1] <- as.numeric(actual.file[[2]])
            } else {
                actual.file <- data.table::fread(file.path(DIR, i), select = 2)
                completed.table1[, count] <- as.numeric(actual.file[[1]])
            }
        }
        close(pb)
        if (Var == "gene.files" || Var == "isoform.files"){
            assign("name.row", rownames(completed.table1), envir = get(envir_link))
            assign("completed.table1", completed.table1, envir = get(envir_link))
        } else if (Var == "gene.files.not.tumor" || Var == "isoform.files.not.tumor"){
            assign("name.row.not.tumor", rownames(completed.table1), envir = get(envir_link))
            assign("completed.table1.not.tumor", completed.table1, envir = get(envir_link))
        }
    }

    gene_finder <- function(completed.table, variable_name) {
        to.google <- ifelse(suppressWarnings(is.na(as.numeric(Name))),
                            paste0("(", toupper(Name), "\\|", ")"),
                            paste0("(", "\\|", Name, ")$"))
        gene.selector.final <- grepl(to.google,
                                     rownames(completed.table), perl=TRUE)

        genes.1 <- completed.table[gene.selector.final, ]
        if (length(genes.1) == 0) {
            message("Gene '", Name, "' not found!!!\n")
            stop()
        }
        genes.transpo1 <- as.matrix(genes.1)
        colnames(genes.transpo1) <- rownames(completed.table)[gene.selector.final]
        assign(paste0(variable_name, toupper(Name)), genes.transpo1, envir = get(envir_link))
    }

    isoform_finder <- function(completed.table, variable_name){
        Isoforms.1 <- completed.table[Name, , drop = FALSE]
        if (length(Isoforms.1) == 0) {
            message("'", Name, "' not found!!!\n")
            stop()
        }
        Isoforms.transpo1 <- t(Isoforms.1)
        Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
        # colnames(Isoforms.transpo1) <- Name

        assign(paste0(variable_name, toupper(Name)), Isoforms.transpo1, envir = get(envir_link))
    }

    OPEN.meth <- function(meth.files.arg){
        #reading files in one file
        pb <- txtProgressBar(min = 0, max = length(meth.files.arg), style = 3)
        count <- 0
        #try to open with lapply(list, function)
        for (file in file.path(DIR, meth.files.arg)){
            count <- count + 1
            if (count > 1){
                actual.file <- data.table::fread(file, sep = "\t",
                                                 showProgress = FALSE, select = 2)
                meth_table[, count] <- as.numeric(actual.file[[1]])
            } else if(count == 1){
                actual.file <- data.table::fread(file, sep = "\t", showProgress = FALSE)
                # actual.file <- actual.file[, c(1, 3, 6, 4, 5, 2), with = FALSE]
                refererence_meth_table <- as.data.frame(actual.file[, Select_ref, with = FALSE])
                rownames(refererence_meth_table) <- as.character(as.matrix(actual.file[, 1]))
                meth_table <- matrix(nrow = nrow(actual.file),
                                     dimnames = list(rownames(refererence_meth_table), patients),
                                                 ncol = length(meth.files.arg))
                meth_table[, 1] <- as.numeric(actual.file[[2]])
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)
        if (tolower(dataBase) == "legacy"){
            colnames(refererence_meth_table) <- refererence_meth_table[1, ]
            refererence_meth_table <- refererence_meth_table[-1, , drop =FALSE]
            meth_table <- meth_table[-1, , drop =FALSE]
        }
        assign("meth_table", meth_table, envir = get(envir_link))
        assign("refererence_meth_table", refererence_meth_table, envir = get(envir_link))
    }

    filter.meth.FUN <- function(meth_table, refererence_meth_table){

        ### preparing data          FASTER WITH MATRIX
        message("Please wait... Filtering data!")
        pb <- txtProgressBar(min = 0, max = 3, style = 3)
        count <- 0
        setTxtProgressBar(pb, count)
        initial.row <- nrow(meth_table)
        # The beta value (β) estimates the methylation level of the CpG locus.
        # beta value (β) = methylated/unmethylated alleles
        #Remove sites containing NA for beta values based in cutoff
        beta.na.filter <- rowMeans(is.na(meth_table)) < cutoffBetaNA
        meth_table <- meth_table[beta.na.filter, ]
        refererence_meth_table <- refererence_meth_table[beta.na.filter, ]

        count <- 1
        setTxtProgressBar(pb, count)

        #Keep sites for which the beta values have standard deviation value higher than cutoff
        beta.sd.filter <- apply(meth_table, 1, sd, na.rm = TRUE) > cutoffBetasd
        meth_table <- meth_table[beta.sd.filter, ]
        refererence_meth_table <- refererence_meth_table[beta.sd.filter, ]

        count <- 2
        setTxtProgressBar(pb, count)

        #Remove CpGs on sex chromosomes
        sex.filter <- refererence_meth_table$Chromosome != "X" & refererence_meth_table$Chromosome != "Y"
        sex.filter[is.na(sex.filter)] <- as.logical("FALSE")
        meth_table <- meth_table[sex.filter, ]
        refererence_meth_table <- refererence_meth_table[sex.filter, ]

        #only with gene symbol info (can be changed after released)
        gene.symbol.filter <- refererence_meth_table$Gene_Symbol != "" & refererence_meth_table$Gene_Symbol != "."
        meth_table <- meth_table[gene.symbol.filter, ]
        refererence_meth_table <- refererence_meth_table[gene.symbol.filter, ]

        count <- 3
        setTxtProgressBar(pb, count)

        #Number of removed lines after filtering
        message(paste("It was filtered", initial.row - nrow(meth_table),
                      "lines from data!", sep = " "))
        close(pb)
        #normalization (3 options) Batch effect: systematic differences across groups of samp
        # The background was subtracted using the methylumi package (method "noob") [1].
        # The signal intensity values were normalized using the SWAN normalization method, as implemented in the minfi package.
        #or lumi R package or minfi Bioconductor packag(SWAN)
        # library(methylumi) ?? wateRmellow??

        if (nrow(meth_table) == 0){
            message("All data were filtered out! Please, ",
                    "choose a different cutoffBetaNA and cutoffBetasd parameters")
            stop()
        }

        #export to global enviroment
        assign("meth_table", meth_table, envir = get(envir_link))
        assign("refererence_meth_table", refererence_meth_table, envir = get(envir_link))
    }

    OPEN.miRNA <- function(mirna.files.arg){
        if (normalization) {
            if (tolower(dataType) == "mirna gene quantification" || tolower(dataType) == "mirna expression quantification"){
                Select <- c(3,4)
                first_col <- "reads_per_million_miRNA_mapped"
            } else if (tolower(dataType) == "mirna isoform quantification" || tolower(dataType) == "isoform expression quantification") {
                # Select <- c(4,5)
                # first_col <- "reads_per_million_miRNA_mapped"
                message("Under development for this dataType!\n")
                stop()
            }

        } else {
            if (tolower(dataType) == "mirna gene quantification" || tolower(dataType) == "mirna expression quantification"){
                Select <- c(2,4)
                first_col <- "read_count"
            } else if (tolower(dataType) == "mirna isoform quantification" || tolower(dataType) == "isoform expression quantification") {
                # Select <- c(3,5)
                # first_col <- "read_count"
                message("Under development for this dataType!\n")
                stop()
            }
        }
        #reading files in one file
        pb <- txtProgressBar(min = 0, max = length(mirna.files.arg), style = 3)
        count <- 0
        #try to open with lapply(list, function)
        for (file in file.path(DIR, mirna.files.arg)){
            count <- count + 1
            if (count > 1){
                actual.file <- data.table::fread(file, sep = "\t",
                                                 showProgress = FALSE, select = Select)
                data.mirna.table[, count] <- as.numeric(actual.file[[1]])
                cross.mapped.mirna.table[, count] <- as.character(actual.file[[2]])
            } else if(count == 1){
                actual.file <- data.table::fread(file, sep = "\t", showProgress = FALSE)
                #table(sapply(actual.file,class))
                data.mirna.table <- matrix(nrow = nrow(actual.file),
                                           dimnames = list(as.character(actual.file[["miRNA_ID"]]),
                                                           manifest$cases),
                                          ncol = length(mirna.files.arg))
                cross.mapped.mirna.table <- matrix(nrow = nrow(actual.file),
                                                   dimnames = list(as.character(actual.file[["miRNA_ID"]]),
                                                                   manifest$cases),
                                                   ncol = length(mirna.files.arg))
                data.mirna.table[, 1] <- as.numeric(actual.file[["read_count"]])
                cross.mapped.mirna.table[, 1] <- as.character(actual.file[["cross-mapped"]])
            }
            setTxtProgressBar(pb, count)
        }
        close(pb)

        assign("data.mirna.table", data.mirna.table, envir = get(envir_link))
        assign("cross.mapped.mirna.table", cross.mapped.mirna.table, envir = get(envir_link))
    }
    # code ####
    #create env if not created yet
    if (missing(env) && tumorData && !onlyFilter){
        assign(paste(toupper(tumor), toupper(dataBase), gsub(" ", "_", tolower(dataType)), "tumor_data", sep = "_"),
               new.env(parent=emptyenv()), envir = .GlobalEnv)
        envir_link <- paste(toupper(tumor), toupper(dataBase), gsub(" ", "_", tolower(dataType)), "tumor_data", sep = "_")
    } else if (missing(env) && !tumorData && !onlyFilter){
        assign(paste(toupper(tumor), toupper(dataBase), gsub(" ", "_", tolower(dataType)), "both_data", sep = "_"),
               new.env(parent=emptyenv()), envir = .GlobalEnv)
        envir_link <- paste(toupper(tumor), toupper(dataBase), gsub(" ", "_", tolower(dataType)), "both_data", sep = "_")
    } else if (missing(env) && onlyFilter){
        message("Please, before using 'onlyFilter' argument, insert the Environment name.")
    } else {
       envir_link <- deparse(substitute(env))
    }

    string_vars <- list(envir_link = get(envir_link))
    attr(envir_link, "name" ) = "Environment created by GDCRtools package, use its name in 'env' argument"

    dir.create(path = file.path(workDir, "GDCRtools", toupper(tumor), "Analyses"),
               showWarnings = FALSE)

    # mutation ####
    if (tolower(dataType) == "mutation"){
        if (tolower(dataBase) == "legacy"){
            #selecting platform file
            DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "mutation_data")
            manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"), select = c("file_name", "platform"))
            # all_maf.files <- dir(path = DIR, pattern = ".maf$")
            if (tolower(Platform) %in% "illumina ga"){
                maf.files <- as.character(manifest[[1]][manifest$platform == "Illumina GA"])
            } else if (tolower(Platform) %in% "illumina hiseq"){
                maf.files <- as.character(manifest[[1]][manifest$platform == "Illumina HiSeq"])
            }

            # maf.files <- intersect(all_maf.files, maf.files)

            #checking for unmatched data
            if (length(maf.files) > 0){
                message("Loading mutation file...")
            } else {
                stop(message(paste0("There is no data from ", tolower(Platform), " in your local storage!")))
            }

            #SomaticCancerAlterations:::.read_maf modified
            # to.remove <- c("center", "NCBI_Build", "Tumor_Validation_Allele1",
            #                "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
            #                "Match_Norm_Validation_Allele2", "Verification_Status" ,
            #                "Validation_Status", "Mutation_Status", "Sequencing_Phase",
            #                "Sequence_Source", "Validation_Method", "Score", "BAM_file",
            #                "Sequencer", "UniProt_Site", "GO_Biological_Process",
            #                "GO_Cellular_Component", "GO_Molecular_Function")

            message("\nIt will be open ", length(maf.files), " file(s)...\n")
            pb <- txtProgressBar(min = 0, max = length(maf.files), style = 3)
            count <- 0
            for (file in maf.files){
                count <- count + 1
                assign(file, read.delim(file = file.path(DIR, file),
                                        comment.char = "#", sep = "\t",
                                        header = TRUE,
                                        stringsAsFactors = FALSE),
                       envir = get(envir_link))
                message("\nThe file is load as '", file, "' in ", envir_link, " Environment.")

                # to.remove2 <- setdiff(colnames(file), to.remove)
                # bbb <<- colnames(file)
                # to.remove2 <<- to.remove2
                # print(head(file))
                # file <- file[, to.remove2]

                if (saveData) {
                    write.table(x = get(file, envir = get(envir_link)),
                                file = file.path(DIR, paste0(file, ".tsv")),
                                row.names = FALSE, quote = FALSE, sep = "\t")
                }

                setTxtProgressBar(pb, count)
            }
            close(pb)

            if (!missing(Name)){

                tmp <- apply(UCS_LEGACY_mutation_both_data[["bcgsc.ca_UCS.IlluminaHiSeq_DNASeq.1.somatic.maf"]][, 1:2], 1, function(x){
                    paste(x, collapse = "\\|")
                })
                gsub("\\s+", "", tmp)

                # insert gene finder here xxxx
            }

        } else if(tolower(dataBase) == "gdc"){
            DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gdc_mutation_data")
            manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"), select = c(11,12))
            manifest$submitter_id <- unlist(lapply(strsplit(x = manifest$submitter_id,
                                                            split = "-", perl = TRUE), "[[", 3))
            if (tolower(workflowType) %in% "all"){
                to.gunzip <- dir(path = DIR, pattern = "maf.gz$")
            } else if (tolower(workflowType) %in% "varscan"){
                to.gunzip <- as.character(manifest[[1]][manifest$submitter_id == "varscan"])
            } else if (tolower(workflowType) %in% "mutect"){
                to.gunzip <- as.character(manifest[[1]][manifest$submitter_id == "mutect"])
            } else if (tolower(workflowType) %in% "muse"){
                to.gunzip <- as.character(manifest[[1]][manifest$submitter_id == "muse"])
            } else if (tolower(workflowType) %in% "somaticsniper"){
                to.gunzip <- as.character(manifest[[1]][manifest$submitter_id == "somaticsniper"])
            }

            tryCatch(R.utils::gunzip(file.path(DIR, to.gunzip), remove = FALSE),
                     error = function(e){})
            import.maf <- gsub(".gz", "", to.gunzip)
            #openning files
            message("\nLoading mutation file...\n")
            maf <- read.delim(file = paste(DIR, import.maf, sep = "/"),
                               comment.char = "#", sep = "\t",
                               header = TRUE,
                               stringsAsFactors = FALSE)#, quote = FALSE)

            assign(import.maf, maf, envir = get(envir_link))
            message("\nThe file is open as '", import.maf, "' in ", envir_link, " environment.")
            if (saveData) {
                write.table(x = get(import.maf, envir = get(envir_link)),
                            file = file.path(DIR, paste0(import.maf, ".tsv")),
                            row.names = FALSE, quote = FALSE, sep = "\t")
            }
        }
    }

    # methylation ####
    if (tolower(dataType) == "methylation"){
        if (!onlyFilter){
            if (tolower(dataBase) == "legacy"){
                DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "methylation_data")
                Select_ref <- c(3:5)
                patient_selector <- "submitter_id"
            } else if (tolower(dataBase) == "gdc"){
                DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gdc_methylation_data")
                Select_ref <- c(3:11)
                patient_selector <- "cases"
            }

            manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"),
                                          select = c("file_name", patient_selector, "platform"))

            tmp <- !grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                          sep = "|"),
                          x = dir(path = DIR))

            meth.files.downloaded <- dir(path = DIR)[tmp]
            manifest <- manifest[manifest$file_name %in% meth.files.downloaded, ]

            #only meth data!!
            if (sum(tmp) > 0){
                if (tolower(Platform) %in% "illumina human methylation 450"){
                    meth.files <- manifest[manifest$platform == "Illumina Human Methylation 450", ]
                } else if (tolower(Platform) %in% "illumina human methylation 27"){
                    meth.files <- manifest[manifest$platform == "Illumina Human Methylation 27", ]
                }

                #checking for unmatched data
                if (nrow(meth.files) > 0){
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

            if(tumorData){
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"),
                                           meth.files[[2]])
                #not related the control samples
                meth.tumor.files <- as.character(meth.files[[1]][only.tumor.tissue])
                patients <- as.character(as.matrix(meth.files[only.tumor.tissue, 2]))
                assign("patients", patients, envir = get(envir_link))
                suppressWarnings(OPEN.meth(meth.tumor.files))
                #saving (should it happends?)
                # write.table(x = meth_table, file = "./methylation_table/meth_table_completed.csv",
                #           row.names = TRUE, quote = FALSE)
                ##filtering
                meth_table <- string_vars[["envir_link"]]$meth_table
                refererence_meth_table <- string_vars[["envir_link"]]$refererence_meth_table
                filter.meth.FUN(meth_table, refererence_meth_table)
                tumor_meth_table_filtered <- string_vars[["envir_link"]]$meth_table
                reference_table_filtered <- string_vars[["envir_link"]]$refererence_meth_table

                if (tolower(dataBase) == "gdc"){
                    reference_table_filtered[, 4] <- unlist(lapply(strsplit(x = reference_table_filtered[, 4],
                                                                       split = ";", perl = TRUE), "[", 1))
                }

                assign("methylation_tumor_filtered", tumor_meth_table_filtered,
                       envir = get(envir_link))
                assign("reference_table_filtered", reference_table_filtered,
                       envir = get(envir_link))

                if (saveData){
                    message("\nSaving data, this could take a while...\n")
                    #saving
                    write.table(x = cbind(refererence_meth_table, tumor_meth_table_filtered),
                                file = file.path(DIR, "methylation_tumor_filtered.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }

                if (!missing(Name)){
                    select.gene <- grepl(pattern = Name,
                                        x = reference_table_filtered$Gene_Symbol)
                    if (sum(select.gene) == 0){
                        message(cat(paste0(Name, " not found...\n\n"),
                                    "Please insert a valid gene symbol!!"))
                    }
                    assign(paste0("reference_table_filtered_selected_", Name),
                           reference_table_filtered[select.gene, ],
                           envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    write.table(x = cbind(reference_table_filtered[select.gene, ],
                                          tumor_meth_table_filtered[select.gene, ]),
                                file = file.path(DIR,
                                                 paste0("methylation_tumor_filtered_selected_",
                                                        Name,
                                                        ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")

                    assign(paste0("methylation_tumor_filtered_selected_", Name),
                           tumor_meth_table_filtered[select.gene, ],
                           envir = get(envir_link))
                    assign("Name", Name, envir = get(envir_link))
                }

            } else {
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                not.tumor.regex <- paste(10:19, collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"), meth.files)
                not.tumor.tissue <- grepl(pattern = paste0("-", not.tumor.regex, "[aAbB]-"), meth.files)
                #not related the control samples

                if (sum(not.tumor.tissue) == 0) {
                    message("\nThere is no 'normal' data available to ",
                            tumor, " tumor!!\n\n")
                    stop()
                }

                ##filtering
                #tumor
                meth.tumor.files <- meth.files[only.tumor.tissue]
                patients <- as.character(as.matrix(meth.files[only.tumor.tissue, 2]))
                patients.tumor <- patients
                suppressWarnings(OPEN.meth(meth.tumor.files))
                meth_table <- string_vars[["envir_link"]]$meth_table
                refererence_meth_table <- string_vars[["envir_link"]]$refererence_meth_table
                filter.meth.FUN(meth_table, refererence_meth_table)
                tumor_meth_table_filtered <- string_vars[["envir_link"]]$meth_table
                reference_table_filtered <- string_vars[["envir_link"]]$refererence_meth_table

                ##filtering
                #not tumor
                meth.not.tumor.files <- string_vars[["envir_link"]]$meth.files[not.tumor.tissue]
                patients <- as.character(as.matrix(meth.files[not.tumor.tissue, 2]))
                patients.not.tumor <- patients
                suppressWarnings(OPEN.meth(meth.tumor.files))
                meth_table <- string_vars[["envir_link"]]$meth_table
                refererence_meth_table <- string_vars[["envir_link"]]$refererence_meth_table
                filter.meth.FUN(meth_table, refererence_meth_table)
                not_tumor_meth_table_filtered <- string_vars[["envir_link"]]$meth_table

                assign("methylation_tumor_filtered", tumor_meth_table_filtered,
                       envir = get(envir_link))
                assign("methylation_not_tumor_filtered", not_tumor_meth_table_filtered,
                       envir = get(envir_link))
                assign("reference_table_filtered", reference_table_filtered,
                       envir = get(envir_link))

                #saving
                if (saveData){
                    message("\nSaving data...\n")
                    write.table(x = cbind(refererence_meth_table, tumor_meth_table_filtered),
                                file = file.path(DIR, "tumor_meth_table_completed_filtered.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                    write.table(x = cbind(refererence_meth_table, not_tumor_meth_table_filtered),
                                file = file.path(DIR, "not_tumor_meth_table_completed_filtered.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }

                if (!missing(Name)){
                    select.gene <- grepl(pattern = Name,
                                         x = reference_table_filtered$Gene_Symbol)
                    if (sum(select.gene) == 0){
                        message(cat(paste0(Name, " not found...\n\n"),
                                    "Please insert a valid gene symbol!!"))
                        suppressWarnings(remove(meth_table, refererence_meth_table,
                                                envir = get(envir_link)))
                        stop()
                    }
                    assign(paste0("reference_table_filtered_selected_", Name),
                           reference_table_filtered[select.gene, ],
                           envir = get(envir_link))

                    assign(paste0("methylation_tumor_filtered_selected_", Name),
                           tumor_meth_table_filtered[select.gene, ],
                           envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    write.table(x = cbind(reference_table_filtered[select.gene, ],
                                          tumor_meth_table_filtered[select.gene, ]),
                                file = file.path(DIR,
                                                 paste0("methylation_tumor_filtered_selected_",
                                                        Name,
                                                        ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")

                    # not tumor
                    assign(paste0("methylation_not_tumor_filtered_selected_", Name),
                           not_tumor_meth_table_filtered[select.gene, ],
                           envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    write.table(x = cbind(reference_table_filtered[select.gene, ],
                                          not_tumor_meth_table_filtered[select.gene, ]),
                                file = file.path(DIR,
                                                 paste0("methylation_not_tumor_filtered_selected_",
                                                        Name,
                                                        ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")

                    assign("Name", Name, envir = get(envir_link))
                    # assign("Name.e", paste0(Name, "_meth"), envir = get(envir_link))
                }
            }
            #removing useless variables
            suppressWarnings(remove(meth_table, refererence_meth_table,
                                    envir = get(envir_link)))
            invisible(gc())
        } else {
            if (!missing(Name)){
                select.gene <- grepl(pattern = Name,
                                     x = string_vars[["envir_link"]]$reference_table_filtered$Gene_Symbol)
                if (sum(select.gene) == 0){
                    message(cat(paste0(Name, " not found...\n\n"),
                                "Please insert a valid gene symbol!!"))
                    stop()
                }
                assign(paste0("reference_table_filtered_selected_", Name),
                       string_vars[["envir_link"]]$reference_table_filtered[select.gene, ],
                       envir = get(envir_link))

                assign(paste0("methylation_tumor_filtered_selected_", Name),
                       string_vars[["envir_link"]]$methylation_tumor_filtered[select.gene, ],
                       envir = get(envir_link))

                message("\nSaving required data...\n")
                #saving
                write.table(x = cbind(string_vars[["envir_link"]]$reference_table_filtered[select.gene, ],
                                      string_vars[["envir_link"]]$methylation_tumor_filtered[select.gene, ]),
                            file = file.path(DIR,
                                             paste0("methylation_tumor_filtered_selected_",
                                                    Name,
                                                    ".tsv")),
                            row.names = TRUE, quote = FALSE, sep = "\t")

                # not tumor
                if (!tumorData){
                    assign(paste0("methylation_not_tumor_filtered_selected_", Name),
                           string_vars[["envir_link"]]$methylation_not_tumor_filtered[select.gene, ],
                           envir = get(envir_link))

                    message("\nSaving required data...\n")
                    #saving
                    write.table(x = cbind(string_vars[["envir_link"]]$reference_table_filtered[select.gene, ],
                                          string_vars[["envir_link"]]$methylation_not_tumor_filtered[select.gene, ]),
                                file = file.path(DIR,
                                                 paste0("methylation_tumor_filtered_selected_",
                                                        Name,
                                                        ".tsv")),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            }
        }
    }

    # gene and isoform ####
    if ("gene" == tolower(dataType) || "isoform" == tolower(dataType)){
        if(!onlyFilter){
            if (tolower(dataType) == "gene"){
                if (tolower(dataBase) == "legacy"){
                    DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gene_data")
                    #prepare data selectors
                    gene.isoform()
                    message(tumor)
                    OPEN.genexpress(files = get(envir_link)$gene.files, Var = "gene.files")
                    patients <- string_vars[["envir_link"]]$gene.files$submitter_id
                    colnames(string_vars[["envir_link"]]$completed.table1) <- patients

                    patients_short <- unname(sapply(patients, function(w){
                        paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                    }))
                    # duplicate fix
                    if (length(patients) != length(unique(patients_short))){
                        coluns_2rename <- numeric()
                        names_2rename <- character()
                        for (Patients in unique(patients_short)){
                            selector <- Patients == patients_short
                            if (sum(selector) > 1) {
                                coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                                names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                            }
                        }
                        colnames(string_vars[["envir_link"]]$completed.table1)[coluns_2rename] <- names_2rename
                        colnames(string_vars[["envir_link"]]$completed.table1)[-coluns_2rename] <- unique(patients_short)
                    } else {
                        colnames(string_vars[["envir_link"]]$completed.table1) <- patients_short
                    }

                    if (tumorData){

                        if (!missing(Name)){
                        # gene.selector.final <- completed.table1[, lapply(.SD,
                        #                                            function(x) grepl(paste("(", "\\|", Name, ")", "$", sep = ""),
                        #                                                                    x, perl=TRUE))]
                            if (normalization){
                                # assign("Name.e", Name, envir = get(envir_link))
                                to.google <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                    paste0("(", Name, "\\|", ")"))
                                gene.selector.final <- grepl(to.google,
                                                             rownames(string_vars[["envir_link"]]$completed.table1), perl=TRUE)
                                genes.1 <- string_vars[["envir_link"]]$completed.table1[gene.selector.final, , drop = FALSE]

                                # genes.1 <- completed.table1[rowSums(gene.selector.final)>0, ]
                                if (nrow(genes.1) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1 <- t(genes.1)
                                # genes.transpo1 <- genes.transpo1[-1, ]
                                # genes.transpo1 <- as.data.frame(genes.transpo1)
                                colnames(genes.transpo1) <- Name

                                assign(paste0("gene_tumor_normalized_selected_", toupper(Name)), genes.transpo1, envir = get(envir_link))

                            } else {
                                to.google <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                    paste0("(", Name, "\\|", ")"))
                                gene.selector.final <- grepl(to.google,
                                                             rownames(string_vars[["envir_link"]]$completed.table1), perl=TRUE)
                                genes.1 <- string_vars[["envir_link"]]$completed.table1[gene.selector.final, , drop = FALSE]

                                # genes.1 <- completed.table1[rowSums(gene.selector.final)>0, ]
                                if (nrow(genes.1) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1 <- t(genes.1)
                                # genes.transpo1 <- genes.transpo1[-1, ]
                                # genes.transpo1 <- as.data.frame(genes.transpo1)
                                colnames(genes.transpo1) <- Name

                                assign(paste0("gene_tumor_not_normalized_selected_", toupper(Name)), genes.transpo1, envir = get(envir_link))

                            }
                        }
                    } else {
                        # #tumor and not tumor data
                        # DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gene_data")
                        # #prepare data selectors
                        # gene.isoform()
                        # message(tumor)
                        # #tumor
                        # OPEN.genexpress(files = string_vars[["envir_link"]]$gene.files, Var = "gene.files")
                        # patients <- string_vars[["envir_link"]]$gene.files$submitter_id
                        # colnames(string_vars[["envir_link"]]$completed.table1) <- patients

                        if (!missing(Name)){
                            if (normalization){
                                # assign("Name.e", Name, envir = get(envir_link))
                                to.google <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                    paste0("(", Name, "\\|", ")"))
                                gene.selector.final <- grepl(to.google,
                                                             rownames(string_vars[["envir_link"]]$completed.table1), perl=TRUE)
                                genes.1 <- string_vars[["envir_link"]]$completed.table1[gene.selector.final, , drop = FALSE]
                                if (nrow(genes.1) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1 <- t(genes.1)
                                colnames(genes.transpo1) <- Name

                                assign(paste0("gene_tumor_normalized_selected_", toupper(Name)), genes.transpo1, envir = get(envir_link))
                            } else {
                                to.google <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                    paste0("(", Name, "\\|", ")"))
                                gene.selector.final <- grepl(to.google,
                                                             rownames(string_vars[["envir_link"]]$completed.table1), perl=TRUE)
                                genes.1 <- string_vars[["envir_link"]]$completed.table1[gene.selector.final, , drop = FALSE]
                                if (nrow(genes.1) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1 <- t(genes.1)
                                colnames(genes.transpo1) <- Name

                                assign(paste0("gene_tumor_not_normalized_selected_", toupper(Name)), genes.transpo1, envir = get(envir_link))

                            }
                        }

                        #not tumor
                        if(length(string_vars[["envir_link"]]$gene.files.not.tumor$file_name) == 0) {
                            message("There is no 'Normal' data in ", tumor,
                                    " tumor folder! Please, rerun after set 'tumorData = TRUE'.",
                                    "\n\n")
                            stop()
                        }
                        OPEN.genexpress(files = string_vars[["envir_link"]]$gene.files.not.tumor, Var = "gene.files.not.tumor")
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$completed.table1.not.tumor
                        patients.not.tumor <- string_vars[["envir_link"]]$gene.files.not.tumor$submitter_id
                        colnames(completed.table1.not.tumor) <- patients.not.tumor

                        patients_short <- unname(sapply(patients.not.tumor, function(w){
                            paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                        }))
                        # duplicate fix
                        if (length(patients.not.tumor) != length(unique(patients_short))){
                            coluns_2rename <- numeric()
                            names_2rename <- character()
                            for (Patients in unique(patients_short)){
                                selector <- Patients == patients_short
                                if (sum(selector) > 1) {
                                    coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                                    names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                                }
                            }
                            colnames(completed.table1.not.tumor)[coluns_2rename] <- names_2rename
                            colnames(completed.table1.not.tumor)[-coluns_2rename] <- unique(patients_short)
                        } else {
                            colnames(completed.table1.not.tumor) <- patients_short
                        }

                        rownames(completed.table1.not.tumor) <- string_vars[["envir_link"]]$name.row.not.tumor

                        if (!missing(Name)){
                            if (normalization) {
                                to.google.not.tumor <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                              paste0("(", Name, "\\|", ")"))
                                gene.selector.final.not.tumor <- grepl(to.google.not.tumor,
                                                                       rownames(string_vars[["envir_link"]]$completed.table1), perl=TRUE)
                                genes.1.not.tumor <- completed.table1.not.tumor[gene.selector.final.not.tumor, , drop = FALSE]

                                if (nrow(genes.1.not.tumor) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1.not.tumor <- t(genes.1.not.tumor)
                                colnames(genes.transpo1.not.tumor) <- Name

                                FILTERED_RESULTS.not.tumor <- genes.transpo1.not.tumor
                                assign(paste0("gene_not_tumor_normalized_selected_", toupper(Name)), genes.transpo1.not.tumor, envir = get(envir_link))
                            } else {
                                to.google.not.tumor <- ifelse(is.numeric(Name), paste0("(", "\\|", Name, ")$"),
                                                              paste0("(", Name, "\\|", ")"))
                                gene.selector.final.not.tumor <- grepl(to.google.not.tumor,
                                                                       rownames(string_vars[["envir_link"]]$completed.table1), perl=TRUE)
                                genes.1.not.tumor <- completed.table1.not.tumor[gene.selector.final.not.tumor, , drop = FALSE]

                                if (nrow(genes.1.not.tumor) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1.not.tumor <- t(genes.1.not.tumor)
                                colnames(genes.transpo1.not.tumor) <- Name

                                FILTERED_RESULTS.not.tumor <- genes.transpo1.not.tumor
                                assign(paste0("gene_not_tumor_not_normalized_selected_", Name), genes.transpo1.not.tumor, envir = get(envir_link))

                            }
                        }
                    }

                } else if (tolower(dataBase) == "gdc"){
                    DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gdc_gene_data")
                    #selecting specific HTSeq data
                    if (tolower(HTSeq) == "counts"){
                        to.gunzip <- dir(path = DIR, pattern = "counts.gz$")
                    } else {
                        to.gunzip <- dir(path = DIR,
                                         pattern = paste0(toupper(HTSeq), ".txt.gz$"))
                    }

                    tmp <- sapply(paste(DIR, to.gunzip, sep = "/"),
                                  function(x){R.utils::gunzip(x, remove = FALSE, overwrite=TRUE)})

                    #open manifest
                    codes <- dir(path = DIR, pattern = ".sdrf$")
                    codigos <- read.table(file.path(DIR, codes[1]), stringsAsFactors = FALSE,
                                          header = TRUE, sep = "\t")
                    codigos <- codigos[, c("file_name", "cases")]
                    codigos <- codigos[codigos$file_name %in% to.gunzip, ]
                    tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"),
                                         collapse = "[aAbB]-|-")
                    #file names to concatenate
                    codigos$file_name <- gsub(".gz", "", codigos$file_name)

                    if (tumorData) {
                        codigos <- codigos[grepl(paste0("-", tumor.regex, "-"),
                                                 x = codigos$cases), ]
                        remove.duplicated <- duplicated(codigos$file_name)
                        codigos <- codigos[!remove.duplicated, ]
                        rownames(codigos) <- codigos$file_name

                        #concatenate and finishing table
                        OPEN.genexpress(files = codigos, Var = "gene.files")
                        rownames(string_vars[["envir_link"]]$completed.table1) <- string_vars[["envir_link"]]$name.row
                        # rownames(string_vars[["envir_link"]]$completed.table1) <- string_vars[["envir_link"]]$completed.table1[, 1]
                        #completed.table1 <- completed.table1[, -1]
                        patients <- codigos$cases
                        colnames(string_vars[["envir_link"]]$completed.table1) <- patients

                        patients_short <- unname(sapply(patients, function(w){
                            paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                        }))
                        # duplicate fix
                        if (length(patients) != length(unique(patients_short))){
                            coluns_2rename <- numeric()
                            names_2rename <- character()
                            for (Patients in unique(patients_short)){
                                selector <- Patients == patients_short
                                if (sum(selector) > 1) {
                                    coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                                    names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                                }
                            }
                            colnames(string_vars[["envir_link"]]$completed.table1)[coluns_2rename] <- names_2rename
                            colnames(string_vars[["envir_link"]]$completed.table1)[-coluns_2rename] <- unique(patients_short)
                        } else {
                            colnames(string_vars[["envir_link"]]$completed.table1) <- patients_short
                        }

                        #saving completed.table
                        write.table(string_vars[["envir_link"]]$completed.table1,
                                    paste0(DIR, "/concatenate_", HTSeq, ".tsv"), sep = "\t")


                        if (!missing(Name)){
                            if (normalization) {
                                genes.1 <- string_vars[["envir_link"]]$completed.table1[Name, , drop = FALSE]

                                # genes.1 <- string_vars[["envir_link"]]$completed.table1[rowSums(gene.selector.final)>0, ]
                                if (nrow(genes.1) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1 <- t(genes.1)
                                # genes.transpo1 <- genes.transpo1[-1, ]
                                # genes.transpo1 <- as.data.frame(genes.transpo1)
                                colnames(genes.transpo1) <- Name

                                FILTERED_RESULTS <- genes.transpo1
                                assign(paste0("gene_tumor_normalized_selected_", toupper(Name)), genes.transpo1, envir = get(envir_link))
                            } else {
                                genes.1 <- string_vars[["envir_link"]]$completed.table1[Name, , drop = FALSE]

                                # genes.1 <- string_vars[["envir_link"]]$completed.table1[rowSums(gene.selector.final)>0, ]
                                if (nrow(genes.1) == 0) {
                                    stop("Gene not found")
                                }
                                genes.transpo1 <- t(genes.1)
                                # genes.transpo1 <- genes.transpo1[-1, ]
                                # genes.transpo1 <- as.data.frame(genes.transpo1)
                                colnames(genes.transpo1) <- Name

                                FILTERED_RESULTS <- genes.transpo1
                                assign(paste0("gene_tumor_not_normalized_selected_", toupper(Name)), genes.transpo1, envir = get(envir_link))
                            }
                        }
                    } else {
                        not.tumor.regex <- paste(10:19, collapse = "[aAbB]-|-")

                        codigos.not.tumor <- codigos[grepl(paste0("-", not.tumor.regex, "-"),
                                                           x = codigos$cases), ]
                        remove.duplicated.not.tumor <- duplicated(codigos.not.tumor$file_name)
                        codigos.not.tumor <- codigos.not.tumor[!remove.duplicated.not.tumor, ]

                        rownames(codigos.not.tumor) <- codigos.not.tumor$file_name

                        #concatenate
                        if(length(codigos.not.tumor$file_name) == 0) {
                            message("There is no 'Normal' data in ", tumor,
                                    " tumor folder! Please, rerun after set 'tumorData = TRUE'.",
                                    "\n\n")
                            stop()
                        }
                        OPEN.genexpress(files = codigos.not.tumor, Var = "gene.files.not.tumor")
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$completed.table1.not.tumor
                        patients.not.tumor <- string_vars[["envir_link"]]$gene.files.not.tumor$submitter_id
                        colnames(completed.table1.not.tumor) <- patients.not.tumor

                        patients_short <- unname(sapply(patients.not.tumor, function(w){
                            paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                        }))
                        # duplicate fix
                        if (length(patients.not.tumor) != length(unique(patients_short))){
                            coluns_2rename <- numeric()
                            names_2rename <- character()
                            for (Patients in unique(patients_short)){
                                selector <- Patients == patients_short
                                if (sum(selector) > 1) {
                                    coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                                    names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                                }
                            }
                            colnames(completed.table1.not.tumor)[coluns_2rename] <- names_2rename
                            colnames(completed.table1.not.tumor)[-coluns_2rename] <- unique(patients_short)
                        } else {
                            colnames(completed.table1.not.tumor) <- patients_short
                        }
                    }
                    #remove no longer need files
                    suppressWarnings(file.remove(file.path(DIR, codigos$file_name)))
                    if (tolower(HTSeq) != "counts"){
                        assign("HTSeq_normalized", toupper(HTSeq), envir = get(envir_link))
                    }
                }

            }

            if (tolower(dataType) == "isoform"){
                DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "isoform_data")
                #prepare data selectors
                gene.isoform()
                message(tumor)
                OPEN.genexpress(files = string_vars[["envir_link"]]$isoform.files, Var = "isoform.files")
                completed.table1.tumor <- string_vars[["envir_link"]]$completed.table1
                patients <- string_vars[["envir_link"]]$isoform.files$submitter_id
                colnames(completed.table1.tumor) <- patients

                patients_short <- unname(sapply(patients, function(w){
                    paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                }))
                # duplicate fix
                if (length(patients) != length(unique(patients_short))){
                    coluns_2rename <- numeric()
                    names_2rename <- character()
                    for (Patients in unique(patients_short)){
                        selector <- Patients == patients_short
                        if (sum(selector) > 1) {
                            coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                            names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                        }
                    }
                    colnames(completed.table1.tumor)[coluns_2rename] <- names_2rename
                    colnames(completed.table1.tumor)[-coluns_2rename] <- unique(patients_short)
                } else {
                    colnames(completed.table1.tumor) <- patients_short
                }


                if (!tumorData) {
                    if (nrow(string_vars[["envir_link"]]$isoform.files.not.tumor) == 0){
                        message("There is no 'Normal' data in ", tumor,
                                " tumor folder! Please, rerun after set 'tumorData = TRUE'.",
                                "\n\n")
                        stop()
                    }
                    OPEN.genexpress(files = string_vars[["envir_link"]]$isoform.files.not.tumor, Var = "isoform.files.not.tumor")
                    patients.not.tumor <- string_vars[["envir_link"]]$isoform.files.not.tumor$submitter_id
                    colnames(string_vars[["envir_link"]]$completed.table1.not.tumor) <- patients.not.tumor

                    patients_short <- unname(sapply(patients.not.tumor, function(w){
                        paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                    }))
                    # duplicate fix
                    if (length(patients.not.tumor) != length(unique(patients_short))){
                        coluns_2rename <- numeric()
                        names_2rename <- character()
                        for (Patients in unique(patients_short)){
                            selector <- Patients == patients_short
                            if (sum(selector) > 1) {
                                coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                                names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                            }
                        }
                        colnames(string_vars[["envir_link"]]$completed.table1.not.tumor)[coluns_2rename] <- names_2rename
                        colnames(string_vars[["envir_link"]]$completed.table1.not.tumor)[-coluns_2rename] <- unique(patients_short)
                    } else {
                        colnames(string_vars[["envir_link"]]$completed.table1.not.tumor) <- patients.not.tumor
                    }

                    if (!missing(Name)){
                        if (normalization){
                            # assign("Name.e", Name, envir = get(envir_link))

                            #tumor
                            Isoforms.1 <- completed.table1.tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms.1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms.transpo1 <- t(Isoforms.1)
                            Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
                            colnames(Isoforms.transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms.transpo1
                            assign(paste0("isoform_tumor_normalized_selected_", toupper(Name)), Isoforms.transpo1, envir = get(envir_link))


                            #not tumor
                            Isoforms.1 <- string_vars[["envir_link"]]$completed.table1.not.tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms.1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms.transpo1 <- t(Isoforms.1)
                            Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
                            colnames(Isoforms.transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms.transpo1
                            assign(paste0("isoform_not_tumor_normalized_selected_", toupper(Name)), Isoforms.transpo1, envir = get(envir_link))

                        } else {
                            #tumor
                            Isoforms.1 <- completed.table1.tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms.1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms.transpo1 <- t(Isoforms.1)
                            Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
                            colnames(Isoforms.transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms.transpo1
                            assign(paste0("isoform_tumor_not_normalized_selected_", toupper(Name)), Isoforms.transpo1, envir = get(envir_link))


                            #not tumor
                            Isoforms.1 <- string_vars[["envir_link"]]$completed.table1.not.tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms.1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms.transpo1 <- t(Isoforms.1)
                            Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
                            colnames(Isoforms.transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms.transpo1
                            assign(paste0("isoform_not_tumor_not_normalized_selected_", toupper(Name)), Isoforms.transpo1, envir = get(envir_link))

                        }
                    }
                } else {
                    if (!missing(Name)){
                        if (normalization){
                            # assign("Name.e", Name, envir = get(envir_link))
                            Isoforms.1 <- completed.table1.tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms.1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms.transpo1 <- t(Isoforms.1)
                            Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
                            colnames(Isoforms.transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms.transpo1
                            assign(paste0("isoform_tumor_normalized_selected_", toupper(Name)), Isoforms.transpo1, envir = get(envir_link))

                        } else {
                            Isoforms.1 <- completed.table1.tumor[Name, , drop = FALSE]

                            if (nrow(Isoforms.1) == 0) {
                                stop("Isoform not found!!!")
                            }
                            Isoforms.transpo1 <- t(Isoforms.1)
                            Isoforms.transpo1 <- as.data.frame(Isoforms.transpo1)
                            colnames(Isoforms.transpo1) <- Name

                            FILTERED_RESULTS <- Isoforms.transpo1
                            assign(paste0("isoform_tumor_not_normalized_selected_", toupper(Name)), Isoforms.transpo1, envir = get(envir_link))

                        }
                    }
                }
            }

            #Getting and export the final results
            if (tolower(dataType) == "gene"){
                if (tumorData){
                    if (normalization || (tolower(HTSeq) != "counts" && !is.null(HTSeq))){
                        assign("gene_tumor_normalized", string_vars[["envir_link"]]$completed.table1, envir = get(envir_link))
                    } else {
                        assign("gene_tumor_not_normalized", string_vars[["envir_link"]]$completed.table1, envir = get(envir_link))
                    }
                    if (saveData){
                        message("\nSaving your data...\n")
                        write.table(string_vars[["envir_link"]]$completed.table1,
                                    paste0(DIR, "/", tumor, "_tumor_data.tsv"), sep = "\t")
                    }
                    assign("patients", patients, envir = get(envir_link))

                } else {
                    #tumor
                    if (normalization || (tolower(HTSeq) != "counts" && !is.null(HTSeq))){
                        assign("gene_tumor_normalized", string_vars[["envir_link"]]$completed.table1, envir = get(envir_link))
                    } else {
                        assign("gene_tumor_not_normalized", string_vars[["envir_link"]]$completed.table1, envir = get(envir_link))

                    }
                    #not tumor
                    row.names(string_vars[["envir_link"]]$completed.table1.not.tumor) <- string_vars[["envir_link"]]$name.row.not.tumor
                    if (normalization || tolower(HTSeq) != "counts"){
                        assign("gene_not_tumor_normalized", string_vars[["envir_link"]]$completed.table1.not.tumor, envir = get(envir_link))
                    } else {
                        assign("gene_not_tumor_not_normalized", string_vars[["envir_link"]]$completed.table1.not.tumor, envir = get(envir_link))
                    }
                    if (saveData){
                        message("\nSaving your data...\n")
                        write.table(string_vars[["envir_link"]]$completed.table1,
                                    paste0(DIR, "/", tumor, "_tumor_data.tsv"), sep = "\t")
                        write.table(string_vars[["envir_link"]]$completed.table1.not.tumor,
                                    paste0(DIR, "/", tumor, "_not_tumor_data.tsv"), sep = "\t")
                    }
                    assign("patients", patients, envir = get(envir_link))
                    assign("patients.not.tumor", patients.not.tumor, envir = get(envir_link))
                }
            } else if (tolower(dataType) == "isoform"){
                if (tumorData){
                    if (normalization){
                        assign("isoform_tumor_normalized", completed.table1.tumor, envir = get(envir_link))
                    } else {
                        assign("isoform_tumor_not_normalized", completed.table1.tumor, envir = get(envir_link))
                    }
                    if (saveData){
                        message("\nSaving your data...\n")
                        write.table(completed.table1.tumor,
                                    paste0(DIR, "/", tumor, "_tumor_data.tsv"), sep = "\t")
                    }
                    assign("patients", patients, envir = get(envir_link))

                } else {
                    #tumor
                    if (normalization){
                        assign("isoform_tumor_normalized", completed.table1.tumor, envir = get(envir_link))
                    } else {
                        assign("isoform_tumor_not_normalized", completed.table1.tumor, envir = get(envir_link))

                    }
                    #not tumor
                    row.names(string_vars[["envir_link"]]$completed.table1.not.tumor) <- string_vars[["envir_link"]]$name.row.not.tumor
                    if (normalization){
                        assign("isoform_not_tumor_normalized", string_vars[["envir_link"]]$completed.table1.not.tumor, envir = get(envir_link))
                    } else {
                        assign("isoform_not_tumor_not_normalized", string_vars[["envir_link"]]$completed.table1.not.tumor, envir = get(envir_link))
                    }
                    if (saveData){
                        message("\nSaving your data...\n")
                        write.table(completed.table1.tumor,
                                    paste0(DIR, "/", tumor, "_tumor_data.tsv"), sep = "\t")
                        write.table(string_vars[["envir_link"]]$completed.table1.not.tumor,
                                    paste0(DIR, "/", tumor, "_not_tumor_data.tsv"), sep = "\t")
                    }
                    assign("patients", patients, envir = get(envir_link))
                    assign("patients.not.tumor", patients.not.tumor, envir = get(envir_link))
                }
            }

            suppressWarnings(remove(genes.normalized, isoforms.normalized,
                                    genes.not.normalized, isoforms.not.normalized,
                                    genes.normalized.not.tumor, isoforms.normalized.not.tumor,
                                    genes.not.normalized.not.tumor, isoforms.not.normalized.not.tumor,
                                    envir = get(envir_link)))
            suppressWarnings(remove(completed.table1, completed.table1.not.tumor,
                                    codigos.not.tumor, isoform.files,
                                    isoform.files.not.tumor, gene.files,
                                    codigos.tumor, gene.files.not.tumor,
                                    envir = get(envir_link)))
            invisible(gc(verbose = FALSE))
            # if (nrow(string_vars[["envir_link"]]$gene.files) == 0) {
            #     remove(gene.files, envir = get(envir_link))
            # } else if (nrow(string_vars[["envir_link"]]$isoform.files) == 0){
            #     remove(isoform.files, envir = get(envir_link))
            # }
        } else {
            if (tolower(dataType) == "gene" && tolower(dataBase) == "legacy"){

                if (tumorData){
                    if (normalization){
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        gene_finder(completed.table1.tumor, variable_name = "gene_tumor_normalized_selected_")
                    } else {
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        gene_finder(completed.table1.tumor, variable_name = "gene_tumor_not_normalized_selected_")
                    }
                } else {
                    #tumor
                    if (normalization){
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        gene_finder(completed.table1.tumor, variable_name = "gene_tumor_normalized_selected_")
                    } else {
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        gene_finder(completed.table1.tumor, variable_name = "gene_tumor_not_normalized_selected_")
                    }
                    #not tumor
                    if (normalization){
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$gene_not_tumor_normalized
                        gene_finder(completed.table1.not.tumor, variable_name = "gene_not_tumor_normalized_selected_")
                    } else {
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$gene_not_tumor_not_normalized
                        gene_finder(completed.table1.not.tumor, variable_name = "gene_not_tumor_not_normalized_selected_")
                    }
                }

            } else if (tolower(dataType) == "gene" && tolower(dataBase) == "gdc") {
                if (tumorData){
                    if (normalization){
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "gene_tumor_normalized_selected_")
                    } else {
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "gene_tumor_not_normalized_selected_")
                    }
                } else {
                    #tumor
                    if (normalization){
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "gene_tumor_normalized_selected_")
                    } else {
                        completed.table1.tumor <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "gene_tumor_not_normalized_selected_")
                    }
                    #not tumor
                    if (normalization){
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$gene_not_tumor_normalized
                        isoform_finder(completed.table1.not.tumor, variable_name = "gene_not_tumor_normalized_selected_")
                    } else {
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$gene_not_tumor_not_normalized
                        isoform_finder(completed.table1.not.tumor, variable_name = "gene_not_tumor_not_normalized_selected_")
                    }
                }

            } else if (tolower(dataType) == "isoform"){

                if (tumorData){
                    if (normalization){
                        completed.table1.tumor <- string_vars[["envir_link"]]$isoform_tumor_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "isoform_tumor_normalized_selected_")
                    } else {
                        completed.table1.tumor <- string_vars[["envir_link"]]$isoform_tumor_not_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "isoform_tumor_not_normalized_selected_")
                    }
                } else {
                    #tumor
                    if (normalization){
                        completed.table1.tumor <- string_vars[["envir_link"]]$isoform_tumor_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "isoform_tumor_normalized_selected_")
                    } else {
                        completed.table1.tumor <- string_vars[["envir_link"]]$isoform_tumor_not_normalized
                        isoform_finder(completed.table1.tumor, variable_name = "isoform_tumor_not_normalized_selected_")
                    }
                    #not tumor
                    if (normalization){
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$isoform_not_tumor_normalized
                        isoform_finder(completed.table1.not.tumor, variable_name = "isoform_not_tumor_normalized_selected_")
                    } else {
                        completed.table1.not.tumor <- string_vars[["envir_link"]]$isoform_not_tumor_not_normalized
                        isoform_finder(completed.table1.not.tumor, variable_name = "isoform_not_tumor_not_normalized_selected_")
                    }
                }
            }
        }
    }

    # clinical ####
    if ("clinical" == tolower(dataType)){
        if (tolower(dataBase) == "legacy"){
            #selecting platform file
            DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "clinical_data")
            patient_file <- dir(path = DIR, pattern = "clinical_patient")

            clinical_df <- data.table::fread(input = file.path(DIR, patient_file))

            clinical_df <- as.data.frame(clinical_df[-c(1,2),])
            rownames(clinical_df) <- clinical_df[, 2]
            clinical_df <- clinical_df[, -2]

            # clinical_df <- clinical_df[, c()] #which colunms to use?

            assign("clinical_df", clinical_df, envir = get(envir_link))

        } else if(tolower(dataBase) == "gdc"){
            DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gdc_biospecimen_data")
            patient_files <- dir(path = DIR, pattern = "xml$")

            pb <- txtProgressBar(min = 0, max = length(patient_files), style = 3)
            count <- 0
            #concatenating all the XML files at once
            for (i in file.path(DIR, patient_files)){
                data_df <- XML::xmlToDataFrame(file.path(DIR, i))
                data_df <- na.omit(reshape2::melt(t(data_df)))
                rownames(data_df) <- data_df[, 1, drop = TRUE]
                data_df <- t(data_df[, 3, drop = FALSE])
                if(count == 0){
                    clinical_df <- as.data.frame(matrix(ncol = ncol(data_df),
                                                        nrow = length(patient_files)))
                    RowNames <- character(length = length(patient_files))
                    colnames(clinical_df) <- colnames(data_df)
                    clinical_df[1, ] <- data_df
                    RowNames[1] <- data_df[, "bcr_patient_barcode"]
                } else if (count <= length(file.path(DIR, patient_files))){
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
    if ("protein" == tolower(dataType)){
        if (!onlyFilter) {
            if (tolower(dataBase) == "gdc") {
                stop(message("\nThrere is no protein expression data in GDC data base!!",
                             "\nPlease use 'legacy' data base"))
            }
            DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "protein_data")

            if (!dir.exists(file.path(workDir, "GDCRtools", toupper(tumor), "mage_data"))) {
                message("Downloading magetab data...")
                suppressWarnings(download_gdc(dataType = "mage",
                                              dataBase = dataBase,
                                              tumor = tumor,
                                              workDir = workDir))
            }

            desing_array_DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "mage_data")

            to.gunzip <- file.path(desing_array_DIR,
                                   paste0("mdanderson.org_", toupper(tumor), ".MDA_RPPA_Core.mage-tab.1.1.0.tar.gz"))
            ## or, if you just want to extract the target file:
            untar(to.gunzip, exdir = desing_array_DIR)

            design_array <- data.table::fread(file.path(desing_array_DIR,
                                                        paste0("mdanderson.org_", toupper(tumor), ".MDA_RPPA_Core.mage-tab.1.1.0"),
                                                        paste0("mdanderson.org_", toupper(tumor), ".MDA_RPPA_Core.array_design.txt")),
                                              select = c(6,7))

            design_array <- unique(design_array)

            rownames(design_array) <- design_array[[2]]

            patient_files <- dir(path = DIR, pattern = "protein_expression")

            pb <- txtProgressBar(min = 0, max = length(patient_files), style = 3)
            count <- 0
            for (i in file.path(DIR, patient_files)){
                count <- count + 1
                if (count == 1){
                    RowNames <- character(length = length(patient_files))
                    patient_id <- character(length = length(patient_files))
                    pr <- data.table::fread(i)
                    uuid <- colnames(pr)[2]
                    protein_table <- matrix(ncol = length(patient_files),
                                            nrow = nrow(pr)-1)
                    colnames(protein_table) <- 1:ncol(protein_table)
                    patient_id[1] <- as.character(design_array[[1]][grep(uuid, design_array[[2]])])
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


            if(tumorData){
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"), colnames(protein_table))
                #not related the control samples
                patients <- as.character(colnames(protein_table)[only.tumor.tissue])
                protein_table <- protein_table[, patients]

                patients_short <- unname(sapply(patients, function(w){
                    paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
                }))
                # duplicate fix
                if (length(patients) != length(unique(patients_short))){
                    coluns_2rename <- numeric()
                    names_2rename <- character()
                    for (Patients in unique(patients_short)){
                        selector <- Patients == patients_short
                        if (sum(selector) > 1) {
                            coluns_2rename <- c(coluns_2rename, grep(Patients, patients_short)[-1])
                            names_2rename <- c(names_2rename, paste0(Patients, seq((sum(selector)-1))))
                        }
                    }
                    colnames(protein_table)[coluns_2rename] <- names_2rename
                    colnames(protein_table)[-coluns_2rename] <- unique(patients_short)
                } else {
                    colnames(protein_table) <- patients_short
                }

                assign("protein_patients", patients, envir = get(envir_link))
                assign("protein_tumor_normalized", protein_table, envir = get(envir_link))

                if (saveData){
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = protein_table,
                                file = file.path(DIR, "protein_tumor_normalized.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            } else {
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                not.tumor.regex <- paste(10:19, collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"), colnames(protein_table))
                not.tumor.tissue <- grepl(pattern = paste0("-", not.tumor.regex, "[aAbB]-"), colnames(protein_table))

                tumor_patients <- as.character(colnames(protein_table)[only.tumor.tissue])
                tumor_protein_table <- protein_table[, tumor_patients]
                assign("tumor_patients", tumor_patients, envir = get(envir_link))
                assign("protein_tumor_normalized", tumor_protein_table, envir = get(envir_link))

                if (sum(not.tumor.tissue) == 0){
                    message("There is no 'Normal' data in ", tumor,
                            " tumor folder! Please, rerun after set 'tumorData = TRUE'.",
                            "\n\n")
                    stop()
                }
                not_tumor_patients <- as.character(colnames(protein_table)[not.tumor.tissue])
                protein_not_tumor <- protein_table[, not_tumor_patients]
                assign("not_tumor_patients", not_tumor_patients, envir = get(envir_link))
                assign("protein_not_tumor_normalized", protein_not_tumor, envir = get(envir_link))

                if (saveData){
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = tumor_protein_table,
                                file = file.path(DIR, "protein_tumor_normalized.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                    write.table(x = protein_not_tumor,
                                file = file.path(DIR, "protein_not_tumor_normalized.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            }

            if (!missing(Name)){
                if (tumorData) {
                    if (normalization){
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein.selector.final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector, perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein.selector.final, ])
                    # if (length(genes.1) == 0) {
                    #     stop("protein not found!!!\n\n")
                    # }
                    # genes.transpo1 <- as.matrix(genes.1)
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_", toupper(Name)),
                           protein_exp_selected, envir = get(envir_link))
                } else {
                    # tumor
                    if (normalization){
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein.selector.final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector, perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein.selector.final, ])
                    # if (length(genes.1) == 0) {
                    #     stop("protein not found!!!\n\n")
                    # }
                    # genes.transpo1 <- as.matrix(genes.1)
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_", toupper(Name)),
                           protein_exp_selected, envir = get(envir_link))

                    # not tumor
                    protein_selector <- rownames(string_vars[["envir_link"]]$protein_not_tumor_normalized)
                    protein.selector.final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector, perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_not_tumor_normalized[protein.selector.final, ])
                    # if (length(genes.1) == 0) {
                    #     stop("protein not found!!!\n\n")
                    # }
                    # genes.transpo1 <- as.matrix(genes.1)
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_not_tumor_normalized_selected_", toupper(Name)),
                           protein_exp_selected, envir = get(envir_link))
                }
            }

        } else {
            if (!missing(Name)){
                if (tumorData) {
                    if (normalization){
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein.selector.final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector, perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein.selector.final, , drop = FALSE])
                    # if (length(genes.1) == 0) {
                    #     stop("protein not found!!!\n\n")
                    # }
                    # genes.transpo1 <- as.matrix(genes.1)
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_", toupper(Name)),
                           protein_exp_selected, envir = get(envir_link))
                } else {
                    # tumor
                    if (normalization){
                        protein_selector <- rownames(string_vars[["envir_link"]]$protein_tumor_normalized)
                    } else {
                        message("There is no not normalized protein data!")
                        stop()
                    }
                    protein.selector.final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector, perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_tumor_normalized[protein.selector.final, , drop = FALSE])
                    # if (length(genes.1) == 0) {
                    #     stop("protein not found!!!\n\n")
                    # }
                    # genes.transpo1 <- as.matrix(genes.1)
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_tumor_normalized_selected_", toupper(Name)),
                           protein_exp_selected, envir = get(envir_link))

                    # not tumor
                    protein_selector <- rownames(string_vars[["envir_link"]]$protein_not_tumor_normalized)
                    protein.selector.final <- grepl(paste0("^", Name, "$"),
                                                    protein_selector, perl = TRUE)

                    protein_exp_selected <- t(string_vars[["envir_link"]]$protein_not_tumor_normalized[protein.selector.final, , drop = FALSE])
                    # if (length(genes.1) == 0) {
                    #     stop("protein not found!!!\n\n")
                    # }
                    # genes.transpo1 <- as.matrix(genes.1)
                    Name <- gsub("-", "_", Name)
                    colnames(protein_exp_selected) <- Name
                    assign(paste0("protein_not_tumor_normalized_selected_", toupper(Name)),
                           protein_exp_selected, envir = get(envir_link))
                }
            }
        }
    }

    # mirna ####
    if ("mirna" %in% strsplit(tolower(dataType), split = " ")[[1]][1] || "isoform expression quantification" %in% tolower(dataType)){
        if (!onlyFilter) {
            if (tolower(dataBase) == "legacy"){
                DIR <- file.path(workDir, "GDCRtools", toupper(tumor), paste0(tolower(dataType), "_data"))

                manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"), select = c("file_name", "cases", "platform"))

                tmp <- !grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                              sep = "|"),
                              x = dir(path = DIR))

                mirna.files.downloaded <- dir(path = DIR)[tmp]
                if (use.hg19.mirbase20){
                    mirna.files.downloaded <- mirna.files.downloaded[grepl("hg19.mirbase20", mirna.files.downloaded, fixed=TRUE)]
                } else {
                    mirna.files.downloaded <- mirna.files.downloaded[!grepl("hg19.mirbase20", mirna.files.downloaded, fixed=TRUE)]
                }
                manifest <- manifest[manifest$file_name %in% mirna.files.downloaded, ]

                #only mirna data!!
                if (sum(tmp) > 0){
                    if (tolower(Platform) %in% "illumina hiseq"){
                        mirna.files <- manifest[manifest$platform == "Illumina HiSeq", ]
                    } else if (tolower(Platform) %in% "illumina ga"){
                        mirna.files <- manifest[manifest$platform == "Illumina GA", ]
                    } else if (tolower(Platform) %in% "h-mirna_8x15kv2") {
                        mirna.files <- manifest[manifest$platform == "H-miRNA_8x15Kv2", ]
                    } else if (tolower(Platform) %in% "h-mirna_8x15kv") {
                        mirna.files <- manifest[manifest$platform == "H-miRNA_8x15Kv", ]
                    }

                } else {
                    stop(message("No miRNA data was downloaded!"))
                }

            } else if (tolower(dataBase) == "gdc"){
                DIR <- file.path(workDir, "GDCRtools", toupper(tumor), paste0("gdc_", tolower(dataType), "_data"))

                manifest <- data.table::fread(file.path(DIR, "manifest.sdrf"), select = c("file_name", "cases"))
                tmp <- !grepl(pattern = paste(".sdrf", "Data_access_time.txt",
                                              sep = "|"),
                              x = dir(path = DIR))

                mirna.files.downloaded <- dir(path = DIR)[tmp]
                manifest <- manifest[manifest$file_name %in% mirna.files.downloaded, ]
                mirna.files <- manifest
            }

            #checking for unmatched data
            if (nrow(mirna.files) > 0){
                message("Reading data...")
            } else {
                message(paste0("Wrong Plataform '", tolower(Platform), "' chosen!"))
                stop()
            }

            if(tumorData){
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"), mirna.files$cases)
                #not related the control samples
                mirna.tumor.files <- as.character(mirna.files[[1]][only.tumor.tissue])
                patients <- as.character(as.matrix(mirna.files[only.tumor.tissue, "cases"]))
                assign("patients", patients, envir = get(envir_link))
                OPEN.miRNA(mirna.tumor.files)
                #saving (should it happends?)
                # write.table(x = data.mirna.table, file = "./mirnaylation_table/mirna_table_completed.csv",
                #           row.names = TRUE, quote = FALSE)
                assign("mirna_tumor_cross_mapped", string_vars[["envir_link"]]$cross.mapped.mirna.table, envir = get(envir_link))

                if (normalization) {
                    assign("mirna_tumor_normalized", string_vars[["envir_link"]]$data.mirna.table, envir = get(envir_link))
                    if (saveData){
                        message("\nSaving your data...\n")
                        #saving
                        write.table(x = string_vars[["envir_link"]]$data.mirna.table,
                                    file = file.path(DIR, "mirna_tumor_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = string_vars[["envir_link"]]$cross.mapped.mirna.table,
                                    file = file.path(DIR, "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                    }
                } else {
                    assign("mirna_tumor_not_normalized", string_vars[["envir_link"]]$data.mirna.table, envir = get(envir_link))
                    if (saveData){
                        message("\nSaving your data...\n")
                        #saving
                        write.table(x = string_vars[["envir_link"]]$data.mirna.table,
                                    file = file.path(DIR, "mirna_tumor_not_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = string_vars[["envir_link"]]$cross.mapped.mirna.table,
                                    file = file.path(DIR, "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                    }
                }
            } else {
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                not.tumor.regex <- paste(10:19, collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "-"), mirna.files$cases)
                not.tumor.tissue <- grepl(pattern = paste0("-", not.tumor.regex, "-"), mirna.files$cases)

                # tumor
                mirna.tumor.files <- as.character(mirna.files[[1]][only.tumor.tissue])
                patients <- as.character(as.matrix(mirna.files[only.tumor.tissue, "cases"]))
                OPEN.miRNA(mirna.tumor.files)
                tumor.data.mirna.table <- string_vars[["envir_link"]]$data.mirna.table
                tumor_cross_mapped_mirna_table <- string_vars[["envir_link"]]$cross.mapped.mirna.table
                patients.tumor <- patients

                if (normalization) {
                    assign("mirna_tumor_normalized", tumor.data.mirna.table, envir = get(envir_link))
                } else {
                    assign("mirna_tumor_not_normalized", tumor.data.mirna.table, envir = get(envir_link))
                }

                suppressWarnings(remove(data.mirna.table,
                                        cross.mapped.mirna.table,
                                        envir = get(envir_link)))

                # not tumor
                mirna.not.tumor.files <- as.character(mirna.files[[1]][not.tumor.tissue])
                if(length(mirna.not.tumor.files) == 0) {
                    message("There is no 'Normal' data in ", tumor,
                            " tumor folder! Please, rerun after set 'tumorData = TRUE'.",
                            "\n\n")
                    stop()
                }
                patients <- as.character(as.matrix(mirna.files[not.tumor.tissue, "cases"]))
                OPEN.miRNA(mirna.not.tumor.files)
                not.tumor.data.mirna.table <- string_vars[["envir_link"]]$data.mirna.table
                not.tumor.cross.mapped.mirna.table <- string_vars[["envir_link"]]$cross.mapped.mirna.table
                patients.not.tumor <- patients

                if (normalization) {
                    assign("mirna_not_tumor_normalized", not.tumor.data.mirna.table, envir = get(envir_link))
                } else {
                    assign("mirna_not_tumor_not_normalized", not.tumor.data.mirna.table, envir = get(envir_link))
                }

                if (normalization) {
                    if (saveData){
                        message("\nSaving your data...\n")
                        write.table(x = tumor.data.mirna.table,
                                    file = file.path(DIR, "mirna_tumor_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = tumor_cross_mapped_mirna_table,
                                    file = file.path(DIR, "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = not.tumor.data.mirna.table,
                                    file = file.path(DIR, "mirna_not_tumor_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = not.tumor.cross.mapped.mirna.table,
                                    file = file.path(DIR, "mirna_not_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                    }
                } else {
                    if (saveData){
                        message("\nSaving your data...\n")
                        write.table(x = tumor.data.mirna.table,
                                    file = file.path(DIR, "mirna_tumor_not_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = tumor_cross_mapped_mirna_table,
                                    file = file.path(DIR, "mirna_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = not.tumor.data.mirna.table,
                                    file = file.path(DIR, "mirna_not_tumor_not_normalized.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                        write.table(x = not.tumor.cross.mapped.mirna.table,
                                    file = file.path(DIR, "mirna_not_tumor_cross_mapped.tsv"),
                                    row.names = TRUE, quote = FALSE, sep = "\t")
                    }
                }
            }

            if (!missing(Name)){
                if (tumorData){
                    if (normalization){
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA.selector.final <- grepl(paste0("^", Name, "$"),
                                                      miRNA_selector, perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])

                        if (sum(miRNA.selector.final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- Name
                        colnames(cross_mapped_selected) <- Name
                        Name <- gsub("-", "_", Name)
                        assign(paste0("mirna_tumor_normalized_selected_", toupper(Name)), miRNA_exp_selected, envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_", toupper(Name)), cross_mapped_selected, envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to separate patients!")
                        stop()
                    }
                } else {
                    # tumor
                    if (normalization){
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA.selector.final <- grepl(paste0("^", Name, "$"),
                                                      miRNA_selector, perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])

                        if (sum(miRNA.selector.final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- Name
                        colnames(cross_mapped_selected) <- Name
                        Name <- gsub("-", "_", Name)
                        assign(paste0("mirna_tumor_normalized_selected_", toupper(Name)), miRNA_exp_selected, envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_", toupper(Name)), cross_mapped_selected, envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to separate patients!")
                        stop()
                    }

                    # not tumor
                    if (normalization){
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_not_tumor_normalized)
                        miRNA.selector.final <- grepl(paste0("^", Name, "$"),
                                                      miRNA_selector, perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA.selector.final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA.selector.final, , drop = FALSE])

                        colnames(miRNA_exp_selected) <- Name
                        colnames(cross_mapped_selected) <- Name
                        Name <- gsub("-", "_", Name)
                        assign(paste0("mirna_not_tumor_normalized_selected_", toupper(Name)), miRNA_exp_selected, envir = get(envir_link))
                        assign(paste0("mirna_not_tumor_cross_mapped_selected_", toupper(Name)), cross_mapped_selected, envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to separate patients!")
                        stop()
                    }
                }
            }

        } else {
            if (!missing(Name)){
                if (tumorData){
                    if (normalization){
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA.selector.final <- grepl(paste0("^", Name, "$"),
                                                      miRNA_selector, perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])

                        if (sum(miRNA.selector.final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- Name
                        colnames(cross_mapped_selected) <- Name
                        Name <- gsub("-", "_", Name)
                        assign(paste0("mirna_tumor_normalized_selected_", toupper(Name)), miRNA_exp_selected, envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_", toupper(Name)), cross_mapped_selected, envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to separate patients!\n")
                        stop()
                    }
                } else {
                    # tumor
                    if (normalization){
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_tumor_normalized)
                        miRNA.selector.final <- grepl(paste0("^", Name, "$"),
                                                      miRNA_selector, perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_tumor_normalized[miRNA.selector.final, , drop = FALSE])

                        if (sum(miRNA.selector.final) == 0) {
                            stop(message("\nmiRNA not found!!!\n\n"))
                        }

                        colnames(miRNA_exp_selected) <- Name
                        colnames(cross_mapped_selected) <- Name
                        Name <- gsub("-", "_", Name)
                        assign(paste0("mirna_tumor_normalized_selected_", toupper(Name)), miRNA_exp_selected, envir = get(envir_link))
                        assign(paste0("mirna_tumor_cross_mapped_selected_", toupper(Name)), cross_mapped_selected, envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to separate patients!\n")
                        stop()
                    }

                    # not tumor
                    if (normalization){
                        miRNA_selector <- rownames(string_vars[["envir_link"]]$mirna_not_tumor_normalized)
                        miRNA.selector.final <- grepl(paste0("^", Name, "$"),
                                                      miRNA_selector, perl = TRUE)

                        miRNA_exp_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA.selector.final, , drop = FALSE])
                        cross_mapped_selected <- t(string_vars[["envir_link"]]$mirna_not_tumor_normalized[miRNA.selector.final, , drop = FALSE])

                        colnames(miRNA_exp_selected) <- Name
                        colnames(cross_mapped_selected) <- Name
                        Name <- gsub("-", "_", Name)
                        assign(paste0("mirna_not_tumor_normalized_selected_", toupper(Name)), miRNA_exp_selected, envir = get(envir_link))
                        assign(paste0("mirna_not_tumor_cross_mapped_selected_", toupper(Name)), cross_mapped_selected, envir = get(envir_link))
                    } else {
                        message("Only normalized data can be used to separate patients!\n")
                        stop()
                    }
                }
            }
        }
        suppressWarnings(remove(data.mirna.table,
                                cross.mapped.mirna.table,
                                envir = get(envir_link)))
    }
        # stop(message("Please use one of these: ",
        #              "\nLegacy - 'miRNA gene quantification' and 'miRNA isoform quantification'",
        #              "\nGDC - 'miRNA expression quantification' and 'isoform expression quantification'"))
        #

    # Exon quantification ####
    if ("exon" == strsplit(tolower(dataType), split = " ")[[1]][1]){
        if (!onlyFilter) {
            if (tolower(dataBase) == "gdc") {
                stop(message("\nThrere is no protein expression data in GDC data base!!",
                             "\nPlease use 'Legacy' data base\n"))
            }
            DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "exon quantification_data")


            if (!dir.exists(file.path(workDir, "GDCRtools", toupper(tumor), "mage_data"))) {
                message("Downloading magetab data...")
                suppressWarnings(download_gdc(dataType = "mage",
                                              dataBase = dataBase,
                                              tumor = tumor,
                                              workDir = workDir))
            }

            desing_array_DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "mage_data")

            if (tolower(Platform) == ""){
                message("Please, insert 'illumina hiseq' or 'illumina ga' for Exon quantification data.")
                stop()
            } else if (tolower(Platform) == "illumina ga"){
                # to.gunzip <- file.path(desing_array_DIR,
                #                        paste0("unc.edu_", toupper(tumor), ".IlluminaHiSeq_RNASeqV2.mage-tab.1.1.0.tar.gz"))
                # ## or, if you just want to extract the target file:
                # untar(to.gunzip, exdir = desing_array_DIR)
                #
                # design_array <- data.table::fread(file.path(desing_array_DIR,
                #                                             paste0("unc.edu_", toupper(tumor), ".IlluminaHiSeq_RNASeqV2.mage-tab.1.1.0"),
                #                                             paste0("unc.edu_", toupper(tumor), ".IlluminaHiSeq_RNASeqV2.1.1.0.sdrf.txt")),
                #                                   select = c(2, 22))
            } else if (tolower(Platform) == "illumina hiseq"){
                to.gunzip <- file.path(desing_array_DIR,
                                       paste0("unc.edu_", toupper(tumor), ".IlluminaHiSeq_RNASeqV2.mage-tab.1.1.0.tar.gz"))
                ## or, if you just want to extract the target file:
                untar(to.gunzip, exdir = desing_array_DIR)

                design_array <- data.table::fread(file.path(desing_array_DIR,
                                                            paste0("unc.edu_", toupper(tumor), ".IlluminaHiSeq_RNASeqV2.mage-tab.1.1.0"),
                                                            paste0("unc.edu_", toupper(tumor), ".IlluminaHiSeq_RNASeqV2.1.1.0.sdrf.txt")),
                                                  select = c(2, 22))
            }

            exon_files <- dir(DIR, pattern = "unc.edu")

            patient_code <- as.character(unlist(design_array[design_array[[2]] %in% exon_files, 1]))

            if (normalization){
                select_this_column <- 4 #RPKM
            } else {
                select_this_column <- 2 #raw counts
            }

            pb <- txtProgressBar(min = 0, max = length(exon_files), style = 3)
            count <- 0
            for (i in file.path(DIR, exon_files)){
                count <- count + 1
                if (count == 1){
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

            if(tumorData){
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"), patient_code)
                #not related the control samples
                patients_tumor <- patient_code[only.tumor.tissue]
                assign("patients", patients_tumor, envir = get(envir_link))
                assign("exon_tumor_", exon_table, envir = get(envir_link))

                if (saveData){
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = exon_table,
                                file = file.path(DIR, "tumor_exon_table.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }
            } else {
                #separing files
                tumor.regex <- paste(formatC(seq(1:9), width=2, flag="0"), collapse = "[aAbB]-|-")
                not.tumor.regex <- paste(10:19, collapse = "[aAbB]-|-")
                only.tumor.tissue <- grepl(pattern = paste0("-", tumor.regex, "[aAbB]-"), patient_code)
                not.tumor.tissue <- grepl(pattern = paste0("-", not.tumor.regex, "[aAbB]-"), patient_code)

                tumor_patients <- as.character(colnames(exon_table)[only.tumor.tissue])
                tumor_protein_table <- exon_table[, tumor_patients]
                assign("tumor_patients", tumor_patients, envir = get(envir_link))
                assign("exon_tumor_", tumor_protein_table, envir = get(envir_link))

                if (sum(not.tumor.tissue) == 0){
                    message("There is no 'Normal' data in ", tumor,
                            " tumor folder! Please, rerun after set 'tumorData = TRUE'.",
                            "\n\n")
                    stop()
                }
                not_tumor_patients <- as.character(colnames(exon_table)[not.tumor.tissue])
                exon_not_tumor <- exon_table[, not_tumor_patients]
                assign("not_tumor_patients", not_tumor_patients, envir = get(envir_link))
                assign("exon_not_tumor_", exon_not_tumor, envir = get(envir_link))

                if (saveData){
                    message("\nSaving your data...\n")
                    #saving
                    write.table(x = tumor_protein_table,
                                file = file.path(DIR, "tumor_protein_table.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                    write.table(x = exon_not_tumor,
                                file = file.path(DIR, "exon_not_tumor.tsv"),
                                row.names = TRUE, quote = FALSE, sep = "\t")
                }

            }

        } else {
            if (normalization){
                protein_selector <- rownames(string_vars[["envir_link"]]$tumor_protein_table)
            } else {
                # completed.table1 <- local(COMPLETED.TABLE.NOT.NORMALIZED, envir = get(envir_link))
                stop(message("Only normalized data can be used to separate patients!"))
            }
            protein.selector.final <- grepl(Name,
                                            protein_selector, perl=TRUE)

            protein_exp_selected <- t(string_vars[["envir_link"]]$tumor_protein_table[protein.selector.final, ])
            # if (length(genes.1) == 0) {
            #     stop("protein not found!!!\n\n")
            # }
            # genes.transpo1 <- as.matrix(genes.1)
            colnames(protein_exp_selected) <- Name
            colnames(cross_mapped_selected) <- Name
            assign(paste0("exon_exp_selected_", Name), protein_exp_selected, envir = get(envir_link))
        }


    }
    # mage ####
    # if ("mage" == tolower(dataType)){
    #     if (tolower(dataBase) == "legacy"){
    #         DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "mage_data")
    #     } #else if (tolower(dataBase) == "gdc"){
    #     #     DIR <- file.path(workDir, "GDCRtools", toupper(tumor), "gdc_protein_data")
    #     # }
    #     patient_file <- dir(path = DIR, pattern = "mage-tab")
    #     count <- 0
    #
    # }

    # end ####
    assign("tumor", tumor, envir = get(envir_link))
    assign("dataBase", dataBase, envir = get(envir_link))
    assign("dataType", dataType, envir = get(envir_link))
    assign("workDir", workDir, envir = get(envir_link))
    # assign("DIR", DIR, envir = get(envir_link))
    if (!missing("Name")) {
        assign("Name", Name, envir = get(envir_link))
        assign("Name.e", Name, envir = get(envir_link))
    }

    message("\nDone!")
}
