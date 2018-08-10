#' From global to GDCRtools environment
#'
#' @param table A matrix or data frame with data used in this package (i.e Gene,
#'   isoform and protein expression, methylation and mutation data)
#' @param dataType
#' @param normalization
#' @param dataBase If the data was downloaded from GDC/Legacy data base, however
#'   not using GDCRtools, please specified which data base. If the data do not
#'   come from GDC/Legacy data base, and it is related with genome version GRCh
#'   37, please insert "legacy" in this argument. Otherwise, if it is related
#'   with genome version GRCh 38, please insert "GDC" in this argument.
#' @param tumor
#' @param tumorData
#' @param env
#' @inheritParams concatenate_files
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' patient <- paste0(paste0("patient_", LETTERS[1:4]))
#' genes <- paste0("gene_", seq(1, 5))
#' # generate a simulated gene expression matrix
#' example_gene_table <- matrix(runif(20, 0.0, 90.5), 5, 4, TRUE, list(genes, patient))
#' # without env created
#' table2GDCRtools(example_gene_table, dataType = "gene", dataBase = "legacy", tumor = "UCS")
#' # with env created
#' table2GDCRtools(example_gene_table, env = "env name without quotes")
#' }
table2GDCRtools <- function(table, dataType,
                            normalization = TRUE,
                            dataBase, tumor,
                            tumorData = TRUE,
                            env){
    if(missing(env)) {
        if (!(paste(tumor, dataBase, dataType, sep = "_") %in% ls(all.names = TRUE, envir = .GlobalEnv)) && tumorData){
            assign(paste(tumor, dataBase, dataType, "tumor_data", sep = "_"),
                   new.env(parent=emptyenv()), envir = .GlobalEnv)
            envir_link <- paste(tumor, dataBase, dataType, "tumor_data", sep = "_")
            message(paste0("It was created", paste(tumor, dataBase, dataType, "tumor_data", sep = "_"), "Environment"))
        } else if (!(paste(tumor, dataBase, dataType, sep = "_") %in% ls(all.names = TRUE, envir = .GlobalEnv))){
            assign(paste(tumor, dataBase, dataType, "both_data", sep = "_"),
                   new.env(parent=emptyenv()), envir = .GlobalEnv)
            envir_link <- paste(tumor, dataBase, dataType, "both_data", sep = "_")
            message(paste0("It was created", paste(tumor, dataBase, dataType, "both_data", sep = "_"), "Environment"))
        }
    }

    envir_link <- deparse(substitute(env))

    new_object_name <- deparse(substitute(table))

    assign(new_object_name, table, envir = get(envir_link))
}
