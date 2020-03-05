#' From global to DOAGDC environment
#'
#' @param table A matrix or data frame with data used in this package (i.e Gene,
#'   isoform and protein expression, methylation and mutation data)
#' @param data_type
#' @param normalization
#' @param data_base If the data was downloaded from GDC/Legacy data base,
#'   however
#'   not using DOAGDC, please specified which data base. If the data do not come
#'   from GDC/Legacy data base, and it is related with genome version GRCh 37,
#'   please insert "legacy" in this argument. Otherwise, if it is related with
#'   genome version GRCh 38, please insert "GDC" in this argument.
#' @param tumor
#' @param tumor_data
#' @param env
#' @inheritParams concatenate_exon
#' @inheritParams download_gdc
#'
#' @return the objects imported are stored inside the determined environment
#'   name.
#' @export
#'
#' @examples
#' patient <- paste0(paste0("patient_", LETTERS[1:4]))
#' genes <- paste0("gene_", seq(1, 5))
#' # generate a simulated gene expression matrix
#' example_gene_table <- matrix(
#'    runif(20, 0.0, 90.5), 5, 4, TRUE,
#'    list(genes, patient)
#' )
#' # without env created
#' table2_doagdc(example_gene_table,
#'    data_type = "gene",
#'    data_base = "legacy",
#'    tumor = "UCS"
#' )
table2_doagdc <- function(table, data_type,
                        normalization = TRUE,
                        data_base, tumor,
                        tumor_data = TRUE,
                        env) {
    if (missing(env)) {
        tmp <- (paste(
            tumor, data_base, data_type,
            sep = "_"
        ) %in% ls(
            all.names = TRUE, envir = .GlobalEnv
        ))

        if (!tmp && tumor_data) {
            assign(paste(tumor, data_base, data_type, "tumor_data", sep = "_"),
                new.env(parent = emptyenv()),
                envir = .GlobalEnv
            )
            envir_link <- paste(tumor, data_base, data_type, "tumor_data",
                sep = "_"
            )
            message(paste0("It was created", paste(tumor, data_base, data_type,
                "tumor_data",
                sep = "_"
            ), "Environment"))
        } else if (!(paste(tumor, data_base, data_type,
            sep = "_"
        ) %in% ls(all.names = TRUE, envir = .GlobalEnv))) {
            assign(paste(tumor, data_base, data_type, "both_data", sep = "_"),
                new.env(parent = emptyenv()),
                envir = .GlobalEnv
            )
            envir_link <- paste(tumor, data_base, data_type, "both_data",
                sep = "_"
            )
            message(paste0("It was created", paste(tumor, data_base, data_type,
                "both_data",
                sep = "_"
            ), "Environment"))
        }
    }

    envir_link <- deparse(substitute(env))

    new_object_name <- deparse(substitute(table))

    assign(new_object_name, table, envir = get(envir_link))
}
