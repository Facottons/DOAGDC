#' Convert beta values to mvalues
#'
#' @param probe_name A character string containing the probe name desired from
#'    the selected \code{'name'} probes.
#' @param env
#' @param save_data Logical value where \code{TRUE} indicates that the
#'    concatenate and filtered matrix should be saved in local storage. The
#'    default is FALSE.
#'
#' @inheritParams concatenate_exon
#'
#' @return the beta values converted in mValues stored inside the determined
#'    environment name.
#' @export
#'
#' @examples
#' library(DOAGDC)
#'
#' download_gdc(
#'    data_type = "methylation",
#'    data_base = "legacy", work_dir = "~/Desktop", tumor = "CHOL",
#'    platform = "Illumina Human Methylation 450"
#' )
#'
#' concatenate_methylation(
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     platform = "Illumina Human Methylation 450",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
#'
#' beta2m_values(
#'    probe_name = "cg16672562",
#'    save_data = TRUE,
#'    env = CHOL_LEGACY_methylation_tumor_data
#' )
beta2m_values <- function(probe_name,
                            env,
                            save_data = FALSE) {

    ## from beta to mValues (already filtered data!) Comparison of Beta-value
    # and M-value methods for quantifying methylation levels by microarray
    # analysis Pan Du1,3*, Xiao Zhang2, Chiang-Ching Huang2,

    if (missing(env)) {
        message(
            "Please, before using this function, ",
            "insert the Environment name."
        )
    } else {
        envir_link <- deparse(substitute(env))
    }
    string_vars <- list(envir_link = get(envir_link))

    message("Converting the data \nPlease wait...")

    name <- string_vars[["envir_link"]]$name
    patients <- string_vars[["envir_link"]]$patients
    work_dir <- string_vars[["envir_link"]]$work_dir
    tumor <- string_vars[["envir_link"]]$tumor

    meth_table <- eval(parse(text = paste0(
        "string_vars[['envir_link']]",
        "$methylation",
        "_tumor_filtered_selected_",
        toupper(name)
    )))

    meth_matrix <- meth_table[probe_name, , drop = FALSE]
    assign(paste0("methylation_tumor_filtered_selected_", toupper(probe_name)),
        t(meth_matrix),
        envir = get(envir_link)
    )
    m_values_vec <- sapply(meth_matrix, function(x) {
        if (is.na(x)) {
            x <- NA
        } else {
            x <- log2(x / (1 - x))
        }
    })
    m_values_matrix <- as.matrix(m_values_vec)
    dimnames(m_values_matrix) <- list(patients, tolower(probe_name))

    # saving
    if (save_data) {
        message("\nSaving your data...\n")
        # saving
        write.csv(
            x = m_values_matrix,
            file = file.path(
                work_dir, "DOAGDC",
                toupper(tumor), "Analyses", name, "Methylation",
                paste0(tolower(probe_name), "_mValues.csv")
            ),
            row.names = TRUE, quote = FALSE
        )
    }

    assign(paste0(
        "methylation_tumor_filtered_selected_", toupper(probe_name),
        "_mValues"
    ),
    m_values_matrix,
    envir = get(envir_link)
    )
    assign("probe_name", probe_name, envir = get(envir_link))
    assign("name_e", probe_name, envir = get(envir_link))

    message("\nDone!")
}
