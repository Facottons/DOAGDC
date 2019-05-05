#' Convert beta values to mvalues
#'
#' @param probeName A character string containing the probe name desired from
#'   the selected \code{'Name'} probes.
#' @param env
#' @param saveData Logical value where \code{TRUE} indicates that the
#'   concatenate and filtered matrix should be saved in local storage. The
#'   default is FALSE.
#'
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#'
#' @return the beta values converted in mValues stored inside the determined
#'   environment name.
#' @export
#'
#' @examples
#' \dontrun{
#' beta2mValues(probeName = "cg16672562", saveData = TRUE, env = "env name without quotes")
#' }
beta2mValues <- function(probeName,
                         env,
                         saveData = FALSE){

    ##from beta to mValues (already filtered data!)
    # Comparison of Beta-value and M-value methods for quantifying methylation
    # levels by microarray analysis Pan Du1,3*, Xiao Zhang2, Chiang-Ching Huang2,

    if (missing(env)) {
        message("Please, before using this function, insert the Environment name.")
    } else {
        envir_link <- deparse(substitute(env))
    }
    string_vars <- list(envir_link = get(envir_link))

    message("Converting the data \nPlease wait...")

    Name <- string_vars[["envir_link"]]$Name
    patients <- string_vars[["envir_link"]]$patients
    workDir <- string_vars[["envir_link"]]$workDir

    meth.table <- eval(parse(text= paste0("string_vars[['envir_link']]$methylation",
                                          "_tumor_filtered_selected_",
                                          toupper(Name))))

    meth.matrix <- meth.table[probeName, , drop = FALSE]
    assign(paste0("methylation_tumor_filtered_selected_", toupper(probeName)),
           t(meth.matrix), envir = get(envir_link))
    mValues.vec <- sapply(meth.matrix, function(x) {
        if(is.na(x)) {
            x <- NA
        } else {
            x <- log2(x/(1 - x))
        }})
    mValues.matrix <- as.matrix(mValues.vec)
    dimnames(mValues.matrix) <- list(patients, tolower(probeName))

    # saving
    if (saveData){
        message("\nSaving your data...\n")
        #saving
        write.csv(x = mValues.matrix,
                  file = file.path(workDir, paste0(tolower(probeName), "_mValues.tsv")),
                  row.names = TRUE, quote = FALSE, sep = "\t")
    }

    assign(paste0("methylation_tumor_filtered_selected_", toupper(probeName), "_mValues"),
           mValues.matrix, envir = get(envir_link))
    assign("probeName", probeName, envir = get(envir_link))
    assign("Name.e", probeName, envir = get(envir_link))

    message("\nDone!")
}
