# .onLoad <- function(libname, pkgname) {
#    # something
# }

# http://stackoverflow.com/questions/15261619/
# sample-gives-different-values-with-same-set-seed
# message(sprintf("On %s I realized %s was...\n%s by the street",
# Sys.Date(), person, action))
.onAttach <- function(lib, pkg) {
    gdc_terms <- paste0("By using this package you accept the terms",
                                        " in 'https://portal.gdc.cancer.gov'")
    if (interactive()) {
        # figlet (-f doom) DOAGDC
        packageStartupMessage("
    DDD   OOO   AA   GGG  DDD   CCC
    D  D O   O A  A G     D  D C
    D  D O   O AAAA G  GG D  D C
    D  D O   O A  A G   G D  D C
    DDD   OOO  A  A  GGG  DDD   CCC    \n\nVersion: ", version, "\n",
        "\n", gdc_terms, "\n")

    } else {
        packageStartupMessage("Package 'DOAGDC' \n\nversion ", version,
                                                        "\n", "\n", gdc_terms)
    }

    packageStartupMessage("For citing this R package in publications, ",
                        "type 'citation(\"DOAGDC\")'.")
    invisible()
}
