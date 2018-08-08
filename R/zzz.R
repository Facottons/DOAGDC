# .onLoad <- function(libname, pkgname)
# {
#     library.dynam("GDCRtools", pkgname, libname)
# }

# http://stackoverflow.com/questions/15261619/sample-gives-different-values-with-same-set-seed
# message(sprintf("On %s I realized %s was...\n%s by the street", Sys.Date(), person, action))
# .onAttach <- function (libname, pkgname){
.onAttach <- function(lib, pkg) {
    # unlockBinding(".GDCRtools", asNamespace("GDCRtools"))
    # version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    #are you human?
    GDCTerms <- "By using this package you accept the terms in 'https://portal.gdc.cancer.gov'"
    if(interactive())
    { # figlet (-f doom) GDCRtools
        packageStartupMessage(
"  ____ ____   ____ ____  _              _
 / ___|  _ \\ / ___|  _ \\| |_ ___   ___ | |___
| |  _| | | | |   | |_) | __/ _ \\ / _ \\| / __|
| |_| | |_| | |___|  _ <| || (_) | (_) | \\__ \\
 \\____|____/ \\____|_| \\_\\\\__\\___/ \\___/|_|___/  \n\nVersion: ", version, "\n",
        "\n", GDCTerms, "\n")

    }
    else
    { packageStartupMessage("Package 'GDCRtools' \n\nversion ", version, "\n", "\n", GDCTerms) }

    packageStartupMessage("For citing this R package in publications, ",
                          "type 'citation(\"GDCRtools\")'.")
    invisible()

    # suppressWarnings(remove(gene.idtype.bods, gene.idtype.list, cpd.simtypes))
}

# packageStartupMessage("
#  ___________  _____ ______ _              _
# |  __ \\  _  \\/  __ \\| ___ \\ |            | |
# | |  \\/ | | || /  \\/| |_/ / |_ ___   ___ | |___
# | | __| | | || |    |    /| __/ _ \\ / _ \\| / __|
# | |_\\ \\ |/ / | \\__/\\| |\\ \\| || (_) | (_) | \\__ \\
#  \\____/___/   \\____/\\_| \\_|\\__\\___/ \\___/|_|___/
# ")
