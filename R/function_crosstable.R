#' Venn diagram of differential expression genes list
#'
#' @param FinalData A character string indicating the name of which diferential
#'   expression package should be used to get the statistical values in the
#'   final list. The default is "EBSeq".
#' @param n.pack A numerical value indicating the number of expression analysis
#'   to be used in venn diagram. It is expected the number 2 or 3. The default
#'   is 3.
#' @param packageNames A character vector indicating the names of at least two
#'   diferential expression packages used in previous steps: "DESeq2", "edgeR",
#'   "DESeq2, or "All".
#' @param pairName
#' @param Width,Height,Res,Unit
#' @param Colors A character vector indicating the colors to be used in the venn
#'   diagram. The default is c('green', 'blue', "red").
#' @param VennDiagram.imagetype A character string indicating the image format
#'   (e.g. "tiff", "png" or "svg"). The default is "png".
#' @param workDir
#' @param env
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#'
#' @return
#' @export
#'
#' @examples
#' CrossThem("DESeq2", 2, packageNames = c("edgeR", "DESeq2"), env = "env name without quotes")
CrossThem <- function(FinalData, n.pack = 3,
                      packageNames,
                      pairName = "G2_over_G1",
                      Width = 2000,
                      Height = 2000,
                      Res = 300,
                      Unit = "px",
                      Colors = c('green', 'blue', "red"),
                      VennDiagram.imagetype = "png",
                      workDir,
                      env) {

    #verifying if the package is already installed
    # to.load <- c("stringi", "stringr", "VennDiagram")

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    dir.create(path = paste0(PATH, "/CrossData_", tolower(FinalData)), showWarnings = FALSE)
    DIR <- paste0(PATH, "/CrossData_", tolower(FinalData))

    if (is.null(FinalData)){
        stop("Please insert which DE method is going to be used in result!! (EBSEq, DESeq2 or edgeR?)")
    }

    if (tolower(FinalData) == "ebseq") {
        DEcross <- string_vars[["envir_link"]]$resultadosDE.EBSeq[[pairName]]
        complete.cross <- string_vars[["envir_link"]]$Results.Completed.EBSeq[[pairName]]
    } else if (tolower(FinalData) == "edger") {
        DEcross <- string_vars[["envir_link"]]$resultadosDE.edgeR[[pairName]]
        complete.cross <- string_vars[["envir_link"]]$Results.Completed.edgeR[[pairName]]
    } else if (tolower(FinalData) == "deseq2") {
        DEcross <- string_vars[["envir_link"]]$resultadosDE.DESeq2[[pairName]]
        complete.cross <- string_vars[["envir_link"]]$Results.Completed.DESeq2[[pairName]]
    }

    if (tolower(packageNames) == "all") {
        packageNames <- c("edgeR", "DESeq2", "EBSeq")
    }

    if (n.pack == 3){
        p1 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE.",
                     packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE.",
                     packageNames[tolower(packageNames) != tolower(FinalData)][1], "[[pairName]]")
        p3 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE.",
                     packageNames[tolower(packageNames) != tolower(FinalData)][2], "[[pairName]]")
        p1completed <- paste0("string_vars[['envir_link']]",
                              "$Results.Completed.",
                              packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2completed <- paste0("string_vars[['envir_link']]",
                              "$Results.Completed.",
                              packageNames[tolower(packageNames) != tolower(FinalData)][1], "[[pairName]]")
        p3completed <- paste0("string_vars[['envir_link']]",
                              "$Results.Completed.",
                              packageNames[tolower(packageNames) != tolower(FinalData)][2], "[[pairName]]")
        DEcross <- DEcross[which(rownames(eval(parse(text=p1)))%in%rownames(eval(parse(text=p2)))), ]
        DEcross <- DEcross[which(rownames(DEcross)%in%rownames(eval(parse(text=p3)))), ]
        complete.cross <- complete.cross[which(rownames(eval(parse(text=p1completed)))%in%rownames(eval(parse(text=p2completed)))), ]
        complete.cross <- complete.cross[which(rownames(complete.cross)%in%rownames(eval(parse(text=p3completed)))), ]
        write.csv(x = DEcross, file = file.path(DIR,
                                               "DE_crossData.csv"))
        write.csv(x = complete.cross, file = file.path(DIR,
                                                      "complete_crossData.csv"))

        tmp <- paste0(tolower(packageNames), c(rep("_Vs_", 2), ""), collapse = "")

        #ploting venn diagram
        VennDiagram::venn.diagram(
            x = list(rownames(eval(parse(text=p1))), rownames(eval(parse(text=p2))),
                     rownames(eval(parse(text=p3)))),
            category.names = c(packageNames[tolower(packageNames) == tolower(FinalData)],
                               packageNames[tolower(packageNames) != tolower(FinalData)][1],
                               packageNames[tolower(packageNames) != tolower(FinalData)][2]),
            filename = file.path(DIR, paste0("venn_diagramm_",
                                             pairName, "_",
                                             tmp,".", tolower(VennDiagram.imagetype))),
            resolution = Res,
            output = TRUE, imagetype = VennDiagram.imagetype, height = Height, width = Width,
            units = Unit,
            compression = "lzw",
            lwd = 2, lty = 'blank',
            fill = Colors,
            cex = 1.2, print.mode = c("raw", "percent"), sigdigs = 2,
            fontface = "bold", fontfamily = "sans",
            cat.cex = 1.6, cat.fontface = "bold", cat.default.pos = "outer",
            cat.pos = c(-40, 40, 135), cat.dist = c(0.07, 0.07, 0.07),
            cat.fontfamily = "sans")
    }
    else if (n.pack == 2) {
        p1 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE.",
                     packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2 <- paste0("string_vars[['envir_link']]",
                     "$resultadosDE.",
                     packageNames[tolower(packageNames) != tolower(FinalData)][1], "[[pairName]]")
        p1completed <- paste0("string_vars[['envir_link']]",
                              "$Results.Completed.",
                              packageNames[tolower(packageNames) == tolower(FinalData)], "[[pairName]]")
        p2completed <- paste0("string_vars[['envir_link']]",
                              "$Results.Completed.",
                              packageNames[tolower(packageNames) != tolower(FinalData)], "[[pairName]]")
        DEcross <- DEcross[which(rownames(eval(parse(text=p1)))%in%rownames(eval(parse(text=p2)))), ]
        complete.cross <- complete.cross[which(rownames(eval(parse(text=p1completed)))%in%rownames(eval(parse(text=p2completed)))), ]
        write.csv(x = DEcross, file = file.path(DIR,
                                               "DE_crossData.csv"))
        write.csv(x = complete.cross, file = file.path(DIR,
                                                      "complete_crossData.csv"))

        tmp <- paste0(tolower(packageNames), c("_Vs_", ""), collapse = "")

        #ploting venn diagram
        VennDiagram::venn.diagram(
            x = list(rownames(eval(parse(text=p1))), rownames(eval(parse(text=p2)))),
            category.names = c(packageNames[tolower(packageNames) == tolower(FinalData)],
                               packageNames[tolower(packageNames) != tolower(FinalData)]),
            filename = file.path(DIR, paste0("venn_diagramm_",
                                             pairName, "_",
                                             tmp,".", tolower(VennDiagram.imagetype))),
            resolution = Res,
            output = TRUE, imagetype= VennDiagram.imagetype, height = Height, width = Width,
            units = Unit,
            compression = "lzw",
            lwd = 2, lty = 'blank',
            fill = Colors,
            cex = 1.2, print.mode = c("raw", "percent"), sigdigs = 2,
            fontface = "bold", fontfamily = "sans",
            cat.cex = 1.6, cat.fontface = "bold", cat.default.pos = "outer",
            cat.pos = c(-35, 30), cat.dist = c(0.055, -0.05),
            cat.fontfamily = "sans")
    }

    tested1 <- vector("list", 1)
    names(tested1) <- pairName
    tested1[[pairName]] <- DEcross

    tested2 <- vector("list", 1)
    names(tested2) <- pairName
    tested2[[pairName]] <- complete.cross

    assign("FinalData", FinalData, envir = get(envir_link))
    assign("Results.Completed.crossed", tested2, envir = get(envir_link))
    assign("resultadosDE.crossed", tested1, envir = get(envir_link))
    if (tolower(FinalData) == "ebseq") {
        string_vars[["envir_link"]]$Tool <- "CrossTable.EBSeq"
    } else if (tolower(FinalData) == "edger") {
        string_vars[["envir_link"]]$Tool <- "CrossTable.edgeR"
    } else if (tolower(FinalData) == "deseq2") {
        string_vars[["envir_link"]]$Tool <- "CrossTable.DESeq2"
    }

    message("Done!\n")
}
