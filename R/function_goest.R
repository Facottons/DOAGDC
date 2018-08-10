#' Perform Gene-Ontology Pathways Enrichment
#'
#' @param condition A character string containing which condition should be
#'   used: "Upregulated", "Downregulated" or "all".
#' @param use_genes_without_cat Logical value where \code{FALSE} indicate that
#'   genes outside the category being tested will be ignored in the calculation
#'   of p-values. The default is "FALSE".
#' @param FDR.cutoff
#' @param Width,Height,Res,Unit,image.format
#' @param Tool A character string indicating which differential expression
#'   analysis tool was last used.
#' @param ID A character string indicating which ID should be used: "HUGO",
#'   "GeneSymbol", "ensembl" , "refGene" or "GeneID". The default is \code{"GeneID"}.
#' @param pairName
#' @param env
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#'
#' @return Enriched terms.
#' @export
#'
#' @examples
#' \dontrun{
#' GOnto(condition = "Upregulated", Tool = "edgeR", env = "env name without quotes")
#' }
GOnto <- function(condition,
                  use_genes_without_cat = TRUE,
                  FDR.cutoff = 0.05,
                  Width = 2000,
                  Height = 2000,
                  Res = 300,
                  Unit = "px",
                  image.format = "png",
                  Tool,
                  ID = "GeneID",
                  pairName = "G2_over_G1",
                  env){

    # local functions ####
    prepar_path_enrich <- function(ID = "GeneID",
                                   pairName = "G2_over_G1",
                                   env,
                                   Tool) {

        #verifying if the package is already installed
        # to.load <- c("RColorBrewer", "GO.db", "goseq",
        #              "methods", "annotate", "org.Hs.eg.db", "DO.db", "reactome.db",
        #              "goseq", "clusterProfiler", "DOSE", "ReactomePA", "igraph",
        #              "pathview", "gage", "gageData")

        # @param File A character string indicating the object name containg the DE
        #   matrix. It must be named \code{Results.Completed.} and it must be inside
        #   the given \code{env}. Use the function \code{table2GDCRtools} to allocate
        #   your table inside the given \code{env}.

        if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

        envir_link <- deparse(substitute(env))
        string_vars <- list(envir_link = get(envir_link))

        if (exists("Name.e", envir = get(envir_link))){
            PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
        } else {
            PATH <- string_vars[["envir_link"]]$PATH
        }

        #Most recent Tool used
        if (missing(Tool)){Tool <- string_vars[["envir_link"]]$Tool}

        Name <- string_vars[["envir_link"]]$Name
        dataBase <- string_vars[["envir_link"]]$dataBase

        if (grepl("crosstable", tolower(Tool))) {
            if (tolower(Tool) == "crosstable.deseq2") {
                DIR <- paste0(PATH, "/CrossData_deseq2")
            } else if (tolower(Tool) == "crosstable.edger") {
                DIR <- paste0(PATH, "/CrossData_edger")
            } else if (tolower(Tool) == "crosstable.ebseq") {
                DIR <- paste0(PATH, "/CrossData_ebseq")
            }
            dir.create(file.path(DIR, "Ontology_Results"), showWarnings = FALSE)
            DIR <- file.path(DIR, "Ontology_Results")

            File <- "resultadosDE.crossed"
            File2 <- "Results.Completed.crossed"

        } else {
            File <- paste("resultadosDE", Tool, sep = ".")
            File2 <- paste("Results.Completed", Tool, sep = ".")

            dir.create(paste0(PATH, "/Ontology_Results_", Tool, "_", toupper(Name)),
                       showWarnings = FALSE)
            DIR <- paste0(PATH, "/Ontology_Results_", Tool, "_", toupper(Name))
        }

        # generate annotation table
        annotation_table <- GDCRtools::annotation_table
        # GO_annotations <- GDCRtools::GO_annotations

        #input DE genes
        resultadosDE <- get(File, envir = string_vars[["envir_link"]])[[pairName]]
        Results.Completed <- get(File, envir = string_vars[["envir_link"]])[[pairName]]

        #pre-GO        gbutils::isNA
        if (tolower(ID) == "geneid") {

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$GeneID <- annotation_table[resultadosDE$ensembl, "GeneID"]
                Results.Completed$GeneID <- annotation_table[Results.Completed$ensembl, "GeneID"]
            }

            #remove any duplicates and NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneID), ]
            resultadosDE <- resultadosDE[!duplicated(resultadosDE$GeneID), ]

            Results.Completed <- Results.Completed[!is.na(Results.Completed$GeneID), ]
            Results.Completed <- Results.Completed[!duplicated(Results.Completed$GeneID), ]
            # rownames(DEGenes.up) <- resultadosDE[which(resultadosDE$FDR < p.cutoff & resultadosDE$log2FC > log2(FC.cutoff)),
            #                                        "GeneID"]
            # rownames(DEGenes.down) <- resultadosDE[which(resultadosDE$FDR < p.cutoff & resultadosDE$log2FC < -log2(FC.cutoff)),
            #                                          "GeneID"]


            #upregulated DE
            DEGenes.up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes.up", DEGenes.up, envir = get(envir_link))

            #downregulated DE
            DEGenes.down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes.down", DEGenes.down, envir = get(envir_link))

            #Up e Down
            DEGenes.all <- resultadosDE
            assign("DEGenes.all", DEGenes.all, envir = get(envir_link))


            gene.vector.up <- as.integer(Results.Completed$GeneID%in%DEGenes.up$GeneID)
            names(gene.vector.up) <- Results.Completed$GeneID
            assign("gene.vector.up", gene.vector.up, envir = get(envir_link))

            gene.vector.down <- as.integer(Results.Completed$GeneID%in%DEGenes.down$GeneID)
            names(gene.vector.down) <- Results.Completed$GeneID
            assign("gene.vector.down", gene.vector.down, envir = get(envir_link))

            gene.vector.all <- as.integer(Results.Completed$GeneID%in%DEGenes.all$GeneID)
            names(gene.vector.all) <- Results.Completed$GeneID
            assign("gene.vector.all", gene.vector.all, envir = get(envir_link))
        } else if (tolower(ID) == "genesymbol"){

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$GeneSymbol <- annotation_table[resultadosDE$ensembl, "UCSC_GeneSymbol"]
                Results.Completed$GeneSymbol <- annotation_table[Results.Completed$ensembl, "UCSC_GeneSymbol"]
            }

            #remove any duplicates and NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneSymbol), ]
            resultadosDE[resultadosDE$GeneSymbol == "SLC35E2", "GeneSymbol"][1] <- "SLC35E2B"
            resultadosDE[duplicated(resultadosDE$GeneSymbol), "GeneSymbol"] <- paste0("?",
                                                                                      1:length(resultadosDE[duplicated(resultadosDE$GeneSymbol),
                                                                                                            "GeneSymbol"]))

            Results.Completed <- Results.Completed[!is.na(Results.Completed$GeneSymbol), ]
            Results.Completed[Results.Completed$GeneSymbol == "SLC35E2", "GeneSymbol"][1] <- "SLC35E2B"
            Results.Completed[duplicated(Results.Completed$GeneSymbol), "GeneSymbol"] <- paste0("?",
                                                                                      1:length(Results.Completed[duplicated(Results.Completed$GeneSymbol),
                                                                                                            "GeneSymbol"]))

            #upregulated DE
            DEGenes.up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes.up", DEGenes.up, envir = get(envir_link))

            #downregulated DE
            DEGenes.down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes.down", DEGenes.down, envir = get(envir_link))

            #Up e Down
            DEGenes.all <- resultadosDE
            assign("DEGenes.all", DEGenes.all, envir = get(envir_link))



            gene.vector.up <- as.integer(Results.Completed$GeneSymbol%in%DEGenes.up$GeneSymbol)
            names(gene.vector.up) <- Results.Completed$GeneSymbol
            assign("gene.vector.up", gene.vector.up, envir = get(envir_link))

            gene.vector.down <- as.integer(Results.Completed$GeneSymbol%in%DEGenes.down$GeneSymbol)
            names(gene.vector.down) <- Results.Completed$GeneSymbol
            assign("gene.vector.down", gene.vector.down, envir = get(envir_link))

            gene.vector.all <- as.integer(Results.Completed$GeneSymbol%in%DEGenes.all$GeneSymbol)
            names(gene.vector.all) <- Results.Completed$GeneSymbol
            assign("gene.vector.all", gene.vector.all, envir = get(envir_link))
        } else if(toupper(ID) == "HUGO"){

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$HUGO <- annotation_table[resultadosDE$ensembl, "HUGO"]
                Results.Completed$HUGO <- annotation_table[Results.Completed$ensembl, "HUGO"]
            } else {
                #remove any duplicates and NA
                resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneID), ]
                resultadosDE <- resultadosDE[!duplicated(resultadosDE$GeneID), ]

                resultadosDE$HUGO <- as.character(annotation_table[match(resultadosDE$GeneID,
                                                                         annotation_table$GeneID), "HUGO"])

                Results.Completed <- Results.Completed[!is.na(Results.Completed$GeneID), ]
                Results.Completed <- Results.Completed[!duplicated(Results.Completed$GeneID), ]


                Results.Completed$HUGO <- as.character(annotation_table[match(Results.Completed$GeneID,
                                                                         annotation_table$GeneID), "HUGO"])
            }

            #remove any NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$HUGO), ]

            Results.Completed <- Results.Completed[!is.na(Results.Completed$HUGO), ]



            #upregulated DE
            DEGenes.up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes.up", DEGenes.up, envir = get(envir_link))

            #downregulated DE
            DEGenes.down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes.down", DEGenes.down, envir = get(envir_link))

            #Up e Down
            DEGenes.all <- resultadosDE
            assign("DEGenes.all", DEGenes.all, envir = get(envir_link))



            gene.vector.up <- as.integer(Results.Completed$HUGO%in%DEGenes.up$HUGO)
            names(gene.vector.up) <- resultadosDE$GeneID
            assign("gene.vector.up", gene.vector.up, envir = get(envir_link))

            gene.vector.down <- as.integer(Results.Completed$HUGO%in%DEGenes.down$HUGO)
            names(gene.vector.down) <- resultadosDE$GeneID
            assign("gene.vector.down", gene.vector.down, envir = get(envir_link))

            gene.vector.all <- as.integer(Results.Completed$HUGO%in%DEGenes.all$HUGO)
            names(gene.vector.all) <- resultadosDE$GeneID
            assign("gene.vector.all", gene.vector.all, envir = get(envir_link))
        } else if (tolower(ID) == "ensembl") {

            if (tolower(dataBase) == "legacy") {
                rownames(annotation_table) <- annotation_table$GeneID
                resultadosDE$ensembl <- annotation_table[resultadosDE$GeneID, "ensembl"]

                Results.Completed$ensembl <- annotation_table[Results.Completed$GeneID, "ensembl"]

                #remove any duplicates and NA
                resultadosDE <- resultadosDE[!is.na(resultadosDE$ensembl), ]
                resultadosDE <- resultadosDE[!duplicated(resultadosDE$ensembl), ]

                Results.Completed <- Results.Completed[!is.na(Results.Completed$ensembl), ]
                Results.Completed <- Results.Completed[!duplicated(Results.Completed$ensembl), ]

                #upregulated DE
                DEGenes.up <- resultadosDE[resultadosDE$FC > 0, ]
                assign("DEGenes.up", DEGenes.up, envir = get(envir_link))

                #downregulated DE
                DEGenes.down <- resultadosDE[resultadosDE$FC < 0, ]
                assign("DEGenes.down", DEGenes.down, envir = get(envir_link))

                #Up e Down
                DEGenes.all <- resultadosDE
                assign("DEGenes.all", DEGenes.all, envir = get(envir_link))



                gene.vector.up <- as.integer(Results.Completed$ensembl%in%DEGenes.up$ensembl)
                names(gene.vector.up) <- resultadosDE$ensembl

                gene.vector.down <- as.integer(Results.Completed$ensembl%in%DEGenes.down$ensembl)
                names(gene.vector.down) <- resultadosDE$ensembl

                gene.vector.all <- as.integer(Results.Completed$ensembl%in%DEGenes.all$ensembl)
                names(gene.vector.all) <- Results.Completed$ensembl
            } else {

                #upregulated DE
                DEGenes.up <- resultadosDE[resultadosDE$FC > 0, ]
                assign("DEGenes.up", DEGenes.up, envir = get(envir_link))

                #downregulated DE
                DEGenes.down <- resultadosDE[resultadosDE$FC < 0, ]
                assign("DEGenes.down", DEGenes.down, envir = get(envir_link))

                #Up e Down
                DEGenes.all <- resultadosDE
                assign("DEGenes.all", DEGenes.all, envir = get(envir_link))


                #remove any duplicates and NA
                gene.vector.up <- as.integer(rownames(Results.Completed)%in%rownames(DEGenes.up))
                names(gene.vector.up) <- rownames(resultadosDE)

                gene.vector.down <- as.integer(rownames(Results.Completed)%in%rownames(DEGenes.down))
                names(gene.vector.down) <- rownames(resultadosDE)

                gene.vector.all <- as.integer(rownames(Results.Completed)%in%rownames(DEGenes.all))
                names(gene.vector.all) <- rownames(resultadosDE)
            }

            #remove any duplicates and NA
            assign("gene.vector.up", gene.vector.up, envir = get(envir_link))

            assign("gene.vector.down", gene.vector.down, envir = get(envir_link))

            assign("gene.vector.all", gene.vector.all, envir = get(envir_link))
        } else if (tolower(ID) == "refgene") {

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$refGene <- annotation_table[resultadosDE$ensembl, "RefSeq"]

                Results.Completed$refGene <- annotation_table[Results.Completed$ensembl, "RefSeq"]

            } else {
                #remove any duplicates and NA
                resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneID), ]
                resultadosDE <- resultadosDE[!duplicated(resultadosDE$GeneID), ]

                Results.Completed <- Results.Completed[!is.na(Results.Completed$GeneID), ]
                Results.Completed <- Results.Completed[!duplicated(Results.Completed$GeneID), ]

                resultadosDE$refGene <- as.character(annotation_table[match(resultadosDE$GeneID,
                                                                            annotation_table$GeneID), "RefSeq"])

                Results.Completed$refGene <- as.character(annotation_table[match(Results.Completed$GeneID,
                                                                            annotation_table$GeneID), "RefSeq"])
            }

            #remove any duplicates and NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$refGene), ]
            resultadosDE <- resultadosDE[!duplicated(resultadosDE$refGene), ]


            Results.Completed <- Results.Completed[!is.na(Results.Completed$refGene), ]
            Results.Completed <- Results.Completed[!duplicated(Results.Completed$refGene), ]



            #upregulated DE
            DEGenes.up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes.up", DEGenes.up, envir = get(envir_link))

            #downregulated DE
            DEGenes.down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes.down", DEGenes.down, envir = get(envir_link))

            #Up e Down
            DEGenes.all <- resultadosDE
            assign("DEGenes.all", DEGenes.all, envir = get(envir_link))



            gene.vector.up <- as.integer(Results.Completed$refGene%in%DEGenes.up$refGene)
            names(gene.vector.up) <- Results.Completed$refGene
            assign("gene.vector.up", gene.vector.up, envir = get(envir_link))

            gene.vector.down <- as.integer(Results.Completed$refGene%in%DEGenes.down$refGene)
            names(gene.vector.down) <- Results.Completed$refGene
            assign("gene.vector.down", gene.vector.down, envir = get(envir_link))

            gene.vector.all <- as.integer(Results.Completed$refGene%in%DEGenes.all$refGene)
            names(gene.vector.all) <- Results.Completed$refGene
            assign("gene.vector.all", gene.vector.all, envir = get(envir_link))
        }

        # assign("resultadosDE", resultadosDE, envir = get(envir_link))
        assign("pairName", pairName, envir = get(envir_link))
    }


    GO.wall.FUN <- function(x, n) {

        if (nrow(x) != 0) {

            if (nrow(x) >= 15){
                PlotNowGO <- x[1:15, ]
            } else{
                PlotNowGO <- x
            }

            PlotNowGO[, 2] <- gsub(" of", "", PlotNowGO[, 2])
            # large <- lengths(strsplit(PlotNowGO[, 2], "\\W+")) > 7
            large <- lengths(strsplit(PlotNowGO[, 2], " ")) > 7
            PlotNowGO[large, 2] <- unname(sapply(PlotNowGO[large, 2],
                                                    function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                                      collapse = " ")}))

            log.10.GO <- -log10(as.numeric(PlotNowGO[, "over_represented_BH"]))
            new.inf <- log.10.GO[order(log.10.GO, decreasing = TRUE)]
            log.10.GO[log.10.GO == "Inf"] <- (new.inf[new.inf != "Inf"][1]+1)

            if (tolower(image.format) == "png") {
                png(filename = paste0(DIR,
                                  "/GO_Output/GraphOutput/GOEnrichPlot_smallest_FDR_",
                                  MainNames_short[n], "_", tolower(condition), "_", ID,
                                  "_", pairName,".png"),
                    width = Width, height = Height, res = Res, units = Unit)
            } else if (tolower(image.format) == "svg") {
                svg(filename = paste0(DIR,
                                  "/GO_Output/GraphOutput/GOEnrichPlot_smallest_FDR_",
                                  MainNames_short[n], "_", tolower(condition), "_", ID,
                                  "_", pairName,".svg"),
                    width = Width, height = Height, onefile = TRUE)
            } else {
                stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
            }

            Count <- PlotNowGO$numDEInCat

            p <- ggplot2::ggplot(PlotNowGO, ggplot2::aes(x = log.10.GO,
                                                         y = forcats::fct_reorder(Description, log.10.GO))) +
                ggplot2::geom_point(ggplot2::aes(size = Count, color = GeneRatio)) +
                ggplot2::scale_colour_gradient(limits=c(min(GeneRatio), 1), low="red", high = "blue") +
                ggplot2::labs(y = "", x = "-log(FDR)", title = MainNames[n]) +
                ggplot2::theme_bw(base_size = 10) +
                ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold",
                                                                    size = 16),
                               axis.text = ggplot2::element_text(face = "bold",
                                                                 color = "#011600", size = 12),
                               title = ggplot2::element_text(face = "bold",
                                                             size = 18),
                               plot.title = ggplot2::element_text(hjust = 0.5))
            print(p)
            dev.off()

            if (tolower(image.format) == "png") {
                png(filename = paste0(DIR,
                                      "/GO_Output/GraphOutput/GOEnrichPlot_smallest_FDR_",
                                      MainNames_short[n], "_", tolower(condition), "_", ID,
                                      "_", pairName,"_2.png"),
                    width = Width*3, height = Height *1.5, res = Res, units = Unit)
            } else if (tolower(image.format) == "svg") {
                svg(filename = paste0(DIR,
                                      "/GO_Output/GraphOutput/GOEnrichPlot_smallest_FDR_",
                                      MainNames_short[n], "_", tolower(condition), "_", ID,
                                      "_", pairName,"_2.svg"),
                    width = Width, height = Height, onefile = TRUE)
            } else {
                stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
            }
            par(mar = c(5, 45, 2, 2), lwd = 2.5)
            barplot(log.10.GO,
                    names = PlotNowGO[, "term"], las = 2, horiz = TRUE,
                    xlab = "-log(FDR)", main = MainNames[n], lwd = 2,
                    cex.lab = 1.7, cex.axis = 2, cex.main=2,
                    axes = FALSE, col = RColorBrewer::brewer.pal(8,"Set1")[n],
                    border = "white",
                    cex.names = 1.5, space = 0.001)

            if (-log10(min(as.numeric(PlotNowGO[, "over_represented_BH"]))) == "Inf"){
                abline(v = 1:floor(0), col = "white", lwd = 4)
            } else {
                abline(v = 1:floor(-log10(min(as.numeric(PlotNowGO[, "over_represented_BH"]), na.rm = TRUE))),
                       col = "white", lwd = 4)
            }

            abline(v = -log10(FDR.cutoff), lwd = 3)

            axis(side = 1, 0:(max(log.10.GO)+2),
                 lwd = 4, cex.axis = 1.2)

            par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 1, 0), new=TRUE)
            plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')

            legend("topright", "-log(FDR.cutoff)", lty = 1, lwd = 4,
                   bty = "n", cex = 1)

            dev.off()

        } else {
            message(paste0("There are no terms with FDR < ", FDR.cutoff, " to plot."))
        }
    }

    # code GO ####
    message("Running GO preparations...")

    prepar_path_enrich(ID = ID,
                       pairName = pairName,
                       env = env,
                       Tool = Tool)

    message("Starting GO estimation...")

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    if (missing(Tool)){Tool <- string_vars[["envir_link"]]$Tool}

    Name <- string_vars[["envir_link"]]$Name
    pairName <- string_vars[["envir_link"]]$pairName
    dataBase <- string_vars[["envir_link"]]$dataBase

    if (grepl("crosstable", tolower(Tool))) {
        if (tolower(Tool) == "crosstable.deseq2") {
            DIR <- paste0(PATH, "/CrossData_deseq2")
        } else if (tolower(Tool) == "crosstable.edger") {
            DIR <- paste0(PATH, "/CrossData_edger")
        } else if (tolower(Tool) == "crosstable.ebseq") {
            DIR <- paste0(PATH, "/CrossData_ebseq")
        }
        dir.create(file.path(DIR, "Ontology_Results"), showWarnings = FALSE)
        DIR <- file.path(DIR, "Ontology_Results")
    } else {
        DIR <- paste0(PATH, "/Ontology_Results_", Tool, "_", toupper(Name))
    }

    dir.create(file.path(DIR, "GO_Output"), showWarnings = FALSE)
    dir.create(file.path(DIR, "GO_Output", "GraphOutput"),
               showWarnings = FALSE)
    dir.create(file.path(DIR, "GO_Output", "TextOutput"),
               showWarnings = FALSE)

    if (tolower(ID) == "geneid"){
        formatted.ID <- "knownGene"
    } else if (tolower(ID) == "genesymbol") {
        formatted.ID <- "geneSymbol"
    } else if (tolower(ID) == "ensembl") {
        formatted.ID <- "ensGene"
    } else if (tolower(ID) == "refgene") {
        formatted.ID <- "refGene"
    }

    if (tolower(image.format) == "png") {
        png(filename = paste0(DIR,
                   "/GO_Output/GraphOutput/ProbabilityWeightingFunction_",
                   tolower(condition), "_", ID, "_", pairName, ".png"),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image.format) == "svg") {
        svg(filename = paste0(DIR,
                          "/GO_Output/GraphOutput/ProbabilityWeightingFunction_",
                          tolower(condition), "_", ID, "_", pairName, ".svg"),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
    }

    if (tolower(dataBase) == "legacy") {
        genome_version <- "hg19"
    } else {
        genome_version <- "hg38"
    }

    #Probability Weighting Function
    if (tolower(condition) == "upregulated") {
        gene.vector.up <- string_vars[["envir_link"]]$gene.vector.up
        suppressPackageStartupMessages(pwf <- goseq::nullp(gene.vector.up ,
                                                           genome_version,
                                                           formatted.ID))
        dev.off()
    } else if (tolower(condition) == "downregulated") {
        gene.vector.down <- string_vars[["envir_link"]]$gene.vector.down
        suppressPackageStartupMessages(pwf <- goseq::nullp(gene.vector.down,
                                                           genome_version,
                                                           formatted.ID))
        dev.off()
    } else if (tolower(condition) == "all") {
        gene.vector.all <- string_vars[["envir_link"]]$gene.vector.all
        suppressPackageStartupMessages(pwf <- goseq::nullp(gene.vector.all,
                                                           genome_version,
                                                           formatted.ID))
        dev.off()
    }

    #GO Analysis Wallenius
    message("Running GO analysis Wallenius\n")
    suppressPackageStartupMessages(GO.wall <- goseq::goseq(pwf, genome_version, formatted.ID,
                            test.cats=c("GO:CC", "GO:BP", "GO:MF"),
                            method = "Wallenius",
                            use_genes_without_cat = use_genes_without_cat))
    # assign("GO.wall", GO.wall, envir = get(envir_link))


    #FDR
    # # Create table for FDR
    # GO.wall_WithFDR <- matrix(ncol = 8, nrow = length(row.names(GO.wall)))
    # # put col names
    # colnames(GO.wall_WithFDR) <- c(colnames(GO.wall), "Benjamini.Hochberg_FDR")
    #
    # #Fill the new table until last column
    # for(w in 1:7){
    #     GO.wall_WithFDR[, w] <- GO.wall[, w]
    # }

    #FDR Add to Go.wall
    GO.wall$over_represented_BH <- p.adjust(GO.wall$over_represented_pvalue,
                                            method = "BH")

    GO.wall$GeneRatio <- GO.wall$numDEInCat/GO.wall$numInCat

    GO.wall <- GO.wall[GO.wall$over_represented_BH < FDR.cutoff, ]

    #Escreve as tabelas
    GO.wall_BH_CC <- GO.wall[GO.wall[, "ontology"] == 'CC', ]
    GO.wall_BH_BP <- GO.wall[GO.wall[, "ontology"] == 'BP', ]
    GO.wall_BH_MF <- GO.wall[GO.wall[, "ontology"] == 'MF', ]

    suppressPackageStartupMessages(GO.wall <- goseq::goseq(pwf, genome_version, formatted.ID,
                                                           test.cats=c("KEGG"),
                                                           method = "Wallenius",
                                                           use_genes_without_cat = use_genes_without_cat))

    GO.wall$over_represented_BH <- p.adjust(GO.wall$over_represented_pvalue,
                                            method = "BH")
    GO.wall$category <- paste0("hsa", GO.wall$category)
    colnames(GO.wall)[1] <- "Pathway"
    GO.wall_BH_KEGG <- GO.wall[GO.wall[, "over_represented_BH"] < FDR.cutoff, ]

    message("Writting tables...\n")
    write.csv(GO.wall_BH_CC, file = paste0(DIR,
                                        "/GO_Output/TextOutput/CC_",
                                        tolower(condition), "_", ID,
                                        "_FDR_", FDR.cutoff,
                                        "_", pairName,".csv"))
    write.csv(GO.wall_BH_BP, file = paste0(DIR,
                                        "/GO_Output/TextOutput/BP_",
                                        tolower(condition), "_", ID,
                                        "_FDR_", FDR.cutoff,
                                        "_", pairName,".csv"))
    write.csv(GO.wall_BH_MF, file = paste0(DIR,
                                        "/GO_Output/TextOutput/MF_",
                                        tolower(condition), "_", ID,
                                        "_FDR_", FDR.cutoff,
                                        "_", pairName,".csv"))
    write.csv(GO.wall_BH_KEGG, file = paste0(DIR,
                                           "/GO_Output/TextOutput/KEGG_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR.cutoff,
                                           "_", pairName,".csv"))

    # Define category order
    MainNames <- c("Cellular component", "Biological Process", "Molecular Funcion", "KEGG")
    MainNames_short <- c("CC", "BP", "MF", "KEGG")


    GO.wall.FUN(x = GO.wall_BH_CC, n = 1)
    GO.wall.FUN(x = GO.wall_BH_BP, n = 2)
    GO.wall.FUN(x = GO.wall_BH_MF, n = 3)
    GO.wall.FUN(x = GO.wall_BH_KEGG, n = 4)
    #################### GO ENRICHMENT
    #################### END
    remove(hg19.knownGene.LENGTH, envir = .GlobalEnv)
    gc()


    # code enrichGO ####

    if (grepl("crosstable", tolower(Tool))) {
        if (tolower(Tool) == "crosstable.deseq2") {
            DIR <- paste0(PATH, "/CrossData_deseq2")
        } else if (tolower(Tool) == "crosstable.edger") {
            DIR <- paste0(PATH, "/CrossData_edger")
        } else if (tolower(Tool) == "crosstable.ebseq") {
            DIR <- paste0(PATH, "/CrossData_ebseq")
        }
        dir.create(file.path(DIR, "Ontology_Results"), showWarnings = FALSE)
        DIR <- file.path(DIR, "Ontology_Results")
    } else {
        DIR <- paste0(PATH, "/Ontology_Results_", Tool, "_", toupper(Name))
    }

    dir.create(file.path(DIR, "enrichGO_Output"), showWarnings = FALSE)
    dir.create(file.path(DIR, "enrichGO_Output", "GraphOutput"),
               showWarnings = FALSE)
    dir.create(file.path(DIR, "enrichGO_Output", "TextOutput"),
               showWarnings = FALSE)

    if (tolower(ID) == "geneid"){
        formatted.ID2 <- "ENTREZID"
    } else if (tolower(ID) == "genesymbol") {
        formatted.ID2 <- "SYMBOL"
    } else if (tolower(ID) == "ensembl") {
        formatted.ID2 <- "ENSEMBL"
    } else if (tolower(ID) == "refgene") {
        formatted.ID2 <- "REFSEQ"
    }

    # clusterProfiler::bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

    org.Hs.eg.db <- org.Hs.eg.db::org.Hs.eg.db

    erichgo <- function(Ont, Condition, n) {
        if (tolower(Condition) == "upregulated") {
            gene_vector <- string_vars[["envir_link"]]$gene.vector.up
        } else if (tolower(Condition) == "downregulated") {
            gene_vector <- string_vars[["envir_link"]]$gene.vector.down
        } else if (tolower(Condition) == "all") {
            gene_vector <- string_vars[["envir_link"]]$gene.vector.all
        }

        enrichGO_obj <- clusterProfiler::enrichGO(gene = names(gene_vector),
                                         OrgDb = org.Hs.eg.db,
                                         keytype = formatted.ID2,
                                         ont = Ont,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         # readable = TRUE,
                                         qvalueCutoff  = 0.05)

        enrichGO_resuts <- enrichGO_obj@result

        enrichGO_resuts <- enrichGO_resuts[order(enrichGO_resuts$p.adjust), ]

        if (nrow(enrichGO_resuts) >= 10){
            PlotNowGO <- enrichGO_resuts[1:10, ]
        } else{
            PlotNowGO <- enrichGO_resuts
        }

        PlotNowGO[, 2] <- gsub(" of", "", PlotNowGO[, 2])
        # large <- lengths(strsplit(PlotNowGO[, 2], "\\W+")) > 7
        large <- lengths(strsplit(PlotNowGO[, 2], " ")) > 7
        PlotNowGO[large, 2] <- unname(sapply(PlotNowGO[large, 2],
                                             function(w){paste(unlist(strsplit(w, " "))[1:7],
                                                               collapse = " ")}))

        log.10.GO <- -log10(as.numeric(PlotNowGO[, "p.adjust"]))
        new.inf <- log.10.GO[order(log.10.GO, decreasing = TRUE)]
        log.10.GO[log.10.GO == "Inf"] <- (new.inf[new.inf != "Inf"][1]+1)

        if (tolower(image.format) == "png") {
            png(filename = paste0(DIR,
                                  "/enrichGO_Output/GraphOutput/enrichGO_smallest_FDR_Ont_", Ont, "_",
                                  tolower(condition), "_", ID,
                                  "_", pairName,".png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = paste0(DIR,
                                  "/enrichGO_Output/GraphOutput/enrichGO_smallest_FDR_", Ont, "_",
                                  tolower(condition), "_", ID,
                                  "_", pairName,".svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }

        ## plot alternative to barplot
        p <- ggplot2::ggplot(PlotNowGO, ggplot2::aes(x = log.10.GO,
                                                y = forcats::fct_reorder(Description, log.10.GO))) +
            ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
            ggplot2::scale_colour_gradient(limits=c(0, FDR.cutoff), low="red", high = "blue") +
            ggplot2::labs(y = "", x = "-log(FDR)", title = MainNames[n]) +
            ggplot2::theme_bw(base_size = 10) +
            ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold",
                                                                size = 16),
                           axis.text = ggplot2::element_text(face = "bold",
                                                             color = "#011600", size = 12),
                           title = ggplot2::element_text(face = "bold",
                                                                size = 18),
                           plot.title = ggplot2::element_text(hjust = 0.5))
        print(p)

        dev.off()

        return(enrichGO_resuts)
    }

    erichgo_CC <- erichgo(Ont = "CC", condition, n = 1)
    erichgo_BP <- erichgo(Ont = "BP", condition, n = 2)
    erichgo_MF <- erichgo(Ont = "MF", condition, n = 3)

    message("Writting tables...\n")
    write.csv(erichgo_CC, file = paste0(DIR,
                                           "/enrichGO_Output/TextOutput/CC_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR.cutoff,
                                           "_", pairName,".csv"))
    write.csv(erichgo_BP, file = paste0(DIR,
                                           "/enrichGO_Output/TextOutput/BP_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR.cutoff,
                                           "_", pairName,".csv"))
    write.csv(erichgo_MF, file = paste0(DIR,
                                           "/enrichGO_Output/TextOutput/MF_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR.cutoff,
                                           "_", pairName,".csv"))

    # ggo <- clusterProfiler::groupGO(gene = names(gene.vector.up),
    #                                 OrgDb = org.Hs.eg.db,
    #                                 keytype = formatted.ID2,
    #                                 ont = "CC",
    #                                 readable = TRUE) #Gene ID will be mapped to gene Symbol
    # cnetplot
    # DOSE::cnetplot(GO2, categorySize="pvalue",
    #                foldChange=resultadosDE$FC[gene.vector.up == 1],
    #                showCategory = 1,
    #                fixed = TRUE)

    # ego3 <- clusterProfiler::gseGO(geneList = names(gene.vector.up)[sort(names(gene.vector.up), decreasing = TRUE)],
    #               OrgDb = org.Hs.eg.db,
    #               ont = "CC",
    #               nPerm = 1000,
    #               minGSSize = 100,
    #               maxGSSize = 500,
    #               pvalueCutoff = 0.05,
    #               keytype = formatted.ID2,
    #               verbose = FALSE)
#
    # kk2 <- clusterProfiler::gseKEGG(geneList = names(gene.vector.up)[sort(names(gene.vector.up), decreasing = TRUE)],
    #                organism = 'hsa',
    #                nPerm = 1000,
    #                minGSSize = 120,
    #                pvalueCutoff = 0.05,
    #                verbose = FALSE)
#
    # clusterProfiler::enrichDAVID(gene = names(gene.vector.up),
    #                      idType = "ENTREZ_GENE_ID",
    #                      listType = "Gene",
    #                      annotation = "KEGG_PATHWAY",
    #                      david.user = "clusterProfiler@hku.hk")

    message("Done!\n")
 }
