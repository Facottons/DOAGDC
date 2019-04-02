#' Perform Gene-Ontology Pathways Enrichment
#'
#' @param condition A character string containing which condition should be
#'   used: "Upregulated", "Downregulated" or "all".
#' @param use_genes_without_cat Logical value where \code{FALSE} indicate that
#'   genes outside the category being tested will be ignored in the calculation
#'   of p-values. The default is "FALSE".
#' @param FDR_cutoff
#' @param Width,Height,Res,Unit,image_format
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
                  FDR_cutoff = 0.05,
                  Width = 2000,
                  Height = 2000,
                  Res = 300,
                  Unit = "px",
                  image_format = "png",
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
        #   matrix. It must be named \code{Results_Completed.} and it must be inside
        #   the given \code{env}. Use the function \code{table2GDCtools} to allocate
        #   your table inside the given \code{env}.

        if(missing(env)){
            stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))
        }

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
            dir.create(file.path(DIR, paste0("Ontology_Results", tolower(groupGen))), showWarnings = FALSE)
            DIR <- file.path(DIR, paste0("Ontology_Results", tolower(groupGen)))

            File <- "resultadosDE_crossed"
            File2 <- "Results_Completed.crossed"

        } else {
            File <- paste("resultadosDE", Tool, sep = ".")
            File2 <- paste("Results_Completed", Tool, sep = ".")

            dir.create(paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_",
                              Tool, "_", toupper(Name)), showWarnings = FALSE)
            DIR <- paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_",
                          Tool, "_", toupper(Name))
        }

        # generate annotation table
        annotation_table <- GDCtools::annotation_table
        # GO_annotations <- GDCtools::GO_annotations

        #input DE genes
        resultadosDE <- get(File, envir = string_vars[["envir_link"]])[[pairName]]
        Results_Completed <- get(File2, envir = string_vars[["envir_link"]])[[pairName]]

        #pre-GO        gbutils::isNA
        if (tolower(ID) == "geneid") {

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$GeneID <- annotation_table[resultadosDE$ensembl, "GeneID"]
                Results_Completed$GeneID <- annotation_table[Results_Completed$ensembl, "GeneID"]


                # If you haven't already installed devtools...
                # install.packages("devtools")

                # Use devtools to install the package
                # devtools::install_github("stephenturner/annotables")
                # library(annotables)
                # head(grch38)

                resultadosDE$ensembl <- gsub(pattern = "\\..*", "", rownames(resultadosDE))
                Results_Completed$ensembl <- gsub(pattern = "\\..*", "", rownames(Results_Completed))
                annotation_table <- GDCtools::annotation_table
                selected <- intersect(annotation_table$ensembl, resultadosDE$ensembl)
                geneID <- unique(annotation_table[annotation_table$ensembl %in% selected, c("ensembl", "GeneID")])
                # Outer_join
                resultadosDE <- merge(x = resultadosDE, y = geneID, by = "ensembl", all = TRUE)
                selected <- intersect(annotation_table$ensembl, Results_Completed$ensembl)
                geneID <- unique(annotation_table[annotation_table$ensembl %in% selected, c("ensembl", "GeneID")])
                # Outer_join
                Results_Completed <- merge(x = Results_Completed, y = geneID, by = "ensembl", all = TRUE)



            }

            #remove any duplicates and NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneID), ]
            resultadosDE <- resultadosDE[!duplicated(resultadosDE$GeneID), ]

            Results_Completed <- Results_Completed[!is.na(Results_Completed$GeneID), ]
            Results_Completed <- Results_Completed[!duplicated(Results_Completed$GeneID), ]


            #upregulated DE
            DEGenes_up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes_up", DEGenes_up, envir = get(envir_link))

            #downregulated DE
            DEGenes_down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes_down", DEGenes_down, envir = get(envir_link))

            #Up e Down
            DEGenes_all <- resultadosDE
            assign("DEGenes_all", DEGenes_all, envir = get(envir_link))


            gene_vector_up <- as.integer(Results_Completed$GeneID%in%DEGenes_up$GeneID)
            names(gene_vector_up) <- Results_Completed$GeneID
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            gene_vector_down <- as.integer(Results_Completed$GeneID%in%DEGenes_down$GeneID)
            names(gene_vector_down) <- Results_Completed$GeneID
            assign("gene_vector_down", gene_vector_down, envir = get(envir_link))

            gene_vector_all <- as.integer(Results_Completed$GeneID%in%DEGenes_all$GeneID)
            names(gene_vector_all) <- Results_Completed$GeneID
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (tolower(ID) == "genesymbol"){

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$GeneSymbol <- annotation_table[resultadosDE$ensembl, "UCSC_GeneSymbol"]
                Results_Completed$GeneSymbol <- annotation_table[Results_Completed$ensembl, "UCSC_GeneSymbol"]
            }

            #remove any duplicates and NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneSymbol), ]
            resultadosDE[resultadosDE$GeneSymbol == "SLC35E2", "GeneSymbol"][1] <- "SLC35E2B"
            resultadosDE[duplicated(resultadosDE$GeneSymbol),
                         "GeneSymbol"] <- paste0("?", 1:length(resultadosDE[duplicated(resultadosDE$GeneSymbol),
                                                                                                            "GeneSymbol"]))

            Results_Completed <- Results_Completed[!is.na(Results_Completed$GeneSymbol), ]
            Results_Completed[Results_Completed$GeneSymbol == "SLC35E2", "GeneSymbol"][1] <- "SLC35E2B"
            Results_Completed[duplicated(Results_Completed$GeneSymbol),
                              "GeneSymbol"] <- paste0("?",
                                                      1:length(Results_Completed[duplicated(Results_Completed$GeneSymbol),
                                                                                                            "GeneSymbol"]))

            #upregulated DE
            DEGenes_up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes_up", DEGenes_up, envir = get(envir_link))

            #downregulated DE
            DEGenes_down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes_down", DEGenes_down, envir = get(envir_link))

            #Up e Down
            DEGenes_all <- resultadosDE
            assign("DEGenes_all", DEGenes_all, envir = get(envir_link))



            gene_vector_up <- as.integer(Results_Completed$GeneSymbol%in%DEGenes_up$GeneSymbol)
            names(gene_vector_up) <- Results_Completed$GeneSymbol
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            gene_vector_down <- as.integer(Results_Completed$GeneSymbol%in%DEGenes_down$GeneSymbol)
            names(gene_vector_down) <- Results_Completed$GeneSymbol
            assign("gene_vector_down", gene_vector_down, envir = get(envir_link))

            gene_vector_all <- as.integer(Results_Completed$GeneSymbol%in%DEGenes_all$GeneSymbol)
            names(gene_vector_all) <- Results_Completed$GeneSymbol
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if(toupper(ID) == "HUGO"){

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$HUGO <- annotation_table[resultadosDE$ensembl, "HUGO"]
                Results_Completed$HUGO <- annotation_table[Results_Completed$ensembl, "HUGO"]
            } else {
                #remove any duplicates and NA
                resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneID), ]
                resultadosDE <- resultadosDE[!duplicated(resultadosDE$GeneID), ]

                resultadosDE$HUGO <- as.character(annotation_table[match(resultadosDE$GeneID,
                                                                         annotation_table$GeneID), "HUGO"])

                Results_Completed <- Results_Completed[!is.na(Results_Completed$GeneID), ]
                Results_Completed <- Results_Completed[!duplicated(Results_Completed$GeneID), ]


                Results_Completed$HUGO <- as.character(annotation_table[match(Results_Completed$GeneID,
                                                                         annotation_table$GeneID), "HUGO"])
            }

            #remove any NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$HUGO), ]

            Results_Completed <- Results_Completed[!is.na(Results_Completed$HUGO), ]



            #upregulated DE
            DEGenes_up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes_up", DEGenes_up, envir = get(envir_link))

            #downregulated DE
            DEGenes_down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes_down", DEGenes_down, envir = get(envir_link))

            #Up e Down
            DEGenes_all <- resultadosDE
            assign("DEGenes_all", DEGenes_all, envir = get(envir_link))



            gene_vector_up <- as.integer(Results_Completed$HUGO%in%DEGenes_up$HUGO)
            names(gene_vector_up) <- resultadosDE$GeneID
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            gene_vector_down <- as.integer(Results_Completed$HUGO%in%DEGenes_down$HUGO)
            names(gene_vector_down) <- resultadosDE$GeneID
            assign("gene_vector_down", gene_vector_down, envir = get(envir_link))

            gene_vector_all <- as.integer(Results_Completed$HUGO%in%DEGenes_all$HUGO)
            names(gene_vector_all) <- resultadosDE$GeneID
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (tolower(ID) == "ensembl") {

            if (tolower(dataBase) == "legacy") {
                rownames(annotation_table) <- annotation_table$GeneID
                resultadosDE$ensembl <- annotation_table[resultadosDE$GeneID, "ensembl"]

                Results_Completed$ensembl <- annotation_table[Results_Completed$GeneID, "ensembl"]

                #remove any duplicates and NA
                resultadosDE <- resultadosDE[!is.na(resultadosDE$ensembl), ]
                resultadosDE <- resultadosDE[!duplicated(resultadosDE$ensembl), ]

                Results_Completed <- Results_Completed[!is.na(Results_Completed$ensembl), ]
                Results_Completed <- Results_Completed[!duplicated(Results_Completed$ensembl), ]

                #upregulated DE
                DEGenes_up <- resultadosDE[resultadosDE$FC > 0, ]
                assign("DEGenes_up", DEGenes_up, envir = get(envir_link))

                #downregulated DE
                DEGenes_down <- resultadosDE[resultadosDE$FC < 0, ]
                assign("DEGenes_down", DEGenes_down, envir = get(envir_link))

                #Up e Down
                DEGenes_all <- resultadosDE
                assign("DEGenes_all", DEGenes_all, envir = get(envir_link))



                gene_vector_up <- as.integer(Results_Completed$ensembl%in%DEGenes_up$ensembl)
                names(gene_vector_up) <- resultadosDE$ensembl

                gene_vector_down <- as.integer(Results_Completed$ensembl%in%DEGenes_down$ensembl)
                names(gene_vector_down) <- resultadosDE$ensembl

                gene_vector_all <- as.integer(Results_Completed$ensembl%in%DEGenes_all$ensembl)
                names(gene_vector_all) <- Results_Completed$ensembl
            } else {

                #upregulated DE
                DEGenes_up <- resultadosDE[resultadosDE$FC > 0, ]
                assign("DEGenes_up", DEGenes_up, envir = get(envir_link))

                #downregulated DE
                DEGenes_down <- resultadosDE[resultadosDE$FC < 0, ]
                assign("DEGenes_down", DEGenes_down, envir = get(envir_link))

                #Up e Down
                DEGenes_all <- resultadosDE
                assign("DEGenes_all", DEGenes_all, envir = get(envir_link))


                #remove any duplicates and NA
                gene_vector_up <- as.integer(rownames(Results_Completed)%in%rownames(DEGenes_up))
                names(gene_vector_up) <- rownames(resultadosDE)

                gene_vector_down <- as.integer(rownames(Results_Completed)%in%rownames(DEGenes_down))
                names(gene_vector_down) <- rownames(resultadosDE)

                gene_vector_all <- as.integer(rownames(Results_Completed)%in%rownames(DEGenes_all))
                names(gene_vector_all) <- rownames(resultadosDE)
            }

            #remove any duplicates and NA
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            assign("gene_vector_down", gene_vector_down, envir = get(envir_link))

            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (tolower(ID) == "refgene") {

            if (tolower(dataBase) == "gdc") {
                rownames(annotation_table) <- annotation_table$ensembl
                resultadosDE$refGene <- annotation_table[resultadosDE$ensembl, "RefSeq"]

                Results_Completed$refGene <- annotation_table[Results_Completed$ensembl, "RefSeq"]

            } else {
                #remove any duplicates and NA
                resultadosDE <- resultadosDE[!is.na(resultadosDE$GeneID), ]
                resultadosDE <- resultadosDE[!duplicated(resultadosDE$GeneID), ]

                Results_Completed <- Results_Completed[!is.na(Results_Completed$GeneID), ]
                Results_Completed <- Results_Completed[!duplicated(Results_Completed$GeneID), ]

                resultadosDE$refGene <- as.character(annotation_table[match(resultadosDE$GeneID,
                                                                            annotation_table$GeneID), "RefSeq"])

                Results_Completed$refGene <- as.character(annotation_table[match(Results_Completed$GeneID,
                                                                            annotation_table$GeneID), "RefSeq"])
            }

            #remove any duplicates and NA
            resultadosDE <- resultadosDE[!is.na(resultadosDE$refGene), ]
            resultadosDE <- resultadosDE[!duplicated(resultadosDE$refGene), ]


            Results_Completed <- Results_Completed[!is.na(Results_Completed$refGene), ]
            Results_Completed <- Results_Completed[!duplicated(Results_Completed$refGene), ]



            #upregulated DE
            DEGenes_up <- resultadosDE[resultadosDE$FC > 0, ]
            assign("DEGenes_up", DEGenes_up, envir = get(envir_link))

            #downregulated DE
            DEGenes_down <- resultadosDE[resultadosDE$FC < 0, ]
            assign("DEGenes_down", DEGenes_down, envir = get(envir_link))

            #Up e Down
            DEGenes_all <- resultadosDE
            assign("DEGenes_all", DEGenes_all, envir = get(envir_link))



            gene_vector_up <- as.integer(Results_Completed$refGene%in%DEGenes_up$refGene)
            names(gene_vector_up) <- Results_Completed$refGene
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            gene_vector_down <- as.integer(Results_Completed$refGene%in%DEGenes_down$refGene)
            names(gene_vector_down) <- Results_Completed$refGene
            assign("gene_vector_down", gene_vector_down, envir = get(envir_link))

            gene_vector_all <- as.integer(Results_Completed$refGene%in%DEGenes_all$refGene)
            names(gene_vector_all) <- Results_Completed$refGene
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        }

        # assign("resultadosDE", resultadosDE, envir = get(envir_link))
        assign("pairName", pairName, envir = get(envir_link))
    }


    GO_wall_FUN <- function(x, n, KEGG = FALSE) {

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
            new_inf <- log.10.GO[order(log.10.GO, decreasing = TRUE)]
            log.10.GO[log.10.GO == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

            longest_word <- max(stringr::str_count(PlotNowGO$term))

            if (longest_word > 50) {
                longest_Width <- Width * 1.2
            } else {
                longest_Width <- Width
            }

            if (KEGG){
                PlotNowGO$term <- PlotNowGO$Pathway

                if (!is.numeric(PlotNowGO$GeneRatio)) {
                    skip <- TRUE
                } else {
                    skip <- FALSE
                }
            } else {
                skip <- FALSE
            }

            if (!skip) {
                if (tolower(image_format) == "png") {
                    png(filename = paste0(DIR,
                                          "/GO_Output/GraphOutput/GOEnrichPlot_smallest_FDR_",
                                          MainNames_short[n], "_", tolower(condition), "_", ID,
                                          "_", pairName,".png"),
                        width = longest_Width, height = Height, res = Res, units = Unit)
                } else if (tolower(image_format) == "svg") {
                    svg(filename = paste0(DIR,
                                          "/GO_Output/GraphOutput/GOEnrichPlot_smallest_FDR_",
                                          MainNames_short[n], "_", tolower(condition), "_", ID,
                                          "_", pairName,".svg"),
                        width = longest_Width, height = Height, onefile = TRUE)
                } else {
                    stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
                }

                Count <- PlotNowGO$numDEInCat
                Gene_Ratio <- round(as.numeric(PlotNowGO$GeneRatio), 2)

                p <- ggplot2::ggplot(PlotNowGO, ggplot2::aes(x = log.10.GO,
                                                             y = forcats::fct_reorder(term, log.10.GO))) +
                    ggplot2::geom_point(ggplot2::aes(size = Count, color = Gene_Ratio)) +
                    ggplot2::scale_colour_gradient(limits=c(0, 1),
                                                   low="red", high = "blue") +
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

            }

        } else {
            message(paste0("There are no terms with FDR < ", FDR_cutoff, " to plot."))
        }
    }

    # code GO ####

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
    dataBase <- string_vars[["envir_link"]]$dataBase
    groupGen <- string_vars[["envir_link"]]$groupGen

    message("Running GO preparations...")

    prepar_path_enrich(ID = ID,
                       pairName = pairName,
                       env = env,
                       Tool = Tool)

    message("Starting GO estimation...")
    if (grepl("crosstable", tolower(Tool))) {
        if (tolower(Tool) == "crosstable.deseq2") {
            DIR <- paste0(PATH, "/CrossData_deseq2")
        } else if (tolower(Tool) == "crosstable.edger") {
            DIR <- paste0(PATH, "/CrossData_edger")
        } else if (tolower(Tool) == "crosstable.ebseq") {
            DIR <- paste0(PATH, "/CrossData_ebseq")
        }
        dir.create(file.path(DIR, paste0("Ontology_Results", tolower(groupGen))), showWarnings = FALSE)
        DIR <- file.path(DIR, paste0("Ontology_Results", tolower(groupGen)))
    } else {
        dir.create(paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_",
                          Tool, "_", toupper(Name)), showWarnings = FALSE)
        DIR <- paste0(PATH, "/Ontology_Results_", tolower(groupGen), "_", Tool, "_", toupper(Name))
    }

    dir.create(file.path(DIR, "GO_Output"), showWarnings = FALSE)
    dir.create(file.path(DIR, "GO_Output", "GraphOutput"),
               showWarnings = FALSE)
    dir.create(file.path(DIR, "GO_Output", "TextOutput"),
               showWarnings = FALSE)

    if (tolower(ID) == "geneid"){
        formatted_ID <- "knownGene"
    } else if (tolower(ID) == "genesymbol") {
        formatted_ID <- "geneSymbol"
    } else if (tolower(ID) == "ensembl") {
        formatted_ID <- "ensGene"
    } else if (tolower(ID) == "refgene") {
        formatted_ID <- "refGene"
    }

    if (tolower(image_format) == "png") {
        png(filename = paste0(DIR,
                   "/GO_Output/GraphOutput/ProbabilityWeightingFunction_",
                   tolower(condition), "_", ID, "_", pairName, ".png"),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = paste0(DIR,
                          "/GO_Output/GraphOutput/ProbabilityWeightingFunction_",
                          tolower(condition), "_", ID, "_", pairName, ".svg"),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }

    if (tolower(dataBase) == "legacy") {
        genome_version <- "hg19"
    } else {
        genome_version <- "hg38"
    }

    #Probability Weighting Function
    if (tolower(condition) == "upregulated") {
        gene_vector_up <- string_vars[["envir_link"]]$gene_vector_up
        suppressPackageStartupMessages(pwf <- goseq::nullp(gene_vector_up ,
                                                           genome_version,
                                                           formatted_ID))
        dev.off()
    } else if (tolower(condition) == "downregulated") {
        gene_vector_down <- string_vars[["envir_link"]]$gene_vector_down
        suppressPackageStartupMessages(pwf <- goseq::nullp(gene_vector_down,
                                                           genome_version,
                                                           formatted_ID))
        dev.off()
    } else if (tolower(condition) == "all") {
        gene_vector_all <- string_vars[["envir_link"]]$gene_vector_all
        suppressPackageStartupMessages(pwf <- goseq::nullp(gene_vector_all,
                                                           genome_version,
                                                           formatted_ID))
        dev.off()
    }

    #GO Analysis Wallenius
    message("Running GO analysis Wallenius\n")
    suppressPackageStartupMessages(GO_wall <- goseq::goseq(pwf, genome_version, formatted_ID,
                            test.cats=c("GO:CC", "GO:BP", "GO:MF"),
                            method = "Wallenius",
                            use_genes_without_cat = use_genes_without_cat))
    # assign("GO_wall", GO_wall, envir = get(envir_link))


    #FDR
    # # Create table for FDR
    # GO_wall_WithFDR <- matrix(ncol = 8, nrow = length(row.names(GO_wall)))
    # # put col names
    # colnames(GO_wall_WithFDR) <- c(colnames(GO_wall), "Benjamini.Hochberg_FDR")
    #
    # #Fill the new table until last column
    # for(w in 1:7){
    #     GO_wall_WithFDR[, w] <- GO_wall[, w]
    # }

    #FDR Add to GO_wall
    GO_wall$over_represented_BH <- p.adjust(GO_wall$over_represented_pvalue,
                                            method = "BH")

    GO_wall$GeneRatio <- GO_wall$numDEInCat/GO_wall$numInCat

    GO_wall <- GO_wall[GO_wall$over_represented_BH < FDR_cutoff, ]

    #Escreve as tabelas
    GO_wall_BH_CC <- GO_wall[GO_wall[, "ontology"] == 'CC', ]
    GO_wall_BH_BP <- GO_wall[GO_wall[, "ontology"] == 'BP', ]
    GO_wall_BH_MF <- GO_wall[GO_wall[, "ontology"] == 'MF', ]

    # Define category order
    MainNames <- c("Cellular component", "Biological Process", "Molecular Funcion", "KEGG")
    MainNames_short <- c("CC", "BP", "MF", "KEGG")

    GO_wall_FUN(x = GO_wall_BH_CC, n = 1)
    GO_wall_FUN(x = GO_wall_BH_BP, n = 2)
    GO_wall_FUN(x = GO_wall_BH_MF, n = 3)

    suppressPackageStartupMessages(GO_wall <- goseq::goseq(pwf, genome_version, formatted_ID,
                                                           test.cats=c("KEGG"),
                                                           method = "Wallenius",
                                                           use_genes_without_cat = use_genes_without_cat))

    GO_wall$over_represented_BH <- p.adjust(GO_wall$over_represented_pvalue,
                                            method = "BH")
    GO_wall$category <- paste0("hsa", GO_wall$category)
    colnames(GO_wall)[1] <- "Pathway"
    GO_wall$GeneRatio <- GO_wall$numDEInCat/GO_wall$numInCat
    GO_wall_BH_KEGG <- GO_wall[GO_wall[, "over_represented_BH"] < FDR_cutoff, ]

    message("Writting tables...\n")
    write.csv(GO_wall_BH_CC, file = paste0(DIR,
                                        "/GO_Output/TextOutput/CC_",
                                        tolower(condition), "_", ID,
                                        "_FDR_", FDR_cutoff,
                                        "_", pairName,".csv"))
    write.csv(GO_wall_BH_BP, file = paste0(DIR,
                                        "/GO_Output/TextOutput/BP_",
                                        tolower(condition), "_", ID,
                                        "_FDR_", FDR_cutoff,
                                        "_", pairName,".csv"))
    write.csv(GO_wall_BH_MF, file = paste0(DIR,
                                        "/GO_Output/TextOutput/MF_",
                                        tolower(condition), "_", ID,
                                        "_FDR_", FDR_cutoff,
                                        "_", pairName,".csv"))
    write.csv(GO_wall_BH_KEGG, file = paste0(DIR,
                                           "/GO_Output/TextOutput/KEGG_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR_cutoff,
                                           "_", pairName,".csv"))

    GO_wall_FUN(x = GO_wall_BH_KEGG, n = 4, KEGG = TRUE)
    #################### GO ENRICHMENT
    #################### END
    remove(hg19.knownGene.LENGTH, envir = .GlobalEnv)
    gc()


    # code enrichGO ####

    # if (grepl("crosstable", tolower(Tool))) {
    #     if (tolower(Tool) == "crosstable.deseq2") {
    #         DIR <- paste0(PATH, "/CrossData_deseq2")
    #     } else if (tolower(Tool) == "crosstable.edger") {
    #         DIR <- paste0(PATH, "/CrossData_edger")
    #     } else if (tolower(Tool) == "crosstable.ebseq") {
    #         DIR <- paste0(PATH, "/CrossData_ebseq")
    #     }
    #     dir.create(file.path(DIR, "Ontology_Results"), showWarnings = FALSE)
    #     DIR <- file.path(DIR, "Ontology_Results")
    # } else {
    #     DIR <- paste0(PATH, "/Ontology_Results_", Tool, "_", toupper(Name))
    # }

    dir.create(file.path(DIR, "enrichGO_Output"), showWarnings = FALSE)
    dir.create(file.path(DIR, "enrichGO_Output", "GraphOutput"),
               showWarnings = FALSE)
    dir.create(file.path(DIR, "enrichGO_Output", "TextOutput"),
               showWarnings = FALSE)

    if (tolower(ID) == "geneid"){
        formatted_ID2 <- "ENTREZID"
    } else if (tolower(ID) == "genesymbol") {
        formatted_ID2 <- "SYMBOL"
    } else if (tolower(ID) == "ensembl") {
        formatted_ID2 <- "ENSEMBL"
    } else if (tolower(ID) == "refgene") {
        formatted_ID2 <- "REFSEQ"
    }

    # clusterProfiler::bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

    erichgo <- function(Ont, Condition, n) {
        if (tolower(Condition) == "upregulated") {
            gene_vector <- string_vars[["envir_link"]]$gene_vector_up
        } else if (tolower(Condition) == "downregulated") {
            gene_vector <- string_vars[["envir_link"]]$gene_vector_down
        } else if (tolower(Condition) == "all") {
            gene_vector <- string_vars[["envir_link"]]$gene_vector_all
        }

        enrichGO_obj <- clusterProfiler::enrichGO(gene = names(gene_vector)[gene_vector != 0],
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                         keytype = formatted_ID2,
                                         ont = Ont,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = FDR_cutoff)

        enrichGO_resuts <- enrichGO_obj@result

        if (nrow(enrichGO_resuts) > 0) {
            enrichGO_resuts <- enrichGO_resuts[order(enrichGO_resuts$p.adjust), ]

            ratios_splited <- unname(sapply(enrichGO_resuts$GeneRatio, function(w){
                unlist(strsplit(w, "\\/"))}))

            num_ratios <- as.numeric(ratios_splited[1, ])/as.numeric(ratios_splited[2, ])


            enrichGO_resuts$GeneRatioNum <- num_ratios


            if (nrow(enrichGO_resuts) >= 15){
                PlotNowGO <- enrichGO_resuts[1:15, ]
            } else{
                PlotNowGO <- enrichGO_resuts
            }

            PlotNowGO[, 2] <- gsub(" of", "", PlotNowGO[, 2])
            # large <- lengths(strsplit(PlotNowGO[, 2], "\\W+")) > 7
            PlotNowGO[, 2] <- gsub("-", " ", PlotNowGO[, 2])
            large <- lengths(strsplit(PlotNowGO[, 2], " ")) > 3
            PlotNowGO[large, 2] <- unname(sapply(PlotNowGO[large, 2],
                                                 function(w){paste(unlist(strsplit(w, " "))[1:4],
                                                                   collapse = " ")}))

            log.10.GO <- -log10(as.numeric(PlotNowGO[, "p.adjust"]))
            new_inf <- log.10.GO[order(log.10.GO, decreasing = TRUE)]
            log.10.GO[log.10.GO == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

            longest_word <- max(stringr::str_count(PlotNowGO$Description))

            if (longest_word > 40) {
                longest_Width <- Width * 1.2
            } else {
                longest_Width <- Width
            }

            if (nrow(enrichGO_resuts) > 0) {
                if (tolower(image_format) == "png") {
                    png(filename = paste0(DIR,
                                          "/enrichGO_Output/GraphOutput/enrichGO_smallest_FDR_Ont_", Ont, "_",
                                          tolower(condition), "_", ID,
                                          "_", pairName,".png"),
                        width = longest_Width, height = Height, res = Res, units = Unit)
                } else if (tolower(image_format) == "svg") {
                    svg(filename = paste0(DIR,
                                          "/enrichGO_Output/GraphOutput/enrichGO_smallest_FDR_", Ont, "_",
                                          tolower(condition), "_", ID,
                                          "_", pairName,".svg"),
                        width = longest_Width, height = Height, onefile = TRUE)
                } else {
                    stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
                }

                # Count <- PlotNowGO$numDEInCat
                Gene_Ratio <- round(as.numeric(PlotNowGO$GeneRatioNum), 2)

                ## plot alternative to barplot
                p <- ggplot2::ggplot(PlotNowGO, ggplot2::aes(x = log.10.GO,
                                                             y = forcats::fct_reorder(Description, log.10.GO))) +
                    ggplot2::geom_point(ggplot2::aes(size = Count, color = Gene_Ratio)) +
                    ggplot2::scale_colour_gradient(limits=c(0, 1),
                                                   low="red", high = "blue") +
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
            }
        }

        return(enrichGO_resuts)
    }

    erichgo_CC <- erichgo(Ont = "CC", condition, n = 1)
    erichgo_BP <- erichgo(Ont = "BP", condition, n = 2)
    erichgo_MF <- erichgo(Ont = "MF", condition, n = 3)

    message("Writting tables...\n")
    write.csv(erichgo_CC, file = paste0(DIR,
                                           "/enrichGO_Output/TextOutput/CC_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR_cutoff,
                                           "_", pairName,".csv"))
    write.csv(erichgo_BP, file = paste0(DIR,
                                           "/enrichGO_Output/TextOutput/BP_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR_cutoff,
                                           "_", pairName,".csv"))
    write.csv(erichgo_MF, file = paste0(DIR,
                                           "/enrichGO_Output/TextOutput/MF_",
                                           tolower(condition), "_", ID,
                                           "_FDR_", FDR_cutoff,
                                           "_", pairName,".csv"))

    # ggo <- clusterProfiler::groupGO(gene = names(gene_vector_up),
    #                                 OrgDb = org.Hs.eg.db,
    #                                 keytype = formatted_ID2,
    #                                 ont = "CC",
    #                                 readable = TRUE) #Gene ID will be mapped to gene Symbol
    # cnetplot
    # DOSE::cnetplot(GO2, categorySize="pvalue",
    #                foldChange=resultadosDE$FC[gene_vector_up == 1],
    #                showCategory = 1,
    #                fixed = TRUE)

    # ego3 <- clusterProfiler::gseGO(geneList = names(gene_vector_up)[sort(names(gene_vector_up), decreasing = TRUE)],
    #               OrgDb = org.Hs.eg.db,
    #               ont = "CC",
    #               nPerm = 1000,
    #               minGSSize = 100,
    #               maxGSSize = 500,
    #               pvalueCutoff = 0.05,
    #               keytype = formatted_ID2,
    #               verbose = FALSE)
#
    # kk2 <- clusterProfiler::gseKEGG(geneList = names(gene_vector_up)[sort(names(gene_vector_up), decreasing = TRUE)],
    #                organism = 'hsa',
    #                nPerm = 1000,
    #                minGSSize = 120,
    #                pvalueCutoff = 0.05,
    #                verbose = FALSE)
#
    # clusterProfiler::enrichDAVID(gene = names(gene_vector_up),
    #                      idType = "ENTREZ_GENE_ID",
    #                      listType = "Gene",
    #                      annotation = "KEGG_PATHWAY",
    #                      david.user = "clusterProfiler@hku.hk")

    message("Done!\n")
 }
