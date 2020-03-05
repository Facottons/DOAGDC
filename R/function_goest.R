#' Perform Gene-Ontology Pathways Enrichment
#'
#' @param condition A character string containing which condition should be
#'    used: "Upregulated", "Downregulated" or "all".
#' @param use_genes_without_cat Logical value where \code{FALSE} indicate that
#'    genes outside the category being tested will be ignored in the calculation
#'    of p-values. The default is "FALSE".
#' @param fdr_cutoff
#' @param width,height,res,unit,image_format
#' @param tool A character string indicating which differential expression
#'    analysis tool was last used.
#' @param id A character string indicating which id should be used: "hugo",
#'    "gene_symbol", "ensembl" , "refGene" or "geneid". The default is
#'    \code{"geneid"}.
#' @param pair_name
#' @param env
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_ebseq
#' @inheritParams concatenate_expression
#'
#' @return Enriched terms.
#' @export
#'
#' @importFrom stringr str_count
#' @importFrom goseq nullp
#' @importFrom goseq goseq
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @examples
#' # data already downloaded using the 'download_gdc' function
#' concatenate_exon("gene",
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     work_dir = "~/Desktop"
#' )
#'
#' # separating gene HIF3A expression data patients in two groups
#' groups_identification_mclust("gene", 2,
#'     name = "HIF3A",
#'     modelName = "E",
#'     env = CHOL_LEGACY_gene_tumor_data,
#'     tumor = "CHOL"
#' )
#'
#' # load not normalized data
#' concatenate_exon("gene",
#'     normalization = FALSE,
#'     name = "HIF3A",
#'     data_base = "legacy",
#'     tumor = "CHOL",
#'     env = CHOL_LEGACY_gene_tumor_data,
#'     work_dir = "~/Desktop"
#' )
#'
#' # start DE analysis
#' # considering concatenate_exon and groups_identification already runned
#' dea_edger(
#'     name = "HIF3A",
#'     group_gen = "mclust",
#'     env = CHOL_LEGACY_gene_tumor_data
#' )
#'
#' gonto(
#'     condition = "Upregulated",
#'     tool = "edgeR", env = CHOL_LEGACY_gene_tumor_data
#' )
gonto <- function(condition,
                use_genes_without_cat = TRUE,
                fdr_cutoff = 0.05,
                width = 2000,
                height = 2000,
                res = 300,
                unit = "px",
                image_format = "png",
                tool,
                id = "geneID",
                pair_name = "G2_over_G1",
                env) {

    # local functions ####
    prepar_path_enrich <- function(id = "geneid",
                                    pair_name = "G2_over_G1",
                                    env,
                                    tool) {

        if (grepl("crosstable", tool)) {
            if (tool == "crosstable_deseq2") {
                dir <- paste0(path, "/CrossData_deseq2")
            } else if (tool == "crosstable_edger") {
                dir <- paste0(path, "/CrossData_edger")
            } else if (tool == "crosstable_ebseq") {
                dir <- paste0(path, "/CrossData_ebseq")
            }

            if (grepl("co_exp", tool)) {
                file <- "resultados_de_crossed_Co"
                dir <- file.path(dir, "co_expression")
            } else {
                file <- "resultados_de_crossed"
            }

            file2 <- "results_completed_crossed"

            dir.create(file.path(dir, paste0(
                "Ontology_Results",
                tolower(group_gen)
            )), showWarnings = FALSE)
            dir <- file.path(dir, paste0(
                "Ontology_Results",
                tolower(group_gen)
            ))
        } else {
            file <- paste("resultados_de", tool, sep = "_")
            file2 <- paste("results_completed", tool, sep = "_")

            dir.create(paste0(
                path, "/Ontology_Results_",
                tolower(group_gen), "_",
                tool, "_", toupper(name)
            ), showWarnings = FALSE)
            dir <- paste0(
                path, "/Ontology_Results_", tolower(group_gen), "_",
                tool, "_", toupper(name)
            )
        }

        # generate annotation table
        annotation_table <- DOAGDC::annotation_table
        colnames(annotation_table) <- tolower(colnames(annotation_table))

        # input DE genes
        resultados_de <- get(file,
            envir = string_vars[["envir_link"]]
        )[[pair_name]]
        results_completed <- get(file2,
            envir = string_vars[["envir_link"]]
        )[[pair_name]]

        # pre-GO        gbutils::isNA
        if (tolower(id) == "geneid") {
            if (tolower(data_base) == "gdc") {
                # resultados_de
                resultados_de$ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(resultados_de)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    resultados_de$ensembl
                )
                tmp <- annotation_table$ensembl %in% selected
                geneid <- unique(annotation_table[tmp, c("ensembl", "geneid")])
                # Outer_join
                resultados_de <- merge(
                    x = resultados_de, y = geneid,
                    by = "ensembl", all = TRUE
                )

                # results_completed
                results_completed_ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(results_completed)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    results_completed_ensembl
                )
                geneid <- unique(annotation_table[tmp, c("ensembl", "geneid")])
                # Outer_join
                results_completed <- merge(
                    x = results_completed, y = geneid,
                    by = "ensembl", all = TRUE
                )
            }

            # remove any duplicates and NA
            resultados_de <- resultados_de[!is.na(resultados_de$geneid), ]
            resultados_de <- resultados_de[!duplicated(resultados_de$geneid), ]

            tmp <- results_completed$geneid
            results_completed <- results_completed[!is.na(tmp), ]
            results_completed <- results_completed[!duplicated(tmp), ]

            # upregulated DE
            de_genes_up <- resultados_de[resultados_de$fc > 0, ]
            assign("de_genes_up", de_genes_up, envir = get(envir_link))

            # downregulated DE
            de_genes_down <- resultados_de[resultados_de$fc < 0, ]
            assign("de_genes_down", de_genes_down, envir = get(envir_link))

            # Up e Down
            de_genes_all <- resultados_de
            assign("de_genes_all", de_genes_all, envir = get(envir_link))

            gene_vector_up <- as.integer(tmp %in% de_genes_up$geneid)
            names(gene_vector_up) <- tmp
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            gene_vector_down <- as.integer(tmp %in% de_genes_down$geneid)
            names(gene_vector_down) <- tmp
            assign("gene_vector_down", gene_vector_down,
                envir = get(envir_link)
            )

            gene_vector_all <- as.integer(tmp %in% de_genes_all$geneid)
            names(gene_vector_all) <- tmp
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (tolower(id) == "genesymbol") {
            if (tolower(data_base) == "gdc") {

                # resultados_de
                resultados_de$ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(resultados_de)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    resultados_de$ensembl
                )
                tmp <- annotation_table$ensembl %in% selected
                gene_symbol <- unique(annotation_table[tmp, c(
                    "ensembl",
                    "ucsc_genesymbol"
                )])
                # Outer_join
                resultados_de <- merge(
                    x = resultados_de, y = gene_symbol,
                    by = "ensembl", all = TRUE
                )

                # results_completed
                results_completed_ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(results_completed)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    results_completed_ensembl
                )
                gene_symbol <- unique(annotation_table[tmp, c(
                    "ensembl",
                    "ucsc_genesymbol"
                )])
                # Outer_join
                results_completed <- merge(
                    x = results_completed,
                    y = gene_symbol, by = "ensembl", all = TRUE
                )
            }

            # remove any duplicates and NA
            resultados_de <- resultados_de[!is.na(resultados_de$gene_symbol), ]
            tmp <- resultados_de$gene_symbol == "SLC35E2"
            resultados_de[tmp, "gene_symbol"][1] <- "SLC35E2B"
            resultados_de[
                duplicated(resultados_de$gene_symbol),
                "gene_symbol"
            ] <- paste0("?", seq_len(length(resultados_de[
                duplicated(resultados_de$gene_symbol),
                "gene_symbol"
            ])))

            tmp <- is.na(results_completed$gene_symbol)
            results_completed <- results_completed[!tmp, ]
            tmp <- results_completed$gene_symbol == "SLC35E2"
            results_completed[tmp, "gene_symbol"][1] <- "SLC35E2B"
            results_completed[
                duplicated(results_completed$gene_symbol),
                "gene_symbol"
            ] <- paste0(
                "?",
                seq(1, length(results_completed[
                    duplicated(results_completed$gene_symbol),
                    "gene_symbol"
                ]))
            )

            # upregulated DE
            de_genes_up <- resultados_de[resultados_de$fc > 0, ]
            assign("de_genes_up", de_genes_up, envir = get(envir_link))

            # downregulated DE
            de_genes_down <- resultados_de[resultados_de$fc < 0, ]
            assign("de_genes_down", de_genes_down, envir = get(envir_link))

            # Up e Down
            de_genes_all <- resultados_de
            assign("de_genes_all", de_genes_all, envir = get(envir_link))

            tmp <- results_completed$gene_symbol %in% de_genes_up$gene_symbol
            gene_vector_up <- as.integer(tmp)
            names(gene_vector_up) <- results_completed$gene_symbol
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            tmp <- results_completed$gene_symbol %in% de_genes_down$gene_symbol
            gene_vector_down <- as.integer(tmp)
            names(gene_vector_down) <- results_completed$gene_symbol
            assign("gene_vector_down", gene_vector_down,
                envir = get(envir_link)
            )

            tmp <- results_completed$gene_symbol %in% de_genes_all$gene_symbol
            gene_vector_all <- as.integer(tmp)
            names(gene_vector_all) <- results_completed$gene_symbol
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (toupper(id) == "HUGO") {
            if (tolower(data_base) == "gdc") {

                # resultados_de
                resultados_de$ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(resultados_de)
                )
                tmp <- annotation_table$ensembl %in% selected

                selected <- intersect(annotation_table$ensembl,
                                                        resultados_de$ensembl)
                hugo <- unique(annotation_table[tmp, c("ensembl", "hugo")])
                # Outer_join
                resultados_de <- merge(
                    x = resultados_de, y = hugo,
                    by = "ensembl", all = TRUE
                )

                # results_completed
                results_completed_ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(results_completed)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    results_completed_ensembl
                )
                hugo <- unique(annotation_table[tmp, c("ensembl", "hugo")])
                # Outer_join
                results_completed <- merge(
                    x = results_completed, y = hugo,
                    by = "ensembl", all = TRUE
                )
            }

            # remove any duplicates and NA
            resultados_de <- resultados_de[!is.na(resultados_de$geneid), ]
            resultados_de <- resultados_de[!duplicated(resultados_de$geneid), ]

            resultados_de$hugo <- as.character(annotation_table[match(
                resultados_de$geneid,
                annotation_table$geneid
            ), "HUGO"])

            tmp <- is.na(results_completed$geneid)
            results_completed <- results_completed[!tmp, ]

            tmp <- duplicated(results_completed$geneid)
            results_completed <- results_completed[!tmp, ]


            results_completed_hugo <- as.character(annotation_table[match(
                results_completed$geneid,
                annotation_table$geneid
            ), "HUGO"])

            # remove any NA
            resultados_de <- resultados_de[!is.na(resultados_de$hugo), ]

            tmp <- is.na(results_completed$hugo)
            results_completed <- results_completed[!tmp, ]



            # upregulated DE
            de_genes_up <- resultados_de[resultados_de$fc > 0, ]
            assign("de_genes_up", de_genes_up, envir = get(envir_link))

            # downregulated DE
            de_genes_down <- resultados_de[resultados_de$fc < 0, ]
            assign("de_genes_down", de_genes_down, envir = get(envir_link))

            # Up e Down
            de_genes_all <- resultados_de
            assign("de_genes_all", de_genes_all, envir = get(envir_link))

            tmp <- results_completed$hugo %in% de_genes_up$hugo
            gene_vector_up <- as.integer(tmp)
            names(gene_vector_up) <- resultados_de$geneid
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            tmp <- results_completed$hugo %in% de_genes_down$hugo
            gene_vector_down <- as.integer(tmp)
            names(gene_vector_down) <- resultados_de$geneid
            assign("gene_vector_down", gene_vector_down,
                                                    envir = get(envir_link))

            tmp <- results_completed$hugo %in% de_genes_all$hugo
            gene_vector_all <- as.integer(tmp)
            names(gene_vector_all) <- resultados_de$geneid
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (tolower(id) == "ensembl") {
            if (tolower(data_base) == "legacy") {
                # resultados_de
                resultados_de$geneid <- gsub(
                    pattern = "\\..*", "",
                    rownames(resultados_de)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    resultados_de$ensembl
                )
                tmp <- annotation_table$ensembl %in% selected
                ensembl <- unique(annotation_table[tmp, c("geneid", "ensembl")])
                # Outer_join
                resultados_de <- merge(
                    x = resultados_de, y = ensembl,
                    by = "geneid", all = TRUE
                )

                # results_completed
                results_completed_geneid <- gsub(
                    pattern = "\\..*", "",
                    rownames(results_completed)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    results_completed_ensembl
                )
                ensembl <- unique(annotation_table[tmp, c("geneid", "ensembl")])
                # Outer_join
                results_completed <- merge(
                    x = results_completed, y = ensembl,
                    by = "geneid", all = TRUE
                )


                # remove any duplicates and NA
                resultados_de <- resultados_de[!is.na(resultados_de$ensembl), ]
                tmp <- duplicated(resultados_de$ensembl)
                resultados_de <- resultados_de[!tmp, ]

                tmp <- is.na(results_completed$ensembl)
                results_completed <- results_completed[!tmp, ]

                tmp <- duplicated(results_completed$ensembl)
                results_completed <- results_completed[!tmp, ]

                # upregulated DE
                de_genes_up <- resultados_de[resultados_de$fc > 0, ]
                assign("de_genes_up", de_genes_up, envir = get(envir_link))

                # downregulated DE
                de_genes_down <- resultados_de[resultados_de$fc < 0, ]
                assign("de_genes_down", de_genes_down, envir = get(envir_link))

                # Up e Down
                de_genes_all <- resultados_de
                assign("de_genes_all", de_genes_all, envir = get(envir_link))


                tmp <- results_completed$ensembl %in% de_genes_up$ensembl
                gene_vector_up <- as.integer(tmp)
                names(gene_vector_up) <- resultados_de$ensembl

                tmp <- results_completed$ensembl %in% de_genes_down$ensembl
                gene_vector_down <- as.integer(tmp)
                names(gene_vector_down) <- resultados_de$ensembl

                tmp <- results_completed$ensembl %in% de_genes_all$ensembl
                gene_vector_all <- as.integer(tmp)
                names(gene_vector_all) <- results_completed_ensembl
            } else {

                # upregulated DE
                de_genes_up <- resultados_de[resultados_de$fc > 0, ]
                assign("de_genes_up", de_genes_up, envir = get(envir_link))

                # downregulated DE
                de_genes_down <- resultados_de[resultados_de$fc < 0, ]
                assign("de_genes_down", de_genes_down, envir = get(envir_link))

                # Up e Down
                de_genes_all <- resultados_de
                assign("de_genes_all", de_genes_all, envir = get(envir_link))


                # remove any duplicates and NA
                tmp <- rownames(results_completed) %in% rownames(de_genes_up)
                gene_vector_up <- as.integer(tmp)
                names(gene_vector_up) <- rownames(resultados_de)

                tmp <- rownames(results_completed) %in% rownames(de_genes_down)
                gene_vector_down <- as.integer(tmp)
                names(gene_vector_down) <- rownames(resultados_de)

                tmp <- rownames(results_completed) %in% rownames(de_genes_all)
                gene_vector_all <- as.integer(tmp)
                names(gene_vector_all) <- rownames(resultados_de)
            }

            # remove any duplicates and NA
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            assign("gene_vector_down", gene_vector_down,
                envir = get(envir_link)
            )

            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        } else if (tolower(id) == "refgene") {
            if (tolower(data_base) == "gdc") {
                # resultados_de
                resultados_de$ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(resultados_de)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    resultados_de$ensembl
                )
                tmp <- annotation_table$ensembl %in% selected
                ref_gene <- unique(annotation_table[, c("ensembl", "refseq")])
                # Outer_join
                resultados_de <- merge(
                    x = resultados_de, y = ref_gene,
                    by = "ensembl", all = TRUE
                )

                # results_completed
                results_completed_ensembl <- gsub(
                    pattern = "\\..*", "",
                    rownames(results_completed)
                )
                selected <- intersect(
                    annotation_table$ensembl,
                    results_completed_ensembl
                )
                ref_gene <- unique(annotation_table[tmp, c(
                    "ensembl",
                    "refseq"
                )])
                # Outer_join
                results_completed <- merge(
                    x = results_completed,
                    y = ref_gene, by = "ensembl", all = TRUE
                )
            }

            # remove any duplicates and NA
            resultados_de <- resultados_de[!is.na(resultados_de$geneid), ]
            resultados_de <- resultados_de[!duplicated(resultados_de$geneid), ]

            tmp <- is.na(results_completed$geneid)
            results_completed <- results_completed[!tmp, ]

            tmp <- duplicated(results_completed$geneid)
            results_completed <- results_completed[!tmp, ]

            resultados_de$ref_gene <- as.character(annotation_table[match(
                resultados_de$geneid,
                annotation_table$geneid
            ), "refseq"])

            results_completed_refgene <- as.character(annotation_table[match(
                results_completed$geneid,
                annotation_table$geneid
            ), "refseq"])


            # remove any duplicates and NA
            resultados_de <- resultados_de[!is.na(resultados_de$ref_gene), ]
            tmp <- duplicated(resultados_de$ref_gene)
            resultados_de <- resultados_de[!tmp, ]


            tmp <- is.na(results_completed$ref_gene)
            results_completed <- results_completed[!tmp, ]
            tmp <- duplicated(results_completed$ref_gene)
            results_completed <- results_completed[!tmp, ]

            # upregulated DE
            de_genes_up <- resultados_de[resultados_de$fc > 0, ]
            assign("de_genes_up", de_genes_up, envir = get(envir_link))

            # downregulated DE
            de_genes_down <- resultados_de[resultados_de$fc < 0, ]
            assign("de_genes_down", de_genes_down, envir = get(envir_link))

            # Up e Down
            de_genes_all <- resultados_de
            assign("de_genes_all", de_genes_all, envir = get(envir_link))

            tmp <- results_completed$ref_gene %in% de_genes_up$ref_gene
            gene_vector_up <- as.integer(tmp)
            names(gene_vector_up) <- results_completed_refgene
            assign("gene_vector_up", gene_vector_up, envir = get(envir_link))

            tmp <- results_completed$ref_gene %in% de_genes_down$ref_gene
            gene_vector_down <- as.integer(tmp)
            names(gene_vector_down) <- results_completed$ref_gene
            assign("gene_vector_down", gene_vector_down,
                envir = get(envir_link)
            )

            tmp <- results_completed$ref_gene %in% de_genes_all$ref_gene
            gene_vector_all <- as.integer(tmp)
            names(gene_vector_all) <- results_completed$ref_gene
            assign("gene_vector_all", gene_vector_all, envir = get(envir_link))
        }

        assign("pair_name", pair_name, envir = get(envir_link))
    }


    go_wall_fun <- function(x, n, kegg = FALSE) {
        if (nrow(x) != 0) {
            if (nrow(x) >= 15) {
                plot_now_go <- x[1:15, ]
            } else {
                plot_now_go <- x
            }

            plot_now_go[, 2] <- gsub(" of", "", plot_now_go[, 2])
            large <- lengths(strsplit(plot_now_go[, 2], " ")) > 7
            plot_now_go[large, 2] <- unname(sapply(
                plot_now_go[large, 2],
                function(w) {
                    paste(unlist(strsplit(w, " "))[1:7],
                        collapse = " "
                    )
                }
            ))

            log_10_go <- -log10(
                as.numeric(plot_now_go[, "over_represented_bh"])
            )
            new_inf <- log_10_go[order(log_10_go, decreasing = TRUE)]
            log_10_go[log_10_go == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

            longest_word <- max(stringr::str_count(plot_now_go$term))

            longest_width <- ifelse(longest_word > 50, width * 1.2, width)

            if (kegg) {
                plot_now_go$term <- plot_now_go$Pathway

                if (!is.numeric(plot_now_go$gene_ratio)) {
                    skip <- TRUE
                } else {
                    skip <- FALSE
                }
            } else {
                skip <- FALSE
            }

            if (!skip) {
                if (tolower(image_format) == "png") {
                    png(
                        filename = paste0(
                            dir,
                            "/GO_Output/GraphOutput/",
                            "GOEnrichPlot_smallest_FDR_",
                            main_names_short[n], "_",
                            tolower(condition), "_", id,
                            "_", pair_name, ".png"
                        ),
                        width = longest_width, height = height,
                        res = res, units = unit
                    )
                } else if (tolower(image_format) == "svg") {
                    svg(
                        filename = paste0(
                            dir,
                            "/GO_Output/GraphOutput/",
                            "GOEnrichPlot_smallest_FDR_",
                            main_names_short[n], "_",
                            tolower(condition), "_", id,
                            "_", pair_name, ".svg"
                        ),
                        width = longest_width, height = height, onefile = TRUE
                    )
                } else {
                    stop(message(
                        "Please, Insert a valid ",
                        "image_format! ('png' or 'svg')"
                    ))
                }

                count <- plot_now_go$numDEInCat
                gene_ratio <- round(as.numeric(plot_now_go$gene_ratio), 2)

                p <- ggplot2::ggplot(plot_now_go, ggplot2::aes(
                    x = log_10_go,
                    y = forcats::fct_reorder(term, log_10_go)
                )) +
                    ggplot2::geom_point(ggplot2::aes(
                        size = count,
                        color = gene_ratio
                    )) +
                    ggplot2::scale_colour_gradient(
                        limits = c(0, 1),
                        low = "red", high = "blue"
                    ) +
                    ggplot2::labs(
                        y = "", x = "-log(FDR)",
                        title = main_names[n]
                    ) +
                    ggplot2::theme_bw(base_size = 10) +
                    ggplot2::theme(
                        axis.title.x = ggplot2::element_text(
                            face = "bold",
                            size = 16
                        ),
                        axis.text = ggplot2::element_text(
                            face = "bold",
                            color = "#011600",
                            size = 12
                        ),
                        title = ggplot2::element_text(
                            face = "bold",
                            size = 18
                        ),
                        plot.title = ggplot2::element_text(hjust = 0.5)
                    )
                print(p)
                dev.off()
            }
        } else {
            message(paste0(
                "There are no terms with FDR < ",
                fdr_cutoff, " to plot."
            ))
        }
    }

    # code GO ####
    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, ",
            "please insert the 'env' name and try again!"
        ))
    }

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    path <- ifelse(
        exists("name_e", envir = get(envir_link)), file.path(
            string_vars[["envir_link"]]$path,
            string_vars[["envir_link"]]$name_e
        ), string_vars[["envir_link"]]$path
    )

    tool <- ifelse(missing(tool), string_vars[["envir_link"]]$tool, tool)

    name <- string_vars[["envir_link"]]$name
    data_base <- string_vars[["envir_link"]]$data_base
    group_gen <- string_vars[["envir_link"]]$group_gen

    message("Running GO preparations...")

    prepar_path_enrich(
        id = id,
        pair_name = pair_name,
        env = env,
        tool = tolower(tool)
    )

    message("Starting GO estimation...")
    if (grepl("crosstable", tolower(tool))) {
        if (tolower(tool) == "crosstable_deseq2") {
            dir <- paste0(path, "/CrossData_deseq2")
        } else if (tolower(tool) == "crosstable_edger") {
            dir <- paste0(path, "/CrossData_edger")
        } else if (tolower(tool) == "crosstable_ebseq") {
            dir <- paste0(path, "/CrossData_ebseq")
        }

        if (grepl("co_exp", tolower(tool))) {
            dir <- file.path(dir, "co_expression")
        }

        dir.create(file.path(dir, paste0(
            "Ontology_Results",
            tolower(group_gen)
        )), showWarnings = FALSE)
        dir <- file.path(dir, paste0("Ontology_Results", tolower(group_gen)))
    } else {
        dir.create(paste0(
            path, "/Ontology_Results_", tolower(group_gen), "_",
            tolower(tool), "_", toupper(name)
        ), showWarnings = FALSE)
        dir <- paste0(
            path, "/Ontology_Results_", tolower(group_gen), "_",
            tolower(tool), "_", toupper(name)
        )
    }

    dir.create(file.path(dir, "GO_Output"), showWarnings = FALSE)
    dir.create(file.path(dir, "GO_Output", "GraphOutput"),
        showWarnings = FALSE
    )
    dir.create(file.path(dir, "GO_Output", "TextOutput"),
        showWarnings = FALSE
    )

    if (tolower(id) == "geneid") {
        formatted_id <- "knownGene"
    } else if (tolower(id) == "genesymbol") {
        formatted_id <- "geneSymbol"
    } else if (tolower(id) == "ensembl") {
        formatted_id <- "ensGene"
    } else if (tolower(id) == "refgene") {
        formatted_id <- "refGene"
    }

    if (tolower(image_format) == "png") {
        png(
            filename = paste0(
                dir,
                "/GO_Output/GraphOutput/ProbabilityWeightingFunction_",
                tolower(condition), "_", id, "_", pair_name, ".png"
            ),
            width = width, height = height, res = res, units = unit
        )
    } else if (tolower(image_format) == "svg") {
        svg(
            filename = paste0(
                dir,
                "/GO_Output/GraphOutput/ProbabilityWeightingFunction_",
                tolower(condition), "_", id, "_", pair_name, ".svg"
            ),
            width = width, height = height, onefile = TRUE
        )
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }

    if (tolower(data_base) == "legacy") {
        genome_version <- "hg19"
    } else {
        genome_version <- "hg38"
    }

    # Probability Weighting Function
    if (tolower(condition) == "upregulated") {
        gene_vector_up <- string_vars[["envir_link"]]$gene_vector_up
        suppressPackageStartupMessages(pwf <- goseq::nullp(
            gene_vector_up,
            genome_version,
            formatted_id
        ))
        dev.off()
    } else if (tolower(condition) == "downregulated") {
        gene_vector_down <- string_vars[["envir_link"]]$gene_vector_down
        suppressPackageStartupMessages(pwf <- goseq::nullp(
            gene_vector_down,
            genome_version,
            formatted_id
        ))
        dev.off()
    } else if (tolower(condition) == "all") {
        gene_vector_all <- string_vars[["envir_link"]]$gene_vector_all
        suppressPackageStartupMessages(pwf <- goseq::nullp(
            gene_vector_all,
            genome_version,
            formatted_id
        ))
        dev.off()
    }

    # GO Analysis Wallenius
    message("Running GO analysis Wallenius\n")
    suppressPackageStartupMessages(go_wall <- goseq::goseq(pwf,
        genome_version, formatted_id,
        test.cats = c("GO:CC", "GO:BP", "GO:MF"),
        method = "Wallenius",
        use_genes_without_cat = use_genes_without_cat
    ))

    # FDR Add to go_wall
    go_wall$over_represented_bh <- p.adjust(go_wall$over_represented_pvalue,
        method = "BH"
    )

    go_wall$gene_ratio <- go_wall$numDEInCat / go_wall$numInCat

    go_wall <- go_wall[go_wall$over_represented_bh < fdr_cutoff, ]

    # Escreve as tabelas
    go_wall_bh_cc <- go_wall[go_wall[, "ontology"] == "CC", ]
    go_wall_bh_bp <- go_wall[go_wall[, "ontology"] == "BP", ]
    go_wall_bh_mf <- go_wall[go_wall[, "ontology"] == "MF", ]

    # Define category order
    main_names <- c(
        "Cellular component", "Biological Process",
        "Molecular Funcion", "KEGG"
    )
    main_names_short <- c("CC", "BP", "MF", "KEGG")

    go_wall_fun(x = go_wall_bh_cc, n = 1)
    go_wall_fun(x = go_wall_bh_bp, n = 2)
    go_wall_fun(x = go_wall_bh_mf, n = 3)

    suppressPackageStartupMessages(go_wall <- goseq::goseq(pwf, genome_version,
        formatted_id,
        test.cats = c("KEGG"),
        method = "Wallenius",
        use_genes_without_cat = use_genes_without_cat
    ))

    go_wall$over_represented_bh <- p.adjust(go_wall$over_represented_pvalue,
        method = "BH"
    )
    go_wall$category <- paste0("hsa", go_wall$category)
    colnames(go_wall)[1] <- "Pathway"
    go_wall$gene_ratio <- go_wall$numDEInCat / go_wall$numInCat
    go_wall_bh_kegg <- go_wall[go_wall[, "over_represented_bh"] < fdr_cutoff, ]

    message("Writting tables...\n")
    write.csv(go_wall_bh_cc, file = paste0(
        dir,
        "/GO_Output/TextOutput/CC_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))
    write.csv(go_wall_bh_bp, file = paste0(
        dir,
        "/GO_Output/TextOutput/BP_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))
    write.csv(go_wall_bh_mf, file = paste0(
        dir,
        "/GO_Output/TextOutput/MF_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))
    write.csv(go_wall_bh_kegg, file = paste0(
        dir,
        "/GO_Output/TextOutput/KEGG_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))

    go_wall_fun(x = go_wall_bh_kegg, n = 4, kegg = TRUE)
    #################### GO ENRICHMENT
    #################### END

    # code enrichGO ####
    dir.create(file.path(dir, "enrichGO_Output"), showWarnings = FALSE)
    dir.create(file.path(dir, "enrichGO_Output", "GraphOutput"),
        showWarnings = FALSE
    )
    dir.create(file.path(dir, "enrichGO_Output", "TextOutput"),
        showWarnings = FALSE
    )

    if (tolower(id) == "geneid") {
        formatted_id2 <- "ENTREZID"
    } else if (tolower(id) == "genesymbol") {
        formatted_id2 <- "SYMBOL"
    } else if (tolower(id) == "ensembl") {
        formatted_id2 <- "ENSEMBL"
    } else if (tolower(id) == "refgene") {
        formatted_id2 <- "REFSEQ"
    }

    erichgo <- function(ont, condition, n) {
        if (tolower(condition) == "upregulated") {
            gene_vector <- string_vars[["envir_link"]]$gene_vector_up
        } else if (tolower(condition) == "downregulated") {
            gene_vector <- string_vars[["envir_link"]]$gene_vector_down
        } else if (tolower(condition) == "all") {
            gene_vector <- string_vars[["envir_link"]]$gene_vector_all
        }

        enrichgo_obj <- clusterProfiler::enrichGO(
            gene = names(gene_vector)[gene_vector != 0],
            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            keyType = formatted_id2,
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = fdr_cutoff
        )

        enrichgo_resuts <- enrichgo_obj@result

        if (nrow(enrichgo_resuts) > 0) {
            tmp <- order(enrichgo_resuts$p.adjust)
            enrichgo_resuts <- enrichgo_resuts[tmp, ]

            ratios_splited <- unname(sapply(
                enrichgo_resuts$GeneRatio,
                function(w) {
                    unlist(strsplit(w, "\\/"))
                }
            ))

            tmp <- as.numeric(ratios_splited[1, ])
            num_ratios <- tmp / as.numeric(ratios_splited[2, ])


            enrichgo_resuts$gene_ratio_num <- num_ratios


            if (nrow(enrichgo_resuts) >= 15) {
                plot_now_go <- enrichgo_resuts[1:15, ]
            } else {
                plot_now_go <- enrichgo_resuts
            }

            plot_now_go[, 2] <- gsub(" of", "", plot_now_go[, 2])
            plot_now_go[, 2] <- gsub("-", " ", plot_now_go[, 2])
            large <- lengths(strsplit(plot_now_go[, 2], " ")) > 3
            plot_now_go[large, 2] <- unname(sapply(
                plot_now_go[large, 2],
                function(w) {
                    paste(unlist(strsplit(w, " "))[1:4],
                        collapse = " "
                    )
                }
            ))

            log_10_go <- -log10(as.numeric(plot_now_go[, "p.adjust"]))
            new_inf <- log_10_go[order(log_10_go, decreasing = TRUE)]
            log_10_go[log_10_go == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

            longest_word <- max(stringr::str_count(plot_now_go$Description))

            if (longest_word > 40) {
                longest_width <- width * 1.2
            } else {
                longest_width <- width
            }

            if (nrow(enrichgo_resuts) > 0) {
                if (tolower(image_format) == "png") {
                    png(
                        filename = paste0(
                            dir,
                            "/enrichGO_Output/GraphOutput/",
                            "enrichGO_smallest_FDR_Ont_", ont, "_",
                            tolower(condition), "_", id,
                            "_", pair_name, ".png"
                        ),
                        width = longest_width, height = height, res = res,
                        units = unit
                    )
                } else if (tolower(image_format) == "svg") {
                    svg(
                        filename = paste0(
                            dir,
                            "/enrichGO_Output/GraphOutput/",
                            "enrichGO_smallest_FDR_", ont, "_",
                            tolower(condition), "_", id,
                            "_", pair_name, ".svg"
                        ),
                        width = longest_width, height = height, onefile = TRUE
                    )
                } else {
                    stop(message(
                        "Please, Insert a valid image_format!",
                        " ('png' or 'svg')"
                    ))
                }

                count <- plot_now_go$numDEInCat

                gene_ratio <- round(as.numeric(plot_now_go$gene_ratio_num), 2)

                ## plot alternative to barplot
                p <- ggplot2::ggplot(plot_now_go, ggplot2::aes(
                    x = log_10_go,
                    y = forcats::fct_reorder(Description, log_10_go)
                )) +
                    ggplot2::geom_point(ggplot2::aes(
                        size = count,
                        color = gene_ratio
                    )) +
                    ggplot2::scale_colour_gradient(
                        limits = c(0, 1),
                        low = "red", high = "blue"
                    ) +
                    ggplot2::labs(
                        y = "", x = "-log(FDR)",
                        title = main_names[n]
                    ) +
                    ggplot2::theme_bw(base_size = 10) +
                    ggplot2::theme(
                        axis.title.x = ggplot2::element_text(
                            face = "bold",
                            size = 16
                        ),
                        axis.text = ggplot2::element_text(
                            face = "bold",
                            color = "#011600",
                            size = 12
                        ),
                        title = ggplot2::element_text(
                            face = "bold",
                            size = 18
                        ),
                        plot.title = ggplot2::element_text(hjust = 0.5)
                    )
                print(p)

                dev.off()
            }
        }

        return(enrichgo_resuts)
    }

    erichgo_cc <- erichgo(ont = "CC", condition, n = 1)
    erichgo_bp <- erichgo(ont = "BP", condition, n = 2)
    erichgo_mf <- erichgo(ont = "MF", condition, n = 3)

    message("Writting tables...\n")
    write.csv(erichgo_cc, file = paste0(
        dir,
        "/enrichGO_Output/TextOutput/CC_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))
    write.csv(erichgo_bp, file = paste0(
        dir,
        "/enrichGO_Output/TextOutput/BP_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))
    write.csv(erichgo_mf, file = paste0(
        dir,
        "/enrichGO_Output/TextOutput/MF_",
        tolower(condition), "_", id,
        "_FDR_", fdr_cutoff,
        "_", pair_name, ".csv"
    ))

    message("Done!\n")
}
