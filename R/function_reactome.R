#' DESEASE-ONTOLOGY and REACTOME ENRICHMENT
#'
#' @param p_cutoff
#' @param fdr_cutoff
#' @param width,height,res,unit,image_format
#' @param tool
#' @param pair_name
#' @param env
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_ebseq
#' @inheritParams gonto
#'
#' @return Enriched terms.
#' @export
#'
#' @importFrom DOSE enrichDO
#' @importFrom clusterProfiler bitr
#' @importFrom stringr str_count
#' @importFrom ReactomePA enrichPathway
#'
#' @examples
#' library(DOAGDC)
#' # data already downloaded using the 'download_gdc' function
#' concatenate_expression("gene",
#'    name = "HIF3A",
#'    data_base = "legacy",
#'    tumor = "CHOL",
#'    work_dir = "~/Desktop"
#' )
#'
#' # separating gene HIF3A expression data patients in two groups
#' groups_identification_mclust("gene", 2,
#'    name = "HIF3A",
#'    modelName = "E",
#'    env = CHOL_LEGACY_gene_tumor_data,
#'    tumor = "CHOL"
#' )
#'
#' # load not normalized data
#' concatenate_expression("gene",
#'    normalization = FALSE,
#'    name = "HIF3A",
#'    data_base = "legacy",
#'    tumor = "CHOL",
#'    env = CHOL_LEGACY_gene_tumor_data,
#'    work_dir = "~/Desktop"
#' )
#'
#' # start DE analysis
#' # considering concatenate_expression and groups_identification already runned
#' dea_edger(
#'    name = "HIF3A",
#'    group_gen = "mclust",
#'    env = CHOL_LEGACY_gene_tumor_data
#' )
#'
#' gonto(
#'    condition = "Upregulated",
#'    tool = "edgeR", env = CHOL_LEGACY_gene_tumor_data
#' )
#'
#' do_react_enrich(tool = "edgeR", env = CHOL_LEGACY_gene_tumor_data)
do_react_enrich <- function(p_cutoff = 0.05,
                            fdr_cutoff = 0.05,
                            width = 8,
                            height = 4,
                            res = 300,
                            unit = "in",
                            image_format = "png",
                            tool = "edgeR",
                            pair_name = "G2_over_G1",
                            env) {
    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please insert the",
            " 'env' name and try again!"
        ))
    }

    ev <- deparse(substitute(env))
    sv <- list(ev = get(ev))

    path <- ifelse(exists("name_e", envir = get(ev)),
                file.path(
                    sv[["ev"]]$path,
                    sv[["ev"]]$name_e
                ), sv[["ev"]]$path
            )

    tool <- ifelse(missing(tool), tolower(sv[["ev"]]$tool),
        tolower(tool)
    )

    name <- sv[["ev"]]$name
    data_base <- sv[["ev"]]$data_base
    group_gen <- sv[["ev"]]$group_gen

    if (grepl("crosstable", tool)) {
        tcga_expression <- sv[["ev"]]$results_completed_crossed
    } else {
        tcga_expression <- eval(parse(text = paste0(
            "sv[[",
            "'ev']]$results_completed_",
            tool
        )))
    }
    tcga_expression <- tcga_expression[[pair_name]]

    if (grepl("crosstable", tool)) {
        de_genes_all <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
    } else {
        de_genes_all <- get(paste("resultados_de", tool, sep = "_"),
            envir = sv[["ev"]]
        )[[pair_name]]
    }

    if (grepl("crosstable", tool)) {
        if (tool == "crosstable.deseq2") {
            dir <- paste0(path, "/CrossData_deseq2")
        } else if (tool == "crosstable.edger") {
            dir <- paste0(path, "/CrossData_edger")
        } else if (tool == "crosstable.ebseq") {
            dir <- paste0(path, "/CrossData_ebseq")
        }
        dir.create(file.path(dir, paste0(
            "Ontology_Results",
            tolower(group_gen)
        )), showWarnings = FALSE)
        dir <- file.path(dir, paste0("Ontology_Results", tolower(group_gen)))
    } else {
        dir.create(paste0(
            path, "/Ontology_Results_", tolower(group_gen), "_",
            tool, "_", toupper(name)
        ), showWarnings = FALSE)
        dir <- paste0(
            path, "/Ontology_Results_", tolower(group_gen), "_", tool,
            "_", toupper(name)
        )
    }

    dir.create(file.path(dir, "REACTOME_Output"), showWarnings = FALSE)
    dir.create(file.path(dir, "DO_Output"), showWarnings = FALSE)

    # "geneid" equal to "entrez gene id"
    if (tolower(data_base) == "gdc") {
        geneid_list_de <- clusterProfiler::bitr(de_genes_all$ensembl,
            fromType = "ENSEMBL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db"
        )$ENTREZID
        geneid_list_universe <- clusterProfiler::bitr(tcga_expression$ensembl,
            fromType = "ENSEMBL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db"
        )$ENTREZID
    } else {
        geneid_list_de <- de_genes_all$geneid
        geneid_list_universe <- tcga_expression$geneid
    }

    # DO ####

    # Perform DO Enrichment analysis
    message("\nPerforming DO enrichment...\n")
    do_enrichment <- DOSE::enrichDO(
        gene = geneid_list_de,
        ont = "DO",
        pvalueCutoff = fdr_cutoff,
        pAdjustMethod = "BH",
        universe = geneid_list_universe,
        minGSSize = 2, # at least a pair
        qvalueCutoff = 1,
        readable = TRUE
    )

    # Makes summary
    do_enriched_summary <- as.data.frame(do_enrichment)

    # Write
    write.csv(do_enriched_summary,
        file = paste0(
            dir,
            "/DO_Output/DO_enrichment_", pair_name, ".csv"
        ),
        row.names = FALSE
    )

    # Get fdr <0.05
    tmp <- do_enriched_summary$p.adjust < fdr_cutoff
    do_enriched_fdr <- do_enriched_summary[which(tmp), ]

    # DO Plot
    message("Plotting DO enrichment...\n")

    if (nrow(do_enriched_fdr) != 0) {
        if (nrow(do_enriched_fdr) >= 10) {
            plot_now_do <- do_enriched_fdr[1:10, ]
        } else {
            plot_now_do <- do_enriched_fdr
        }

        plot_now_do[, 2] <- gsub(" of", "", plot_now_do[, 2])
        large <- lengths(strsplit(plot_now_do[, 2], " ")) > 7
        plot_now_do[large, 2] <- unname(sapply(
            plot_now_do[large, 2],
            function(w) {
                paste(unlist(strsplit(w, " "))[1:7],
                    collapse = " "
                )
            }
        ))

        log_10_do <- -log10(as.numeric(plot_now_do[, "p.adjust"]))
        new_inf <- log_10_do[order(log_10_do, decreasing = TRUE)]
        log_10_do[log_10_do == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        longest_word <- max(stringr::str_count(plot_now_do$Description))

        gene_ratio1 <- unname(sapply(
            plot_now_do$GeneRatio,
            function(w) {
                unlist(strsplit(w, "/"))[1]
            }
        ))
        gene_ratio2 <- unname(sapply(
            plot_now_do$GeneRatio,
            function(w) {
                unlist(strsplit(w, "/"))[2]
            }
        ))

        gene_ratio <- round(
            as.numeric(gene_ratio1) / as.numeric(gene_ratio2), 2
        )

        if (longest_word > 50) {
            longest_width <- width * 1.5
        } else {
            longest_width <- width
        }

        if (tolower(image_format) == "png") {
            png(
                filename = paste0(
                    dir,
                    "/DO_Output/DO_EnrichPlot_10first_", pair_name,
                    ".png"
                ),
                width = longest_width, height = height, res = res,
                units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(
                    dir,
                    "/DO_Output/DO_EnrichPlot_10first_", pair_name,
                    ".svg"
                ),
                width = longest_width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format!",
                " ('png' or 'svg')"
            ))
        }

        p <- ggplot2::ggplot(plot_now_do, ggplot2::aes(
            x = log_10_do,
            y = forcats::fct_reorder(Description, log_10_do)
        )) +
            ggplot2::geom_point(ggplot2::aes(
                size = Count,
                color = gene_ratio
            )) +
            ggplot2::scale_colour_gradient(
                limits = c(0, 1), low = "red",
                high = "blue"
            ) +
            ggplot2::labs(y = "", x = "-log(FDR)") +
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
    } else {
        message(cat("There is nothing to show in Disease Ontology\n"))
    }

    # REACTOME ####
    # perform enrichment in REACTOME
    message("Performing reactome enrichment...\n")
    reactome_enriched <- ReactomePA::enrichPathway(geneid_list_de,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        universe = geneid_list_universe,
        minGSSize = 2,
        readable = TRUE
    )
    # Makes summary
    reactome_enriched_summary <- as.data.frame(reactome_enriched)

    # Write
    write.csv(reactome_enriched_summary,
        file = paste0(
            dir,
            "/REACTOME_Output/REAC_enrichment_", pair_name, ".csv"
        ),
        row.names = FALSE
    )

    # Get fdr <0.05
    tmp <- reactome_enriched_summary$p.adjust < fdr_cutoff
    reactome_enriched_fdr <- reactome_enriched_summary[which(tmp), ]

    tmp <- order(reactome_enriched_fdr$p.adjust, decreasing = FALSE)
    reactome_enriched_fdr <- reactome_enriched_fdr[tmp, ]

    # REACTOME Plot
    message("Plotting Reactome enrichment\n")
    if (nrow(reactome_enriched_fdr) != 0) {
        if (nrow(reactome_enriched_fdr) >= 10) {
            plot_now_react <- reactome_enriched_fdr[1:10, ]
        } else {
            plot_now_react <- reactome_enriched_fdr
        }

        plot_now_react[, 2] <- gsub(" of", "", plot_now_react[, 2])
        large <- lengths(strsplit(plot_now_react[, 2], " ")) > 7
        plot_now_react[large, 2] <- unname(sapply(
            plot_now_react[large, 2],
            function(w) {
                paste(unlist(strsplit(w, " "))[1:7],
                    collapse = " "
                )
            }
        ))

        log_10_react <- -log10(as.numeric(plot_now_react[, "p.adjust"]))
        new_inf <- log_10_react[order(log_10_react, decreasing = TRUE)]
        log_10_react[log_10_react == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

        longest_word <- max(stringr::str_count(plot_now_react$Description))

        gene_ratio1 <- unname(sapply(
            plot_now_react$GeneRatio,
            function(w) {
                unlist(strsplit(w, "/"))[1]
            }
        ))
        gene_ratio2 <- unname(sapply(
            plot_now_react$GeneRatio,
            function(w) {
                unlist(strsplit(w, "/"))[2]
            }
        ))

        gene_ratio <- round(
            as.numeric(gene_ratio1) / as.numeric(gene_ratio2), 2
        )

        if (longest_word > 50) {
            longest_width <- width * 1.5
        } else {
            longest_width <- width
        }

        if (tolower(image_format) == "png") {
            png(
                filename = paste0(
                    dir,
                    "/REACTOME_Output/REACT_EnrichPlot_10first_",
                    pair_name, ".png"
                ),
                width = longest_width, height = height, res = res,
                units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(
                    dir,
                    "/REACTOME_Output/REACT_EnrichPlot_10first_",
                    pair_name, ".svg"
                ),
                width = longest_width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format!",
                " ('png' or 'svg')"
            ))
        }

        p <- ggplot2::ggplot(plot_now_react, ggplot2::aes(
            x = log_10_react,
            y = forcats::fct_reorder(Description, log_10_react)
        )) +
            ggplot2::geom_point(ggplot2::aes(
                size = Count,
                color = gene_ratio
            )) +
            ggplot2::scale_colour_gradient(
                limits = c(0, 1), low = "red",
                high = "blue"
            ) +
            ggplot2::labs(y = "", x = "-log(FDR)") +
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
    } else {
        message(cat("There is nothing to show from Reactome Enrichment"))
    }

    # Write
    write.csv(reactome_enriched_summary,
        file = paste0(
            dir,
            "/REACTOME_Output/REAC_enrichment_under_cutoff_",
            pair_name, ".csv"
        ),
        row.names = FALSE
    )

    gc()
    message("Done!\n")
}
