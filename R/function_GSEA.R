#' Gene Set Enrichment Analysis
#'
#' @param fdr_cutoff
#' @param width,height,res,unit,image_format
#' @param tool
#' @param id
#' @param pair_name
#' @param env
#' @inheritParams dea_ebseq
#' @inheritParams groups_identification_mclust
#' @inheritParams gonto
#' @inheritParams concatenate_exon
#'
#' @return Enriched terms.
#' @export
#'
#' @importFrom AnnotationDbi keytypes
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom clusterProfiler bitr
#' @importFrom ReactomePA gsePathway
#' @importFrom stringr str_count
#' @importFrom forcats fct_reorder
#'
#' @examples
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
#' gsea(tool = "edgeR", env = CHOL_LEGACY_gene_tumor_data)
gsea <- function(fdr_cutoff = 0.05,
                    width = 10,
                    height = 3,
                    res = 500,
                    unit = "in",
                    image_format = "png",
                    tool = "edgeR",
                    id = "geneid",
                    pair_name = "G2_over_G1",
                    env) {
    message("Performing Gene Set Enrichment Analysis - over Reactome...\n")
    # Looks for enrichment by fold direction
    # List of input id types

    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please",
            " insert the 'env' name and try again!"
        ))
    }

    swp <- function(x) {
        suppressPackageStartupMessages(x)
    }

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    path <- ifelse(exists("name_e", envir = get(envir_link)),
        file.path(string_vars[["envir_link"]]$path,
        string_vars[["envir_link"]]$name_e),
        string_vars[["envir_link"]]$path
    )

    tool <- ifelse(missing(tool), tolower(string_vars[["envir_link"]]$tool),
                                                                tolower(tool))

    name <- string_vars[["envir_link"]]$name
    group_gen <- string_vars[["envir_link"]]$group_gen

    tcga_expression <- ifelse(grepl("crosstable", tool),
                        string_vars[["envir_link"]]$results_completed_crossed,
                        eval(parse(text = paste0(
                                "string_vars",
                                "[['envir_link']]$results_completed_",
                                tool
                            )))
                        )

    tcga_expression <- tcga_expression[[pair_name]]

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
            path, "/Ontology_Results_", tolower(group_gen), "_",
            tool, "_", toupper(name)
        )
    }

    dir.create(file.path(dir, "GSEA_Output"), showWarnings = FALSE)

    swp(input_types <- as.matrix(
        x = AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
    ))
    write.table(input_types, paste0(
        dir,
        "/GSEA_Output/input_types",
        id, "_", pair_name, ".txt"
    ), row.names = FALSE)

    gene_set_vector <- tcga_expression$fc

    if (tolower(id) == "geneid") {
        # It demands FOld change in one vector decreasing ordered
        names(gene_set_vector) <- tcga_expression$geneid
        gene_set_vector <- gene_set_vector[order(gene_set_vector,
            decreasing = TRUE
        )]
    } else if (tolower(id) == "genesymbol") {
        # It demands FOld change in one vector decreasing ordered
        names(gene_set_vector) <- tcga_expression$gene_symbol
        gene_set_vector <- gene_set_vector[order(gene_set_vector,
            decreasing = TRUE
        )]

        # Again, change the names for id using HUGO names
        ids <- clusterProfiler::bitr(names(gene_set_vector),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db"
        )
        gene_set_vector <- gene_set_vector[ids$ALIAS]
        names(gene_set_vector) <- ids$ENTREZID
    } else if (tolower(id) == "ensembl") {
        names(gene_set_vector) <- rownames(tcga_expression)
        gene_set_vector <- gene_set_vector[order(gene_set_vector,
            decreasing = TRUE
        )]

        # Again, change the names for id using HUGO names
        ids <- clusterProfiler::bitr(names(gene_set_vector),
            fromType = "ENSEMBL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db"
        )
        gene_set_vector <- gene_set_vector[ids$ENSEMBL]
        names(gene_set_vector) <- ids$ENTREZID
    } else if (tolower(id) == "refgene") {
        names(gene_set_vector) <- tcga_expression$Ensembl
        gene_set_vector <- gene_set_vector[order(gene_set_vector,
            decreasing = TRUE
        )]

        # Again, change the names for id using HUGO names
        ids <- clusterProfiler::bitr(names(gene_set_vector),
            fromType = "REFSEQ",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db"
        )
        gene_set_vector <- gene_set_vector[ids$ENSEMBL]
        names(gene_set_vector) <- ids$ENTREZID
    }

    # Perform GSEA
    gsea_reactome <- ReactomePA::gsePathway(gene_set_vector,
        pvalueCutoff = fdr_cutoff, nPerm = 10000,
        pAdjustMethod = "BH", verbose = TRUE
    )

    # Get summary
    gsea_reactome_summary_fdr <- as.data.frame(gsea_reactome)

    # Write
    write.csv(gsea_reactome_summary_fdr,
        file = paste0(
            dir,
            "/GSEA_Output/",
            "REAC.GSEA.enrichment",
            id, "_", pair_name, ".csv"
        ),
        row.names = FALSE
    )
    assign("gene_set_vector", gene_set_vector, envir = get(envir_link))


    # REACTOME-GSEA Plot
    if (nrow(gsea_reactome_summary_fdr) != 0) {
        if (nrow(gsea_reactome_summary_fdr) >= 15) {
            plot_now_react <- gsea_reactome_summary_fdr[1:15, ]
        } else {
            plot_now_react <- gsea_reactome_summary_fdr
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

        # is it better to use the manually adjusted?
        log_10_gsea <- -log10(as.numeric(plot_now_react[, "p.adjust"]))
        new_inf <- log_10_gsea[order(log_10_gsea, decreasing = TRUE)]
        log_10_gsea[log_10_gsea == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        longest_word <- max(stringr::str_count(plot_now_react$Description))


        gene_ratio <- unname(sapply(
            plot_now_react$leading_edge,
            function(w) {
                unlist(strsplit(w, ", "))[1]
            }
        ))

        gene_ratio <- round(
            as.numeric(gsub("tags=|%", "", gene_ratio)) / 100, 2
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
                    "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", id,
                    "_", pair_name, ".png"
                ),
                width = longest_width, height = height, res = res,
                units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(
                    dir,
                    "/GSEA_Output/REACT-GSEA_EnrichPlot_10first", id,
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

        count <- ceiling(plot_now_react$setSize * gene_ratio)


        p <- ggplot2::ggplot(plot_now_react, ggplot2::aes(
            x = log_10_gsea,
            y = forcats::fct_reorder(Description, log_10_gsea)
        )) +
            ggplot2::geom_point(ggplot2::aes(
                size = count,
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
        message(
            "There is nothing to show in Gene Set Enrichment",
            " Analysis - over Reactome"
        )
    }
    gc()
    message("Done!\n")
}
