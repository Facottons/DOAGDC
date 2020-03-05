#' Generate plot for Principal Component Analysis
#'
#' @param tool
#' @param name
#' @param work_dir
#' @param pair_name
#' @param width,height,res,unit,image_format
#' @param env
#' @inheritParams gonto
#' @inheritParams concatenate_exon
#' @inheritParams download_gdc
#' @inheritParams dea_ebseq
#' @inheritParams groups_identification_mclust
#'
#' @return the PCAs plots.
#' @export
#'
#' @importFrom ggbiplot ggbiplot
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' library(DOAGDC)
#'
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
#' pca_analysis("edgeR", "HIF3A", "~/Desktop",
#'    pair_name = "G2_over_G1",
#'    env = CHOL_LEGACY_gene_tumor_data
#' )
pca_analysis <- function(tool,
                        name, work_dir,
                        pair_name = "G2_over_G1",
                        width = 4,
                        height = 2,
                        res = 300,
                        unit = "in",
                        image_format = "png",
                        env) {

    # local function ####
    pca_local <- function(size, file) {
        if (tolower(size) == "all") {
            genetop <- size
            para_heatmaps <- normalized_expression[match(
                rownames(resultados_de),
                rownames(normalized_expression)
            ), ]

            ir_pca <- prcomp(t(log2(para_heatmaps + 1)),
                center = TRUE,
                scale. = TRUE
            )

            if (tolower(image_format) == "png") {
                png(
                    filename = paste0(dir, file, ".png"),
                    width = width, height = height, res = res, units = unit
                )
            } else if (tolower(image_format) == "svg") {
                svg(
                    filename = paste0(dir, file, ".svg"),
                    width = width, height = height, onefile = TRUE
                )
            } else {
                stop(message(
                    "Please, Insert a valid image_format!",
                    " ('png' or 'svg')"
                ))
            }
            a <- ggbiplot::ggbiplot(ir_pca,
                choices = 1:2,
                obs.scale = 1,
                var.scale = 1,
                groups = cond_heatmap,
                ellipse = TRUE,
                circle = TRUE,
                alpha = 0.5,
                var.axes = FALSE
            ) +
                scale_color_manual(
                    name = "Groups",
                    values = color_pallete[unique(cond_heatmap)]
                ) +
                theme(axis.title.y = element_text(
                    size = rel(0.6),
                    angle = 90
                )) +
                theme(axis.title.x = element_text(
                    size = rel(0.6),
                    angle = 00
                )) +
                theme(legend.title = element_text(size = 4, face = "bold")) +
                theme(legend.text = element_text(size = 4, face = "bold"))
            print(a)
            dev.off()
        } else if (is.numeric(size)) {
            genetop <- size
            tmp <- order(resultados_de$fdr, -abs(resultados_de$fc))
            resultados_plot <- resultados_de[tmp, ]
            resultados_plot <- resultados_plot[1:genetop, ]

            # Separacao para os heatmaps
            tmp <- match(rownames(resultados_plot),
                                            rownames(normalized_expression))
            para_heatmaps <- normalized_expression[tmp, ]
            para_heatmaps <- na.exclude(para_heatmaps)

            ir_pca <- prcomp(t(log2(para_heatmaps + 1)),
                center = TRUE,
                scale. = TRUE
            )

            if (tolower(image_format) == "png") {
                png(
                    filename = paste0(dir, file, "2.png"),
                    width = width, height = height, res = res, units = unit
                )
            } else if (tolower(image_format) == "svg") {
                svg(
                    filename = paste0(dir, file, "2.svg"),
                    width = width, height = height, onefile = TRUE
                )
            } else {
                stop(message(
                    "Please, Insert a valid image_format!",
                    " ('png' or 'svg')"
                ))
            }
            a <- ggbiplot::ggbiplot(ir_pca,
                choices = 1:2,
                obs.scale = 1,
                var.scale = 1,
                groups = cond_heatmap,
                ellipse = TRUE,
                circle = TRUE,
                alpha = 0.5,
                var.axes = FALSE
            ) +
                scale_color_manual(
                    name = "Groups",
                    values = color_pallete[unique(cond_heatmap)]
                ) +
                theme(axis.title.y = element_text(
                    size = rel(0.6),
                    angle = 90
                )) +
                theme(axis.title.x = element_text(
                    size = rel(0.6),
                    angle = 00
                )) +
                theme(legend.title = element_text(size = 4, face = "bold")) +
                theme(legend.text = element_text(size = 4, face = "bold"))
            print(a)
            dev.off()
        }
    }

    # code ####
    name <- gsub("-", "_", name)

    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please insert the ",
            "'env' name and try again!"
        ))
    }

    ev <- deparse(substitute(env))
    sv <- list(ev = get(ev))

    work_dir <- ifelse(missing("work_dir"), sv[["ev"]]$work_dir, work_dir)

    path <- ifelse(exists("name_e", envir = get(ev)),
        file.path(
            sv[["ev"]]$path,
            sv[["ev"]]$name_e
        ), sv[["ev"]]$path)

    group_gen <- sv[["ev"]]$group_gen

    # creating the dir to outputs
    if (tolower(tool) == "ebseq") {
        dir <- paste0(
            path, "/EBSeq_Results.",
            tolower(group_gen), "_", toupper(name)
        )
        resultados_de <- sv[["ev"]]$resultados_de_ebseq[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_ebseq
    } else if (tolower(tool) == "edger") {
        dir <- paste0(
            path, "/edgeR_Results.",
            tolower(group_gen), "_", toupper(name)
        )
        resultados_de <- sv[["ev"]]$resultados_de_edger[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_edger
    } else if (tolower(tool) == "deseq2") {
        dir <- paste0(
            path, "/DESeq2_Results.",
            tolower(group_gen), "_", toupper(name)
        )
        resultados_de <- sv[["ev"]]$resultados_de_deseq2[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_deseq2
    } else if (tolower(tool) == "crosstable.deseq2") {
        dir <- paste0(path, "/CrossData_deseq2")
        resultados_de <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_deseq2
    } else if (tolower(tool) == "crosstable.edger") {
        dir <- paste0(path, "/CrossData_edger")
        resultados_de <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_edger
    } else if (tolower(tool) == "crosstable.ebseq") {
        dir <- paste0(path, "/CrossData_ebseq")
        resultados_de <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_ebseq
    } else {
        stop(message(
            "Please, insert a valid tool name!",
            " ('EBSeq', 'DESeq2' or 'edgeR')"
        ))
    }
    dir.create(file.path(dir, "PCA_Plots"), showWarnings = FALSE)

    cond_heatmap <- eval(parse(text = paste0(
        "sv",
        "[['ev']]$cond_heatmap"
    )))

    patients_stay <- unlist(strsplit(gsub("G", "", pair_name), "_over_"))

    # keep the desired group pair
    tmp <- cond_heatmap %in% patients_stay
    normalized_expression <- normalized_expression[, tmp]

    cond_heatmap <- droplevels(cond_heatmap[cond_heatmap %in% patients_stay])

    color_pallete <- RColorBrewer::brewer.pal(8, "Set2")

    if (nrow(resultados_de) != 0) {
        pca_local(size = "All", file = paste0(
            "/PCA_Plots/",
            "PCA_GENETOP=all_", pair_name
        ))
        pca_local(
            size = nrow(resultados_de) %/% 4,
            file = paste0(
                "/PCA_Plots/PCA_GENETOP=",
                nrow(resultados_de) %/% 4,
                "_", pair_name
            )
        )
        pca_local(
            size = nrow(resultados_de) %/% 2,
            file = paste0(
                "/PCA_Plots/PCA_GENETOP=",
                nrow(resultados_de) %/% 2,
                "_", pair_name
            )
        )
        pca_local(
            size = (3 * nrow(resultados_de)) %/% 4,
            file = paste0(
                "/PCA_Plots/PCA_GENETOP=",
                (3 * nrow(resultados_de)) %/% 4,
                "_", pair_name
            )
        )
    } else {
        message("There's no DE gene in your query!!")
    }
}
