#' Generate plot for Principal Component Analysis
#'
#' @param Tool
#' @param dataType
#' @param Name
#' @param workDir
#' @param pairName
#' @param Width,Height,Res,Unit,image.format
#' @param env
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#' @inheritParams GOnto
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' PCA_Analysis("EBSeq", "gene", "HIF3A", pairName = "G2_over_G1", env = "env name without quotes")
#' }
PCA_Analysis <- function(Tool,
                         dataType, Name, workDir,
                         pairName = "G2_over_G1",
                         Width = 4,
                         Height = 2,
                         Res = 300,
                         Unit = "in",
                         image.format = "png",
                         env) {

    # local function ####
    PCA_Anal_local <- function(SIZE, FILE){
        if (tolower(SIZE) == "all"){

            GENETOP <- SIZE
            ParaHeatmaps <- NormalizedExpression[match(rownames(resultadosDE), rownames(NormalizedExpression)), ]
            # colnames(ParaHeatmaps) <- colnames(NormalizedExpression)

            ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                             center = TRUE,
                             scale. = TRUE)

            if (tolower(image.format) == "png") {
                png(filename = paste0(DIR, FILE, ".png"),
                    width = Width, height = Height, res = Res, units = Unit)
            } else if (tolower(image.format) == "svg") {
                svg(filename = paste0(DIR, FILE, ".svg"),
                    width = Width, height = Height, onefile = TRUE)
            } else {
                stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
            }
            a <- ggbiplot::ggbiplot(ir.pca, choices = 1:2,
                                    obs.scale = 1,
                                    var.scale = 1,
                                    groups = condHeatmap,
                                    ellipse = TRUE,
                                    circle = TRUE,
                                    alpha = 0.5,
                                    var.axes = FALSE) +
                scale_color_discrete(name = 'Groups') +
                theme(axis.title.y = element_text(size = rel(0.6), angle = 90)) +
                theme(axis.title.x = element_text(size = rel(0.6), angle = 00)) +
                theme(legend.title = element_text(size=4, face = "bold")) +
                theme(legend.text = element_text(size=4, face = "bold"))
            print(a)
            #  theme(legend.direction = 'horizontal', legend.position = 'top')
            dev.off()

        } else if(is.numeric(SIZE)) {

            GENETOP <- SIZE
            resultadosDE2 <- resultadosDE[order(resultadosDE$FC), ]
            resultadosDE_Down <- resultadosDE2[1:GENETOP, ]
            #row.names(resultadosDE_Down) <- resultadosDE2[1:GENETOP, "gene"]
            resultadosDE3 <- resultadosDE[order(resultadosDE$FC, decreasing = TRUE), ]
            resultadosDE_Up <- resultadosDE3[1:GENETOP, ]
            #row.names(resultadosDE_Up) <- resultadosDE3[1:GENETOP, "gene"]
            resultadosPlot <- rbind(resultadosDE_Down, resultadosDE_Up)

            # Separacao para os heatmaps
            ParaHeatmaps <- NormalizedExpression[match(rownames(resultadosPlot), rownames(NormalizedExpression)), ]
            ParaHeatmaps <- na.exclude(ParaHeatmaps)

            ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                             center = TRUE,
                             scale. = TRUE)

            if (tolower(image.format) == "png") {
                png(filename = paste0(DIR, FILE, "2.png"),
                    width = Width, height = Height, res = Res, units = Unit)
            } else if (tolower(image.format) == "svg") {
                svg(filename = paste0(DIR, FILE, "2.svg"),
                    width = Width, height = Height, onefile = TRUE)
            } else {
                stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
            }
            a <- ggbiplot::ggbiplot(ir.pca, choices = 1:2,
                                    obs.scale = 1,
                                    var.scale = 1,
                                    groups = condHeatmap,
                                    ellipse = TRUE,
                                    circle = TRUE,
                                    alpha = 0.5,
                                    var.axes = FALSE) +
                scale_color_discrete(name = 'Groups') +
                theme(axis.title.y = element_text(size = rel(0.6), angle = 90)) +
                theme(axis.title.x = element_text(size = rel(0.6), angle = 00)) +
                theme(legend.title = element_text(size=4, face = "bold")) +
                theme(legend.text = element_text(size=4, face = "bold"))
            print(a)
            #  theme(legend.direction = 'horizontal', legend.position = 'top')
            dev.off()
        }
    }

    # code ####
    dataType <- gsub(" ", "_", dataType)
    Name <- gsub("-", "_", Name)

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (missing("workDir")){
        workDir <- string_vars[["envir_link"]]$workDir
    }

    # assign("PATH", file.path(workDir, "GDCRtools", toupper(string_vars[["envir_link"]]$tumor),
    #                          "Analyses"), envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
        # dir.create(PATH, showWarnings = FALSE)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    #creating the dir to outputs
    if (tolower(Tool) == "ebseq") {
        DIR <- paste0(PATH, "/EBSeq_Results.",
                      tolower(dataType), "_", toupper(Name))
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.EBSeq[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.EBSeq
    } else if (tolower(Tool) == "edger") {
        DIR <- paste0(PATH, "/edgeR_Results.",
                      tolower(dataType), "_", toupper(Name))
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.edgeR[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.edgeR
    } else if (tolower(Tool) == "deseq2") {
        DIR <- paste0(PATH, "/DESeq2_Results.",
                      tolower(dataType), "_", toupper(Name))
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.DESeq2[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.DESeq2
    } else if (tolower(Tool) == "crosstable.deseq2") {
        DIR <- paste0(PATH, "/CrossData_deseq2")
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.crossed[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.DESeq2
    } else if (tolower(Tool) == "crosstable.edger") {
        DIR <- paste0(PATH, "/CrossData_edger")
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.crossed[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.edgeR
    } else if (tolower(Tool) == "crosstable.ebseq") {
        DIR <- paste0(PATH, "/CrossData_ebseq")
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.crossed[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.EBSeq
    } else {
        stop(message("Please, insert a valid Tool name! ('EBSeq', 'DESeq2' or 'edgeR')"))
    }
    dir.create(file.path(DIR, "PCA_Plots"), showWarnings = FALSE)

    condHeatmap <- string_vars[["envir_link"]]$condHeatmap

    # cores <- ifelse(ParaGrafico == "Low", "#E41A1C","#377EB8")

    if (nrow(resultadosDE) != 0 ){
        PCA_Anal_local(SIZE = "All", FILE = paste0("/PCA_Plots/PCA_GENETOP=all_",
                                                   pairName))
        PCA_Anal_local(SIZE = nrow(resultadosDE)%/%4,
                 FILE = paste0("/PCA_Plots/PCA_GENETOP=", nrow(resultadosDE)%/%4,
                               "_", pairName))
        PCA_Anal_local(SIZE = nrow(resultadosDE)%/%2,
                 FILE = paste0("/PCA_Plots/PCA_GENETOP=", nrow(resultadosDE)%/%2,
                               "_", pairName))
        PCA_Anal_local(SIZE = (3*nrow(resultadosDE))%/%4,
                 FILE = paste0("/PCA_Plots/PCA_GENETOP=", (3*nrow(resultadosDE))%/%4,
                               "_", pairName))
    } else {
        message("There's no DE gene in your query!!")
    }
}
