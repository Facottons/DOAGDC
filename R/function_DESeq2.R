#' Run DESeq2 gene Differential Expression Analysis (DEA).
#'
#' @param Name
#' @param coreNumber A numeric value indicating how many CPU cores should be
#'   used in the analysis. The default value is 2.
#' @param test A character string indicating which test should be used:
#'   \code{"LRT", "wald" or "Default Test"}. The default is \code{"Default
#'   Test"}.
#' @param clinical_pair
#' @param groupGen
#' @param FC_cutoff
#' @param workDir
#' @param tumor
#' @param FDR_cutoff
#' @param Width,Height,Res,Unit,image_format
#' @param env
#' @param cooksCutoff Cooks distance remove outliers from the analysis; More
#'   details in DESeq2 \link{results} page. The default is \code{FALSE}.
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#'
#' @return A matrix with DE genes in row and statistical values in columns.
#' @export
#'
#' @examples
#' \dontrun{
#' #considering concatenate_files and groups_identification already runned
#' dea_DESeq2(Name = "HIF3A", test = "LRT", env = "env name without quotes")
#' }
dea_DESeq2 <- function(Name,
                       coreNumber = 2,
                       test = "Default Test",
                       groupGen,
                       clinical_pair,
                       FC_cutoff = 2,
                       workDir,
                       tumor,
                       FDR_cutoff = 0.05,
                       Width = 2000,
                       Height = 1500,
                       Res = 300,
                       Unit = "px",
                       image_format = "png",
                       env,
                       cooksCutoff = FALSE){

    # DESeq2 normalization does not account for gene length, and there are sound reasons for making
    # that choice when using the data for statistical hypothesis testing. Visualization based on
    # regularized log transformation should not be biased based on gene length.  However, gene
    # expression in RNA-seq does have a gene length bias "built-in"; this is a function of the "count"
    # nature of RNA-seq and not due to any software processing of the data.


    # #verifying if the package is already installed
    # to.load <- c("DESeq2", "gage", "BiocParallel", "IHW")

    #local functions ####

    deseq_plots <- function(dds) {
        if (tolower(image_format) == "png") {
            png(filename = paste0(DIR, "/plotDispEsts.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = paste0(DIR, "/plotDispEsts.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        DESeq2::plotDispEsts(dds, main = paste0("Test = ", toupper(test)),
                             las = 1, ylab = "Dispersion",
                             xlab = "Mean of normalized counts")
        dev.off()

        if (tolower(image_format) == "png") {
            png(filename = paste0(DIR,
                                  "/deseq2_MAplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = paste0(DIR,
                                  "/deseq2_MAplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        DESeq2::plotMA(dds, ylim = c(-log2(FC_cutoff)-1, log2(FC_cutoff)+1), main='DESeq2',
                       las = 1, ylab = "Log2FC", xlab = "Mean of normalized counts")
        abline(h = log2(FC_cutoff), col = "dodgerblue", lwd = 2)
        abline(h = -log2(FC_cutoff), col = "dodgerblue", lwd = 2)
        dev.off()

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "transformation_effect.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "transformation_effect.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(4,4,1,2))
        px <- DESeq2::counts(dds)[,1] / DESeq2::sizeFactors(dds)[1]
        ord <- order(px)
        ord <- ord[px[ord]<150]
        ord <- ord[seq(1, length(ord), length=50)]
        last <- ord[length(ord)]
        vstcol <- c('blue', 'black')
        matplot(px[ord], cbind(SummarizedExperiment::assay(varTrans)[, 1], log2(px))[ord, ], type="l", lty=1,
                col=vstcol, xlab='n', ylab='f(n)', las = 1)
        legend('bottomright', legend = c(expression('variance stabilizing transformation'),
                                         expression(log[2](n/s[1]))), fill=vstcol,
               cex = 0.8)
        dev.off()
    }

    volcano <- function(Results_Completed, Pairs){

        if (Pairs > 0) {
            comb_name <- names(Results_Completed)[Pairs]
            Results_Completed_local <- Results_Completed[[Pairs]]
        } else {
            comb_name <- "G2_over_G1"
            Results_Completed_local <- Results_Completed
        }

        # START VOLCANO PLOT
        Results_Completed_local$Colour = rgb(100, 100, 100, 50, maxColorValue = 255)

        # Set new column values to appropriate colours
        Results_Completed_local$Colour[Results_Completed_local$log2FC >= log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff] <- rgb(222, 22, 22, 50, maxColorValue = 255)
        Results_Completed_local$Colour[Results_Completed_local$log2FC < -log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff] <- rgb(56, 50, 237, 50, maxColorValue = 255)

        ####Volcano Plot
        axislimits_x <- ceiling(max(c(-min(Results_Completed_local$log2FC, na.rm = TRUE) - 1,
                                      max(Results_Completed_local$log2FC, na.rm = TRUE) + 1)))

        log_10_FDR <- -log10(Results_Completed_local$FDR)
        new_inf <- log_10_FDR[order(log_10_FDR, decreasing = TRUE)]
        log_10_FDR[log_10_FDR == "Inf"] <- (new_inf[new_inf != "Inf"][1]+1)

        axislimits_y <- ceiling(max(log_10_FDR, na.rm = TRUE))+1

        message("Start volcano plot...")
        #Volcano Plot
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, paste0("VolcanoPlot_Basic_",
                                     comb_name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, paste0("VolcanoPlot_Basic_",
                                     comb_name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(4,6,3,2), mgp = c(2,.7,0), tck = -0.01)
        plot(Results_Completed_local$log2FC, log_10_FDR, axes = FALSE,
             xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
             xlab = expression('log'[2]*'(FC)'),
             # xlab = bquote(.("") ~ 'log'[2]*.('(FC)')),
             ylab = "",
             # main = "Volcano Plot",
             cex.lab = 1.5, cex.main = 2,
             cex.sub = 2,
             pch = 16, col = Results_Completed_local$Colour, cex = 2.5, las = 1)
        title(ylab = "-log(FDR)", line = 4, cex.lab = 1.5, family = "Calibri Light")
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5, las = 1)
        abline(v = log2(FC_cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(v = -log2(FC_cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(h = -log10(FDR_cutoff), col = "black", lwd = 4, lty = 3)
        text(x = -axislimits_x + 0.3, y = axislimits_y/10, labels =
                 length(Results_Completed_local$Colour[Results_Completed_local$log2FC <= -log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff]),
             cex = 1, col = "blue")
        # text(x = -axislimits_x + 0.3, y = ((axislimits_y/10)+1.1), labels = "DOWN",
        #      cex = 0.8, col = "blue")
        text(x = axislimits_x - 0.4, y = axislimits_y/10, labels =
                 length(Results_Completed_local$Colour[Results_Completed_local$log2FC >= log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff]),
             cex = 1, col = "red")
        # text(x = axislimits_x - 0.4, y = ((axislimits_y/10)+1.1), labels = "UP",
        #      cex = 0.8, col = "red")
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

        legend("topright", border = FALSE, bty = "n",
               legend=c(paste0("-log(", FDR_cutoff, ") = ", round(-log10(FDR_cutoff), 2)),
                        paste0("\u00B1", " log", "\u2082","(", FC_cutoff, ") = ",log2(FC_cutoff)),
                        "UP", "DOWN"), pt.cex = c(0, 0, 1.8, 1.8),
               lty = c(3, 6, 0, 0), pch = c(-1, -1, 16, 16), cex = c(0.8, 0.8, 0.8, 0.8), lwd = c(2, 2, 5, 5),
               col = c("Black", "Black", "red3", "blue3"))
        dev.off()
    }

    # Code ####

    BiocParallel::register(BiocParallel::SnowParam(coreNumber))

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))

    string_vars <- list(envir_link = get(envir_link))

    # dataType <- string_vars[["envir_link"]]$dataType
    # dataType <- gsub(" ", "_", dataType)
    Name <- gsub("-", "_", Name)

    if (missing("workDir")){
        workDir <- string_vars[["envir_link"]]$workDir
    }

    dataBase <- string_vars[["envir_link"]]$dataBase

    assign("PATH", file.path(workDir, "GDCtools",
                             toupper(string_vars[["envir_link"]]$tumor), "Analyses"), envir = get(envir_link))

    assign("groupGen", groupGen, envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
        dir.create(PATH, showWarnings = FALSE)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    dir.create(paste0(PATH, "/DESeq2_Results.", tolower(groupGen), "_", toupper(Name)),
               showWarnings = FALSE)
    DIR <- paste0(PATH, "/DESeq2_Results.", tolower(groupGen), "_", toupper(Name))
    # dir.create(file.path(DIR,
    #                   "PCA_Plots"), showWarnings = FALSE)

    # Cooks distance remove outliers from the analysis, it looks to see how much
    # each sample contributes to a genes overall value fold change, with samples
    # that cause extreme effects removed.

    #Preparing groups data
    if (tolower(groupGen) == "mclust") {
        # Grupos_DESeq2 <- local(GROUPS, envir = get(envir_link))
        Grupos_DESeq2 <- as.data.frame(string_vars[["envir_link"]]$GROUPS)
        Grupos_DESeq2[, "Selected_classification"] <- factor(Grupos_DESeq2[, "Selected_classification"])
        colnames(Grupos_DESeq2) <- c("condition", "type")
        Grupos_DESeq2[, 2] <- c("paired-end")
    } else if (tolower(groupGen) == "clinical") {
        Grupos_DESeq2 <- as.data.frame(string_vars[["envir_link"]]$clinical_groups_clinical[[clinical_pair]])
        Grupos_DESeq2[, "Selected_classification"] <- factor(Grupos_DESeq2[, "Selected_classification"])
        colnames(Grupos_DESeq2) <- c("condition", "type")
        Grupos_DESeq2[, 2] <- c("paired-end")
    } else if (tolower(groupGen) == "coxhr") {
        Grupos_DESeq2 <- as.data.frame(string_vars[["envir_link"]]$clinical_groups)
        Grupos_DESeq2$classification <- gsub("low", "1", Grupos_DESeq2$classification)
        Grupos_DESeq2$classification <- gsub("high", "2", Grupos_DESeq2$classification)
        Grupos_DESeq2[, "classification"] <- factor(Grupos_DESeq2[, "classification"])
        colnames(Grupos_DESeq2) <- c("type", "condition")
        Grupos_DESeq2[, 1] <- c("paired-end")
        Grupos_DESeq2 <- Grupos_DESeq2[, c(2,1)]
    } else {
        stop(message("Please insert a valid 'groupGen' value!! ('mlcust', 'coxHR' or 'clinical')"))
    }


    # check patient in common
    tmp <- rownames(Grupos_DESeq2) %in% colnames(string_vars[["envir_link"]]$gene_tumor_not_normalized)

    Grupos_DESeq2 <- Grupos_DESeq2[tmp, ]

    assign("condHeatmap", Grupos_DESeq2[, 1], envir = get(envir_link))

    #selecting specifics patients
    completed_matrix <- string_vars[["envir_link"]]$gene_tumor_not_normalized[, rownames(Grupos_DESeq2)]

    # colData <- factor(Grupos_DESeq2$condition, levels = c("High", "Low"))
    message("Filtering data ...")
    #remove residual data exclude rows with all zeroes (genes with no counts)
    filter_rows <- rowSums(completed_matrix) >= 10
    completed_matrix <- completed_matrix[filter_rows, ]
    message("It was filtered = ", as.numeric(table(filter_rows)[1]), " genes!")

    completed_matrix <- round(completed_matrix)
    storage.mode(completed_matrix) = "integer"

    dataBase <- string_vars[["envir_link"]]$dataBase
    #Importing data into DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = completed_matrix,
        colData = Grupos_DESeq2,
        design = ~condition)


    # separate conditoins for more than 1 pair
    group2_number <- max(as.numeric(levels(Grupos_DESeq2[, 1])))
    combinations <- combn(1:group2_number, 2)
    tested <- vector("list", group2_number)
    resultadosDE <- vector("list", group2_number)
    Results_Completed <- vector("list", group2_number)
    combinations_names <-  apply(combinations, 2, function(x) {
        paste0("G", x[2], "_over_", "G", x[1])
    })
    names(tested) <- combinations_names
    names(resultadosDE) <- combinations_names
    names(Results_Completed) <- combinations_names
    if (group2_number > 2) {

        message("Performing differential expressed analysis\n")
        for (Pairs in 1:ncol(combinations)) {
            # dds$condition <- factor(dds$condition,
            #                         levels = c(as.character(combinations[1, Pairs]),
            #                                    as.character(combinations[2, Pairs])))
            # dds$condition <- droplevels(dds$condition)
            if (tolower(test) == "lrt") {
                dds <- DESeq2::DESeq(dds, test = "LRT",
                                     reduced =~ 1, parallel = TRUE)
            } else if (tolower(test) == "wald") {
                dds <- DESeq2::DESeq(dds, test = "Wald", parallel = TRUE)
            } else {
                dds <- DESeq2::estimateSizeFactors(dds)
                dds <- DESeq2::estimateDispersions(dds)
                dds <- DESeq2::nbinomWaldTest(dds)
            }

            NormalizedExpression <- DESeq2::counts(dds, normalized=TRUE)
            varTrans <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)

            deseq_plots(dds)


            res <- DESeq2::results(dds, parallel = TRUE, addMLE=FALSE, cooksCutoff = cooksCutoff,
                                   contrast = c("condition", as.character(combinations[2, Pairs]),
                                                as.character(combinations[1, Pairs])))#, alpha = FDR_cutoff) alfa is the FDR

            # resLFC <- ::lfcShrink(dds, coef=paste0("condition_",
            #                                      as.character(combinations[2, Pairs]), "_vs_",
            #                                      as.character(combinations[1, Pairs])), type="apeglm")

            if (tolower(image_format) == "png") {
                png(filename = file.path(DIR, "MA_plot1.png"),
                    width = Width, height = Height, res = Res, units = Unit)
            } else if (tolower(image_format) == "svg") {
                svg(filename = file.path(DIR, "MA_plot1.svg"),
                    width = Width, height = Height, onefile = TRUE)
            } else {
                stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
            }
            DESeq2::plotMA(res, main = "DESeq2", ylim = c(-log2(FC_cutoff)-1, log2(FC_cutoff)+1), las = 1, ylab = "Log2FC",
                           xlab = "Mean of normalized counts")
            abline(h = log2(FC_cutoff), col = "dodgerblue", lwd = 2)
            abline(h = -log2(FC_cutoff), col = "dodgerblue", lwd = 2)
            dev.off()

            Results_Completed_local <- as.data.frame(res)

            if (tolower(dataBase) == "legacy") {
                GeneSymbol <- strsplit(row.names(Results_Completed_local), split = "\\|")
                GeneSymbol <- as.data.frame(GeneSymbol)
                GeneSymbol <- t(GeneSymbol)

                Results_Completed_local$GeneSymbol <- GeneSymbol[, 1]
                Results_Completed_local$GeneID <- GeneSymbol[, 2]

                colnames(Results_Completed_local)[c(2, 5, 6)] <- c("log2FC", "Pvalue", "FDR")

                Results_Completed_local$FC <- 2**(Results_Completed_local$log2FC)

                Results_Completed_local <- Results_Completed_local[, c(7,8,5,6,9,2,1,3,4)]
            } else {
                colnames(Results_Completed_local)[c(2, 5, 6)] <- c("log2FC", "Pvalue", "FDR")

                Results_Completed_local$FC <- 2**(Results_Completed_local$log2FC)

                Results_Completed_local <- Results_Completed_local[, c(5,6,7,2,1,3,4)]
            }

            Results_Completed[[Pairs]] <- Results_Completed_local

            #filtering threashold
            res_filtered <- Results_Completed_local[Results_Completed_local$FDR < FDR_cutoff, ]
            res_filtered <- res_filtered[abs(res_filtered$log2FC) > log2(FC_cutoff), ]
            res_filtered <- na.exclude(res_filtered)
            resultadosDE[[Pairs]] <- as.data.frame(res_filtered)

            write.csv(resultadosDE[[Pairs]], file = paste0(DIR, "/ResultsDE_", combinations_names, ".csv"))
            write.csv(x = NormalizedExpression, file = paste0(DIR,
                                                              "/Normalized_Expression_", tolower(test), ".csv"))

            suppressWarnings(volcano(Results_Completed, Pairs))
        }
    } else {
        message("Performing differential expressed analysis\n")
        # dds$condition <- factor(dds$condition, levels = c("G1","G2"))
        if (tolower(test) == "lrt" || tolower(test) == "wald") {
            dds <- DESeq2::DESeq(dds, test = toupper(test),
                                 reduced =~ 1, parallel = TRUE)
        } else {
            dds <- DESeq2::estimateSizeFactors(dds)
            dds <- DESeq2::estimateDispersions(dds)
            dds <- DESeq2::nbinomWaldTest(dds)
        }

        # save normalized counts
        NormalizedExpression <- DESeq2::counts(dds, normalized=TRUE)

        # rlogTransformation() stabilizes the variance across the range of counts, so genes
        # have a nearly equal effect on the distances and in the PCA plot for example.
        # the rlog is not biased towards long genes (high count genes).
        ## rlogTrans <- rlogTransformation(dds, blind = TRUE)
        # varianceStabilizingTransformation is a faster choice
        varTrans <-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)

        deseq_plots(dds)

        res <- DESeq2::results(dds, parallel = TRUE, addMLE=FALSE, cooksCutoff = cooksCutoff,
                               contrast = c("condition", "2", "1"))#, alpha = FDR_cutoff) alfa is the FDR
                               #, filterFun=ihw)

        # resLFC <- lfcShrink(dds, coef="condition_G2_vs_G1", type="apeglm")

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "MA_plot1.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "MA_plot1.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        # MA plot. Points which fall out of the window are plotted as open triangles
        # pointing either up or down
        DESeq2::plotMA(res, main = "DESeq2", ylim = c(-log2(FC_cutoff)-1, log2(FC_cutoff)+1), las = 1, ylab = "Log2FC",
                       xlab = "Mean of normalized counts")
        abline(h = log2(FC_cutoff), col = "dodgerblue", lwd = 2)
        abline(h = -log2(FC_cutoff), col = "dodgerblue", lwd = 2)
        dev.off()

        Results_Completed_local <- as.data.frame(res)

        if (tolower(dataBase) == "legacy") {
            GeneSymbol <- strsplit(row.names(Results_Completed_local), split = "\\|")
            GeneSymbol <- as.data.frame(GeneSymbol)
            GeneSymbol <- t(GeneSymbol)

            Results_Completed_local$GeneSymbol <- GeneSymbol[, 1]
            Results_Completed_local$GeneID <- GeneSymbol[, 2]

            colnames(Results_Completed_local)[c(2, 5, 6)] <- c("log2FC", "Pvalue", "FDR")

            Results_Completed_local$FC <- 2**(Results_Completed_local$log2FC)
            Results_Completed_local$FC <- ifelse(Results_Completed_local$FC < 1,
                                                 (-1/Results_Completed_local$FC),
                                                 Results_Completed_local$FC)


            Results_Completed_local <- Results_Completed_local[, c(7,8,5,6,9,2,1,3,4)]
        } else {
            colnames(Results_Completed_local)[c(2, 5, 6)] <- c("log2FC", "Pvalue", "FDR")

            Results_Completed_local$FC <- 2**(Results_Completed_local$log2FC)

            Results_Completed_local$FC <- ifelse(Results_Completed_local$FC < 1,
                                                 (-1/Results_Completed_local$FC),
                                                 Results_Completed_local$FC)

            Results_Completed_local <- Results_Completed_local[, c(5,6,7,2,1,3,4)]
        }

        Results_Completed[[1]] <- Results_Completed_local

        #filtering threashold
        res_filtered <- Results_Completed_local[Results_Completed_local$FDR < FDR_cutoff, ]
        res_filtered <- res_filtered[abs(res_filtered$log2FC) > log2(FC_cutoff), ]
        res_filtered <- na.exclude(res_filtered)
        resultadosDE[[1]] <- as.data.frame(res_filtered)

        write.csv(resultadosDE[[1]], file = paste0(DIR, "/ResultsDE_", combinations_names, ".csv"))
        write.csv(x = NormalizedExpression, file = paste0(DIR,
                                                          "/Normalized_Expression_", tolower(test), ".csv"))

        suppressWarnings(volcano(Results_Completed[[1]], 0))
    }

    # #for cross tables later
    # resultadosDE <- tableDE
    # colnames(resultadosDE)[1] <- "FDR"
    #
    # resultadosDE.edgeR <<- resultadosDE
    # resultadosDE <<- resultadosDE
    #
    assign("NormalizedExpression.DESeq2", NormalizedExpression, envir = get(envir_link))
    assign("Results_Completed.DESeq2", Results_Completed, envir = get(envir_link))
    assign("resultadosDE.DESeq2", resultadosDE, envir = get(envir_link))
    # assign("tested", tested, envir = get(envir_link))
    #colnames(tested)[1:2] <- c("log2FC", "log2CPM")

    assign("Tool", "DESeq2", envir = get(envir_link))
    # PCA_Analysis(Tool = "DESeq2", dataType = "gene", Name = Name, env = env)

    remove(first_time, envir = .GlobalEnv)

    message("Done!")
}
