#' Run EBSeq gene Differential Expression Analysis (DEA).
#'
#' @param dataType Type of data. It could be \code{"gene"} or \code{"isoform"}.
#' @param workDir
#' @param Name
#' @param env
#' @param tumor
#' @param groupGen A character string indicating the groups generation function:
#'   \itemize{\item{"mclust" - }{\code{groups_identification_mclust()};}
#'   \item{"CoxHR" - }{\code{groups_identification_coxHR()};} \item{"clinical" -
#'   }{\code{groups_identification_clinical()}.} }
#' @param clinical_pair A character string containing one of the groups selected
#'   after statistical analysis runned in \code{check_clinical_terms()}
#'   function.
#' @param pairName A character string indicating which condition name should be
#'   used. When there are only two groups the default is \code{"G2_over_G1"}.
#' @param rounds Numerical value indicating the number of iterations. It is
#'   recommended to check the Alpha and Beta convergence plots in output and
#'   adjust this value until the hyper-parameter estimations converged. The
#'   default is \code{7}.
#' @param normType A character string indicating which EBSeq normalization
#'   factors type should be used in the analysis "QuantileNorm" or "MedianNorm".
#'   The default is \code{"all"}.
#' @param EBTest.Qtrm,EBTest.QtrmCut Numerical value. It is removed from the
#'   analysis genes with EBTest.Qtrm th quantile < = EBTest.QtrmCut. More
#'   details in EBSeq \link{EBTest} page. The default is \code{EBTest.Qtrm =
#'   0.75} and \code{EBTest.QtrmCut = 10}.
#' @param p.cutoff Numerical value indicating the maximum value for p-values.
#'   The default is \code{0.05}.
#' @param FDR.cutoff Numerical value indicating the maximum value for FDR
#'   values. The default is \code{0.05}.
#' @param FC.cutoff Numerical value indicating the maximum value for Fold Change
#'   (FC). The default is \code{2}.
#' @param Width,Height,Res,Unit,image.format
#' @param Bullard.quantile Numerical value indicating the quantile for the
#'   Bullard's normalization. The default is \code{0.75}.
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#'
#' @return A matrix with DE genes in row and statistical values in columns.
#' @export
#'
#' @examples
#' \dontrun{
#' #considering concatenate_files and groups_identification already runned
#' dea_EBSeq("gene", pairName = "G2_over_G1", rounds = 7, Name = "HIF3A", env = "env_name")
#' }
dea_EBSeq <- function(dataType, workDir, Name, env, tumor,
                      groupGen,
                      clinical_pair,
                      pairName = "G2_over_G1",
                      rounds = 7,
                      normType = "QuantileNorm",
                      EBTest.Qtrm = 0.75,
                      EBTest.QtrmCut = 10,
                      p.cutoff = 0.05,
                      FDR.cutoff = 0.05,
                      FC.cutoff = 2,
                      Width = 2000,
                      Height = 1500,
                      Res = 300,
                      Unit = "px",
                      image.format = "png",
                      Bullard.quantile = 0.75) {

    #verifying if the package is already installed
    # to.load <- c("RColorBrewer", "ggbiplot", "devtools", "methods", "pheatmap", "ggplot2")

    #local functions ####
    ebseq_plots <- function(test, Pairs){

        if (Pairs > 1) {
            comb_name <- names(tested)[Pairs]
            # tested <- tested[[Pairs]]
        } else {
            comb_name <- pairName#"G2_over_G1"
        }

        #QQP
        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "QQplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "QQplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        # par(mfrow = c(ceiling(group2_number/2), 2))
        par(mfrow = c(1, 2))
        EBSeq::QQP(test)
        dev.off()
        #DenNHist
        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "DenNHistplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "DenNHistplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 2))
        EBSeq::DenNHist(test)
        dev.off()

        ######PlotPostVsRawFC plot
        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "PostVsRawFCplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "PostVsRawFCplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        EBSeq::PlotPostVsRawFC(test, FC)
        dev.off()

        #Get the iters and generate a table    iters = convergência??
        # ITERS <- matrix(ncol = length(table(test$Beta))+2, nrow = rounds)
        ITERS <- cbind(test$Alpha, test$Beta,
                       test$P)
        rownames(ITERS) <- rownames(test$Alpha)
        colnames(ITERS) <- c("Alpha", "Beta_Ng1", "P")
        write.csv(ITERS, file = paste0(DIR, "/ITERS.results_", comb_name, ".csv"))

        # Print the plot
        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "Convergence.plot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "Convergence.plot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        par(mfrow = c(1,3))
        plot(test$Alpha, frame=FALSE, xaxp  = c(1, 10, n=1), xlab="Rounds",
             ylab = "Alpha", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5,
             cex.lab = 1.5)
        lines(test$Alpha, col = "#FDB462", lwd = 2)
        plot(test$Beta, frame = FALSE, xaxp = c(1, 10, n = 1), xlab = "Rounds",
             ylab = "Beta", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5,
             cex.lab = 1.5)
        lines(test$Beta, col = "#FDB462", lwd = 2)
        plot(test$P, frame = FALSE, xaxp = c(1, 10, n = 1), xlab = "Rounds",
             ylab = "P", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5, cex.lab = 1.5)
        lines(test$P, col = "#FDB462", lwd = 2)
        dev.off()

        #Generate table with values of PPEE and PPDE
        EBSeq.RESULTS <- as.data.frame(EBSeq::GetPPMat(test))
        if (tolower(dataBase) == "legacy") {
            GeneSymbol <- strsplit(row.names(EBSeq.RESULTS), split = "\\|")
            GeneSymbol <- as.data.frame(GeneSymbol)
            GeneSymbol <- t(GeneSymbol)
            EBSeq.RESULTS[, 2] <- GeneSymbol[, 1]
            colnames(EBSeq.RESULTS)[1:2] <- c("FDR", "GeneSymbol")
            EBSeq.RESULTS$GeneID <- GeneSymbol[, 2]
            EBSeq.RESULTS$PostFC <- FC$PostFC
            EBSeq.RESULTS$log2FC <- log2(EBSeq.RESULTS$PostFC)

            EBSeq.RESULTS$FC <- ifelse(EBSeq.RESULTS$PostFC < 1,
                                       -1/EBSeq.RESULTS$PostFC,
                                       EBSeq.RESULTS$PostFC)

            # EBSeq.RESULTS <- na.omit(EBSeq.RESULTS)
            EBSeq.RESULTS <- EBSeq.RESULTS[, c(2, 3, 1, 6, 5, 4)]
        } else {
            colnames(EBSeq.RESULTS)[1] <- c("FDR")
            EBSeq.RESULTS$PostFC <- FC$PostFC
            EBSeq.RESULTS$log2FC <- log2(EBSeq.RESULTS$PostFC)

            EBSeq.RESULTS$FC <- ifelse(EBSeq.RESULTS$PostFC < 1,
                                       -1/EBSeq.RESULTS$PostFC,
                                       EBSeq.RESULTS$PostFC)

            # EBSeq.RESULTS <- na.omit(EBSeq.RESULTS)
            EBSeq.RESULTS <- EBSeq.RESULTS[, c(1, 4, 3)]
        }

        if (exists("Results.Completed.EBSeq", envir = get(envir_link))){
            Results.Completed <- string_vars[["envir_link"]]$Results.Completed.EBSeq
            resultadosDE <- string_vars[["envir_link"]]$resultadosDE.EBSeq
        }

        Results.Completed[[comb_name]] <- EBSeq.RESULTS
        assign("Results.Completed.EBSeq", Results.Completed, envir = get(envir_link))

        #for DE
        ResultadosDE <- EBSeq.RESULTS[abs(EBSeq.RESULTS$FC) > FC.cutoff, ]
        ResultadosDE <- ResultadosDE[ResultadosDE$FDR < FDR.cutoff, ]
        ResultadosDE <- ResultadosDE[order(ResultadosDE$FDR, -ResultadosDE$FC), ]
        write.csv(ResultadosDE, file = paste0(DIR, "/ResultsDE_", comb_name, ".csv"))
        resultadosDE[[comb_name]] <- ResultadosDE
        assign("resultadosDE.EBSeq", resultadosDE, envir = get(envir_link))
    }

    volcano <- function(Results.Completed, Pairs){

        comb_name <- Pairs
        EBSeq.RESULTS <- Results.Completed[[Pairs]]

        # START VOLCANO PLOT
        EBSeq.RESULTS$Colour = rgb(100, 100, 100, 50, maxColorValue = 255)

        # Set new column values to appropriate colours
        EBSeq.RESULTS$Colour[EBSeq.RESULTS$log2FC >= log2(FC.cutoff) & EBSeq.RESULTS$FDR <= FDR.cutoff] <- rgb(222, 22, 22, 50, maxColorValue = 255)
        EBSeq.RESULTS$Colour[EBSeq.RESULTS$log2FC <= -log2(FC.cutoff) & EBSeq.RESULTS$FDR <= FDR.cutoff] <- rgb(56, 50, 237, 50, maxColorValue = 255)

        ####Volcano Plot
        axislimits_x <- ceiling(max(c(-min(EBSeq.RESULTS$log2FC, na.rm = TRUE) - 1,
                                      max(EBSeq.RESULTS$log2FC, na.rm = TRUE) + 1)))

        log.10.FDR <- -log10(EBSeq.RESULTS$FDR)
        new.inf <- log.10.FDR[order(log.10.FDR, decreasing = TRUE)]
        log.10.FDR[log.10.FDR == "Inf"] <- (new.inf[new.inf != "Inf"][1]+1)

        axislimits_y <- ceiling(max(log.10.FDR, na.rm = TRUE))+1

        message("Start volcano plot...")
        #Volcano Plot
        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, paste0("VolcanoPlot_Basic_",
                                     comb_name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, paste0("VolcanoPlot_Basic_",
                                     comb_name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        par(mar = c(4,6,3,2), mgp = c(2,.7,0), tck = -0.01)
        plot(EBSeq.RESULTS$log2FC, log.10.FDR, axes = FALSE,
             xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
             xlab = expression('log'[2]*'(FC)'),
             # xlab = bquote(.("") ~ 'log'[2]*.('(FC)')),
             ylab = "",
             # main = "Volcano Plot",
             cex.lab = 1.5, cex.main = 2,
             cex.sub = 2,
             pch = 16, col = EBSeq.RESULTS$Colour, cex = 2.5, las = 1)
        title(ylab = "-log(FDR)", line = 4, cex.lab = 1.5, family = "Calibri Light")
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5, las = 1)
        abline(v = log2(FC.cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(v = -log2(FC.cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(h = -log10(FDR.cutoff), col = "black", lwd = 4, lty = 3)
        text(x = -axislimits_x + 0.3, y = axislimits_y/10, labels =
                 length(EBSeq.RESULTS$Colour[EBSeq.RESULTS$log2FC <= -log2(FC.cutoff) & EBSeq.RESULTS$FDR <= FDR.cutoff]),
             cex = 1, col = "blue")
        # text(x = -axislimits_x + 0.3, y = ((axislimits_y/10)+1.1), labels = "DOWN",
        #      cex = 0.8, col = "blue")
        text(x = axislimits_x - 0.4, y = axislimits_y/10, labels =
                 length(EBSeq.RESULTS$Colour[EBSeq.RESULTS$log2FC >= log2(FC.cutoff) & EBSeq.RESULTS$FDR <= FDR.cutoff]),
             cex = 1, col = "red")
        # text(x = axislimits_x - 0.4, y = ((axislimits_y/10)+1.1), labels = "UP",
        #      cex = 0.8, col = "red")
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        legend("topright", border = FALSE, bty = "n",
               legend=c(paste0("-log(", FDR.cutoff, ") = ", round(-log10(FDR.cutoff), 2)),
                        paste0("\u00B1", " log", "\u2082","(", FC.cutoff, ") = ",log2(FC.cutoff)),
                        "UP", "DOWN"), pt.cex = c(0, 0, 1.8, 1.8),
               lty = c(3, 6, 0, 0), pch = c(-1, -1, 16, 16), cex = c(0.8, 0.8, 0.8, 0.8), lwd = c(2, 2, 5, 5),
               col = c("Black", "Black", "red3", "blue3"))
        dev.off()
    }

    # Code ####
    dataType <- gsub(" ", "_", dataType)
    Name <- gsub("-", "_", Name)

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))

    string_vars <- list(envir_link = get(envir_link))
    if (missing("workDir")){
        workDir <- string_vars[["envir_link"]]$workDir
    }

    dataBase <- string_vars[["envir_link"]]$dataBase

    assign("PATH", file.path(workDir, "GDCRtools", toupper(string_vars[["envir_link"]]$tumor), "Analyses"), envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
        dir.create(PATH, showWarnings = FALSE)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    #creating the dir to outputs
    dir.create(paste0(PATH, "/EBSeq_Results.", tolower(dataType), "_", toupper(Name)),
               showWarnings = FALSE)
    DIR <- paste0(PATH, "/EBSeq_Results.", tolower(dataType), "_", toupper(Name))
    # dir.create(file.patMh(DIR, "PCA_Plots"), showWarnings = FALSE)

    Levels <- gsub("[Gg]", "", unlist(strsplit(pairName, "_over_")), perl = TRUE)

    # For groups
    if (tolower(groupGen) == "mclust") {
        Grupos.EBSeq <- as.data.frame(string_vars[["envir_link"]]$GROUPS)
        group2_number <- max(Grupos.EBSeq[, "Selected_classification"])
        Grupos.EBSeq[, "Selected_classification"] <- factor(Grupos.EBSeq[, "Selected_classification"], levels = Levels)
        colnames(Grupos.EBSeq) <- c("condition", "type")
        Grupos.EBSeq[, 2] <- c("paired-end")
        Grupos.EBSeq <- na.omit(Grupos.EBSeq)
    } else if (tolower(groupGen) == "clinical") {
        Grupos.EBSeq <- as.data.frame(string_vars[["envir_link"]]$clinical_groups_clinical[[clinical_pair]])
        group2_number <- max(Grupos.EBSeq[, "Selected_classification"])
        Grupos.EBSeq[, "Selected_classification"] <- factor(Grupos.EBSeq[, "Selected_classification"])
        colnames(Grupos.EBSeq) <- c("condition", "type")
        Grupos.EBSeq[, 2] <- c("paired-end")
    } else if (tolower(groupGen) == "coxhr") {
        Grupos.EBSeq <- as.data.frame(string_vars[["envir_link"]]$clinical_groups)
        Grupos.EBSeq$classification <- gsub("low", "1", Grupos.EBSeq$classification)
        Grupos.EBSeq$classification <- gsub("high", "2", Grupos.EBSeq$classification)
        group2_number <- max(as.numeric(Grupos.EBSeq[, "classification"]))
        Grupos.EBSeq[, "classification"] <- factor(Grupos.EBSeq[, "classification"])
        colnames(Grupos.EBSeq) <- c("type", "condition")
        Grupos.EBSeq[, 1] <- c("paired-end")
        Grupos.EBSeq <- Grupos.EBSeq[, c(2,1)]
    } else {
        stop(message("Please insert a valid 'groupGen' value!! ('mclust', 'coxHR' or 'clinical')"))
    }

    assign("condHeatmap", Grupos.EBSeq[, 1], envir = get(envir_link))

    #selecting specifics patients
    completed.matrix <- string_vars[["envir_link"]]$gene_tumor_not_normalized[, rownames(Grupos.EBSeq)]

    # message("Starting EBSeq analysis...")
    if (tolower(normType) == "quantilenorm") {
        Sizes <- EBSeq::QuantileNorm(completed.matrix, Bullard.quantile)
    } else if (tolower(normType) == "mediannorm") {
        Sizes <- EBSeq::MedianNorm(completed.matrix)
    } else {
        stop(message("Please insert a valid 'normType' argument!!"))
    }

    if (!exists("Results.Completed.EBSeq", envir = get(envir_link))) {
        combinations <- combn(1:group2_number, 2)
        tested <- vector("list", group2_number)
        resultadosDE <- vector("list", group2_number)
        Results.Completed <- vector("list", group2_number)
        combinations_names <-  apply(combinations, 2, function(x) {
            paste0("G", x[2], "_over_", "G", x[1])
        })
        names(tested) <- combinations_names
        names(resultadosDE) <- combinations_names
        names(Results.Completed) <- combinations_names
    }

    message("Starting EBSeq analysis...")
    # if (group2_number == 2){
        EBOut <- EBSeq::EBTest(Data = completed.matrix,
                        Conditions = Grupos.EBSeq[, "condition"],
                        #filtrando ruído
                        Qtrm = EBTest.Qtrm,
                        QtrmCut = EBTest.QtrmCut,
                        sizeFactors = Sizes,
                        maxround = rounds)
        #generating Fold Change
        FC <- EBSeq::PostFC(EBOut)
        # str(GeneFC)$Direction

        ebseq_plots(EBOut, 0)
        suppressWarnings(volcano(string_vars[["envir_link"]]$Results.Completed.EBSeq, pairName))
    #} else if (group2_number > 2){

    #     PosParti <- EBSeq::GetPatterns(Grupos.EBSeq[, "condition"])
    #     # MultiSize <- EBSeq::MedianNorm(completed.matrix)
    #
    #     MultiOut <- EBSeq::EBMultiTest(Data = completed.matrix,
    #                             NgVector = NULL,
    #                             Conditions = Grupos.EBSeq[, "condition"],
    #                             Qtrm = 0.75,
    #                             QtrmCut = 10,
    #                             AllParti = PosParti,
    #                             sizeFactors = Sizes,
    #                             maxround = rounds)
    #     # EB.RESULTS.partial <- MultiOut
    #     # assign("EB.RESULTS.partial", EB.RESULTS.partial, envir = get(envir_link))
    #     #FoldChanges
    #     FC <- EBSeq::GetMultiFC(MultiOut)
    #     # assign("FC", FC, envir = get(envir_link))
    #     MultiPP <- EBSeq::GetMultiPP(MultiOut)
    #
    #
    #     # EB.RESULTS <- MultiPP
    #     # assign("EB.RESULTS", EB.RESULTS, envir = get(envir_link))
    #
    #
    #     count <- 1
    #     for (Pairs in seq(ncol(combinations), 1)) {
    #         Results.Completed[[count]] <- cbind(FC$FCMat[, Pairs], FC$Log2FCMat[, Pairs])
    #         ebseq_plots(MultiOut, Pairs)# ebseq_plots(, )
    #         volcano(Results.Completed[[Pairs]], Pairs)
    #
    #         count <- count + 1
    #     }
    # } else {
    #     stop("Please, insert valid data type!!")
    # }

    #A disadvantage to and serious risk of using fold change (for DE) is that it is biased [2]
    #and may miss differentially expressed genes with large differences (B-A) but small ratios (A/B),
    #leading to a high miss rate at high intensities.
    #source: wiki

    #Generatin' normalized data matrix
    NormalizedExpression <- EBSeq::GetNormalizedMat(completed.matrix, Sizes)
    write.csv(NormalizedExpression, file = paste0(DIR,
                                                  "/Normalized_Expression.csv"))

    assign("NormalizedExpression.EBSeq", NormalizedExpression, envir = get(envir_link))
    assign("Tool", "EBSeq", envir = get(envir_link))

    message("Done!")
}
