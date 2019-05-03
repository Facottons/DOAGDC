#' Run EBSeq gene Differential Expression Analysis (DEA).
#'
#' @param Name
#' @param workDir
#' @param env
#' @param tumor
#' @param groupGen A character string indicating the groups generation function:
#'   \itemize{\item{"mclust" - }{\code{groups_identification_mclust()};}
#'   \item{"CoxHR" - }{\code{groups_identification_coxHR()};} \item{"clinical" -
#'   }{\code{groups_identification_clinical()}.} }
#' @param clinical_pair A character string containing one of the group pairs
#'   selected after statistical analysis runned in \code{clinical_terms()}
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
#' @param EBTest_Qtrm,EBTest_QtrmCut Numerical value. It is removed from the
#'   analysis genes with EBTest_Qtrm th quantile < = EBTest_QtrmCut. More
#'   details in EBSeq \link{EBTest} page. The default is \code{EBTest_Qtrm =
#'   0.75} and \code{EBTest_QtrmCut = 10}.
#' @param p_cutoff Numerical value indicating the maximum value for p-values.
#'   The default is \code{0.05}.
#' @param FDR_cutoff Numerical value indicating the maximum value for FDR
#'   values. The default is \code{0.05}.
#' @param FC_cutoff Numerical value indicating the maximum value for Fold Change
#'   (FC). The default is \code{2}.
#' @param Width,Height,Res,Unit,image_format
#' @param Bullard_quantile Numerical value indicating the quantile for the
#'   Bullard's normalization. The default is \code{0.75}.
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#'
#' @return A matrix with DE genes in row and statistical values in columns.
#'
#' @examples
#' \dontrun{
#' #considering concatenate_files and groups_identification already runned
#' dea_EBSeq(pairName = "G2_over_G1", rounds = 7, Name = "HIF3A", env = "env_name")
#' }
dea_EBSeq <- function(Name, workDir, env, tumor,
                      groupGen,
                      clinical_pair,
                      pairName = "G2_over_G1",
                      rounds = 7,
                      normType = "QuantileNorm",
                      EBTest_Qtrm = 0.75,
                      EBTest_QtrmCut = 10,
                      p_cutoff = 0.05,
                      FDR_cutoff = 0.05,
                      FC_cutoff = 2,
                      Width = 2000,
                      Height = 1500,
                      Res = 300,
                      Unit = "px",
                      image_format = "png",
                      Bullard_quantile = 0.75) {

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
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "QQplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "QQplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        # par(mfrow = c(ceiling(group2_number/2), 2))
        par(mfrow = c(1, 2))
        EBSeq::QQP(test)
        dev.off()
        #DenNHist
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "DenNHistplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "DenNHistplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 2))
        EBSeq::DenNHist(test)
        dev.off()

        ######PlotPostVsRawFC plot
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "PostVsRawFCplot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "PostVsRawFCplot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
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
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "Convergence.plot.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "Convergence.plot.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1,3))
        plot(test$Alpha, frame=FALSE, xlab="Rounds",
             ylab = "Alpha", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5,
             cex.lab = 1.5)
        lines(test$Alpha, col = "#FDB462", lwd = 2)
        plot(test$Beta, frame = FALSE, xlab = "Rounds",
             ylab = "Beta", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5,
             cex.lab = 1.5)
        lines(test$Beta, col = "#FDB462", lwd = 2)
        plot(test$P, frame = FALSE, xlab = "Rounds",
             ylab = "P", cex = 1.1, col = "#FECB92", pch = 16, cex.axis = 1.5, cex.lab = 1.5)
        lines(test$P, col = "#FDB462", lwd = 2)
        dev.off()

        #Generate table with values of PPEE and PPDE
        EBSeq_RESULTS <- as.data.frame(EBSeq::GetPPMat(test))
        if (tolower(dataBase) == "legacy") {
            GeneSymbol <- strsplit(row.names(EBSeq_RESULTS), split = "\\|")
            GeneSymbol <- as.data.frame(GeneSymbol)
            GeneSymbol <- t(GeneSymbol)
            EBSeq_RESULTS[, 2] <- GeneSymbol[, 1]
            colnames(EBSeq_RESULTS)[1:2] <- c("FDR", "GeneSymbol")
            EBSeq_RESULTS$GeneID <- GeneSymbol[, 2]
            EBSeq_RESULTS$PostFC <- FC$PostFC
            EBSeq_RESULTS$log2FC <- log2(EBSeq_RESULTS$PostFC)

            EBSeq_RESULTS$FC <- ifelse(EBSeq_RESULTS$PostFC < 1,
                                       -1/EBSeq_RESULTS$PostFC,
                                       EBSeq_RESULTS$PostFC)

            # EBSeq_RESULTS <- na.omit(EBSeq_RESULTS)
            EBSeq_RESULTS <- EBSeq_RESULTS[, c(2, 3, 1, 6, 5, 4)]
        } else {
            colnames(EBSeq_RESULTS)[1] <- c("FDR")
            EBSeq_RESULTS$PostFC <- FC$PostFC
            EBSeq_RESULTS$log2FC <- log2(EBSeq_RESULTS$PostFC)

            EBSeq_RESULTS$FC <- ifelse(EBSeq_RESULTS$PostFC < 1,
                                       -1/EBSeq_RESULTS$PostFC,
                                       EBSeq_RESULTS$PostFC)

            # EBSeq_RESULTS <- na.omit(EBSeq_RESULTS)
            EBSeq_RESULTS <- EBSeq_RESULTS[, c(1, 4, 3)]
        }

        if (exists("Results_Completed_EBSeq", envir = get(envir_link))){
            Results_Completed <- string_vars[["envir_link"]]$Results_Completed_EBSeq
            resultadosDE <- string_vars[["envir_link"]]$resultadosDE.EBSeq
        }

        Results_Completed[[comb_name]] <- EBSeq_RESULTS
        assign("Results_Completed_EBSeq", Results_Completed, envir = get(envir_link))

        #for DE
        ResultadosDE <- EBSeq_RESULTS[abs(EBSeq_RESULTS$FC) > FC_cutoff, ]
        ResultadosDE <- ResultadosDE[ResultadosDE$FDR < FDR_cutoff, ]
        ResultadosDE <- ResultadosDE[order(ResultadosDE$FDR, -ResultadosDE$FC), ]
        write.csv(ResultadosDE, file = paste0(DIR, "/ResultsDE_", comb_name, ".csv"))
        resultadosDE[[comb_name]] <- ResultadosDE
        assign("resultadosDE.EBSeq", resultadosDE, envir = get(envir_link))
    }

    volcano <- function(Results_Completed, Pairs){

        comb_name <- Pairs
        EBSeq_RESULTS <- Results_Completed[[Pairs]]

        # START VOLCANO PLOT
        EBSeq_RESULTS$Colour = rgb(100, 100, 100, 50, maxColorValue = 255)

        # Set new column values to appropriate colours
        EBSeq_RESULTS$Colour[EBSeq_RESULTS$log2FC >= log2(FC_cutoff) & EBSeq_RESULTS$FDR <= FDR_cutoff] <- rgb(222, 22, 22, 50, maxColorValue = 255)
        EBSeq_RESULTS$Colour[EBSeq_RESULTS$log2FC <= -log2(FC_cutoff) & EBSeq_RESULTS$FDR <= FDR_cutoff] <- rgb(56, 50, 237, 50, maxColorValue = 255)

        ####Volcano Plot
        axislimits_x <- ceiling(max(c(-min(EBSeq_RESULTS$log2FC, na.rm = TRUE) - 1,
                                      max(EBSeq_RESULTS$log2FC, na.rm = TRUE) + 1)))

        log_10_FDR <- -log10(EBSeq_RESULTS$FDR)
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
        plot(EBSeq_RESULTS$log2FC, log_10_FDR, axes = FALSE,
             xlim = c(-axislimits_x, axislimits_x), ylim = c(0, axislimits_y),
             xlab = expression('log'[2]*'(FC)'),
             # xlab = bquote(.("") ~ 'log'[2]*.('(FC)')),
             ylab = "",
             # main = "Volcano Plot",
             cex.lab = 1.5, cex.main = 2,
             cex.sub = 2,
             pch = 16, col = EBSeq_RESULTS$Colour, cex = 2.5, las = 1)
        title(ylab = "-log(FDR)", line = 4, cex.lab = 1.5, family = "Calibri Light")
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5, las = 1)
        abline(v = log2(FC_cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(v = -log2(FC_cutoff), col = "black", lty = 6, cex = 0.8,
               lwd = 4)
        abline(h = -log10(FDR_cutoff), col = "black", lwd = 4, lty = 3)
        text(x = -axislimits_x + 0.3, y = axislimits_y/10, labels =
                 length(EBSeq_RESULTS$Colour[EBSeq_RESULTS$log2FC <= -log2(FC_cutoff) & EBSeq_RESULTS$FDR <= FDR_cutoff]),
             cex = 1, col = "blue")
        # text(x = -axislimits_x + 0.3, y = ((axislimits_y/10)+1.1), labels = "DOWN",
        #      cex = 0.8, col = "blue")
        text(x = axislimits_x - 0.4, y = axislimits_y/10, labels =
                 length(EBSeq_RESULTS$Colour[EBSeq_RESULTS$log2FC >= log2(FC_cutoff) & EBSeq_RESULTS$FDR <= FDR_cutoff]),
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

    assign("PATH", file.path(workDir, "DOAGDC", toupper(string_vars[["envir_link"]]$tumor), "Analyses"),
           envir = get(envir_link))

    assign("groupGen", groupGen, envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
        dir.create(PATH, showWarnings = FALSE)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    #creating the dir to outputs
    dir.create(paste0(PATH, "/EBSeq_Results.", tolower(groupGen), "_", toupper(Name)),
               showWarnings = FALSE)
    DIR <- paste0(PATH, "/EBSeq_Results.", tolower(groupGen), "_", toupper(Name))
    # dir.create(file.patMh(DIR, "PCA_Plots"), showWarnings = FALSE)

    # level[1] over level[2]
    Levels <- gsub("[Gg]", "", unlist(strsplit(pairName, "_over_")), perl = TRUE)

    # For groups
    if (tolower(groupGen) == "mclust") {
        Grupos_EBSeq <- as.data.frame(string_vars[["envir_link"]]$GROUPS)
        group2_number <- max(Grupos_EBSeq[, "Selected_classification"])
        Grupos_EBSeq[, "Selected_classification"] <- factor(Grupos_EBSeq[, "Selected_classification"],
                                                            levels = Levels)
        colnames(Grupos_EBSeq) <- c("condition", "type")
        Grupos_EBSeq[, 2] <- c("paired-end")
        Grupos_EBSeq <- na.omit(Grupos_EBSeq)
    } else if (tolower(groupGen) == "clinical") {
        Grupos_EBSeq <- as.data.frame(string_vars[["envir_link"]]$clinical_groups_clinical[[clinical_pair]])
        group2_number <- max(Grupos_EBSeq[, "Selected_classification"])
        Grupos_EBSeq[, "Selected_classification"] <- factor(Grupos_EBSeq[, "Selected_classification"])
        colnames(Grupos_EBSeq) <- c("condition", "type")
        Grupos_EBSeq[, 2] <- c("paired-end")
    } else if (tolower(groupGen) == "coxhr") {
        Grupos_EBSeq <- as.data.frame(string_vars[["envir_link"]]$clinical_groups)
        Grupos_EBSeq$classification <- gsub("low", "1", Grupos_EBSeq$classification)
        Grupos_EBSeq$classification <- gsub("high", "2", Grupos_EBSeq$classification)
        group2_number <- max(as.numeric(Grupos_EBSeq[, "classification"]))
        Grupos_EBSeq[, "classification"] <- factor(Grupos_EBSeq[, "classification"])
        colnames(Grupos_EBSeq) <- c("type", "condition")
        Grupos_EBSeq[, 1] <- c("paired-end")
        Grupos_EBSeq <- Grupos_EBSeq[, c(2,1)]
    } else {
        stop(message("Please insert a valid 'groupGen' value!! ('mclust', 'coxHR' or 'clinical')"))
    }

    # check patient in common
    tmp <- rownames(Grupos_EBSeq) %in% colnames(string_vars[["envir_link"]]$gene_tumor_not_normalized)

    Grupos_EBSeq <- Grupos_EBSeq[tmp, ]

    assign("condHeatmapq", Grupos_EBSeq[, 1], envir = get(envir_link))

    #selecting specifics patients
    completed_matrix <- string_vars[["envir_link"]]$gene_tumor_not_normalized[, rownames(Grupos_EBSeq)]

    # message("Starting EBSeq analysis...")
    if (tolower(normType) == "quantilenorm") {
        Sizes <- EBSeq::QuantileNorm(completed_matrix, Bullard_quantile)
    } else if (tolower(normType) == "mediannorm") {
        Sizes <- EBSeq::MedianNorm(completed_matrix)
    } else {
        stop(message("Please insert a valid 'normType' argument!!"))
    }

    if (!exists("Results_Completed_EBSeq", envir = get(envir_link))) {
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
    }

    message("Starting EBSeq analysis...")
    # if (group2_number == 2){
        EBOut <- EBSeq::EBTest(Data = completed_matrix,
                        Conditions = Grupos_EBSeq[, "condition"],
                        #filtrando ruído
                        Qtrm = EBTest_Qtrm,
                        QtrmCut = EBTest_QtrmCut,
                        sizeFactors = Sizes,
                        maxround = rounds)
        #generating Fold Change
        FC <- EBSeq::PostFC(EBOut)
        # str(GeneFC)$Direction

        ebseq_plots(EBOut, 0)
        suppressWarnings(volcano(string_vars[["envir_link"]]$Results_Completed_EBSeq, pairName))
    #} else if (group2_number > 2){

    #     PosParti <- EBSeq::GetPatterns(Grupos_EBSeq[, "condition"])
    #     # MultiSize <- EBSeq::MedianNorm(completed_matrix)
    #
    #     MultiOut <- EBSeq::EBMultiTest(Data = completed_matrix,
    #                             NgVector = NULL,
    #                             Conditions = Grupos_EBSeq[, "condition"],
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
    #         Results_Completed[[count]] <- cbind(FC$FCMat[, Pairs], FC$Log2FCMat[, Pairs])
    #         ebseq_plots(MultiOut, Pairs)# ebseq_plots(, )
    #         volcano(Results_Completed[[Pairs]], Pairs)
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
    NormalizedExpression <- EBSeq::GetNormalizedMat(completed_matrix, Sizes)
    write.csv(NormalizedExpression, file = paste0(DIR,
                                                  "/Normalized_Expression.csv"))

    assign("NormalizedExpression.EBSeq", NormalizedExpression, envir = get(envir_link))
    assign("Tool", "EBSeq", envir = get(envir_link))

    message("Done!")
}
