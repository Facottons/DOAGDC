#' Run edgeR gene Differential Expression Analysis (DEA).
#'
#' @param Name
#' @param Method A character string indicating which method should be used:
#'   \code{"exacttest"} or \code{"glmlrt"}. The default is \code{"exacttest"}.
#' @param clinical_pair
#' @param FC_cutoff
#' @param workDir
#' @param env
#' @param FDR_cutoff
#' @param Width,Height,Res,Unit,image_format
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#'
#' @return A matrix with DE genes in row and statistical values in columns.
#'
#' @import edgeR
#' @export
#'
#' @examples
#' \dontrun{
#' #considering concatenate_files and groups_identification already runned
#' dea_edgeR(Name = "HIF3A", env = "env name without quotes")
#' }
dea_edgeR <- function(Name,
                    Method = "exacttest",
                    clinical_pair,
                    groupGen,
                    FC_cutoff = 2,
                    workDir, env,
                    FDR_cutoff = 0.05,
                    Width = 2000,
                    Height = 1500,
                    Res = 300,
                    Unit = "px",
                    image_format = "png"){

    #local functions ####
    plotBCV.modified <- function(y, xlab = "Average logCPM", ylab = "Biological coefficient of variation",
                 pch = 16, cex = 0.2, col.common = "red", col.trend = "blue",
                 col.tagwise = "black", ...) {
        if (!is(y, "DGEList"))
            stop("y must be a DGEList.")
        A <- y$AveLogCPM
        if (is.null(A))
            A <- aveLogCPM(y$counts, offset = getOffset(y))
        disp <- getDispersion(y)
        if (is.null(disp))
            stop("No dispersions to plot")
        if (attr(disp, "type") == "common")
            disp <- rep(disp, length = length(A))
        plot(A, sqrt(disp), xlab = xlab, ylab = ylab, type = "n",
             ...)
        labels <- cols <- lty <- pt <- NULL
        if (!is.null(y$tagwise.dispersion)) {
            points(A, sqrt(y$tagwise.dispersion), pch = pch, cex = cex,
                   col = col.tagwise)
            labels <- c(labels, "Tagwise")
            cols <- c(cols, col.tagwise)
            lty <- c(lty, -1)
            pt <- c(pt, pch)
        }
        if (!is.null(y$common.dispersion)) {
            abline(h = sqrt(y$common.dispersion), col = col.common,
                   lwd = 2)
            labels <- c(labels, "Common")
            cols <- c(cols, col.common)
            lty <- c(lty, 1)
            pt <- c(pt, -1)
        }
        if (!is.null(y$trended.dispersion)) {
            o <- order(A)
            lines(A[o], sqrt(y$trended.dispersion)[o], col = col.trend,
                  lwd = 2)
            labels <- c(labels, "Trend")
            cols <- c(cols, col.trend)
            lty <- c(lty, 1)
            pt <- c(pt, -1)
        }
        #better without border and reduzed cex, only here
        legend("topright", legend = labels, lty = lty, pch = pt,
               cex = 0.8, lwd = 2, col = cols, bty = "n")
        invisible()
    }

    exact_groups_fix <- function(tested, Pairs){

        if (Pairs > 0) {
            comb_name <- names(tested)[Pairs]
            tested_local <- tested[[Pairs]]
        } else {
            comb_name <- "G2_over_G1"
            tested_local <- tested
        }

        tableDE <- edgeR::topTags(tested_local, n = nrow(tested_local$table),
                                  adjust.method="BH", sort.by="PValue")$table

        tableDE$FC <- 2**tableDE[, "logFC"]
        tableDE <- tableDE[, c(5, 1, 2 ,3 ,4)]
        colnames(tableDE)[c(2, 5)] <- c("log2FC", "FDR")

        # tableDE <- tableDE[order(rownames(tableDE)), ]
        if (tolower(dataBase) == "legacy") {
            GeneSymbol <- strsplit(row.names(tableDE), split = "\\|")
            GeneSymbol <- as.data.frame(GeneSymbol)
            GeneSymbol <- t(GeneSymbol)

            tableDE[, "GeneSymbol"] <- GeneSymbol[, 1]
            tableDE[, "GeneID"] <- GeneSymbol[, 2]
            tableDE <- tableDE[, c(6, 7, 4, 5, 1, 2, 3)]
        } else {
            tableDE <- tableDE[, c(4, 5, 1, 2, 3)]
        }

        Results_Completed_local <- tableDE
        tableDE <- tableDE[tableDE$FDR < FDR_cutoff, ]
        tableDE <- tableDE[abs(tableDE$log2FC) > log2(FC_cutoff), ]

        write.csv(x = tableDE, file = paste0(DIR, "/ResultsDE_exactTest_",
                                             comb_name, ".csv"),
                  row.names = TRUE)

        if (Pairs > 0) {
            string_vars[["envir_link"]]$resultadosDE.edgeR[[Pairs]] <- tableDE
            string_vars[["envir_link"]]$Results_Completed_edgeR[[Pairs]] <- Results_Completed_local
            # assign("Results_Completed_edgeR", Results_Completed, envir = get(envir_link))
            # assign("resultadosDE.edgeR", resultadosDE, envir = get(envir_link))
        } else {
            Results_Completed <- vector("list", 1)
            resultadosDE <- vector("list", 1)
            names(Results_Completed) <- comb_name
            names(resultadosDE) <- comb_name

            Results_Completed[[1]] <- Results_Completed_local
            resultadosDE[[1]] <- tableDE
            assign("Results_Completed_edgeR", Results_Completed, envir = get(envir_link))
            assign("resultadosDE.edgeR", resultadosDE, envir = get(envir_link))
        }

        #plots
        message("Plotting...")
        de <- edgeR::decideTestsDGE(tested_local)
        detags <- rownames(dge)[as.logical(de)]
        assign(paste0("detags_", comb_name), detags, envir = get(envir_link))

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, paste0("exactTest_",
                                     comb_name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, paste0("exactTest_",
                                                 comb_name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(5.1, 4.1, 2, 2.1))
        edgeR::plotSmear(tested_local, de.tags = detags, las = 1, ylab = "log2FC")
        abline(h = c(-log2(FC_cutoff), log2(FC_cutoff)), col = "blue", cex = 2)
        legend("topright", border = FALSE, bty = "n",
               legend=c("p.adj > FDR_cutoff", "p.adj < FDR_cutoff", "log2FC_cutoff"),
               lty=c(0,0,1), pch=c(20, 20, -1), cex=0.8, lwd = c(1, 1, 2),
               col=c("Black", "Red", "blue"))
        dev.off()
    }

    glmlrt_groups_fix <- function(tested, Pairs){

        if (Pairs > 0) {
            comb_name <- names(tested)[Pairs]
            tested_local <- tested[[Pairs]]
        } else {
            comb_name <- "G2_over_G1"
            tested_local <- tested
        }

        tableDE <- edgeR::topTags(tested_local, n = nrow(tested_local$table),
                                  adjust.method="BH", sort.by="PValue")$table
        # tableDE <- cbind(tested_local$table, FDR = p.adjust(tested_local$table$PValue, "fdr"))
        colnames(tableDE)[c(1,2,5)] <- c("log2FC", "log2CPM", "FDR")
        tableDE[, "FC"] <- 2**tableDE[, "log2FC"]
        tableDE <- tableDE[order(rownames(tableDE)), ]

        if (tolower(dataBase) == "legacy") {
            GeneSymbol <- strsplit(row.names(tableDE), split = "\\|")
            GeneSymbol <- as.data.frame(GeneSymbol)
            GeneSymbol <- t(GeneSymbol)

            tableDE[, "GeneSymbol"] <- GeneSymbol[, 1]
            tableDE[, "GeneID"] <- GeneSymbol[, 2]
            tableDE <- tableDE[, c(7, 8, 4, 5, 6, 1, 3, 2)]
        } else {
            tableDE <- tableDE[, c(4, 5, 1, 2, 3)]
        }

        Results_Completed_local <- tableDE
        tableDE <- tableDE[tableDE$FDR < FDR_cutoff, ]
        tableDE <- tableDE[abs(tableDE$log2FC) > log2(FC_cutoff), ]
        write.csv(x = tableDE, file = paste0(DIR, "/ResultsDE_glmLRT_",
                                             comb_name, ".csv"),
                  row.names = FALSE)

        if (Pairs > 0) {
            string_vars[["envir_link"]]$resultadosDE.edgeR[[Pairs]] <- tableDE
            string_vars[["envir_link"]]$Results_Completed_edgeR[[Pairs]] <- Results_Completed_local
        } else {
            Results_Completed <- vector("list", 1)
            resultadosDE <- vector("list", 1)
            names(Results_Completed) <- comb_name
            names(resultadosDE) <- comb_name

            Results_Completed[[1]] <- Results_Completed_local
            resultadosDE[[1]] <- tableDE
            assign("Results_Completed_edgeR", Results_Completed, envir = get(envir_link))
            assign("resultadosDE.edgeR", resultadosDE, envir = get(envir_link))
        }

        #plots
        message("Plotting...")
        de <- edgeR::decideTestsDGE(tested_local)
        detags <- rownames(dge)[as.logical(de)]
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, paste0("glmLRT_",
                                     comb_name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, paste0("glmLRT_",
                                                 comb_name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        edgeR::plotSmear(tested_local, de.tags = detags, las = 1, ylab = "log2FC")
        abline(h = c(-log2(FC_cutoff), log2(FC_cutoff)), col = "blue", cex = 2)
        legend("topright", border = FALSE, bty = "n",
               legend=c("p.adj > FDR_cutoff", "p.adj < FDR_cutoff", "logFC_cutoff"),
               lty=c(0,0,1), pch=c(20, 20, -1), cex=0.8, lwd = c(1, 1, 2),
               col=c("Black", "Red", "blue"))
        dev.off()
    }

    volcano <- function(Results_Completed, Pairs){

        if (Pairs > 0) {
            comb_name <- names(Results_Completed)[Pairs]
            Results_Completed_local <- Results_Completed[[Pairs]]
        } else {
            comb_name <- "G2_over_G1"
            Results_Completed_local <- Results_Completed[[1]]
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
        log_10_FDR[log_10_FDR == "Inf"] <- (new_inf[new_inf != "Inf"][1] + 1)

        axislimits_y <- ceiling(max(log_10_FDR, na.rm = TRUE)) + 1

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
        text(x = axislimits_x - 0.4, y = axislimits_y/10, labels = length(Results_Completed_local$Colour[Results_Completed_local$log2FC >= log2(FC_cutoff) & Results_Completed_local$FDR <= FDR_cutoff]),
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
    dir.create(paste0(PATH, "/edgeR_Results.",
                      tolower(groupGen), "_", toupper(Name)), showWarnings = FALSE)
    DIR <- paste0(PATH, "/edgeR_Results.", tolower(groupGen), "_", toupper(Name))
    # dir.create(file.path(DIR, "PCA_Plots"), showWarnings = FALSE)

    # for the groups
    if (tolower(groupGen) == "mclust") {
        #groups from mclust
        Grupos_edgeR <- as.data.frame(string_vars[["envir_link"]]$GROUPS)
        Grupos_edgeR[, "Selected_classification"] <- factor(Grupos_edgeR[, "Selected_classification"])
        colnames(Grupos_edgeR) <- c("conditions", "type")
        Grupos_edgeR[, 2] <- c("paired-end")
    } else if (tolower(groupGen) == "clinical") {
        stop()

        Grupos_edgeR <- as.data.frame(string_vars[["envir_link"]]$clinical_groups_clinical[[clinical_pair]])
        Grupos_edgeR[, "Selected_classification"] <- factor(Grupos_edgeR[, "Selected_classification"])
        colnames(Grupos_edgeR) <- c("conditions", "type")
        Grupos_edgeR[, 2] <- c("paired-end")
    } else if (tolower(groupGen) == "coxhr") {
        Grupos_edgeR <- as.data.frame(string_vars[["envir_link"]]$clinical_groups)
        Grupos_edgeR$classification <- gsub("low", "1", Grupos_edgeR$classification)
        Grupos_edgeR$classification <- gsub("high", "2", Grupos_edgeR$classification)
        Grupos_edgeR[, "classification"] <- factor(Grupos_edgeR[, "classification"])
        colnames(Grupos_edgeR) <- c("type", "conditions")
        Grupos_edgeR[, 1] <- c("paired-end")
        Grupos_edgeR <- Grupos_edgeR[, c(2,1)]
    } else {
        stop(message("Please insert a valid 'groupGen' value!! ('mlcust', coxHR or 'clinical')"))
    }

    # check patient in common
    tmp <- rownames(Grupos_edgeR) %in% colnames(string_vars[["envir_link"]]$gene_tumor_not_normalized)

    Grupos_edgeR <- Grupos_edgeR[tmp, ]

    assign("condHeatmap", Grupos_edgeR[, 1], envir = get(envir_link))

    #selecting specifics patients
    completed_matrix <- string_vars[["envir_link"]]$gene_tumor_not_normalized[, rownames(Grupos_edgeR)]

    message("Filtering data ...")
    #remove residual data exclude rows with less than 10
    filter_rows <- rowSums(completed_matrix) >= 10
    completed_matrix <- completed_matrix[filter_rows, ]
    message("It was filtered = ", as.numeric(table(filter_rows)[1]), " genes!")

    # round numbers to whole counts
    completed_matrix <- round(completed_matrix, 0)

    # Cria o objeto DGEList
    dge <- edgeR::DGEList(counts = completed_matrix,
                      group = Grupos_edgeR[, 1])

    #Tdge$countso apply TMM normalization, it is convenient to create a DGEList
    #Robinson, M.D. and Oshlack, A. (2010). A scaling normalization method for
    #differential expres- sion analysis of RNA-seq data. Genome Biology 11, R25.
    #The calcNormFactors function normalizes for RNA composition by finding a
    #set of scaling factors for the library sizes that minimize the log-fold
    #changes between the samples for most genes. The default method for
    #computing these scale factors uses a trimmed mean of Mvalues (TMM) between
    #each pair of samples [26]. We call the product of the original library size
    #and the scaling factor the effective library size. The effective library
    #size replaces the original library size in all downsteam analyses. TMM
    #normalization is applied to this dataset to account for compositional
    #difference between the libraries.

    # binomial models are fitted and dispersion estimates are obtained, we can proceed with testing
    # procedures for determining differential expression using the exact test.
    # The GLM likelihood ratio test is based on the idea of fitting negative binomial GLMs
    # with the Cox-Reid dispersion estimates.
    design <- model.matrix(~Grupos_edgeR[, 1])
    colnames(design) <- paste0("G", seq(1, max(as.numeric(Grupos_edgeR[, 1]))))
    rownames(design) <- colnames(dge$counts)

    #PS: CPM = counts-per-million
    # The logCPM values can optionally be converted to RPKM or FPKM by subtracting log2 of gene length, see rpkm()
    if (tolower(Method) == "exacttest") {
        dge <- edgeR::calcNormFactors(dge)

        NormalizedExpression_TMM <- dge$counts

        # For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
        # likelihood (CR) method in estimating dispersions
        # To estimate common dispersion, trended dispersions and tagwise dispersions in one run:
        message("Estimating common dispersion...")
        dge <- edgeR::estimateCommonDisp(dge) #edgeR::estimateDisp(dge)
        group2_number <- max(as.numeric(levels(Grupos_edgeR[, 1])))
        # [1]  1  3  6 10 15 21 28 36 45 combinations
        if (group2_number > 2) {

            message("There are more than two group combinations, this may take a while...\n")

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

            assign("Results_Completed_edgeR", Results_Completed, envir = get(envir_link))
            assign("resultadosDE.edgeR", resultadosDE, envir = get(envir_link))

            count <- 0
            pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
            for (Pairs in seq(1, ncol(combinations))) {
                count <- count + 1
                tested[[Pairs]] <- edgeR::exactTest(dge, pair = combinations[, Pairs])
                exact_groups_fix(tested, Pairs)
                suppressWarnings(volcano(string_vars[["envir_link"]]$Results_Completed_edgeR, Pairs))

                setTxtProgressBar(pb, count)
            }
            close(pb)


        } else {
            tested <- edgeR::exactTest(dge)
            exact_groups_fix(tested, 0)
            suppressWarnings(volcano(string_vars[["envir_link"]]$Results_Completed_edgeR, 0))
        }
    } else if (tolower(Method) == "glmlrt") {
        dge <- edgeR::calcNormFactors(dge)

        NormalizedExpression_TMM <- dge$counts

        message("Estimating common dispersion...")
        # aDGEList <- DGEList(counts = TOC, group = tumorType)
        dge <- edgeR::estimateGLMCommonDisp(dge, design)
        aDGEList <- edgeR::estimateGLMTagwiseDisp(dge, design)
        aGlmFit <- edgeR::glmFit(aDGEList, design, dispersion = aDGEList$tagwise.dispersion,
                                 prior.count.total = 0)
        #To perform likelihood ratio tests:
        # fits the negative binomial GLM for each tag and produces an object of class DGEGLM with
        # some new components
        group2_number <- max(as.numeric(levels(Grupos_edgeR[, 1])))
        if (group2_number > 2) {
            message("There are more than two group combinations, this may take a while...\n")

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

            assign("Results_Completed_edgeR", Results_Completed, envir = get(envir_link))
            assign("resultadosDE.edgeR", resultadosDE, envir = get(envir_link))

            count <- 0
            pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
            for (Pairs in seq(1, ncol(combinations))) {
                if (combinations[ 1, Pairs] == 1) {
                    tested[[Pairs]] <- edgeR::glmLRT(aGlmFit, coef = Pairs + 1)
                } else {
                    count <- count + 1
                    vec <- numeric(ncol(combinations))
                    vec[combinations[ 1, Pairs]] <- -1
                    vec[combinations[ 2, Pairs]] <- 1

                    tested[[Pairs]] <- edgeR::glmLRT(aGlmFit, contrast = vec)
                }

                glmlrt_groups_fix(tested, Pairs)
                suppressWarnings(volcano(string_vars[["envir_link"]]$Results_Completed_edgeR, Pairs))

                setTxtProgressBar(pb, count)
            }
            close(pb)

        } else {
            tested <- edgeR::glmLRT(aGlmFit, coef = 2)
            glmlrt_groups_fix(tested, 0)
            suppressWarnings(volcano(string_vars[["envir_link"]]$Results_Completed_edgeR, 0))
        }
    } else {
        stop("Please insert a valid Method!!! ('exactTest' or 'glmLRT')")
    }

    #Data exploration
    # if (tolower(image_format) == "png") {
    #     png(filename = file.path(DIR, paste0(tolower(Method), "_MDS_PLOT.png")),
    #         width = Width, height = Height, res = Res, units = Unit)
    # } else if (tolower(image_format) == "svg") {
    #     svg(filename = file.path(DIR, paste0(tolower(Method), "_MDS_PLOT.svg")),
    #         width = Width, height = Height, onefile = TRUE)
    # } else {
    #     stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    # }
    # # png(file = paste0("./edgeR_Results.", tolower(dataType), "_", toupper(Name), "/MDS_PLOT.png"),
    # limma::plotMDS(dge, labels = 1:ncol(dge$samples), pch = 2, method = "logFC", las = 1)
    # # limma::plotMDS(dge)
    # dev.off()

    #Estimating the dispersion
    # dge <- edgeR::calcNormFactors(dge)
    # dge <- edgeR::estimateDisp(dge, design, robust=TRUE)
    # # dge$common.dispersion
    # if (tolower(image_format) == "png") {
    #     png(filename = file.path(DIR, paste0(tolower(Method), "_BCV_dispersion_PLOT.png")),
    #         width = Width, height = Height, res = Res, units = Unit)
    # } else if (tolower(image_format) == "svg") {
    #     svg(filename = file.path(DIR, paste0(tolower(Method), "_BCV_dispersion_PLOT.svg")),
    #         width = Width, height = Height, onefile = TRUE)
    # } else {
    #     stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    # }
    # plotBCV.modified(dge, las = 1)
    # dev.off()

    write.csv(x = NormalizedExpression_TMM, file = file.path(DIR, "TMM_Normalized_Expression.csv"))
    assign("NormalizedExpression.edgeR", NormalizedExpression_TMM, envir = get(envir_link))

    # assign("edgeR_results", tested, envir = get(envir_link))

    assign("Tool", "edgeR", envir = get(envir_link))


    message("Done!\n")
}
