#' Co-Expression Analyses
#'
#' Powered by WGCNA
#'
#' @param Data Used for external non-log expression data. Matrix or data frame.
#'   This \code{Data} must have patients/sample code as \code{colnames} and
#'   genes as \code{rownames}.
#' @param Name
#' @param dataBase
#' @param workDir
#' @param tumor
#' @param normalization
#' @param tumorData
#' @param traitData A character string where "default" indicates that the trait
#'   data to be used is provide by the output of \code{check_clinical_terms}
#'   function. Use "custom" for data inserted by the user after the
#'   \code{table2DOAGDC} function and with the table object named as
#'   \code{trait_data}. This object must have patients/sample code as
#'   \code{rownames} and the trait categories as \code{colnames}. By default
#'   trait data is not used.
#' @param networkType A character string indicating which network type should be
#'   used. WGCNA allows: "unsigned", "signed", and "signed hybrid". The default
#'   is "unsigned".
#' @param minModuleSize Numerical value specifying the minimum cluster size. The
#'   default is 15.
#' @param max_softpower Numerical value indicating the maximum soft thresholding
#'   power for which the scale free topology fit indices are to be calculated.
#'   The default is 20.
#' @param nthreads Numerical value indicating how many threads to allow. The
#'   number of threads should not be more than the number of actual
#'   processors/cores. The default is 1.
#' @param MEDissThres Numerical value specifying the maximum dendrogram cut
#'   height for module merging qualified by dissimilarity (i.e., 1-correlation).
#'   The default is 0.25.
#' @param Width,Height,Res,Unit,image_format
#' @param env
#' @param saveCheckpoints Logical value where TRUE indicates that an external
#'   representation of analysis' objects will be saved to file in disk. The
#'   default is FALSE.
#' @param loadCheckpoint Logical value where TRUE indicates that the saved
#'   checkpoint will be loaded and the analysis is going to continue from that
#'   point. The default is FALSE.
#' @param pearsonCutoff Numerical value specifying the minimum Pearson
#'   correlation value. The default is 0.5.
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#' @inheritParams GOnto
#'
#' @return a list of genes co-expressed and store it inside the determined
#'   environment name for further use.
#'
#'
#' @examples
#' \dontrun{
#' co_expression(Data, workDir, tumor, env = "env name without quotes")
#' }
co_expression <- function(Data = NULL,
                          Name,
                          dataBase,
                          workDir, tumor,
                          normalization = TRUE,
                          tumorData = TRUE,
                          traitData = NULL,
                          networkType = "unsigned",
                          minModuleSize = 15,
                          max_softpower = 20,
                          nthreads = 1,
                          MEDissThres = 0.25,
                          Width = 2000,
                          Height = 1500,
                          Res = 300,
                          Unit = "px",
                          image_format = "png",
                          env,
                          saveCheckpoints = FALSE,
                          loadCheckpoint = NULL,
                          pearsonCutoff = 0.5){

    # library(WGCNA, km2gcn, sleuth)

    WGCNA::enableWGCNAThreads(nThreads = nthreads)

    if (missing(env)) {stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))

    string_vars <- list(envir_link = get(envir_link))
    if (missing("workDir")) {
        workDir <- string_vars[["envir_link"]]$workDir
    }

    assign("PATH", file.path(workDir, "DOAGDC", toupper(string_vars[["envir_link"]]$tumor), "Analyses"), envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))) {
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    #creating the dir to outputs
    dir.create(paste0(PATH, "/co_expression_",
                      tolower(networkType), "_", toupper(tumor)), showWarnings = FALSE)
    DIR <- paste0(PATH, "/co_expression_",
                  tolower(networkType), "_", toupper(tumor))

    message("Filtering and data normalization...")

    # checkpoint 1 ####
    if (is.null(loadCheckpoint)) {
        if (is.null(Data)) {
            if (normalization) {
                if (tumorData) {
                    datExpr <- string_vars[["envir_link"]]$gene_tumor_normalized
                } else {
                    datExpr <- string_vars[["envir_link"]]$gene_not_tumor_normalized
                }
            } else {
                if (tumorData) {
                    datExpr <- string_vars[["envir_link"]]$gene_tumor_not_normalized
                } else {
                    datExpr <- string_vars[["envir_link"]]$gene_not_tumor_not_normalized
                }
            }
        } else {
            datExpr <- Data
        }

        datExpr <- log2(datExpr + 1)
        datExpr0 <- t(datExpr)

        # Step 1: Filtering genes and seeking for outliers ####
        # lines with NA and zero-variance
        gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)

        if (!gsg$allOK)
        {
            if (sum(!gsg$goodGenes) > 0)
                dynamicTreeCut::printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
            if (sum(!gsg$goodSamples) > 0)
                dynamicTreeCut::printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
            # Remove the offending genes and samples from the data:
            datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
        }
        sampleTree <- hclust(dist(datExpr0), method = "average")

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "sampleClustering.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "sampleClustering.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(0,4,2,0), cex = 0.6)
        plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
             cex.axis = 1.5, cex.main = 2, las = 1)
        dev.off()

        message("Your plot was saved in ",
                file.path(DIR, "sampleClustering."), image_format,
                ". Please, check this plot in order to insert the cutHeight value.")
        cutHeight <- as.numeric(readline(prompt = "Please, insert the cutHeight value: "))

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "sampleClustering2.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "sampleClustering2.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(0,4,2,0), cex = 0.6)
        plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
             cex.axis = 1.5, cex.main = 2, las = 1)
        abline(h = cutHeight, col = "red")
        dev.off()

        if (tolower(cutHeight) != "none" && is.numeric(cutHeight)) {
            # Determine cluster under the line
            clust <- WGCNA::cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
            # clust 1 contains the samples we want to keep.
            keepSamples <- (clust == 1)
            datExpr <- datExpr0[keepSamples, ]

            # clustering to detect outliers
            sampleTree <- hclust(dist(datExpr), method = "average")
        } else {
            datExpr <- datExpr0
        }

        remove(datExpr0)
        gc(verbose = FALSE)

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "sampleClustering_after_cut.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "sampleClustering_after_cut.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(0,4,2,0), cex = 0.6)
        plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
             cex.axis = 1.5, cex.main = 2, las = 1)
        dev.off()

        # for trait data - part1 ####
        if (!is.null(traitData)) {
            if (tolower(traitData) == "default") {
                traitData <- string_vars[["envir_link"]]$trait_data
            } else if (tolower(traitData) == "custom") {
                traitData <- string_vars[["envir_link"]]$trait_data
            }

            traitRows = match(rownames(datExpr), rownames(traitData))
            # traitData$days_to_death <- gsub(pattern = "\\[Not Applicable\\]", replacement = "NA", x = traitData$days_to_death, perl = TRUE)

            traitData[, 1:ncol(traitData)] <- sapply(colnames(traitData), function(col){
                as.numeric(gsub(pattern = "\\[Not Available\\]", replacement = "NA", x = traitData[, col], perl = TRUE))
            })

            traitData[, 1:ncol(traitData)] <- sapply(colnames(traitData), function(col){
                as.numeric(gsub(pattern = "\\[Not Applicable\\]", replacement = "NA", x = traitData[, col], perl = TRUE))
            })

            traitData[, 1:ncol(traitData)] <- sapply(colnames(traitData), function(col){
                as.numeric(traitData[, col])
            })

            datTraits <- as.matrix(traitData[traitRows, ])
            sampleTree2 <- hclust(dist(datExpr), method = "average")
            # white means low, red high e grey NA
            if (tolower(networkType) == "signed") {
                traitColors <- WGCNA::numbers2colors(datTraits, signed = TRUE, naColor = "grey")
            } else {
                traitColors <- WGCNA::numbers2colors(datTraits, signed = FALSE, naColor = "grey")
            }
            # par(mar = c(0,20,0,0))
            WGCNA::plotDendroAndColors(sampleTree2, traitColors,
                                       groupLabels = names(datTraits),
                                       main = "Sample dendrogram and trait heatmap")
        }

        # Step 2: Network construction and module detection ####
        powers = c(c(1:10), seq(from = 12, to = max_softpower, by = 2))

        sft = WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = networkType, verbose = 5)

        # Plot the results:
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "connectivity.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "connectivity.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1,2))
        cex1 = 0.9
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
             main = paste("Scale independence"))
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels = powers, cex = cex1, col = "red")

        # this line corresponds to using an R^2 cut-off of h
        abline(h = 0.9, col = "red")

        # Mean connectivity as a function of the soft-thresholding power

        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type = "n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1,col = "red")

        dev.off()

        # We choose the power XX, which is the lowest power for which the scale-free topology fit index reaches 0.90.
        # nGenes <- ncol(datExpr)
        nSamples <- nrow(datExpr)

        softPower <- as.numeric(readline(prompt = "Please, insert the soft threshold value: "))

        message("Calculating adjacency...")

        adjacency <- WGCNA::adjacency(datExpr, power = softPower, type = networkType)
        k <- WGCNA::softConnectivity(datExpr, power = softPower, type = networkType)

        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "softConnectivity.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "softConnectivity.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1,2))
        hist(k)
        WGCNA::scaleFreePlot(k, main = "Check scale free topology\n")
        dev.off()

        # Step 3 - TOPOLOGICAL OVERLAP MATRIX ####
        message("starting 'TOPOLOGICAL OVERLAP MATRIX'(TOM).\nIt takes several minutes...")
        TOM <- WGCNA::TOMsimilarity(adjacency, TOMType = networkType)
        if (saveCheckpoints) {
            message("Saving chechpoint 1, please wait...")
            save.image(file.path(DIR,"checkpoint_TOMsimilarity.RData"))
        }
    } else {
        message("Loading chechpoint 1, please wait...")
        load(file.path(DIR,"checkpoint_TOMsimilarity.RData"))
    }

    # continue from checkpoint 1 area####
    remove(adjacency)
    gc(verbose = FALSE)
    dissTOM <- 1 - TOM

    # Call the hierarchical clustering function TO SEE THE TOM
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "gene_treee1.png"),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "gene_treee1.svg"),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    # Plot the resulting clustering tree (dendrogram)
    #(12,9)
    plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
         labels = FALSE, hang = 0.04)
    dev.off()

    gc(verbose = FALSE)

    # NOW IT IS NECESSARY TO PERFORM BRANCH CUTTING
    # We like large modules, so we set the minimum module size relatively high:
    # Module identification using dynamic tree cut:
    dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
    remove(dissTOM)
    gc(verbose = FALSE)

    # LISTS THE SIZE OF THE MODULES
    # Convert numeric lables into colors
    dynamicColors <- WGCNA::labels2colors(dynamicMods)
    #table(dynamicColors)
    gc(verbose = FALSE)
    # Plot the dendrogram and colors underneath
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "dendroandcolors.png"),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "dendroandcolors.svg"),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()

    # Step 4 -Calculate eigengenes ####
    MEList <- WGCNA::moduleEigengenes(datExpr, colors = dynamicColors)
    MEs <- MEList$eigengenes

    # Calculate dissimilarity of module eigengenes
    MEDiss = 1 - cor(MEs, use = 'pairwise.complete.obs')
    MEDiss["MEgrey" == rownames(MEDiss), ] <- 0
    METree = hclust(dist(MEDiss), method = "average")

    # Error NA/NaN/Inf in foreign function call
    # Cluster module eigengenes
    # MEDiss <- 1-cor(MEs)
    # METree <- hclust(as.dist(MEDiss), method = "average")

    # Plot the result
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "MEtree.png"),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "MEtree.svg"),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    plot(METree, main = "Clustering of module eigengenes",
         xlab = "", sub = "")

    # Plot the cut line into the dendrogram
    abline(h = MEDissThres, col = "red")
    # Call an automatic merging function
    merge <- WGCNA::mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors <- merge$colors
    # Eigengenes of the new merged modules:
    mergedMEs <- merge$newMEs
    dev.off()

    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, paste0("genedendro_", MEDissThres, ".png")),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, paste0("genedendro_", MEDissThres, ".svg")),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()


    moduleColors <- mergedColors
    # numerical labels matching color and module
    colorOrder <- c("grey", WGCNA::standardColors(50))

    MEs <- mergedMEs

    # Step 5 - blockwise ####
    # bwnet <- WGCNA::blockwiseModules(datExpr, maxBlockSize = 5000,
    #                          power = softPower, TOMType = networkType, minModuleSize = minModuleSize,
    #                          reassignThreshold = 0, mergeCutHeight = MEDissThres,
    #                          numericLabels = TRUE,
    #                          saveTOMs = FALSE, #or TRUE
    #                          saveTOMFileBase = file.path(DIR,"TOM-blockwise"),
    #                          verbose = 3)
    #
    # moduleLabels <- match(moduleColors, colorOrder)-1
    #
    # # blockwise label again module
    # bwLabels <- WGCNA::matchLabels(bwnet$colors, moduleLabels)
    # # convert labels in colors
    # bwModuleColors <- WGCNA::labels2colors(bwLabels)
    #
    # # Plot dendogram and module1's colors for 1 block
    # if (tolower(image_format) == "png") {
    #     png(filename = file.path(DIR, "Gene_dendrogram_and_module_colors.png"),
    #         width = Width, height = Height, res = Res, units = Unit)
    # } else if (tolower(image_format) == "svg") {
    #     svg(filename = file.path(DIR, "Gene_dendrogram_and_module_colors.svg"),
    #         width = Width, height = Height, onefile = TRUE)
    # } else {
    #     stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    # }
    # WGCNA::plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
    #                     "Module colors", main = "Gene dendrogram and module colors in block 1",
    #                     dendroLabels = FALSE, hang = 0.03,
    #                     addGuide = TRUE, guideHang = 0.05)
    # dev.off()
    #
    # if (tolower(image_format) == "png") {
    #     png(filename = file.path(DIR, "Single_block_gene_dendrogram.png"),
    #         width = Width, height = Height, res = Res, units = Unit)
    # } else if (tolower(image_format) == "svg") {
    #     svg(filename = file.path(DIR, "Single_block_gene_dendrogram.svg"),
    #         width = Width, height = Height, onefile = TRUE)
    # } else {
    #     stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    # }
    # WGCNA::plotDendroAndColors(geneTree,
    #                     cbind(moduleColors, bwModuleColors),
    #                     c("Single block", "2 blocks"),
    #                     main = "Single block gene dendrogram and module colors",
    #                     dendroLabels = FALSE, hang = 0.03,
    #                     addGuide = TRUE, guideHang = 0.05)
    # dev.off()

    # for trait data - part2 ####
    if (!is.null(traitData)) {
        singleBlockMEs <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
        blockwiseMEs <- WGCNA::moduleEigengenes(datExpr, bwModuleColors)$eigengenes

        single2blockwise <- match(names(singleBlockMEs), names(blockwiseMEs))
        signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)

        gc(verbose = FALSE)

        # Recalculate MEs with labels colors
        MEs0 <- singleBlockMEs
        MEs <- orderMEs(MEs0)
        moduleTraitCor <- cor(MEs, datTraits, use = "p")
        moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)

        # Show cor and their p-values
        textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "")
        dim(textMatrix) <- dim(moduleTraitCor)
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "Module_trait_relationships.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "Module_trait_relationships.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        # par(mar = c(8, 4.5, 1, 1))
        # Show cor in heatmap plot
        WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                       xLabels = names(datTraits),
                       yLabels = names(MEs),
                       ySymbols = names(MEs),
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = 0.15,
                       zlim = c(-1,1),
                       cex.lab.y = 0.3,
                       cex.lab.x = 0.8,
                       cex.lab = 0.5, xLabelsAngle = 25,
                       naColor = "black",
                       main = paste("Module-trait relationships"))
        dev.off()
    }

    # genes <- colnames(datExpr)

    # Step 6 - exporting network ####
    # Select module probes
    probes <- colnames(datExpr)
    all_module <- as.data.frame(cbind(probes, mergedColors))
    write.table(all_module, file = file.path(DIR, 'gene_modules.tsv'), col.names = c("gene_id", "module"),
                row.names = FALSE, sep = "\t", quote = FALSE)

    fileNameAll <- paste("LocusIDs_Allcolors",".txt", sep = "")
    write.table(all_module, file = paste(DIR, fileNameAll, sep = "/"),
                col.names = FALSE, quote = FALSE, row.names = FALSE)


    if (tolower(dataBase) == "gdc") {
        module <- all_module$mergedColors[grep(toupper(Name), all_module$probes)]
    } else if (tolower(dataBase) == "legacy") {
        if (suppressWarnings(is.na(as.numeric(Name)))) {
            # Symbol before |
            module <- all_module$mergedColors[grep(paste0(toupper(Name), "\\|"), all_module$probes)]
        } else {
            # Entrez
            module <- all_module$mergedColors[grep(paste0("\\|", toupper(Name)), all_module$probes)]
        }
    }


    inModule <- is.finite(match(moduleColors, module)) #(moduleColors == module)
    modProbes <- probes[inModule]
    # Select the corresponding Topological Overlap
    modTOM <- TOM[inModule, inModule]
    remove(TOM, gsg)
    gc(verbose = FALSE)
    dimnames(modTOM) <- list(modProbes, modProbes)

    if (tolower(networkType) == "unsigned") {
        Threshold <- abs((0.5*pearsonCutoff))
    } else {
        Threshold <- abs(0.5 + (0.5*pearsonCutoff))
    }

    # export to cytoscape
    WGCNA::exportNetworkToCytoscape(modTOM,
                                   edgeFile = paste0(DIR, "/CytoscapeInput_edges_",
                                                     paste(module, collapse = "_"), "_with_threshold_",
                                                     Threshold**softPower, ".txt"),
                                   nodeFile = paste0(DIR, "/CytoscapeInput_nodes_",
                                                     paste(module, collapse = "_"), "_with_threshold_",
                                                     Threshold**softPower, ".txt"),
                                   weighted = TRUE,
                                   threshold = Threshold**softPower,
                                   nodeNames = modProbes,
                                   nodeAttr = moduleColors[inModule])


    # Nodes
    fileNodes <- read.table(paste0(DIR, "/CytoscapeInput_nodes_",
                                   paste(module, collapse = "_"), "_with_threshold_",
                                   Threshold**softPower, ".txt"),
                            header = TRUE, stringsAsFactors = FALSE)


    if (suppressWarnings(is.na(as.numeric(Name)))) {
        # Symbol before |
        selected_from <- grep(paste0(toupper(Name), "\\|"), fileNodes$fromNode)
        selected_to <- grep(paste0(toupper(Name), "\\|"), fileNodes$toNode)
    } else {
        # Entrez
        selected_from <- grep(paste0("\\|", toupper(Name)), fileNodes$fromNode)
        selected_to <- grep(paste0("\\|", toupper(Name)), fileNodes$toNode)
    }

    fileNodes <- fileNodes[c(selected_from, selected_to), ]


    write.table(fileNodes, paste0(DIR, "/", toupper(Name), "_nodes_",
                                  paste(module, collapse = "_"), "_with_threshold_",
                                  Threshold**softPower, ".txt"),
                quote = FALSE, row.names = FALSE)



    # Edges
    fileEdges <- read.table(paste0(DIR, "/CytoscapeInput_edges_",
                                   paste(module, collapse = "_"), "_with_threshold_",
                                   Threshold**softPower, ".txt"),
                            header = TRUE, stringsAsFactors = FALSE)

    fileEdges <- fileEdges[, 1:3]

    possible_cor <- seq(from = -1, to = 1, by = 1e-3)
    threshold = abs(0.5 + (0.5*possible_cor))**softPower
    threshold <- round(threshold, 3)
    names(threshold) <- possible_cor
    pearson_values <- sapply(as.numeric(fileEdges$weight), function(x){
        unique(round(as.numeric(names(threshold)[threshold == round(x, 3)]), 3))[1]
    })

    fileEdges <- dplyr::mutate(fileEdges, weightPearson = pearson_values)

    # fileEdges[, 1] <- unname(sapply(fileEdges[, 1], function(w){
    #     paste0(unlist(strsplit(w, "\\|"))[2])}))
    #
    # fileEdges[, 2] <- unname(sapply(fileEdges[, 2], function(w){
    #     paste0(unlist(strsplit(w, "\\|"))[2])}))

    write.table(fileEdges, paste0(DIR, "/CytoscapeInput_edges_",
                                  paste(module, collapse = "_"), "_with_threshold_",
                                  Threshold**softPower, ".txt"),
                quote = FALSE, row.names = FALSE)


    if (suppressWarnings(is.na(as.numeric(Name)))) {
        # Symbol before |
        selected_from <- grep(paste0(toupper(Name), "\\|"), fileEdges$fromNode)
        selected_to <- grep(paste0(toupper(Name), "\\|"), fileEdges$toNode)
    } else {
        # Entrez
        selected_from <- grep(paste0("\\|", toupper(Name)), fileEdges$fromNode)
        selected_to <- grep(paste0("\\|", toupper(Name)), fileEdges$toNode)
    }

    fileEdges <- fileEdges[c(selected_from, selected_to), ]


    write.table(fileEdges, paste0(DIR, "/", toupper(Name), "_edges_",
                                  paste(module, collapse = "_"), "_with_threshold_",
                                  Threshold**softPower, ".txt"),
                quote = FALSE, row.names = FALSE)


    # put together node and edge
    coexp_genes <- unique(unlist(fileEdges[, 1:2]))


    # export gene list    CHUNTA TUTO XXXX
    assign("coexp_genes", coexp_genes, envir = get(envir_link))


    gc(verbose = FALSE)
    message("\nDone!")
}
