#' Co-Expression Analyses
#'
#' Powered by WGCNA
#'
#' @param data Used for external non-log expression data. Matrix or data frame.
#'    This \code{data} must have patients/sample code as \code{colnames} and
#'    genes as \code{rownames}.
#' @param name
#' @param data_base
#' @param work_dir
#' @param tumor
#' @param normalization
#' @param tumor_data
#' @param trait_data A character string where "default" indicates that the trait
#'    data to be used is provide by the output of \code{check_clinical_terms}
#'    function. Use "custom" for data inserted by the user after the
#'    \code{table2DOAGDC} function and with the table object named as
#'    \code{trait_data}. This object must have patients/sample code as
#'    \code{rownames} and the trait categories as \code{colnames}. By default
#'    trait data is not used.
#' @param network_type A character string indicating which network type should
#'    be used. WGCNA allows: "unsigned", "signed", and "signed hybrid". The
#'    default is "unsigned".
#' @param min_module_size Numerical value specifying the minimum cluster size.
#'    The default is 15.
#' @param max_softpower Numerical value indicating the maximum soft
#'    thresholding power for which the scale free topology fit indices are to be
#'    calculated. The default is 20.
#' @param nthreads Numerical value indicating how many threads to allow. The
#'    number of threads should not be more than the number of actual
#'    processors/cores. The default is 1.
#' @param me_diss_thres Numerical value specifying the maximum dendrogram cut
#'    height for module merging qualified by dissimilarity (i.e.,
#'    1-correlation). The default is 0.25.
#' @param width,height,res,unit,image_format
#' @param env
#' @param save_checkpoints Logical value where TRUE indicates that an external
#'    representation of analysis' objects will be saved to file in disk. The
#'    default is FALSE.
#' @param load_checkpoint Logical value where TRUE indicates that the saved
#'    checkpoint will be loaded and the analysis is going to continue from that
#'    point. The default is FALSE.
#' @param pearson_cutoff Numerical value specifying the minimum Pearson
#'    correlation value. The default is 0.5.
#' @inheritParams concatenate_exon
#' @inheritParams download_gdc
#' @inheritParams groups_identification_mclust
#'
#' @return a list of genes co-expressed and store it inside the determined
#'    environment name for further use.
#' @export
#'
#' @importFrom dynamicTreeCut printFlush
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom dplyr mutate
#'
co_expression <- function(data = NULL,
                        name,
                        data_base,
                        work_dir, tumor,
                        normalization = TRUE,
                        tumor_data = TRUE,
                        trait_data = NULL,
                        network_type = "unsigned",
                        min_module_size = 15,
                        max_softpower = 20,
                        nthreads = 1,
                        me_diss_thres = 0.25,
                        width = 2000,
                        height = 1500,
                        res = 300,
                        unit = "px",
                        image_format = "png",
                        env,
                        save_checkpoints = FALSE,
                        load_checkpoint = NULL,
                        pearson_cutoff = 0.5) {

    WGCNA::enableWGCNAThreads(nThreads = nthreads)

    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please insert",
            " the 'env' name and try again!"
        ))
    }

    ev <- deparse(substitute(env))
    sv <- list(ev = get(ev))

    work_dir <- ifelse(
        missing("work_dir"), sv[["ev"]]$work_dir, work_dir
    )

    assign("path", file.path(
        work_dir, "DOAGDC",
        toupper(sv[["ev"]]$tumor),
        "Analyses"
    ), envir = get(ev))

    path <- ifelse(exists("name_e", envir = get(ev)), file.path(
            sv[["ev"]]$path,
            sv[["ev"]]$name_e
        ), sv[["ev"]]$path)

    # creating the dir to outputs
    dir.create(paste0(
        path, "/co_expression_",
        tolower(network_type), "_", toupper(tumor)
    ),
    showWarnings = FALSE
    )
    dir <- paste0(
        path, "/co_expression_",
        tolower(network_type), "_", toupper(tumor)
    )

    message("Filtering and data normalization...")

    # checkpoint 1 ####
    if (is.null(load_checkpoint)) {
        if (is.null(data)) {
            if (normalization) {
                if (tumor_data) {
                    dat_expr <- sv[["ev"]]$gene_tumor_normalized
                } else {
                    dat_expr <- sv[["ev"]]$gene_not_tumor_normalized
                }
            } else {
                if (tumor_data) {
                    dat_expr <- sv[["ev"]]$gene_tumor_not_normalized
                } else {
                    dat_expr <- sv[["ev"]]$gene_not_tumor_not_normalized
                }
            }
        } else {
            dat_expr <- data
        }

        dat_expr <- log2(dat_expr + 1)
        dat_expr0 <- t(dat_expr)

        # Step 1: Filtering genes and seeking for outliers ####
        # lines with NA and zero-variance
        gsg <- WGCNA::goodSamplesGenes(dat_expr0, verbose = 3)

        if (!gsg$allOK) {
            if (sum(!gsg$goodGenes) > 0) {
                dynamicTreeCut::printFlush(paste(
                    "Removing genes:",
                    paste(names(dat_expr0)[!gsg$goodGenes],
                        collapse = ", "
                    )
                ))
            }
            if (sum(!gsg$goodSamples) > 0) {
                dynamicTreeCut::printFlush(paste(
                    "Removing samples:",
                    paste(rownames(dat_expr0)[!gsg$goodSamples],
                        collapse = ", "
                    )
                ))
            }
            # Remove the offending genes and samples from the data:
            dat_expr0 <- dat_expr0[gsg$goodSamples, gsg$goodGenes]
        }
        sample_tree <- hclust(dist(dat_expr0), method = "average")

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "sampleClustering.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "sampleClustering.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(0, 4, 2, 0), cex = 0.6)
        plot(sample_tree,
            main = "Sample clustering to detect outliers",
            sub = "", xlab = "", cex.lab = 1.5,
            cex.axis = 1.5, cex.main = 2, las = 1
        )
        dev.off()

        message(
            "Your plot was saved in ",
            file.path(dir, "sampleClustering."), image_format,
            ". Please, check this plot in order to insert the ",
            "cutheight value."
        )
        tmp <- "Insert the cutheight value: "
        cutheight <- as.numeric(readline(prompt = tmp))

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "sampleClustering2.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "sampleClustering2.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(0, 4, 2, 0), cex = 0.6)
        plot(sample_tree,
            main = "Sample clustering to detect outliers",
            sub = "", xlab = "", cex.lab = 1.5,
            cex.axis = 1.5, cex.main = 2, las = 1
        )
        abline(h = cutheight, col = "red")
        dev.off()

        if (tolower(cutheight) != "none" && is.numeric(cutheight)) {
            # Determine cluster under the line
            clust <- WGCNA::cutreeStatic(sample_tree,
                cutheight = cutheight,
                minSize = 10
            )
            # clust 1 contains the samples we want to keep.
            keep_samples <- (clust == 1)
            dat_expr <- dat_expr0[keep_samples, ]

            # clustering to detect outliers
            sample_tree <- hclust(dist(dat_expr), method = "average")
        } else {
            dat_expr <- dat_expr0
        }

        remove(dat_expr0)
        gc(verbose = FALSE)

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "sampleClustering_after_cut.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "sampleClustering_after_cut.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            tmp <- "Please, Insert a valid image_format! ('png' or 'svg')"
            stop(message(tmp))
        }
        par(mar = c(0, 4, 2, 0), cex = 0.6)
        plot(sample_tree,
            main = "Sample clustering to detect outliers",
            sub = "", xlab = "", cex.lab = 1.5,
            cex.axis = 1.5, cex.main = 2, las = 1
        )
        dev.off()

        # for trait data - part1 ####
        if (!is.null(trait_data)) {
            trait_data <- sv[["ev"]]$trait_data

            trait_rows <- match(rownames(dat_expr), rownames(trait_data))

            trait_data[, seq(1, ncol(trait_data))] <- sapply(
                colnames(trait_data),
                function(col) {
                    as.numeric(gsub(
                        pattern = "\\[Not Available\\]",
                        replacement = "NA", x = trait_data[, col],
                        perl = TRUE
                    ))
                }
            )

            trait_data[, seq(1, ncol(trait_data))] <- sapply(
                colnames(trait_data),
                function(col) {
                    as.numeric(gsub(
                        pattern = "\\[Not Applicable\\]",
                        replacement = "NA", x = trait_data[, col],
                        perl = TRUE
                    ))
                }
            )

            trait_data[, seq(1, ncol(trait_data))] <- sapply(
                colnames(trait_data),
                function(col) {
                    as.numeric(trait_data[, col])
                }
            )

            dat_traits <- as.matrix(trait_data[trait_rows, ])
            sample_tree2 <- hclust(dist(dat_expr), method = "average")

            trait_colors <- WGCNA::numbers2colors(dat_traits,
                signed = ifelse(
                    tolower(network_type) == "signed", TRUE, FALSE
                ),
                naColor = "grey"
            )

            WGCNA::plotDendroAndcolors(sample_tree2, trait_colors,
                groupLabels = names(dat_traits),
                main = paste0(
                    "Sample dendrogram ",
                    "and trait heatmap"
                )
            )
        }

        # Step 2: Network construction and module detection ####
        powers <- c(c(1:10), seq(from = 12, to = max_softpower, by = 2))

        sft <- WGCNA::pickSoftThreshold(dat_expr,
            powerVector = powers,
            network_type = network_type,
            verbose = 5
        )

        # Plot the results:
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "connectivity.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "connectivity.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 2))
        cex1 <- 0.9
        # Scale-free topology fit index as a function of the soft-thresholding
        # power
        tmpx <- sft$fitIndices[, 1]
        tmpy <- tmpy
        plot(tmpx, tmpy,
            xlab = "Soft Threshold (power)",
            ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
            main = paste("Scale independence")
        )
        text(tmpx, tmpy, labels = powers, cex = cex1, col = "red")

        # this line corresponds to using an R^2 cut-off of h
        abline(h = 0.9, col = "red")

        # Mean connectivity as a function of the soft-thresholding power
        plot(tmpx, sft$fitIndices[, 5],
            xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
            type = "n",
            main = paste("Mean connectivity")
        )
        text(tmpx, sft$fitIndices[, 5],
            labels = powers,
            cex = cex1, col = "red"
        )

        dev.off()

        # We choose the power XX, which is the lowest power for which the
        # scale-free topology fit index reaches 0.90. nGenes <- ncol(dat_expr)
        n_samples <- nrow(dat_expr)

        soft_power <- as.numeric(readline(prompt = paste0(
            "Insert the soft ",
            "threshold value: "
        )))

        message("Calculating adjacency...")

        adjacency <- WGCNA::adjacency(dat_expr,
            power = soft_power,
            type = network_type
        )
        k <- WGCNA::softConnectivity(dat_expr,
            power = soft_power,
            type = network_type
        )

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "softConnectivity.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "softConnectivity.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mfrow = c(1, 2))
        hist(k)
        WGCNA::scaleFreePlot(k, main = "Check scale free topology\n")
        dev.off()

        # Step 3 - TOPOLOGICAL OVERLAP MATRIX ####
        message(
            "starting 'TOPOLOGICAL OVERLAP MATRIX'(TOM).\n",
            "It takes several minutes..."
        )
        tom <- WGCNA::TOMsimilarity(adjacency, TOMType = network_type)
        if (save_checkpoints) {
            message("Saving chechpoint 1, please wait...")
            save.image(file.path(dir, "checkpoint_TOMsimilarity.RData"))
        }
    } else {
        message("Loading chechpoint 1, please wait...")
        load(file.path(dir, "checkpoint_TOMsimilarity.RData"))
    }

    # continue from checkpoint 1 area####
    remove(adjacency)
    gc(verbose = FALSE)
    diss_tom <- 1 - tom

    # Call the hierarchical clustering function TO SEE THE TOM
    gene_tree <- hclust(as.dist(diss_tom), method = "average")
    if (tolower(image_format) == "png") {
        png(
            filename = file.path(dir, "gene_tree1.png"),
            width = width, height = height, res = res, units = unit
        )
    } else if (tolower(image_format) == "svg") {
        svg(
            filename = file.path(dir, "gene_tree1.svg"),
            width = width, height = height, onefile = TRUE
        )
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    # Plot the resulting clustering tree (dendrogram)
    # (12,9)
    plot(gene_tree,
        xlab = "", sub = "",
        main = "Gene clustering on TOM-based dissimilarity",
        labels = FALSE, hang = 0.04
    )
    dev.off()

    gc(verbose = FALSE)

    # NOW IT IS NECESSARY TO PERFORM BRANCH CUTTING
    # We like large modules, so we set the minimum module size relatively high:
    # Module identification using dynamic tree cut:
    dynamic_mods <- dynamicTreeCut::cutreeDynamic(
        dendro = gene_tree,
        distM = diss_tom,
        deepSplit = 2, pamRespectsDendro = FALSE,
        minClusterSize = min_module_size
    )
    remove(diss_tom)
    gc(verbose = FALSE)

    # LISTS THE SIZE OF THE MODULES
    # Convert numeric lables into colors
    dynamic_colors <- WGCNA::labels2colors(dynamic_mods)
    gc(verbose = FALSE)
    # Plot the dendrogram and colors underneath
    if (tolower(image_format) == "png") {
        png(
            filename = file.path(dir, "dendroandcolors.png"),
            width = width, height = height, res = res, units = unit
        )
    } else if (tolower(image_format) == "svg") {
        svg(
            filename = file.path(dir, "dendroandcolors.svg"),
            width = width, height = height, onefile = TRUE
        )
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    WGCNA::plotDendroAndcolors(gene_tree, dynamic_colors, "Dynamic Tree Cut",
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05,
        main = "Gene dendrogram and module colors"
    )
    dev.off()

    # Step 4 -Calculate eigengenes ####
    me_list <- WGCNA::moduleEigengenes(dat_expr, colors = dynamic_colors)
    mes <- me_list$eigengenes

    # Calculate dissimilarity of module eigengenes
    me_diss <- 1 - cor(mes, use = "pairwise.complete.obs")
    me_diss["MEgrey" == rownames(me_diss), ] <- 0
    me_tree <- hclust(dist(me_diss), method = "average")

    # Error NA/NaN/Inf in foreign function call
    # Cluster module eigengenes

    # Plot the result
    if (tolower(image_format) == "png") {
        png(
            filename = file.path(dir, "MEtree.png"),
            width = width, height = height, res = res, units = unit
        )
    } else if (tolower(image_format) == "svg") {
        svg(
            filename = file.path(dir, "MEtree.svg"),
            width = width, height = height, onefile = TRUE
        )
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    plot(me_tree,
        main = "Clustering of module eigengenes",
        xlab = "", sub = ""
    )

    # Plot the cut line into the dendrogram
    abline(h = me_diss_thres, col = "red")
    # Call an automatic merging function
    merge <- WGCNA::mergeCloseModules(dat_expr, dynamic_colors,
        cutheight = me_diss_thres, verbose = 3
    )
    # The merged module colors
    merged_colors <- merge$colors
    # Eigengenes of the new merged modules:
    merged_mes <- merge$newMEs
    dev.off()

    if (tolower(image_format) == "png") {
        png(
            filename = file.path(dir, paste0(
                "genedendro_",
                me_diss_thres, ".png"
            )),
            width = width, height = height, res = res, units = unit
        )
    } else if (tolower(image_format) == "svg") {
        svg(
            filename = file.path(dir, paste0(
                "genedendro_",
                me_diss_thres, ".svg"
            )),
            width = width, height = height, onefile = TRUE
        )
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    WGCNA::plotDendroAndcolors(gene_tree, cbind(dynamic_colors, merged_colors),
        c("Dynamic Tree Cut", "Merged dynamic"),
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05
    )
    dev.off()


    module_colors <- merged_colors
    # numerical labels matching color and module
    color_order <- c("grey", WGCNA::standardcolors(50))

    mes <- merged_mes

    # for trait data - part2 ####
    if (!is.null(trait_data)) {

        # Step 5 - blockwise ####
        bwnet <- WGCNA::blockwiseModules(dat_expr,
            maxBlockSize = 5000,
            power = soft_power, TOMType = network_type,
            min_module_size = min_module_size,
            reassignThreshold = 0,
            mergeCutheight = me_diss_thres,
            numericLabels = TRUE,
            saveTOMs = FALSE, # or TRUE
            saveTOMFileBase = file.path(
                dir,
                "TOM-blockwise"
            ),
            verbose = 3
        )

        module_labels <- match(module_colors, color_order) - 1

        # blockwise label again module
        bw_labels <- WGCNA::matchLabels(bwnet$colors, module_labels)
        # convert labels in colors
        bw_module_colors <- WGCNA::labels2colors(bw_labels)

        # Plot dendogram and module1's colors for 1 block
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(
                    dir,
                    "Gene_dendrogram_and_module_colors.png"
                ),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(
                    dir,
                    "Gene_dendrogram_and_module_colors.svg"
                ),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        WGCNA::plotDendroAndcolors(bwnet$dendrograms[[1]],
            bw_module_colors[bwnet$blockGenes[[1]]],
            "Module colors",
            main = "Gene dendrogram and module colors in block 1",
            dendroLabels = FALSE, hang = 0.03,
            addGuide = TRUE, guideHang = 0.05
        )
        dev.off()

        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "Single_block_gene_dendrogram.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "Single_block_gene_dendrogram.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        WGCNA::plotDendroAndcolors(gene_tree,
            cbind(module_colors, bw_module_colors),
            c("Single block", "2 blocks"),
            main = paste0(
                "Single block gene dendrogram and",
                " module colors"
            ),
            dendroLabels = FALSE, hang = 0.03,
            addGuide = TRUE, guideHang = 0.05
        )
        dev.off()

        single_block_mes <- WGCNA::moduleEigengenes(
            dat_expr,
            module_colors
        )$eigengenes

        blockwise_mes <- WGCNA::moduleEigengenes(
            dat_expr,
            bw_module_colors
        )$eigengenes

        single2blockwise <- match(names(single_block_mes),
                                        names(blockwise_mes))
        signif(diag(cor(blockwise_mes[, single2blockwise],
                                                single_block_mes)), 3)

        gc(verbose = FALSE)

        # Recalculate MEs with labels colors
        mes0 <- single_block_mes
        mes <- WGCNA::orderMEs(mes0)
        module_trait_cor <- cor(mes, dat_traits, use = "p")
        module_trait_pvalue <- WGCNA::corPvalueStudent(module_trait_cor,
                                                                    n_samples)

        # Show cor and their p-values
        text_matrix <- paste(signif(module_trait_cor, 2), "\n(",
            signif(module_trait_pvalue, 1), ")",
            sep = ""
        )
        dim(text_matrix) <- dim(module_trait_cor)
        if (tolower(image_format) == "png") {
            png(
                filename = file.path(dir, "Module_trait_relationships.png"),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = file.path(dir, "Module_trait_relationships.svg"),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message("Insert a valid image_format! ('png' or 'svg')"))
        }
        # Show cor in heatmap plot
        WGCNA::labeledHeatmap(
            Matrix = module_trait_cor,
            xLabels = names(dat_traits),
            yLabels = names(mes),
            ySymbols = names(mes),
            colorLabels = FALSE,
            colors = WGCNA::blueWhiteRed(50),
            text_matrix = text_matrix,
            setStdMargins = FALSE,
            cex.text = 0.15,
            zlim = c(-1, 1),
            cex.lab.y = 0.3,
            cex.lab.x = 0.8,
            cex.lab = 0.5, xLabelsAngle = 25,
            naColor = "black",
            main = paste("Module-trait relationships")
        )
        dev.off()
    }

    # Step 6 - exporting network ####
    # Select module probes
    probes <- colnames(dat_expr)
    all_module <- as.data.frame(cbind(probes, merged_colors))
    write.table(all_module,
        file = file.path(dir, "gene_modules.tsv"),
        col.names = c("geneid", "module"),
        row.names = FALSE, sep = "\t", quote = FALSE
    )

    file_name_all <- paste("LocusIDs_Allcolors", ".txt", sep = "")
    write.table(all_module,
        file = paste(dir, file_name_all, sep = "/"),
        col.names = FALSE, quote = FALSE, row.names = FALSE
    )

    module <- all_module$merged_colors[grep(
        ifelse(tolower(data_base) == "gdc", toupper(name),
            ifelse(
                suppressWarnings(is.na(as.numeric(name))), paste0(
                    toupper(name),
                    "\\|"
                ), paste0(
                    "\\|",
                    toupper(name)
                )
            )
        ),
        all_module$probes
    )]

    in_module <- is.finite(match(module_colors, module))
    mod_probes <- probes[in_module]
    # Select the corresponding Topological Overlap
    mod_tom <- tom[in_module, in_module]
    remove(tom, gsg)
    gc(verbose = FALSE)
    dimnames(mod_tom) <- list(mod_probes, mod_probes)

    threshold <- ifelse(
        tolower(network_type) == "unsigned",
        abs((0.5 * pearson_cutoff)),
        abs(0.5 + (0.5 * pearson_cutoff))
    )

    # export to cytoscape
    WGCNA::exportNetworkToCytoscape(mod_tom,
        edgeFile = paste0(
            dir,
            "/CytoscapeInput_edges_",
            paste(module,
                collapse = "_"
            ),
            "_with_threshold_",
            threshold**soft_power,
            ".txt"
        ),
        nodeFile = paste0(
            dir,
            "/CytoscapeInput_nodes_",
            paste(module,
                collapse = "_"
            ),
            "_with_threshold_",
            threshold**soft_power,
            ".txt"
        ),
        weighted = TRUE,
        threshold = threshold**soft_power,
        nodeNames = mod_probes,
        nodeAttr = module_colors[in_module]
    )


    # Nodes
    file_nodes <- read.table(paste0(
        dir, "/CytoscapeInput_nodes_",
        paste(module, collapse = "_"),
        "_with_threshold_",
        threshold**soft_power, ".txt"
    ),
    header = TRUE, stringsAsFactors = FALSE
    )

    selected_from <- grep(ifelse(
        suppressWarnings(is.na(as.numeric(name))),
        paste0(toupper(name), "\\|"),
        paste0("\\|", toupper(name))
    ), file_nodes$fromNode)

    selected_to <- grep(ifelse(
        suppressWarnings(is.na(as.numeric(name))),
        paste0(toupper(name), "\\|"),
        paste0("\\|", toupper(name))
    ), file_nodes$toNode)

    file_nodes <- file_nodes[c(selected_from, selected_to), ]

    write.table(file_nodes, paste0(
        dir, "/", toupper(name), "_nodes_",
        paste(module, collapse = "_"),
        "_with_threshold_",
        threshold**soft_power, ".txt"
    ),
    quote = FALSE, row.names = FALSE
    )

    # Edges
    file_edges <- read.table(paste0(
        dir, "/CytoscapeInput_edges_",
        paste(module, collapse = "_"),
        "_with_threshold_",
        threshold**soft_power, ".txt"
    ),
    header = TRUE, stringsAsFactors = FALSE
    )

    file_edges <- file_edges[, 1:3]

    possible_cor <- seq(from = -1, to = 1, by = 1e-3)
    threshold <- abs(0.5 + (0.5 * possible_cor))**soft_power
    threshold <- round(threshold, 3)
    names(threshold) <- possible_cor
    pearson_values <- sapply(as.numeric(file_edges$weight), function(x) {
        unique(round(
            as.numeric(names(threshold)[threshold == round(x, 3)]),
            3
        ))[1]
    })

    file_edges <- dplyr::mutate(file_edges, weightPearson = pearson_values)

    write.table(file_edges, paste0(
        dir, "/CytoscapeInput_edges_",
        paste(module, collapse = "_"),
        "_with_threshold_",
        threshold**soft_power, ".txt"
    ),
    quote = FALSE, row.names = FALSE
    )

    selected_from <- grep(ifelse(
        suppressWarnings(is.na(as.numeric(name))),
        paste0(toupper(name), "\\|"),
        paste0("\\|", toupper(name))
    ), file_edges$fromNode)

    selected_to <- grep(ifelse(
        suppressWarnings(is.na(as.numeric(name))),
        paste0(toupper(name), "\\|"),
        paste0("\\|", toupper(name))
    ), file_edges$toNode)

    file_edges <- file_edges[c(selected_from, selected_to), ]

    write.table(file_edges, paste0(
        dir, "/", toupper(name), "_edges_",
        paste(module, collapse = "_"),
        "_with_threshold_",
        threshold**soft_power, ".txt"
    ),
    quote = FALSE, row.names = FALSE
    )


    # put together node and edge
    coexp_genes <- unique(unlist(file_edges[, 1:2]))


    # export gene list  XXXX
    assign("coexp_genes", coexp_genes, envir = get(ev))


    gc(verbose = FALSE)
    message("\nDone!")
}
