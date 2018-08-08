GO_topGO <- function(condition,
                     FDR.cutoff = 0.05,
                     Width = 2000,
                     Height = 2000,
                     Res = 300,
                     Unit = "px",
                     image.format = "png",
                     Tool,
                     ID = "GeneID",
                     env) {

    # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

    # code ####
    message("Starting GO estimation...")

    if(missing(env)){stop(message("The 'env' argument is missing, please insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        PATH <- file.path(string_vars[["envir_link"]]$PATH, string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    if (missing(Tool)){Tool <- string_vars[["envir_link"]]$Tool}

    Name <- string_vars[["envir_link"]]$Name
    pairName <- string_vars[["envir_link"]]$pairName
    dataBase <- string_vars[["envir_link"]]$dataBase

    DIR <- paste0(PATH, "/topGo_Results_", Tool, "_", toupper(Name))
    dir.create(file.path(DIR, "GO_Output"), showWarnings = FALSE)
    dir.create(file.path(DIR, "GO_Output", "GraphOutput"),
               showWarnings = FALSE)
    dir.create(file.path(DIR, "GO_Output", "TextOutput"),
               showWarnings = FALSE)

    resultadosDE <- get(File, envir = string_vars[["envir_link"]])[[pairName]]


    # aaa <- read.csv('/home/fabio/Desktop/mart_export.txt', stringsAsFactors = FALSE)
    # aaa <- aaa[order(aaa$GO.domain, decreasing = TRUE), ]
    # aaa <- aaa[seq(1, (nrow(aaa) - table(aaa$GO.domain)[1])), ]
    # GO_annotations <- vector('list', 3)
    # names(GO_annotations) <- names(table(aaa$GO.domain))
    # GO_annotations[[1]] <- aaa[aaa$GO.domain == "biological_process", -3]
    # GO_annotations[[2]] <- aaa[aaa$GO.domain == "cellular_component", -3]
    # GO_annotations[[3]] <- aaa[aaa$GO.domain == "molecular_function", -3]


    GO_annotations <- GDCRtools::GO_annotations
    gene_annotations <- GDCRtools::annotation_table

    count <- 0
    for (Ontology in c("BP", "CC", "MF")) {
        count <- count + 1
        if (tolower(dataBase) == "legacy") {
            rownames(gene_annotations) <- gene_annotations$HUGO
            GO_annotations[[count]]$GeneID <- gene_annotations[GO_annotations[[count]]$HGNC.symbol, "GeneID"]
            GO_annotations[[count]] <- GO_annotations[[count]][!is.na(GO_annotations[[count]]$GeneID), ]
            genesOfInterest <- resultadosDE$GeneID
        } else if (tolower(dataBase) == "gdc") {
            rownames(gene_annotations) <- gene_annotations$ensembl
            GO_annotations[[count]]$GeneID <- gene_annotations[GO_annotations[[count]]$Gene.stable.ID, "GeneID"]
            GO_annotations[[count]] <- GO_annotations[[count]][!is.na(GO_annotations[[count]]$GeneID), ]
            genesOfInterest <- resultadosDE$GeneID
        }

        geneID2GO <- GO_annotations[[count]][, c(1, 5)]

        geneID2GO <- unstack(geneID2GO)

        unstack <- as.data.frame(matrix(nrow = length(geneID2GO), ncol = 2))
        unstack[, 1] <- names(geneID2GO)
        for (ids in seq(1, nrow(unstack))) {
            unstack[ids , 2] <- paste0(aaa[[ids]], collapse = ",")
        }

        geneID2GO <- unstack

        geneUniverse <- geneID2GO[, 1, drop = TRUE]
        geneList <- factor(as.character(as.integer(geneUniverse %in% genesOfInterest)), levels = c("0", "1"))
        names(geneList) <- geneUniverse
        sampleGOdata <- methods::new("topGOdata", description="My project",
                        ontology=Ontology, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

        resultFisher <- topGO::runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
        resultKS <- topGO::runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
        resultKS.elim <- topGO::runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
        resultTopgo <- topGO::runTest(sampleGOdata, algorithm="weight01", statistic="fisher")
        resultParentchild <- topGO::runTest(sampleGOdata, algorithm="parentchild", statistic="fisher")


        allRes <- topGO::GenTable(sampleGOdata, classicFisher = resultFisher,
                                  classicKS = resultKS, elimFisher = resultKS.elim,
                                  topgoFisher = resultTopgo,
                                  parentchildFisher = resultParentchild,
                                  orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = 10)

        topGO::printGraph(sampleGOdata, resultTopgo, firstSigNodes = numsignif,
                          fn.prefix = paste(output_file,"Topgo", sep="_"),
                          useInfo = "all", pdfSW = TRUE)

        # print out the genes that are annotated with the significantly enriched GO terms:
        myterms <- allRes$GO.ID
        mygenes <- topGO::genesInTerm(myGOdata, myterms)
        for (i in 1:length(myterms))
        {
            myterm <- myterms[i]
            mygenesforterm <- mygenes[myterm][[1]]
            myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
            mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
            mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
            print(paste("Term",myterm,"genes:",mygenesforterm2))
        }

        pValue.classic <- topGO::score(resultKS)
        pValue.elim <- topGO::score(resultKS.elim)[names(pValue.classic)]
        gstat <- topGO::termStat(sampleGOdata, names(pValue.classic))
        gSize <- gstat$Annotated / max(gstat$Annotated) * 4
        gCol <- colMap(gstat$Significant)
        plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
             pch = 19, cex = gSize, col = gCol)


        sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
        cbind(termStat(sampleGOdata, sel.go),
              elim = pValue.elim[sel.go],
              classic = pValue.classic[sel.go])


        topGO::showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
    }
}
