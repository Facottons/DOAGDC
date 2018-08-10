#' Separate patients in groups
#'
#' \code{groups_identification_mclust} is a function designed to separate patients in
#' groups, powered by mclust.
#'
#' @param dataType
#' @param group.number Numerical value indicating how many groups should be
#'   generated.
#' @param modelName A character string indicating which mclust modelName will be
#'   used. For more details please check \link{mclustModelNames} help file.
#' @param uncertaintyCutoff Numerical value indicating which uncertainty value
#'   for the separation should be tolerated. Patients over this threshold will
#'   be removed from the analysis.
#' @param reRunPlots Logical value where TRUE indicate that the function should
#'   run the step of group generation using the \code{uncertaintyCutoff}
#'   parameter for filtering the data. The default is \code{FALSE}.
#' @param nBreaks Numerical value giving the number of cells for the \code{hist}
#'   bars.
#' @param Width,Height,Res,Unit Graphical parameters. See \link{par} for more
#'   details. As default \code{nBreaks = 55, Width = 2000, Height = 1500, Res =
#'   300 and Unit = "px"}.
#' @param image.format A character string indicating which image format will be
#'   used. It could be "png" or "svg". The only unit available in "svg" is
#'   inches ('in'). The default is "png".
#' @param saveData
#' @param env
#' @param tumor
#' @param dataBase
#' @param workDir
#' @param Name
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#'
#' @export
#'
#' @import mclust
#'
#' @examples
#' \dontrun{
#' #separating isoform uc002peh.2 expression data patients in two groups
#' groups_identification_mclust("isoform", 2, Name = "uc002peh.2", modelName = "E")
#' }
groups_identification_mclust <- function(dataType,
                     group.number = "",
                     modelName = NULL,
                     uncertaintyCutoff = 0.05,
                     reRunPlots = FALSE,
                     nBreaks = 55,
                     Width = 2000,
                     Height = 1500,
                     Res = 300,
                     Unit = "px",
                     image.format = "png",
                     saveData = TRUE, env, tumor, dataBase,
                     workDir = "~/Desktop",
                     Name){

    # # based on criterea
    # set.seed(69)
    # x <- rnorm(500)
    # breaks <- pretty(x,10)
    # values <- hist(x, col = col, breaks = breaks)
    # col <- ifelse(values$mids <= -1, "red", "blue")
    # hist(x, col = col, breaks = breaks)

    # to.load <- c("mclust", "MineICA")


    #local functions ####
    #plotMix adapted
    PLOT <- function (mc, data, nbBreaks, traceDensity = TRUE, title = "",
                      xlim, ylim, ...)
    {
        if (length(data) <= 5000)
            shapiro.pval = shapiro.test(data)$p.value
        p <- seq(min(data), to = max(data), length = 1000)
        d <- mclust::cdens(modelName = mc$modelName, data = p, parameters = mc$parameters)
        title <- paste0(mc$G, "-component Gaussian Mixture Model")

        temp1 <- data[mc$classification == 1]

        par(mar = c(5,5,2,5))
        h_freq <- hist(as.numeric(data),
                       breaks = nbBreaks, freq = TRUE,
                       axes = FALSE, xlab = "", ylab = "", main = "")

        color <- vector(length = length(h_freq$mids))
        possible_colors <- RColorBrewer::brewer.pal(8, "Set2")
        last_delimiter <- 0
        classification_number <- length(unique(mc$classification))
        for (groups in seq(1, classification_number)){
            delimiter <- max(data[mc$classification == groups])
            if (groups == 1){
                color[h_freq$mids <= delimiter] <- possible_colors[groups]
            } else if(groups == max(unique(mc$classification))){
                size <- sum(h_freq$mids <= delimiter)
                # na_size <- sum(is.na(color[h_freq$mids <= delimiter]))
                color[(last_delimiter+1):length(h_freq$mids)] <- possible_colors[groups]
            } else {
                size <- sum(h_freq$mids <= delimiter)
                # na_size <- sum(is.na(color[h_freq$mids <= delimiter]))
                color[(last_delimiter+1):size] <- possible_colors[groups]
            }
            last_delimiter <- sum(h_freq$mids <= delimiter)
        }

        suppressWarnings(rug(as.numeric(data)))
        axis(4, las =1)
        mtext("Frequency", side=4, line=2.5, cex.lab=1,las=3)
        par(new=TRUE)

        h_dens <- hist(as.numeric(data),
                       breaks = nbBreaks, col = color, ylab = "",
                       freq = FALSE, las = 1, xlab = expression('Log'[2]*'(expression+1)'),
                       main = (paste(title, "\n", if (length(data) <= 5000){
                           paste("shapiro_pval = ", signif(shapiro.pval, digits = 3),
                                 sep = "")})))
        mtext("Density", side=2, line=2.5, cex.lab=1, las=3)

        if (traceDensity){
            par(new=TRUE)
            plot(MCLUST.RESULT.FINAL, what = "density", axes = FALSE, col = scales::alpha('red',.5),
                 ann = FALSE, main = "", lwd = 3.5)
        }
        # par(xpd = TRUE)
        # text(x = par("usr")[2]+1.8, y = mean(par("usr")[3:4]), "p logrank", srt = 270, cex=1.1)

        if (mc[["G"]] == 2){
            legend("topright", legend=c("Low", "High"), title = "Groups",
                   fill=unique(color), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
        } else {
            legend("topright", legend=paste0("G", seq(1, mc[["G"]])), title = "Groups",
                   fill=unique(color), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
        }
    }

    # Code ####

    dataType <- gsub(" ", "_", dataType)
    Name <- gsub("-", "_", Name)

    if (missing(env) && tumorData && !onlyFilter){
        assign(paste(toupper(tumor), toupper(dataBase),
                     gsub(" ", "_", dataType), "tumor_data", sep = "_"),
               new.env(parent=emptyenv()), envir = .GlobalEnv)
        envir_link <- paste(toupper(tumor), toupper(dataBase),
                            gsub(" ", "_", tolower(dataType)), "tumor_data", sep = "_")
        attr(envir_link, "name" ) = "Environment created by GDCRtools package, use its name in 'env' argument"
    } else if (missing(env) && !tumorData && !onlyFilter){
        assign(paste(toupper(tumor), toupper(dataBase),
                     gsub(" ", "_", tolower(dataType)), "both_data", sep = "_"),
               new.env(parent=emptyenv()), envir = .GlobalEnv)
        envir_link <- paste(toupper(tumor), toupper(dataBase),
                            gsub(" ", "_", tolower(dataType)), "both_data", sep = "_")
        attr(envir_link, "name" ) = "Environment created by GDCRtools package, use its name in 'env' argument"
    } else if (missing(env) && onlyFilter){
        message("Please, before using 'onlyFilter' argument, insert the Environment name.")
    } else {
        envir_link <- deparse(substitute(env))
    }

    string_vars <- list(envir_link = get(envir_link))


    if (missing("workDir")){
        workDir <- string_vars[["envir_link"]]$workDir
    }

    dir.create(path = file.path(workDir, "GDCRtools", toupper(tumor), "Analyses"),
               showWarnings = FALSE)

    PATH <- file.path(workDir, "GDCRtools", toupper(tumor), "Analyses")
    assign("PATH", PATH, envir = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))){
        dir.create(path = file.path(PATH, string_vars[["envir_link"]]$Name.e),
                   showWarnings = FALSE)
        dir.create(path = file.path(PATH, string_vars[["envir_link"]]$Name.e,
                                    paste0("Mclust_Results.", tolower(dataType),
                                           "_", toupper(Name))),
                   showWarnings = FALSE)
        DIR <- file.path(PATH, string_vars[["envir_link"]]$Name.e,
                         paste0("Mclust_Results.", tolower(dataType),
                                "_", toupper(Name)))
    } else {
        dir.create(path = file.path(PATH, paste0("Mclust_Results.",
                                                 tolower(dataType), "_",
                                                 toupper(Name))),
                   showWarnings = FALSE)
        DIR <- file.path(PATH, paste0("Mclust_Results.", tolower(dataType),
                                      "_", toupper(Name)))
    }

    if (tolower(dataType) == "methylation") {
        # COMPLETED.TABLE.NORMALIZED <- eval(parse(text= paste0('string_vars[["envir_link"]]$',
        #                                                       strsplit(tolower(dataType), "_")[[1]][1],
        #                                                       '_tumor_normalized')))

        FILTERED_RESULTS <- eval(parse(text= paste0("string_vars[['envir_link']]$",
                                                    "methylation_",
                                                    "tumor_filtered_selected_",
                                                    toupper(Name))))
    } else {
        # COMPLETED.TABLE.NORMALIZED <- eval(parse(text= paste0('string_vars[["envir_link"]]$',
        #                                                       strsplit(tolower(dataType), "_")[[1]][1],
        #                                                       '_tumor_normalized')))

        FILTERED_RESULTS <- eval(parse(text= paste0("string_vars[['envir_link']]$",
                                                    strsplit(tolower(dataType), "_")[[1]][1],
                                                    "_tumor_normalized_selected_",
                                                    toupper(Name))))
    }

    if (!reRunPlots){
        unique.col.row <- unique(rownames(FILTERED_RESULTS))
        # COMPLETED.TABLE.NORMALIZED <- COMPLETED.TABLE.NORMALIZED[, unique.col.row]
        FILTERED_RESULTS <- FILTERED_RESULTS[unique.col.row, , drop = FALSE]
        #starting
        message("Running Model-Based Clustering...")
        # MCLUST.RESULT <- mclust:::Mclust(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1))

        if (is.numeric(group.number) && is.character(modelName)) {
            MCLUST.RESULT.FINAL <-  mclust::Mclust(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1),
                                           G = group.number, modelNames = modelName)
        } else if (is.character(modelName)) {
            MCLUST.RESULT.FINAL <-  mclust::Mclust(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1),
                                           modelNames = modelName)
        } else if (is.numeric(group.number)){
            MCLUST.RESULT.FINAL <-  mclust::Mclust(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1),
                                           G = group.number)
        } else if (group.number == "" && is.null(modelName)){
            # group.number <- MCLUST.RESULT[["G"]]
            # modelName <- MCLUST.RESULT[["modelName"]]
            MCLUST.RESULT.FINAL <-  mclust::Mclust(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1))
                                           # G = group.number, modelNames = modelName)

        }
            # stop("Invalid arguments!! Please insert valid arguments.")
        ICL <- mclustICL(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1))
        BIC <- mclustBIC(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1))
        if (MCLUST.RESULT.FINAL$G > 1) {
            LRT <- mclustBootstrapLRT(log2(as.numeric(FILTERED_RESULTS[, 1]) + 1),
                                      modelName = MCLUST.RESULT.FINAL[["modelName"]])
        } else {
            message("Using ", Name, " does not generate two or more groups.",
                    " Please, try another 'Name'!")
            stop()
        }
        assign("MCLUST.RESULT.FINAL", MCLUST.RESULT.FINAL, envir = get(envir_link))

        # PLots
        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "mixture_log2_expression_%01d.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "mixture_log2_expression_%01d.svg"),
                width = Width, height = Height, onefile = FALSE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        par(mar = c(5, 6, 4, 3.5))
        PLOT(mc = MCLUST.RESULT.FINAL, data = log2(as.numeric(FILTERED_RESULTS[, 1]) + 1),
             nbBreaks = nBreaks, las = 1,
             traceDensity = TRUE,
             cex.main = 1.4,
             cex.lab = 1,
             cex.axis = 1.1)
        plot(BIC, las = 1)
        plot(ICL, las = 1)
        # plot(MCLUST.RESULT.FINAL, what = "BIC",
        #      xlab = expression('Log'[2]*'(expression+1)'), las = 1)
        plot(MCLUST.RESULT.FINAL, what = "classification",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1)
        plot(MCLUST.RESULT.FINAL, what = "density",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1, ylab = "Density")
        plot(MCLUST.RESULT.FINAL, what = "uncertain",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1)
        dev.off()

        #preparando output final
        GROUPS <- as.data.frame(matrix(ncol = 2, nrow = length(MCLUST.RESULT.FINAL$classification)))
        # if (tolower(dataType) == "gene")  {
        #     colnames(GROUPS) <- c("Selected_Gene_classification","Selected_Gene_uncertainty")
        # } else if (tolower(dataType) == "isoform") {
        #     colnames(GROUPS) <- c("Selected_Isoform_classification","Selected_Isoform_uncertainty")
        # } else {
        #     stop("Invalid arguments!! Please insert valid arguments.")
        # }
        colnames(GROUPS) <- c("Selected_classification","Selected_uncertainty")
        rownames(GROUPS) <- rownames(FILTERED_RESULTS)
        GROUPS[, 1] <- MCLUST.RESULT.FINAL$classification
        if (MCLUST.RESULT.FINAL[["G"]] == 2){
            GROUPS[, 1] <- sapply(GROUPS[, 1], FUN = function(x){ifelse(x == 1, "Low", "High")})
        }
        GROUPS[, 2] <- MCLUST.RESULT.FINAL$uncertainty
        write.csv(GROUPS, file = file.path(DIR, "groups_classifition.csv"),
                  quote = TRUE, row.names = TRUE)
        assign("GROUPS", GROUPS, envir = get(envir_link))

        #Uncertainty possible cutoffs1
        Uncertainty.range <- seq(0, 1, 0.0001)

        remaining.patients <- sapply(Uncertainty.range,
                                     function(pnow){length(which(GROUPS[, "Selected_uncertainty"] < pnow))})

        remaning.example.cutoff <- remaining.patients[which(Uncertainty.range == 0.05)]

        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "uncertainty_possible_cutoffs1.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "uncertainty_possible_cutoffs1.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        plot(Uncertainty.range, remaining.patients, main = "Uncertainty possible cutoffs",
             xlab="Uncertainty cutoff", ylab="Remaining Patients", lwd=2,
             axes=FALSE, cex.lab=1.1, col="green4", type = "l")

        axis(side=2, seq(min(remaining.patients), (max(remaining.patients)+30), 60),
             lwd=1, las = 1, cex.axis=1.2)
        axis(side=1, seq(0, 1, 0.05), lwd = 1, las = 1,
             cex.axis = 1.2)

        #example cutoff = 0.05
        abline(v=0.05, lwd=2)
        text(x=0.05+0.01, y=10,
             labels="0.05",
             cex=1.2, pos=4, srt=90)

        abline(v = Uncertainty.range[which(remaining.patients == max(remaining.patients))[1]], lwd = 2)
        text(x = 1.008*Uncertainty.range[which(remaining.patients == max(remaining.patients))[1]], y = 10,
             labels = paste0("Plateau cut = ",
                             Uncertainty.range[which(remaining.patients == max(remaining.patients))[1]]),
             cex = 1.2, pos = 4, srt = 90)

        abline(h = remaning.example.cutoff, lwd = 2, lty = 2)
        text(x = 0.9, y = remaning.example.cutoff,
             labels = paste0(remaning.example.cutoff," (",
                             scales::percent(remaning.example.cutoff/max(remaining.patients)), ")"),
             cex = 1.2, pos = 1, srt = 0)

        dev.off()

        #Uncertainty possible cutoffs2
        Uncertainty.range <- seq(0, 1, 0.01)

        remaining.patients <- sapply(Uncertainty.range,
                                     function(pnow){length(which(GROUPS[, "Selected_uncertainty"] < pnow))})

        to.plot <- data.frame(xaxis = Uncertainty.range[2:length(Uncertainty.range)],
                              yaxix = sapply(2:length(remaining.patients),
                                             function(n){remaining.patients[n]-remaining.patients[n-1]}))

        if (tolower(image.format) == "png") {
            png(filename = file.path(DIR, "uncertainty_possible_cutoffs2.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = file.path(DIR, "uncertainty_possible_cutoffs2.svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        plot(x = to.plot$xaxis,
             y = to.plot$yaxix,
             ylim = c(0,10), axes = FALSE, cex.lab = 1.1, type = "l",
             col = "blue2", lwd = 3, main = "Uncertainty possible cutoffs",
             xlab = "Uncertainty cutoff", ylab = "\U0394Patients")

        legend("topright", bty = "n", legend = "Cutoff = 0.05",
               lwd = 1, lty = 3, cex = 0.9)

        axis(side=2, lwd=1, las = 1, cex.axis=1.5)
        axis(side=1, lwd=1, las = 1, cex.axis=1.5)

        #example cutoff = 0.05
        abline(v=0.05, lwd=1, lty = 3)

        dev.off()

        assign("MCLUST.GROUP.NUMBER", group.number, envir = get(envir_link))
        sink(paste0(DIR, "/mclust_summary.txt"))
        print(summary(MCLUST.RESULT.FINAL))
        cat("\n\n\n")
        print(summary(BIC))
        cat("\n\n\n")
        print(summary(ICL))
        if (MCLUST.RESULT.FINAL$G > 1) {
            cat("\n\n\n")
            print(LRT)
        }
        sink()

    } else if (reRunPlots) {

        before <- nrow(string_vars[["envir_link"]]$GROUPS)
        GROUPS <- subset(x = string_vars[["envir_link"]]$GROUPS,
                         string_vars[["envir_link"]]$GROUPS[, 2] < uncertaintyCutoff)
        after <- nrow(GROUPS)
        message(paste0("It was removed ", before - after, " patients!"))

        FILTERED_RESULTS.filtered <- FILTERED_RESULTS[rownames(GROUPS), ]

        MCLUST.RESULT <- string_vars[["envir_link"]]$MCLUST.RESULT.FINAL

        group.number <- MCLUST.RESULT[["G"]]
        modelName <- MCLUST.RESULT[["modelName"]]
        MCLUST.RESULT.FINAL <-  mclust::Mclust(log2(FILTERED_RESULTS.filtered + 1),
                                       G = group.number, modelNames = modelName)

        if (tolower(image.format) == "png") {
            png(filename = paste0(DIR,
                              "/mixture_log2_expression_", uncertaintyCutoff, "_filtered_%01d.png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image.format) == "svg") {
            svg(filename = paste0(DIR,
                              "/mixture_log2_expression_", uncertaintyCutoff, "_filtered_%01d.svg"),
                width = Width, height = Height, onefile = FALSE)
        } else {
            stop(message("Please, Insert a valid image.format! ('png' or 'svg')"))
        }
        par(mar = c(5, 6, 4, 3.5))
        PLOT(mc = MCLUST.RESULT.FINAL, data = log2(FILTERED_RESULTS.filtered + 1),
             nbBreaks = nBreaks, las = 1,
             traceDensity = TRUE,
             cex.main = 1.4,
             cex.lab = 1,
             cex.axis = 1.1)
        plot(MCLUST.RESULT.FINAL, what = "BIC",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1)
        plot(MCLUST.RESULT.FINAL, what = "classification",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1)
        plot(MCLUST.RESULT.FINAL, what = "density",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1, ylab = "Density")
        plot(MCLUST.RESULT.FINAL, what = "uncertain",
             xlab = expression('Log'[2]*'(expression+1)'), las = 1)

        dev.off()

        GROUPS[, 1] <- MCLUST.RESULT.FINAL$classification
        GROUPS[, 2] <- MCLUST.RESULT.FINAL$uncertainty
        if (group.number == 2){
            GROUPS[, 1] <- sapply(GROUPS[, 1], FUN = function(x){ifelse(x == 1, "Low", "High")})
        }

        if (saveData){
            #saving
            write.csv(GROUPS, file = paste0(DIR, "/groups_classifition_",
                                                uncertaintyCutoff, "_filtered.csv"),
                      quote = TRUE, row.names = TRUE)

            # write.csv(COMPLETED.TABLE.NORMALIZED,
            #           file = paste0(DIR, "/mclust_Normalized_Expression.csv"))

            assign("GROUPS", GROUPS, envir = get(envir_link))
        }

        assign("MCLUST.RESULT.FINAL.filtered", MCLUST.RESULT.FINAL, envir = get(envir_link))
    }

    message("Done!")

}
