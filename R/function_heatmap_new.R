#' Draw a heatmap
#'
#' @param Tool
#' @param FC_cutoff
#' @param Name
#' @param Method The agglomeration method to be used: \code{"euclidean",
#'    "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson",
#'    "correlation", "abscorrelation", "spearman" or "kendall"}. The default is
#'    \code{"euclidean"}. More details about \code{Method} argument in
#'    \code{amap} package \link{Dist} \code{method} argument.
#' @param pairName A character string indicating the pair name to be used. When
#'    there are only two groups the default is \code{"G2_over_G1"}
#' @param RawValues A logical value. If \code{"TRUE"} the expression values are
#'    going to be converted to Z-Score before draw the heat map.
#' @param Width,Height,Res,Unit,image_format
#' @param env
#' @param ScaleMethod A character string indicating which method of scale should
#'    be used: \code{"none"}, \code{"row"}, \code{"column"}. The default is
#'    \code{"row"}. More details about \code{ScaleMethod} argument in
#'    \link{heatmap} \code{scale} argument.
#' @param outerMargins A numerical vector of the form c(bottom, left, top,
#'    right) giving the outer margins measured in lines of text. The default is
#'    no outer margins, i.e \code{"c(0,0,0,0)"}. Note: This argument is only
#'    used 'when \code{"labRow"} and \code{"labCol"} are set.
#' @param cexCol,cexRow A numerical value giving the amount by which
#'    \code{"labRow"} and \code{"labCol"} should be modified relative to the
#'    default size.
#' @param degree The \code{"labCol"} rotation in degrees. The default value is
#'    45 degrees.
#' @param labRow,labCol A logical value. If \code{"TRUE"} it is displayed row
#'    names (\code{"labRow"}) and col names (\code{"labCol"}). The default is
#'    \code{"FALSE"} for both, in order to kepp the plot clean.
#' @inheritParams groups_identification_mclust
#' @inheritParams dea_EBSeq
#' @inheritParams GOnto
#'
#' @return a heat map image.
#' @export
#'
#' @importFrom gplots colorpanel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom amap Dist
#'
#' @examples
#' \dontrun{
#' draw_heatmap("EBSeq", Name = "HIF3A", env = "env name without quotes")
#' }
draw_heatmap <- function(Tool, FC_cutoff = 2,
                        Name,
                        Method = "euclidean",
                        pairName = "G2_over_G1",
                        RawValues = FALSE,
                        Width = 6,
                        Height = 6,
                        Res = 300,
                        Unit = "in",
                        image_format = "svg",
                        env,
                        ScaleMethod = "row",
                        outerMargins = c(0,0,0,0),
                        cexCol= 0,
                        cexRow = 0,
                        degree = 45,
                        labRow = NULL,
                        labCol = NULL) {

    #verifying if the package is already installed

    #local function ####
    heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.1,
                        #cexRow = 0.2 + 1/log10(max(nr,2)),
                        #cexCol = 0.2 + 1/log10(max(nc,2)),
                        cexRow,
                        cexCol,
                        degree = 65,
                        scaleRangeMin,
                        scaleRangeMax,
                        cex.main = 1,
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        NumColSideColors = 1,
                        NumRowSideColors = 1,
                        KeyValueName,...){

        # created by obigriffith in
        # "https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/master/Analysis/DifferentialExpression/R/heatmap.3.R"
        # # pulled from here, and then tweaked slightly: http://www.biostars.org/p/18211/

        invalid <- function (x) {
            if (missing(x) || is.null(x) || length(x) == 0)
                return(TRUE)
            if (is.list(x))
                return(all(sapply(x, invalid)))
            else if (is.vector(x))
                return(all(is.na(x)))
            else return(FALSE)
        }

        x <- as.matrix(x)
        scale01 <- function(x, low = min(x), high = max(x)) {
            x <- (x - low)/(high - low)
            x
        }

        retval <- list()

        scale <- if (symm && missing(scale))
            "none"
        else match.arg(scale)

        dendrogram <- match.arg(dendrogram)

        trace <- match.arg(trace)

        density.info <- match.arg(density.info)

        if (length(col) == 1 && is.character(col))
            col <- get(col, mode = "function")

        if (!missing(breaks) && (scale != "none"))
            warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                    "specified can produce unpredictable results.",
                    "Please consider using only one or the other.")

        if (is.null(Rowv) || is.na(Rowv))
            Rowv <- FALSE

        if (is.null(Colv) || is.na(Colv))
            Colv <- FALSE
        else if (Colv == "Rowv" && !isTRUE(Rowv))
            Colv <- FALSE

        if (length(di <- dim(x)) != 2 || !is.numeric(x))
            stop("`x' must be a numeric matrix")

        nr <- di[1]
        nc <- di[2]

        if (nr <= 1 || nc <= 1)
            stop("`x' must have at least 2 rows and 2 columns")

        if (!is.numeric(margins) || length(margins) != 2)
            stop("`margins' must be a numeric vector of length 2")

        if (missing(cellnote))
            cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))

        if (!inherits(Rowv, "dendrogram")) {
            if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
                if (is.logical(Colv) && (Colv))
                    dendrogram <- "column"
                else dedrogram <- "none"

                warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                        dendrogram, "'. Omitting row dendogram.")
            }
        }

        if (!inherits(Colv, "dendrogram")) {
            if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
                if (is.logical(Rowv) && (Rowv))
                    dendrogram <- "row"
                else dendrogram <- "none"

                warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                        dendrogram, "'. Omitting column dendogram.")
            }
        }

        if (inherits(Rowv, "dendrogram")) {
            ddr <- Rowv
            rowInd <- order.dendrogram(ddr)
        }
        else if (is.integer(Rowv)) {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
            rowInd <- order.dendrogram(ddr)
            if (nr != length(rowInd))
                stop("row dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Rowv)) {
            Rowv <- rowMeans(x, na.rm = na.rm)
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
            rowInd <- order.dendrogram(ddr)
            if (nr != length(rowInd))
                stop("row dendrogram ordering gave index of wrong length")
        }
        else {
            rowInd <- nr:1
        }

        if (inherits(Colv, "dendrogram")) {
            ddc <- Colv
            colInd <- order.dendrogram(ddc)
        }
        else if (identical(Colv, "Rowv")) {
            if (nr != nc)
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            if (exists("ddr")) {
                ddc <- ddr
                colInd <- order.dendrogram(ddc)
            }
            else colInd <- rowInd
        }
        else if (is.integer(Colv)) {
            hcc <- hclustfun(distfun(if (symm)
                x
                else t(x)))
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, Colv)
            colInd <- order.dendrogram(ddc)
            if (nc != length(colInd))
                stop("column dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Colv)) {
            Colv <- colMeans(x, na.rm = na.rm)
            hcc <- hclustfun(distfun(if (symm)
                x
                else t(x)))
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, Colv)
            colInd <- order.dendrogram(ddc)
            if (nc != length(colInd))
                stop("column dendrogram ordering gave index of wrong length")
        }
        else {
            colInd <- 1:nc
        }

        retval$rowInd <- rowInd
        retval$colInd <- colInd
        retval$call <- match.call()

        x <- x[rowInd, colInd]  # rearrange matrix according to dendrograms
        x.unscaled <- x

        cellnote <- cellnote[rowInd, colInd]  # also rearrange the cellnotes

        # get labels
        if (is.null(labRow))
            labRow <- if (is.null(rownames(x)))
                (1:nr)[rowInd]
        else rownames(x)
        else labRow <- labRow[rowInd]
        if (is.null(labCol))
            labCol <- if (is.null(colnames(x)))
                (1:nc)[colInd]
        else colnames(x)
        else labCol <- labCol[colInd]

        ## do scaling of matrix according to Z-scores
        if (scale == "row") {
            retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
            x <- sweep(x, 1, rm)
            retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
            x <- sweep(x, 1, sx, "/")
        }
        else if (scale == "column") {
            retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
            x <- sweep(x, 2, rm)
            retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
            x <- sweep(x, 2, sx, "/")
        }

        # number of breaks
        if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
            if (missing(col) || is.function(col))
                breaks <- 16
            else breaks <- length(col) + 1
        }

        # set breakpoints
        if (length(breaks) == 1) {
            if (missing(scaleRangeMin))
                scaleRangeMin = min(x, na.rm=na.rm)

            if (missing(scaleRangeMax))
                scaleRangeMax = max(x, na.rm=na.rm)


            if (!symbreaks) {
                breaks <- seq(scaleRangeMin, scaleRangeMax, length=breaks);
            } else {
                #extreme <- max(abs(x), na.rm = TRUE)
                extreme = max(abs(c(scaleRangeMin,scaleRangeMax)), na.rm=na.rm)
                breaks <- seq(-extreme, extreme, length = breaks)
            }
        }

        nbr <- length(breaks)
        ncol <- length(breaks) - 1

        if (class(col) == "function")
            col <- col(ncol)

        min.breaks <- min(breaks)
        max.breaks <- max(breaks)

        # adjust for out-of-range given break settings
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks

        # layout height
        if (missing(lhei) || is.null(lhei))
            lhei <- c(keysize, 4)

        # layout width
        if (missing(lwid) || is.null(lwid))
            lwid <- c(keysize, 4)

        # define the layout
        if (missing(lmat) || is.null(lmat)) {
            lmat <- rbind(4:3, 2:1)

            if (!missing(ColSideColors)) {
                if (!is.character(ColSideColors) || ncol(ColSideColors) != nc)
                    stop("'ColSideColors' must be a matrix of ncol(x) ", nc, " columns")
                lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
                #lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
                side_height = min(side.height.fraction*nrow(ColSideColors), 1);
                lhei=c(lhei[1], side_height, lhei[2])
            }

            if (!missing(RowSideColors)) {
                if (!is.character(RowSideColors) || nrow(RowSideColors) != nr)
                    stop("'RowSideColors' must be a matrix of nrow(x) ", nr, " rows.  It currently has ", nrow(RowSideColors), " rows.")
                lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
                side_width = min(side.height.fraction*ncol(RowSideColors), 1);
                lwid <- c(lwid[1], side_width, lwid[2])
            }
            lmat[is.na(lmat)] <- 0
        }

        if (length(lhei) != nrow(lmat))
            stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
        if (length(lwid) != ncol(lmat))
            stop("lwid must have length = ncol(lmat) =", ncol(lmat))

        op <- par(no.readonly = TRUE)
        on.exit(par(op))

        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

        # Draw the colorbars for the annotations:
        if (!missing(RowSideColors)) {
            if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
            } else {
                par(mar = c(margins[1], 0, 0, 0.5))
                rsc = t(RowSideColors[rowInd, , drop=FALSE])
                rsc.colors = matrix()
                rsc.names = names(table(rsc))
                rsc.i = 1
                for (rsc.name in rsc.names) {
                    rsc.colors[rsc.i] = rsc.name
                    rsc[rsc == rsc.name] = rsc.i
                    rsc.i = rsc.i + 1
                }
                # print(rsc)
                rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
                image(1:nrow(rsc), 1:ncol(rsc), rsc,
                    col = as.vector(rsc.colors), axes = FALSE,
                    xlab="", ylab="")

                # add labels
                if (length(colnames(RowSideColors)) > 0) {
                    axis(1, 1:ncol(RowSideColors),
                                        labels=colnames(RowSideColors), las=2,
                                        tick=FALSE, xlab="", ylab="")
                }
            }
        }

        if (!missing(ColSideColors)) {

            if (!is.matrix(ColSideColors)){
                par(mar = c(0.5, 0, 0, margins[2]))
                image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
            } else {
                par(mar = c(0.5, 0, 0, margins[2]))
                csc = ColSideColors[, colInd, drop=FALSE]
                csc.colors = matrix()
                csc.names = names(table(csc))
                csc.i = 1
                for (csc.name in csc.names) {
                    csc.colors[csc.i] = csc.name
                    csc[csc == csc.name] = csc.i
                    csc.i = csc.i + 1
                }
                csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
                image(1:nrow(t(csc)), 1:ncol(t(csc)), t(csc),
                                    col = as.vector(csc.colors), axes = FALSE,
                                    xlab="", ylab="")

                # add labels
                if (length(rownames(ColSideColors)) > 0) {
                    axis(2, 1:(nrow(ColSideColors)),
                                    labels=rownames(ColSideColors), las = 2,
                                    tick = FALSE)
                }
            }
        }

        par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)
        cellnote <- t(cellnote)
        if (revC) {
            iy <- nr:1
            if (exists("ddr"))
                ddr <- rev(ddr)
            x <- x[, iy]
            cellnote <- cellnote[, iy]
        }
        else iy <- 1:nr

        # draw the central heatmap
        image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr),
                                            axes = FALSE, xlab = "", ylab = "",
                                            col = col, breaks = breaks, ...)

        # store the matrix drawn
        retval$carpet <- x

        # store the dendrograms
        if (exists("ddr"))
            retval$rowDendrogram <- ddr
        if (exists("ddc"))
            retval$colDendrogram <- ddc

        # store the breaks
        retval$breaks <- breaks

        # store the colormap used
        retval$col <- col

        # specially color in the na values
        if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
            mmat <- ifelse(is.na(x), 1, NA)
            image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
                                                    col = na.color, add = TRUE)
        }

        # X-axis column labels
        if (degree == 90){
            axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0)
            # , cex.axis = cexCol
        } else {
            tck <- axis(1, 1:nc, labels = FALSE, las = 2, line = -0.5,
                                                                    tick = 0)
            # ,cex.axis = cexCol
            labels <- labCol
            text(tck, par("usr")[3], labels=labels, srt=degree,
                xpd=NA, adj=c(1,0.8), cex=cexCol)
        }

        # X-axis title
        if (!is.null(xlab))
            mtext(xlab, side = 1, line = margins[1] - 1.25)

        # Y-axis row labeling
        axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0)
        # ,cex.axis = cexRow

        # Y-axis title
        if (!is.null(ylab))
            mtext(ylab, side = 4, line = margins[2] - 1.25)

        if (!missing(add.expr))
            eval(substitute(add.expr))
        if (!missing(colsep))
            for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
        if (!missing(rowsep))
            for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)


        min.scale <- min(breaks)
        max.scale <- max(breaks)
        x.scaled <- scale01(t(x), min.scale, max.scale)

        # column trace
        if (trace %in% c("both", "column")) {
            retval$vline <- vline
            vline.vals <- scale01(vline, min.scale, max.scale)
            for (i in colInd) {
                if (!is.null(vline)) {
                    abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
                }
                xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
                xv <- c(xv[1], xv)
                yv <- 1:length(xv) - 0.5
                lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
            }
        }

        # row trace
        if (trace %in% c("both", "row")) {
            retval$hline <- hline
            hline.vals <- scale01(hline, min.scale, max.scale)
            for (i in rowInd) {
                if (!is.null(hline)) {
                    abline(h = i + hline, col = linecol, lty = 2)
                }
                yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
                yv <- rev(c(yv[1], yv))
                xv <- length(yv):1 - 0.5
                lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
            }
        }

        # add cell labels
        if (!missing(cellnote))
            text(x = c(row(cellnote)), y = c(col(cellnote)),
                                        labels = c(cellnote), col = notecol,
                                        cex = notecex)

        # Plot the row dendrogram
        par(mar = c(margins[1], 0, 0, 0))
        if (dendrogram %in% c("both", "row")) {
            plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
        }
        else plot.new()

        # Plot the column dendrogram
        par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
        if (dendrogram %in% c("both", "column")) {
            plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
        }
        else plot.new()

        if (!is.null(main))
            title(main, cex.main=cex.main) #cex.main = 1.5 * op[["cex.main"]])


        # Add the Color Chart
        if (key) {
            par(mar = c(5, 4, 2, 1), cex = 0.75)
            tmpbreaks <- breaks
            if (symkey) {
                max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
                min.raw <- -max.raw
                tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
                tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
            }
            else {
                min.raw <- min(c(x,breaks), na.rm = TRUE)
                max.raw <- max(c(x,breaks), na.rm = TRUE)
            }

            z <- seq(min.raw, max.raw, length = length(col))
            image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                xaxt = "n", yaxt = "n")
            par(usr = c(0, 1, 0, 1))
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            axis(1, at = xv, labels = lv)
            if (scale == "row")
                mtext(side = 1, "Row Z-Score", line = 2, cex = 0.7)
            else if (scale == "column")
                mtext(side = 1, "Column Z-Score", line = 2, cex = 0.7)
            else mtext(side = 1, KeyValueName, line = 2, cex = 0.7)
            if (density.info == "density") {
                dens <- density(x, adjust = densadj, na.rm = TRUE)
                omit <- dens$x < min(breaks) | dens$x > max(breaks)
                dens$x <- dens$x[-omit]
                dens$y <- dens$y[-omit]
                dens$x <- scale01(dens$x, min.raw, max.raw)
                lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                    lwd = 1)
                axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
                title("Color Key\nand Density Plot")
                par(cex = 0.3)
                mtext(side = 2, "Density", line = 2)
            }
            else if (density.info == "histogram") {
                h <- hist(x, plot = FALSE, breaks = breaks)
                hx <- scale01(breaks, min.raw, max.raw)
                hy <- c(h$counts, h$counts[length(h$counts)])
                lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                    col = denscol)
                axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
                title("Color Key\nand Histogram")
                par(cex = 0.3)
                mtext(side = 2, "Count", line = 2)
            }
            else title("Color Key")
        }
        else plot.new()

        retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], high = retval$breaks[-1], color = retval$col)

        invisible(retval)
    }

    final.heatmap <- function(dist_Method, Scale_Method){

        # is the labRow desired?
        if (!is.null(labRow) && tolower(labRow) == "rownames") {
            labRow <- rownames(ParaHeatmaps)
        } else if(is.null(labRow)){
            labRow <-  ""
        }

        # is the labCol desired?
        if (!is.null(labCol) && tolower(labCol) == "colnames") {
            labCol <- colnames(ParaHeatmaps)
        } else if(is.null(labCol)){
            labCol <-  ""
        }

        if (tolower(image_format) == "png") {
            png(filename = paste0(DIR,
                                "/Heatmaps/", dist_Method, "_",
                                method, "_", Scale_Method,
                                "_FC_cutoff=", FC_cutoff, "_",
                                pairName, ".png"),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = paste0(DIR,
                                "/Heatmaps/", dist_Method, "_",
                                method, "_", Scale_Method,
                                "_FC_cutoff=", FC_cutoff, "_",
                                pairName, ".svg"),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format!",
                                                        " ('png' or 'svg')"))
        }
        par(oma = outerMargins)
        heatmap.3(DF,
                scale=Scale_Method,
                dendrogram="both",
                margins=c(3, 9),
                Rowv=row.dendro,
                keysize = 1.1,
                Colv=col.dendro,
                labRow = labRow,
                labCol = labCol,
                symbreaks=FALSE,
                key=TRUE, symkey=FALSE,
                ColSideColors = t(colour.groups),
                density.info="none",
                trace="none",
                cexCol= cexCol,
                cexRow = cexRow,
                degree = degree, col=RampaDeCor,
                KeyValueName = color_key_name)
        legend("topright", legend=paste0("G", levels(condHeatmap)),
                        title = "Groups",
                        fill=possible_colors[as.numeric(levels(condHeatmap))],
                        border=FALSE,
                        bty="n", y.intersp = 0.7, cex=0.7)
        dev.off()

        patient.order <- as.matrix(patient.order)
        patient.order <- cbind(1:nrow(patient.order), patient.order)
        colnames(patient.order) <- c("Dendrogram Order", "Patient")
        write.csv(patient.order, paste0(DIR, "/Heatmaps/Patient_order_",
                                        dist_Method, "_",
                                        method, "_", Scale_Method,
                                        "_FC_cutoff=", FC_cutoff, "_",
                                        pairName, ".csv"),
                row.names = FALSE)

        row.order <- as.matrix(row.order)
        row.order <- cbind(1:nrow(resultadosDE), row.order)
        colnames(row.order) <- c("Dendrogram Order", "Row")
        write.csv(row.order, paste0(DIR, "/Heatmaps/Row_order_",
                                    dist_Method, "_",
                                    method, "_", Scale_Method,
                                    "_FC_cutoff=", FC_cutoff, "_",
                                    pairName, ".csv"),
                row.names = FALSE)
    }

    #code ####
    Name <- gsub("-", "_", Name)

    if (missing(env)) {stop(message("The 'env' argument is missing, please",
                                    " insert the 'env' name and try again!"))}

    envir_link <- deparse(substitute(env))
    string_vars <- list(envir_link = get(envir_link))

    if (exists("Name.e", envir = get(envir_link))) {
        PATH <- file.path(string_vars[["envir_link"]]$PATH,
                                            string_vars[["envir_link"]]$Name.e)
    } else {
        PATH <- string_vars[["envir_link"]]$PATH
    }

    groupGen <- string_vars[["envir_link"]]$groupGen

    #creating the dir to outputs
    if (tolower(Tool) == "ebseq") {
        DIR <- paste0(PATH, "/EBSeq_Results.",
                    tolower(groupGen), "_", toupper(Name))
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.EBSeq[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.EBSeq
    } else if (tolower(Tool) == "edger") {
        DIR <- paste0(PATH, "/edgeR_Results.",
                    tolower(groupGen), "_", toupper(Name))
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.edgeR[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.edgeR
    } else if (tolower(Tool) == "deseq2") {
        DIR <- paste0(PATH, "/DESeq2_Results.",
                    tolower(groupGen), "_", toupper(Name))
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE.DESeq2[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.DESeq2
    } else if (tolower(Tool) == "crosstable.deseq2") {
        DIR <- paste0(PATH, "/CrossData_deseq2")
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE_crossed[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.DESeq2
    } else if (tolower(Tool) == "crosstable.edger") {
        DIR <- paste0(PATH, "/CrossData_edger")
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE_crossed[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.edgeR
    } else if (tolower(Tool) == "crosstable.ebseq") {
        DIR <- paste0(PATH, "/CrossData_ebseq")
        resultadosDE <- string_vars[["envir_link"]]$resultadosDE_crossed[[pairName]]
        NormalizedExpression <- string_vars[["envir_link"]]$NormalizedExpression.EBSeq
    } else {
        stop(message("Please, insert a valid Tool name!",
                                            " ('EBSeq', 'DESeq2' or 'edgeR')"))
    }

    dir.create(file.path(DIR, "Heatmaps"), showWarnings = FALSE)

    condHeatmap <- eval(parse(text= paste0("string_vars[['envir_link']]",
                                                            "$condHeatmap")))

    patients_stay <- unlist(strsplit(gsub("G", "", pairName), "_over_"))

    # keep the desired group pair
    NormalizedExpression <- NormalizedExpression[, condHeatmap %in% patients_stay]

    condHeatmap <- droplevels(condHeatmap[condHeatmap %in% patients_stay])

    # Create table just with DE and with FC_cutoff
    resultadosDEUP <- resultadosDE[resultadosDE$log2FC > log2(FC_cutoff), ]
    resultadosDEDOWN <- resultadosDE[resultadosDE$log2FC < log2(FC_cutoff), ]
    resultadosDE <- rbind(resultadosDEUP, resultadosDEDOWN)

    #Select only the DE results to use in heatmaps (next step) removing the EE
    ParaHeatmaps <- as.matrix(NormalizedExpression[match(rownames(resultadosDE),
                                                        rownames(NormalizedExpression)), ])
    colnames(ParaHeatmaps) <- 1:length(colnames(ParaHeatmaps))

    #Preparing color ramp
    RampaDeCor <- gplots::colorpanel(512,"blue","white","red")

    #Naming the heatmap cols
    possible_colors <- RColorBrewer::brewer.pal(8, "Set2")
    colour.groups <- possible_colors[condHeatmap]
    colour.groups <- as.matrix(colour.groups)
    rownames(colour.groups) <- 1:nrow(colour.groups)
    colnames(colour.groups) <- "Groups"


    if (tolower(ScaleMethod) == "none" && RawValues) {
        DF <- log2(ParaHeatmaps + 1)
        color_key_name <- "Log2(Expression Values + 1)"
    } else if (tolower(ScaleMethod) == "none") {
        DF <- scale(log2(ParaHeatmaps + 1))
        color_key_name <- "Z-Score"
    } else {
        DF <- log2(ParaHeatmaps + 1)
    }


    #dist(c,Scale_Method="euclidian")
    distMatrix.col <- amap::Dist(t(DF), method = Method)
    distMatrix.row <- amap::Dist(DF, method = Method)
    #complete
    dendro.euc.complete.col <- hclust(d = distMatrix.col, method = "complete")
    #plot(dendro.euc.complete.col)
    dendro.euc.complete.row <- hclust(d = distMatrix.row, method = "complete")
    #plot(dendro.euc.complete.row)
    row.dendro <- as.dendrogram(dendro.euc.complete.row)
    col.dendro <- as.dendrogram(dendro.euc.complete.col)
    patient.order <- colnames(NormalizedExpression)[dendro.euc.complete.col$order]
    row.order <- rownames(resultadosDE)[dendro.euc.complete.row$order]
    method <- "complete"
    if (length(Method) > 1){
        for (Methods in Method){
            final.heatmap(dist_Method = Methods, Scale_Method = ScaleMethod)
        }
    } else {
        final.heatmap(dist_Method = Method, Scale_Method = ScaleMethod)
    }

    #average
    dendro.euc.average.col <- hclust(d = distMatrix.col, method = "average")
    #plot(dendro.euc.average.col)
    dendro.euc.average.row <- hclust(d = distMatrix.row, method = "average")
    #plot(dendro.euc.average.row)
    row.dendro <- as.dendrogram(dendro.euc.average.row)
    col.dendro <- as.dendrogram(dendro.euc.average.col)
    patient.order <- colnames(NormalizedExpression)[dendro.euc.average.col$order]
    row.order <- rownames(resultadosDE)[dendro.euc.average.row$order]
    method <- "average"
    if (length(Method) > 1){
        for (Methods in Method){
            final.heatmap(dist_Method = Methods, Scale_Method = ScaleMethod)
        }
    } else {
        final.heatmap(dist_Method = Method, Scale_Method = ScaleMethod)
    }

    message("\nDone!\n")

}
