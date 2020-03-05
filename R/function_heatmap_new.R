#' Draw a heatmap
#'
#' @param tool
#' @param fc_cutoff
#' @param name
#' @param method The agglomeration method to be used: \code{"euclidean",
#'    "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson",
#'    "correlation", "abscorrelation", "spearman" or "kendall"}. The default is
#'    \code{"euclidean"}. More details about \code{method} argument in
#'    \code{amap} package \link{Dist} \code{method} argument.
#' @param pair_name A character string indicating the pair name to be used. When
#'    there are only two groups the default is \code{"G2_over_G1"}
#' @param raw_values A logical value. If \code{"TRUE"} the expression values are
#'    going to be converted to Z-Score before draw the heat map.
#' @param width,height,res,unit,image_format
#' @param env
#' @param scale_method A character string indicating which method of scale
#'    should
#'    be used: \code{"none"}, \code{"row"}, \code{"column"}. The default is
#'    \code{"row"}. More details about \code{scale_method} argument in
#'    \link{heatmap} \code{scale} argument.
#' @param outer_margins A numerical vector of the form c(bottom, left, top,
#'    right) giving the outer margins measured in lines of text. The default is
#'    no outer margins, i.e \code{"c(0,0,0,0)"}. Note: This argument is only
#'    used 'when \code{"lab_row"} and \code{"lab_col"} are set.
#' @param cex_col,cex_row A numerical value giving the amount by which
#'    \code{"lab_row"} and \code{"lab_col"} should be modified relative to the
#'    default size.
#' @param degree The \code{"lab_col"} rotation in degrees. The default value is
#'    45 degrees.
#' @param lab_row,lab_col A logical value. If \code{"TRUE"} it is displayed row
#'    names (\code{"lab_row"}) and col names (\code{"lab_col"}). The default is
#'    \code{"FALSE"} for both, in order to kepp the plot clean.
#' @inheritParams dea_ebseq
#' @inheritParams concatenate_exon
#' @inheritParams groups_identification_mclust
#'
#' @return a heat map image.
#' @export
#'
#' @importFrom gplots colorpanel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom amap Dist
#'
#' @examples
#' library(DOAGDC)
#'
#' # data already downloaded using the 'download_gdc' function
#' concatenate_expression("gene",
#'    name = "HIF3A",
#'    data_base = "legacy",
#'    tumor = "CHOL",
#'    work_dir = "~/Desktop"
#' )
#'
#' # separating gene HIF3A expression data patients in two groups
#' groups_identification_mclust("gene", 2,
#'    name = "HIF3A",
#'    modelName = "E",
#'    env = CHOL_LEGACY_gene_tumor_data,
#'    tumor = "CHOL"
#' )
#'
#' # load not normalized data
#' concatenate_expression("gene",
#'    normalization = FALSE,
#'    name = "HIF3A",
#'    data_base = "legacy",
#'    tumor = "CHOL",
#'    env = CHOL_LEGACY_gene_tumor_data,
#'    work_dir = "~/Desktop"
#' )
#'
#' # start DE analysis
#' # considering concatenate_expression and groups_identification already runned
#' dea_edger(
#'    name = "HIF3A",
#'    group_gen = "mclust",
#'    env = CHOL_LEGACY_gene_tumor_data
#' )
#'
#' draw_heatmap("edgeR", name = "HIF3A", env = CHOL_LEGACY_gene_tumor_data)
draw_heatmap <- function(tool, fc_cutoff = 2,
                        name,
                        method = "euclidean",
                        pair_name = "G2_over_G1",
                        raw_values = FALSE,
                        width = 6,
                        height = 6,
                        res = 300,
                        unit = "in",
                        image_format = "svg",
                        env,
                        scale_method = "row",
                        outer_margins = c(0, 0, 0, 0),
                        cex_col = 0,
                        cex_row = 0,
                        degree = 45,
                        lab_row = NULL,
                        lab_col = NULL) {

    # verifying if the package is already installed

    # local function ####
    heatmap_3 <- function(x,
                    rowv = TRUE, colv = if (symm) "Rowv" else TRUE,
                    distfun = dist,
                    hclustfun = hclust,
                    dendrogram = c("both", "row", "column", "none"),
                    symm = FALSE,
                    scale = c("none", "row", "column"),
                    na.rm = TRUE,
                    rev_c = identical(colv, "Rowv"),
                    add_expr,
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
                    na_color = par("bg"),
                    trace = c("none", "column", "row", "both"),
                    tracecol = "cyan",
                    hline = median(breaks),
                    vline = median(breaks),
                    linecol = tracecol,
                    margins = c(5, 5),
                    col_sidecolors,
                    row_sidecolors,
                    side_height_fraction = 0.1,
                    cex_row,
                    cex_col,
                    degree = 65,
                    scale_range_min,
                    scale_range_max,
                    cex_main = 1,
                    lab_row = NULL,
                    lab_col = NULL,
                    key = TRUE,
                    keysize = 1.5,
                    density_info = c("none", "histogram", "density"),
                    denscol = tracecol,
                    symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                    densadj = 0.25,
                    main = NULL,
                    xlab = NULL,
                    ylab = NULL,
                    lmat = NULL,
                    lhei = NULL,
                    lwid = NULL,
                    num_col_sidecolors = 1,
                    num_row_sidecolors = 1,
                    key_value_name, ...) {

        # created by obigriffith in
        # "https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/
        # master/Analysis/DifferentialExpression/R/heatmap_3.R"
        # # pulled from here, and then tweaked slightly:
        # http://www.biostars.org/p/18211/

        invalid <- function(x) {
            if (missing(x) || is.null(x) || length(x) == 0) {
                    return(TRUE)
            }
            if (is.list(x)) {
                    return(all(sapply(x, invalid)))
            } else if (is.vector(x)) {
                    return(all(is.na(x)))
            } else {
                    return(FALSE)
            }
        }

        x <- as.matrix(x)
        scale01 <- function(x, low = min(x), high = max(x)) {
            x <- (x - low) / (high - low)
            x
        }

        retval <- list()

        scale <- ifelse(symm && missing(scale), "none", match.arg(scale))

        dendrogram <- match.arg(dendrogram)

        trace <- match.arg(trace)

        density_info <- match.arg(density_info)

        if (length(col) == 1 && is.character(col)) {
                col <- get(col, mode = "function")
        }

        if (!missing(breaks) && (scale != "none")) {
            warning(
                "Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can produce unpredictable results.",
                "Please consider using only one or the other."
            )
        }

        if (is.null(rowv) || is.na(rowv)) {
            rowv <- FALSE
        }

        colv <- FALSE

        if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
            stop("`x' must be a numeric matrix")
        }

        nr <- di[1]
        nc <- di[2]

        if (nr <= 1 || nc <= 1) {
            stop("`x' must have at least 2 rows and 2 columns")
        }

        if (!is.numeric(margins) || length(margins) != 2) {
            stop("`margins' must be a numeric vector of length 2")
        }

        if (missing(cellnote)) {
            cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
        }

        if (!inherits(rowv, "dendrogram")) {
            tmp <- ((!isTRUE(rowv)) || (is.null(rowv)))
            if (tmp && (dendrogram %in% c("both", "row"))) {
                if (is.logical(colv) && (colv)) {
                    dendrogram <- "column"
                } else {
                    dendrogram <- "none"
                }

                warning(
                    "Discrepancy: Rowv is FALSE, while dendrogram is '",
                    dendrogram, "'. Omitting row dendogram."
                )
            }
        }

        if (!inherits(colv, "dendrogram")) {
            tmp <- ((!isTRUE(colv)) || (is.null(colv)))
            if (tmp && (dendrogram %in% c("both", "column"))) {
                if (is.logical(rowv) && (rowv)) {
                    dendrogram <- "row"
                } else {
                    dendrogram <- "none"
                }

                warning(
                    "Discrepancy: Colv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting column dendogram."
                )
            }
        }

        if (inherits(rowv, "dendrogram")) {
            ddr <- rowv
            row_ind <- order.dendrogram(ddr)
        }
        else if (is.integer(rowv)) {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, rowv)
            row_ind <- order.dendrogram(ddr)
            if (nr != length(row_ind)) {
                stop("row dendrogram ordering gave index of wrong length")
            }
        }
        else if (isTRUE(rowv)) {
            rowv <- rowMeans(x, na.rm = na.rm)
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, rowv)
            row_ind <- order.dendrogram(ddr)
            if (nr != length(row_ind)) {
                stop("row dendrogram ordering gave index of wrong length")
            }
        }
        else {
            row_ind <- nr:1
        }

        if (inherits(colv, "dendrogram")) {
            ddc <- colv
            col_ind <- order.dendrogram(ddc)
        }
        else if (identical(colv, "Rowv")) {
            if (nr != nc) {
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            }
            if (exists("ddr")) {
                ddc <- ddr
                col_ind <- order.dendrogram(ddc)
            }
            else {
                col_ind <- row_ind
            }
        }
        else if (is.integer(colv)) {
            hcc <- hclustfun(distfun(if (symm) {
                x
            } else {
                t(x)
            }))
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, colv)
            col_ind <- order.dendrogram(ddc)
            if (nc != length(col_ind)) {
                stop("column dendrogram ordering gave index of wrong length")
            }
        }
        else if (isTRUE(colv)) {
            colv <- col_means(x, na.rm = na.rm)
            hcc <- hclustfun(distfun(if (symm) {
                x
            } else {
                t(x)
            }))
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, colv)
            col_ind <- order.dendrogram(ddc)
            if (nc != length(col_ind)) {
                stop("column dendrogram ordering gave index of wrong length")
            }
        }
        else {
            col_ind <- 1:nc
        }

        retval$row_ind <- row_ind
        retval$col_ind <- col_ind
        retval$call <- match.call()

        x <- x[row_ind, col_ind] # rearrange matrix according to dendrograms

        cellnote <- cellnote[row_ind, col_ind] # also rearrange the cellnotes

        # get labels
        if (is.null(lab_row)) {
            lab_row <- if (is.null(rownames(x))) {
                    (1:nr)[row_ind]
                } else {
                    rownames(x)
                }
        } else {
            lab_row <- lab_row[row_ind]
        }
        if (is.null(lab_col)) {
            lab_col <- if (is.null(colnames(x))) {
                    (1:nc)[col_ind]
                } else {
                    colnames(x)
                }
        } else {
            lab_col <- lab_col[col_ind]
        }

        ## do scaling of matrix according to Z-scores
        if (scale == "row") {
            retval$row_means <- rm <- rowMeans(x, na.rm = na.rm)
            x <- sweep(x, 1, rm)
            retval$row_sds <- sx <- apply(x, 1, sd, na.rm = na.rm)
            x <- sweep(x, 1, sx, "/")
        }
        else if (scale == "column") {
            retval$col_means <- rm <- col_means(x, na.rm = na.rm)
            x <- sweep(x, 2, rm)
            retval$col_sds <- sx <- apply(x, 2, sd, na.rm = na.rm)
            x <- sweep(x, 2, sx, "/")
        }

        # number of breaks
        if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
            breaks <- ifelse(
                missing(col) || is.function(col), 16, length(col) + 1
            )
        }

        # set breakpoints
        if (length(breaks) == 1) {
            if (missing(scale_range_min)) {
                scale_range_min <- min(x, na.rm = na.rm)
            }

            if (missing(scale_range_max)) {
                scale_range_max <- max(x, na.rm = na.rm)
            }


            if (!symbreaks) {
                breaks <- seq(scale_range_min, scale_range_max, length = breaks)
            } else {
                extreme <- max(abs(c(scale_range_min, scale_range_max)),
                                                            na.rm = na.rm)
                breaks <- seq(-extreme, extreme, length = breaks)
            }
        }

        ncol <- length(breaks) - 1

        if (class(col) == "function") {
            col <- col(ncol)
        }

        min.breaks <- min(breaks)
        max.breaks <- max(breaks)

        # adjust for out-of-range given break settings
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks

        # layout height
        if (missing(lhei) || is.null(lhei)) {
            lhei <- c(keysize, 4)
        }

        # layout width
        if (missing(lwid) || is.null(lwid)) {
            lwid <- c(keysize, 4)
        }

        # define the layout
        if (missing(lmat) || is.null(lmat)) {
            lmat <- rbind(4:3, 2:1)

            if (!missing(col_sidecolors)) {
                tmp <- !is.character(col_sidecolors)
                if (tmp || ncol(col_sidecolors) != nc) {
                    stop("'col_sidecolors' must be a matrix of ncol(x) ",
                                                                nc, " columns")
                }
                lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
                side_height <- min(
                    side_height_fraction * nrow(col_sidecolors), 1
                )
                lhei <- c(lhei[1], side_height, lhei[2])
            }

            if (!missing(row_sidecolors)) {
                tmp <- !is.character(row_sidecolors)
                if (tmp || nrow(row_sidecolors) != nr) {
                    stop("'row_sidecolors' must be a matrix of nrow(x) ",
                        nr, " rows.  It currently has ",
                        nrow(row_sidecolors), " rows.")
                }
                lmat <- cbind(
                    lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1),
                    lmat[, 2] + 1
                )
                side_width <- min(
                    side_height_fraction * ncol(row_sidecolors), 1
                )
                lwid <- c(lwid[1], side_width, lwid[2])
            }
            lmat[is.na(lmat)] <- 0
        }

        if (length(lhei) != nrow(lmat)) {
            stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
        }
        if (length(lwid) != ncol(lmat)) {
            stop("lwid must have length = ncol(lmat) =", ncol(lmat))
        }

        op <- par(no.readonly = TRUE)
        on.exit(par(op))

        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

        # Draw the colorbars for the annotations:
        if (!missing(row_sidecolors)) {
            if (!is.matrix(row_sidecolors)) {
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = row_sidecolors[row_ind], axes = FALSE)
            } else {
                par(mar = c(margins[1], 0, 0, 0.5))
                rsc <- t(row_sidecolors[row_ind, , drop = FALSE])
                rsc_colors <- matrix()
                rsc_names <- names(table(rsc))
                rsc_i <- 1
                for (rsc.name in rsc_names) {
                    rsc_colors[rsc_i] <- rsc.name
                    rsc[rsc == rsc.name] <- rsc_i
                    rsc_i <- rsc_i + 1
                }
                rsc <- matrix(as.numeric(rsc), nrow = dim(rsc)[1])
                image(seq_len(nrow(rsc)), seq_len(ncol(rsc)), rsc,
                    col = as.vector(rsc_colors), axes = FALSE,
                    xlab = "", ylab = ""
                )

                # add labels
                if (length(colnames(row_sidecolors)) > 0) {
                    axis(1, seq_len(ncol(row_sidecolors)),
                        labels = colnames(row_sidecolors), las = 2,
                        tick = FALSE, xlab = "", ylab = ""
                    )
                }
            }
        }

        if (!missing(col_sidecolors)) {
            if (!is.matrix(col_sidecolors)) {
                par(mar = c(0.5, 0, 0, margins[2]))
                image(cbind(1:nc), col = col_sidecolors[col_ind], axes = FALSE)
            } else {
                par(mar = c(0.5, 0, 0, margins[2]))
                csc <- col_sidecolors[, col_ind, drop = FALSE]
                csc_colors <- matrix()
                csc_names <- names(table(csc))
                csc_i <- 1
                for (csc.name in csc_names) {
                    csc_colors[csc_i] <- csc.name
                    csc[csc == csc.name] <- csc_i
                    csc_i <- csc_i + 1
                }
                csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
                image(seq(1, nrow(t(csc))), seq(1, ncol(t(csc))), t(csc),
                    col = as.vector(csc_colors), axes = FALSE,
                    xlab = "", ylab = ""
                )

                # add labels
                if (length(rownames(col_sidecolors)) > 0) {
                    axis(2, 1:(nrow(col_sidecolors)),
                        labels = rownames(col_sidecolors), las = 2,
                        tick = FALSE
                    )
                }
            }
        }

        par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)
        cellnote <- t(cellnote)
        if (rev_c) {
            iy <- nr:1
            if (exists("ddr")) {
                ddr <- rev(ddr)
            }
            x <- x[, iy]
            cellnote <- cellnote[, iy]
        }
        else {
            iy <- 1:nr
        }

        # draw the central heatmap
        image(1:nc, 1:nr, x,
            xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr),
            axes = FALSE, xlab = "", ylab = "",
            col = col, breaks = breaks, ...
        )

        # store the matrix drawn
        retval$carpet <- x

        # store the dendrograms
        if (exists("ddr")) {
            retval$row_dendrogram <- ddr
        }
        if (exists("ddc")) {
            retval$col_dendrogram <- ddc
        }

        # store the breaks
        retval$breaks <- breaks

        # store the colormap used
        retval$col <- col

        # specially color in the na values
        if (!invalid(na_color) & any(is.na(x))) { # load library(gplots)
            mmat <- ifelse(is.na(x), 1, NA)
            image(1:nc, 1:nr, mmat,
                axes = FALSE, xlab = "", ylab = "",
                col = na_color, add = TRUE
            )
        }

        # X-axis column labels
        if (degree == 90) {
            axis(1, 1:nc, labels = lab_col, las = 2, line = -0.5, tick = 0)
        } else {
            tck <- axis(1, 1:nc,
                labels = FALSE, las = 2, line = -0.5,
                tick = 0
            )
            labels <- lab_col
            text(tck, par("usr")[3],
                labels = labels, srt = degree,
                xpd = NA, adj = c(1, 0.8), cex = cex_col
            )
        }

        # X-axis title
        if (!is.null(xlab)) {
            mtext(xlab, side = 1, line = margins[1] - 1.25)
        }

        # Y-axis row labeling
        axis(4, iy, labels = lab_row, las = 2, line = -0.5, tick = 0)

        # Y-axis title
        if (!is.null(ylab)) {
            mtext(ylab, side = 4, line = margins[2] - 1.25)
        }

        if (!missing(add_expr)) {
            eval(substitute(add_expr))
        }
        if (!missing(colsep)) {
            for (csep in colsep) rect(
                    xleft = csep + 0.5,
                    ybottom = rep(0, length(csep)),
                    xright = csep + 0.5 + sepwidth[1],
                    ytop = rep(ncol(x) + 1, csep),
                    lty = 1,
                    lwd = 1,
                    col = sepcolor, border = sepcolor
                )
        }
        if (!missing(rowsep)) {
            for (rsep in rowsep) rect(
                    xleft = 0,
                    ybottom = (ncol(x) + 1 - rsep) - 0.5,
                    xright = nrow(x) + 1,
                    ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
                    lty = 1,
                    lwd = 1,
                    col = sepcolor,
                    border = sepcolor
                )
        }


        min.scale <- min(breaks)
        max.scale <- max(breaks)
        x_scaled <- scale01(t(x), min.scale, max.scale)

        # column trace
        if (trace %in% c("both", "column")) {
            retval$vline <- vline
            vline_vals <- scale01(vline, min.scale, max.scale)
            for (i in col_ind) {
                if (!is.null(vline)) {
                    abline(v = i - 0.5 + vline_vals, col = linecol, lty = 2)
                }
                xv <- rep(i, nrow(x_scaled)) + x_scaled[, i] - 0.5
                xv <- c(xv[1], xv)
                yv <- seq_len(length(xv) - 0.5)
                lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
            }
        }

        # row trace
        if (trace %in% c("both", "row")) {
            retval$hline <- hline
            for (i in row_ind) {
                if (!is.null(hline)) {
                    abline(h = i + hline, col = linecol, lty = 2)
                }
                yv <- rep(i, ncol(x_scaled)) + x_scaled[i, ] - 0.5
                yv <- rev(c(yv[1], yv))
                xv <- seq(length(yv), 1) - 0.5
                lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
            }
        }

        # add cell labels
        if (!missing(cellnote)) {
            text(
                x = c(row(cellnote)), y = c(col(cellnote)),
                labels = c(cellnote), col = notecol,
                cex = notecex
            )
        }

        # Plot the row dendrogram
        par(mar = c(margins[1], 0, 0, 0))
        if (dendrogram %in% c("both", "row")) {
            plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
        }
        else {
            plot.new()
        }

        # Plot the column dendrogram
        par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
        if (dendrogram %in% c("both", "column")) {
            plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
        }
        else {
            plot.new()
        }

        if (!is.null(main)) {
            title(main, cex_main = cex_main)
        }

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
                min.raw <- min(c(x, breaks), na.rm = TRUE)
                max.raw <- max(c(x, breaks), na.rm = TRUE)
            }

            z <- seq(min.raw, max.raw, length = length(col))
            image(
                z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                xaxt = "n", yaxt = "n"
            )
            par(usr = c(0, 1, 0, 1))
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            axis(1, at = xv, labels = lv)
            if (scale == "row") {
                mtext(side = 1, "Row Z-Score", line = 2, cex = 0.7)
            } else if (scale == "column") {
                mtext(side = 1, "Column Z-Score", line = 2, cex = 0.7)
            } else {
                mtext(side = 1, key_value_name, line = 2, cex = 0.7)
            }
            if (density_info == "density") {
                dens <- density(x, adjust = densadj, na.rm = TRUE)
                omit <- dens$x < min(breaks) | dens$x > max(breaks)
                dens$x <- dens$x[-omit]
                dens$y <- dens$y[-omit]
                dens$x <- scale01(dens$x, min.raw, max.raw)
                lines(dens$x, dens$y / max(dens$y) * 0.95,
                    col = denscol,
                    lwd = 1
                )
                axis(2, at = pretty(dens$y) / max(dens$y) * 0.95,
                                                                pretty(dens$y))
                title("Color Key\nand Density Plot")
                par(cex = 0.3)
                mtext(side = 2, "Density", line = 2)
            }
            else if (density_info == "histogram") {
                h <- hist(x, plot = FALSE, breaks = breaks)
                hx <- scale01(breaks, min.raw, max.raw)
                hy <- c(h$counts, h$counts[length(h$counts)])
                lines(hx, hy / max(hy) * 0.95,
                    lwd = 1, type = "s",
                    col = denscol
                )
                axis(2, at = pretty(hy) / max(hy) * 0.95, pretty(hy))
                title("Color Key\nand Histogram")
                par(cex = 0.3)
                mtext(side = 2, "Count", line = 2)
            }
            else {
                title("Color Key")
            }
        }
        else {
            plot.new()
        }

        retval$color_table <- data.frame(
            low = retval$breaks[-length(retval$breaks)],
            high = retval$breaks[-1], color = retval$col
        )

        invisible(retval)
    }

    final_heatmap <- function(dist_method, scale_method) {

        # is the lab_row desired?
        if (!is.null(lab_row) && tolower(lab_row) == "rownames") {
            lab_row <- rownames(para_heatmaps)
        } else if (is.null(lab_row)) {
            lab_row <- ""
        }

        # is the lab_col desired?
        if (!is.null(lab_col) && tolower(lab_col) == "colnames") {
            lab_col <- colnames(para_heatmaps)
        } else if (is.null(lab_col)) {
            lab_col <- ""
        }

        if (tolower(image_format) == "png") {
            png(
                filename = paste0(
                    dir,
                    "/Heatmaps/", dist_method, "_",
                    method, "_", scale_method,
                    "_fc_cutoff=", fc_cutoff, "_",
                    pair_name, ".png"
                ),
                width = width, height = height, res = res, units = unit
            )
        } else if (tolower(image_format) == "svg") {
            svg(
                filename = paste0(
                    dir,
                    "/Heatmaps/", dist_method, "_",
                    method, "_", scale_method,
                    "_fc_cutoff=", fc_cutoff, "_",
                    pair_name, ".svg"
                ),
                width = width, height = height, onefile = TRUE
            )
        } else {
            stop(message(
                "Please, Insert a valid image_format!",
                " ('png' or 'svg')"
            ))
        }
        par(oma = outer_margins)
        heatmap_3(df,
            scale = scale_method,
            dendrogram = "both",
            margins = c(3, 9),
            rowv = row_dendro,
            keysize = 1.1,
            colv = col_dendro,
            lab_row = lab_row,
            lab_col = lab_col,
            symbreaks = FALSE,
            key = TRUE, symkey = FALSE,
            col_sidecolors = t(colour_groups),
            density_info = "none",
            trace = "none",
            cex_col = cex_col,
            cex_row = cex_row,
            degree = degree, col = rampa_de_cor,
            key_value_name = color_key_name
        )
        legend("topright",
            legend = paste0("G", levels(cond_heatmap)),
            title = "Groups",
            fill = possible_colors[as.numeric(levels(cond_heatmap))],
            border = FALSE,
            bty = "n", y.intersp = 0.7, cex = 0.7
        )
        dev.off()

        patient_order <- as.matrix(patient_order)
        patient_order <- cbind(seq_len(nrow(patient_order)), patient_order)
        colnames(patient_order) <- c("Dendrogram Order", "Patient")
        write.csv(patient_order, paste0(
            dir, "/Heatmaps/Patient_order_",
            dist_method, "_",
            method, "_", scale_method,
            "_fc_cutoff=", fc_cutoff, "_",
            pair_name, ".csv"
        ),
        row.names = FALSE
        )

        row_order <- as.matrix(row_order)
        row_order <- cbind(seq(1, nrow(resultados_de)), row_order)
        colnames(row_order) <- c("Dendrogram Order", "Row")
        write.csv(row_order, paste0(
            dir, "/Heatmaps/Row_order_",
            dist_method, "_",
            method, "_", scale_method,
            "_fc_cutoff=", fc_cutoff, "_",
            pair_name, ".csv"
        ),
        row.names = FALSE
        )
    }

    # code ####
    name <- gsub("-", "_", name)

    if (missing(env)) {
        stop(message(
            "The 'env' argument is missing, please",
            " insert the 'env' name and try again!"
        ))
    }

    ev <- deparse(substitute(env))
    sv <- list(ev = get(ev))

    path <- ifelse(exists("name_e", envir = get(ev)), file.path(
        sv[["ev"]]$path,
        sv[["ev"]]$name_e
    ), sv[["ev"]]$path)

    group_gen <- sv[["ev"]]$group_gen

    # creating the dir to outputs
    if (tolower(tool) == "ebseq") {
        dir <- paste0(
            path, "/EBSeq_Results.",
            tolower(group_gen), "_", toupper(name)
        )
        resultados_de <- sv[["ev"]]$resultados_de_ebseq[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_ebseq
    } else if (tolower(tool) == "edger") {
        dir <- paste0(
            path, "/edgeR_Results.",
            tolower(group_gen), "_", toupper(name)
        )
        resultados_de <- sv[["ev"]]$resultados_de_edger[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_edger
    } else if (tolower(tool) == "deseq2") {
        dir <- paste0(
            path, "/DESeq2_Results.",
            tolower(group_gen), "_", toupper(name)
        )
        resultados_de <- sv[["ev"]]$resultados_de_deseq2[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_deseq2
    } else if (tolower(tool) == "crosstable.deseq2") {
        dir <- paste0(path, "/CrossData_deseq2")
        resultados_de <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_deseq2
    } else if (tolower(tool) == "crosstable.edger") {
        dir <- paste0(path, "/CrossData_edger")
        resultados_de <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_edger
    } else if (tolower(tool) == "crosstable.ebseq") {
        dir <- paste0(path, "/CrossData_ebseq")
        resultados_de <- sv[["ev"]]$resultados_de_crossed[[pair_name]]
        normalized_expression <- sv[["ev"]]$normalized_expression_ebseq
    } else {
        stop(message(
            "Please, insert a valid tool name!",
            " ('EBSeq', 'DESeq2' or 'edgeR')"
        ))
    }

    dir.create(file.path(dir, "Heatmaps"), showWarnings = FALSE)

    cond_heatmap <- eval(parse(text = paste0(
        "sv[['ev']]",
        "$cond_heatmap"
    )))

    patients_stay <- unlist(strsplit(gsub("G", "", pair_name), "_over_"))

    # keep the desired group pair
    tmp <- cond_heatmap %in% patients_stay
    normalized_expression <- normalized_expression[, tmp]

    cond_heatmap <- droplevels(cond_heatmap[cond_heatmap %in% patients_stay])

    # Create table just with DE and with fc_cutoff
    resultados_de_up <- resultados_de[resultados_de$log2FC > log2(fc_cutoff), ]
    tmp <- resultados_de$log2FC < log2(fc_cutoff)
    resultados_de_down <- resultados_de[tmp, ]
    resultados_de <- rbind(resultados_de_up, resultados_de_down)

    # Select only the DE results to use in heatmaps (next step) removing the EE
    para_heatmaps <- as.matrix(normalized_expression[match(
        rownames(resultados_de),
        rownames(normalized_expression)
    ), ])
    colnames(para_heatmaps) <- seq(1, length(colnames(para_heatmaps)))

    # Preparing color ramp
    rampa_de_cor <- gplots::colorpanel(512, "blue", "white", "red")

    # Naming the heatmap cols
    possible_colors <- RColorBrewer::brewer.pal(8, "Set2")
    colour_groups <- possible_colors[cond_heatmap]
    colour_groups <- as.matrix(colour_groups)
    rownames(colour_groups) <- seq_len(nrow(colour_groups))
    colnames(colour_groups) <- "Groups"


    if (tolower(scale_method) == "none" && raw_values) {
        df <- log2(para_heatmaps + 1)
        color_key_name <- "Log2(Expression Values + 1)"
    } else if (tolower(scale_method) == "none") {
        df <- scale(log2(para_heatmaps + 1))
        color_key_name <- "Z-Score"
    } else {
        df <- log2(para_heatmaps + 1)
    }

    dist_matrix_col <- amap::Dist(t(df), method = method)
    dist_matrix_row <- amap::Dist(df, method = method)
    # complete
    dendro_euc_complete_col <- hclust(d = dist_matrix_col, method = "complete")
    dendro_euc_complete_row <- hclust(d = dist_matrix_row, method = "complete")
    row_dendro <- as.dendrogram(dendro_euc_complete_row)
    col_dendro <- as.dendrogram(dendro_euc_complete_col)
    tmp <- dendro_euc_complete_col$order
    patient_order <- colnames(normalized_expression)[tmp]
    row_order <- rownames(resultados_de)[dendro_euc_complete_row$order]
    method <- "complete"
    if (length(method) > 1) {
        for (Methods in method) {
            final_heatmap(dist_method = Methods, scale_method = scale_method)
        }
    } else {
        final_heatmap(dist_method = method, scale_method = scale_method)
    }

    # average
    dendro_euc_average_col <- hclust(d = dist_matrix_col, method = "average")
    dendro_euc_average_row <- hclust(d = dist_matrix_row, method = "average")
    row_dendro <- as.dendrogram(dendro_euc_average_row)
    col_dendro <- as.dendrogram(dendro_euc_average_col)
    tmp <- dendro_euc_average_col$order
    patient_order <- colnames(normalized_expression)[tmp]
    row_order <- rownames(resultados_de)[dendro_euc_average_row$order]
    method <- "average"
    if (length(method) > 1) {
        for (Methods in method) {
            final_heatmap(dist_method = Methods, scale_method = scale_method)
        }
    } else {
        final_heatmap(dist_method = method, scale_method = scale_method)
    }

    message("\nDone!\n")
}
