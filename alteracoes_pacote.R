ggplot(datalong[datalong$Uniprot.ID%in%c("P12109", "P06756", "P02144"),],
       aes(x=factor(condition), y=value, fill=Uniprot.ID)) +
    geom_boxplot(outlier.size=1.3) +
    xlab("condition") +
    ylab("Log2(Intensity)") +
    stat_boxplot(geom ='errorbar') +
    theme_bw() +
    scale_fill_manual(values=c("red","lightblue", "lightgreen")) +
    theme(axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size = rel(1.5)),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank())

facetsdf <- as.data.frame(mget(facetvars, envir = env)) # or inherits = TRUE

importFrom(magrittr,"%>%")
#' @importFrom packagename functionname Then you run
#' devtools::document()

cores <- ifelse(ParaGrafico == "Low", "#E41A1C","#377EB8")


qlf <- glmQLFTest(fit, coef=2)
go <- goana(qlf, species="Mm")
topGO(go, sort="up")
keg <- kegga(qlf, species="Mm")
topKEGG(keg, sort="up")


# ao rodar deseq2
# ############################################################################################## DO SITE
# # plot(attr(res,'filterNumRej'),type='b', ylab='number of rejections')
#
# #see this plot xlab to wald and not wald test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# stat <- res$stat
# maxCooks <- apply(assays(dds)[['cooks']], 1, max)
# png(file = paste0(DIR,
#                   "/Wald statistic_vs_maximum Cooks distance per gene.png"),
#     width = 2000, height = 1500, res = 300)
# plot(rank(stat[!is.na(stat)]), maxCooks[!is.na(stat)], xlab='rank of Wald statistic',
#      ylab='maximum Cooks distance per gene', las = 1,
#      ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
# m <- ncol(dds)
# p <- 3
# abline(h=qf(.99, p, m - p))
# dev.off()
#
# use <- res$pvalue < p.cutoff
# # table(use)
# h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
# h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
# colori <- c('do not pass'='khaki', 'pass'='powderblue')
# png(file = paste0(DIR,
#                   "/pass_vs_not pass barplot.png"),
#     width = 2000, height = 1500, res = 300)
# barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
#         col = colori, space = 0, xlab = "p-value", ylab='frequency')
# text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
#      adj = c(0.5,1.7), xpd=NA)
# legend('topright', fill=rev(colori), legend=rev(names(colori)), bty = "n")
# dev.off()
#
#
# # A new Bioconductor package, IHW , is now available that implements the
# # method of Independent Hypothesis Weighting [11].
# resIHW <- DESeq2::results(dds, filterFun = ihw, alpha = p.cutoff,
#                   contrast = c("condition", "High", "Low"))
#
# png(file = paste0(DIR,
#                   "/MA_plot.png"), width = 2000, height = 1500, res = 300)
# DESeq2::plotMA(resIHW, main="DESeq2", ylim = c(-2, 2))
# abline(h = log2(FC.cutoff), col = "dodgerblue", lwd = 2)
# abline(h = -log2(FC.cutoff), col = "dodgerblue", lwd = 2)
# dev.off()
#
# #compare dds vs varTrans heatmap
# library('RColorBrewer')
# library('gplots')
# select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
# hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
# #worst (linear)
# heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
#           Rowv = FALSE, Colv = FALSE, scale='none',
#           dendrogram='none', trace='none', margin=c(10,6))
# #better (log)
# heatmap.2(assay(varTrans)[select,], col = hmcol,
#           Rowv = FALSE, Colv = FALSE, scale='none',
#           dendrogram='none', trace='none', margin=c(10, 6))
#
# #more plots from ==> https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
# ############################################################################################## DO SITE
#
