
fragSizeDist <- function (bamFiles, bamFiles.labels, index = bamFiles, ylim = NULL, 
          logYlim = NULL, color = "dark blue") 
{
  opar <- par(c("fig", "mar"))
  on.exit(par(opar))
  pe <- mapply(testPairedEndBam, bamFiles, index)
  if (any(!pe)) {
    stop(paste(bamFiles[!pe], collapse = ", "), "is not Paired-End file.")
  }
  summaryFunction <- function(seqname, seqlength, bamFile, 
                              ind, ...) {
    param <- ScanBamParam(what = c("isize"), which = GRanges(seqname, 
                                                             IRanges(1, seqlength)), flag = scanBamFlag(isSecondaryAlignment = FALSE, 
                                                                                                        isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE))
    table(abs(unlist(sapply(scanBam(bamFile, index = ind, 
                                    ..., param = param), `[[`, "isize"), use.names = FALSE)))
  }
  idxstats <- unique(do.call(rbind, mapply(function(.ele, 
                                                    .ind) idxstatsBam(.ele, index = .ind)[, c("seqnames", 
                                                                                              "seqlength")], bamFiles, index, SIMPLIFY = FALSE)))
idxstats <- idxstats[-which(idxstats$seqnames == "*"),]
  seqnames <- as.character(idxstats[, "seqnames"])
  seqlen <- as.numeric(idxstats[, "seqlength"])
  fragment.len <- mapply(function(bamFile, ind) summaryFunction(seqname = seqnames, 
                                                                seqlength = seqlen, bamFile, ind), bamFiles, index, 
                         SIMPLIFY = FALSE)
  names(fragment.len) <- bamFiles.labels
  minor.ticks.axis <- function(ax, n = 9, t.ratio = 0.5, mn, 
                               mx, ...) {
    lims <- par("usr")
    lims <- if (ax %in% c(1, 3)) 
      lims[1:2]
    else lims[3:4]
    major.ticks <- pretty(lims, n = 5)
    if (missing(mn)) 
      mn <- min(major.ticks)
    if (missing(mx)) 
      mx <- max(major.ticks)
    major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= 
                                 mx]
    labels <- sapply(major.ticks, function(i) as.expression(bquote(10^.(i))))
    axis(ax, at = major.ticks, labels = labels, las = ifelse(ax %in% 
                                                               c(2, 4), 2, 1), ...)
    n <- n + 2
    minors <- log10(pretty(10^major.ticks[1:2], n)) - major.ticks[1]
    minors <- minors[-c(1, n)]
    minor.ticks = c(outer(minors, major.ticks, `+`))
    minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < 
                                 mx]
    axis(ax, at = minor.ticks, tcl = par("tcl") * t.ratio, 
         labels = FALSE)
  }
  null <- mapply(function(frag.len, frag.name) {
    x <- 1:1010
    frag.len <- frag.len[match(x, names(frag.len))]
    frag.len[is.na(frag.len)] <- 0
    y <- frag.len/sum(frag.len)
    y <- as.numeric(y)
    names(y) <- x
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(x, y * 10^3, main = paste(frag.name, "fragment sizes"), col=color,bty="l",
         xlim = c(0, 1010), ylim = ylim, xlab = "Fragment length (bp)", 
         ylab = expression(Normalized ~ read ~ density ~ 
                             x ~ 10^-3), type = "l")
    par(fig = c(0.4, 0.95, 0.4, 0.95), new = TRUE)
    plot(x, log10(y), xlim = c(0, 1010), ylim = logYlim, 
         col = color, xlab = "Fragment length (bp)", ylab = "Norm. read density",bty="l",
         type = "l", yaxt = "n")
    minor.ticks.axis(2)
    par(opar)
  }, fragment.len, names(fragment.len))
  return(invisible(fragment.len))
}

