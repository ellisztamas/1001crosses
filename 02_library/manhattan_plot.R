#!/usr/bin/env Rscript


# Function to generate Manhattan plots
# From https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R

library(lattice)
library(tidyverse)

#' X-axis for a Manhattan plot
#'
#' Convert columns of chromosome labels and SNP positions into a single vector
#' that can be used as an x-axis on a plot.
#'
#' @inheritParams manhattan_plot
#'
#' @return The input dataframe, with an additional column added, showing the
#' cumulative SNP position over each chromosome.
#' @author Tom Ellis from code by Daniel Roeffs
add_base_pair_positions <- function(gwas_data){
  # Get the maximum base-pair position on each chromosome
  data_cum <- gwas_data %>%
    group_by(chr) %>%
    summarise(max_bp = max(ps), .groups = "drop_last") %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    dplyr::select(chr, bp_add)
  # Cumulative bp position for each SNP along the whole chromosome
  gwas_data <- gwas_data %>%
    inner_join(data_cum, by = "chr") %>%
    mutate(ps_cum = ps + bp_add)

  gwas_data
}


manhattan_plot<-function(chr, pos, pvalue,
                         sig.level=NA,
                         annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2,
                         xlab="Chromosome", ylab=expression(-log[10](p-value)), title=NULL,
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")

  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }

  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;

  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }

  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5,
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }

  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    }
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }

  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

  if (length(ann.settings)>1) {
    lcols<-trellis.par.get("superpose.symbol")$col
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch,
                              col=lcols[(i-2) %% length(lcols) +1 ],
                              fill=lfills[(i-2) %% length(lfills) +1 ],
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label,
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label,
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]],
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)

  # reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places),
      pos=round(genpos,thin.pos.places),
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()

  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr),
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }

  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) {
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }

  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings,
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab,main=title,
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}


# QQ-plots for GWAS -log10(p-values)
# Code by Pieter Clauw
plot.qq <- function(gwas.pval, ci = 0.95){
  # get number of tests
  nTests = length(gwas.pval)
  qq.dat <- tibble(
    observed = -log10(sort(gwas.pval)),
    expected = -log10(ppoints(nTests)),
    cLower = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nTests), shape2 = rev(seq(nTests)))),
    cUpper = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nTests), shape2 = rev(seq(nTests)))))
  # make QQ-plot
  qq.plt <- ggplot(qq.dat, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cUpper, ymin = cLower), fill = "grey30", alpha = 0.5) +
    geom_step(color = wes_palettes['Darjeeling1'][[1]][1], size = 1.1, direction = "vh") +
    geom_segment(data = . %>% filter(expected == max(expected)),
                 aes(x = 0, xend = expected, y = 0, yend = expected),
                 size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
    labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
         y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
    theme_minimal()

  return(qq.plt)
}

