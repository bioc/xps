#==============================================================================#
# AffyRNAdeg.R: quality control functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AffyRNAdeg:  
# summaryAffyRNAdeg: 
# plotAffyRNAdeg: 
#==============================================================================#


"AffyRNAdeg" <-
function (xps.data,
          treename = "*",
          qualopt  = "raw",
          log.it   = TRUE)
{
   if (debug.xps()) print("------AffyRNAdeg------")

   if (is(xps.data, "QualTreeSet")) {
      if (!log.it ) stop("currently only log2 is implemented");

      lm.stats <- function(x) {
         index <- 0:(length(x) - 1);
         res   <- summary(lm(x ~ index))$coefficients[2, c(1, 4)];
         return(res);
      }#lm.stats

      rnadeg <- xpsRNAdeg(xps.data, treename=treename, qualopt=qualopt);

      mns <- rnadeg$mns;
      mn  <- mns[, 1];
      mns <- sweep(mns, 1, mn);
      mns <- mns/(rnadeg$ses);

      stats <- apply(mns, 1, lm.stats);

      return(c(rnadeg, slope=list(stats[1,]), pvalue=list(stats[2,])));
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("QualTreeSet")));
   }#if
}#AffyRNAdeg

#------------------------------------------------------------------------------#

"summaryAffyRNAdeg" <-
function (rna.deg, signif.digits = 3)
{
  table <- rbind(signif(rna.deg$slope,  signif.digits),
                 signif(rna.deg$pvalue, signif.digits));

  colnames(table) <- rna.deg$sample.names;
  rownames(table) <- c("slope", "pvalue");

  return(table);
}#summaryAffyRNAdeg

#------------------------------------------------------------------------------#

"plotAffyRNAdeg" <-
function (rna.deg,
          transform  = "shift.scale",
          cols       = NULL,
          add.legend = FALSE,
          ...)
{
   if (debug.xps()) print("------plotAffyRNAdeg------")

   if (!is.element(transform,c("shift.scale","shift.only","none"))) {
      stop("transform must be one of <shift.scale, shift.only, none>");
   }#if

   ylab <- "Mean Intensity";

   mns <- rna.deg$mns;
   if (transform == "shift.scale") {
      sds <- rna.deg$ses;
      mn  <- mns[, 1];
      mns <- sweep(mns, 1, mn);
      mns <- mns/(sds);
      mns <- sweep(mns, 1, 1:nrow(mns), "+");

      ylab <- paste(ylab, ": shifted and scaled");
   } else if (transform == "shift.only") {
      mn  <- mns[, 1];
      mns <- sweep(mns, 1, mn);
      mns <- sweep(mns, 1, 1:nrow(mns), "+");

      ylab <- paste(ylab, ": shifted");
   }#if

   if (is.null(cols)) {
      cols <- c("blue3", "blue2", "blue1", "steelblue3", "steelblue2", "steelblue1",
                "lightblue3", "lightblue2", "lightblue1", "gray60",
                "red3", "red2", "red1", "orange3", "orange2", "orange1",
                "yellow3", "yellow2", "yellow1", "black");
      cols <-rep(cols, (floor(nrow(mns)/length(cols)) + 1));
      lty  <- unlist(lapply(1:5,function(x)rep(x,20)));
      lty  <- rep(lty, (floor(nrow(mns)/(5*length(cols))) + 1));
   }#if

   plot(-2, -1,
        pch = "",
        xlim = range(-1, ncol(mns)),
        ylim = range(min(as.vector(mns)) - 1, max(as.vector(mns)) + 1),
        xlab = "5' <-----> 3'\n Probe Number ",
        ylab = ylab,
        axes = FALSE,
        main = "RNA degradation plot",
        ...)
   axis(1);
   axis(2);

   for (i in 1:nrow(mns)) {
      lines(0:((ncol(mns) - 1)), mns[i, ], col=cols[i]);
   }#for

   if (add.legend) {
      if (is.numeric(add.legend)) {
         n <- min(add.legend, nrow(mns));
      } else {
         n <- nrow(mns);
      }#if

      legend(x      = "topleft",
             legend = namePart(rna.deg$sample.names)[1:n],
             lty    = lty[1:n],
             pt.bg  = "white",
             col    = cols[1:n],
             cex    = 0.6
            );
   }#if
}#plotAffyRNAdeg

#------------------------------------------------------------------------------#
