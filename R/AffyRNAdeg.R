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
          transform = "shift.scale",
          cols      = NULL,
          ...)
{
#   if (debug.xps()) print("------plotAffyRNAdeg------")

   if (!is.element(transform,c("shift.scale","shift.only","neither"))) {
      stop("transform must be one of <shift.scale, shift.only, neither>");
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

   if (is.null(cols)) cols = rep(4, nrow(mns));

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
}#plotAffyRNAdeg

#------------------------------------------------------------------------------#
