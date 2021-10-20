####################################################################################
# QQ plot function of -log10 scale (with/without marginal histogram) used by
#     our project.
# Coding Partially inheritted from function 'qq' (R package: qqman [1])
# [1] Turner, S.D. qqman: an R package for visualizing GWAS results using Q-Q and 
#     manhattan plots. biorXiv DOI: 10.1101/005165 (2014).
####################################################################################

qqMargin <- function (pvector, x_breaks, y_breaks, x_lim, y_lim, margin = 0, ...) 
{
  require(ggplot2)
  require(ggExtra)
  if (!is.numeric(pvector)) 
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                       is.finite(pvector) & pvector < 1 & pvector > 0]
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  data1 = cbind.data.frame(o=o,e=e)
  
  if(!margin){
  g = ggplot(data1,aes(x=e,y=o)) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2, colour = "red") +
    xlab(expression(Expected ~ ~-log[10](italic(p)))) +
    ylab(expression(Observed ~ ~-log[10](italic(p)))) +
    scale_x_continuous(breaks = x_breaks, limits = x_lim) +
    scale_y_continuous(breaks = y_breaks, limits = y_lim)
  g
  }
  else{
  def_args <- list(type="histogram", size=6, bins=100)
  dotargs <- list(...)
  tryCatch(do.call("ggMarginal", c(list(p = g), def_args[!names(def_args) %in% 
                                                            names(dotargs)], dotargs)), warn = stop)
  }
}


