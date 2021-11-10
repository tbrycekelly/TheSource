#' @title Add Bootstrapped Trendline
#' @author Thomas Bryce Kelly
#' @description Add the maximum likelihood trendline to a figured based on a bootstrap estimation.
#' @keywords Statistics
#' @export
add.boot.trendline = function(model, col = 'black', lty = 2, lwd = 1, ...) {
  warning('Depreciated function, please use add.regess.trendline instead.')
  abline(a = median(model$models$b), b = median(model$models$m), col = col, lty = lty, lwd = lwd, ...)
}

#' @title Add Bootstrapped Trendline
#' @author Thomas Bryce Kelly
#' @description Add the maximum likelihood trendline to a figured based on a bootstrap estimation.
#' @keywords Statistics
#' @export
add.regress.trendline = function(model, col = 'black', lty = 2, lwd = 1, ...) {
  abline(a = median(model$models$b), b = median(model$models$m), col = col, lty = lty, lwd = lwd, ...)
}



#' @title Add Linear Model Confidence Intervals
#' @author Thomas Bryce Kelly
#' @description A helper function to plot the confidence intervals determined from a base::lm model.
#' @keywords Statistics
#' @export
add.lm.conf = function(model, x = NULL, name = NULL, col = '#50505030',
                       level = 0.95, ...) {

  ## Apply defaults if not provided
  if (is.null(x)) {
    x = model$model[[2]]
  }
  if (is.null(name)) {
    name = names(model$model)[2]
  }

  ## Plot
  dat = data.frame(a = c(1:length(x)))
  dat[[name]] = x
  pred = predict(model, interval='confidence', newdata = dat, level = level)
  polygon(x = c(x, rev(x)), y = c(pred[,"lwr"], rev(pred[,"upr"])), border = NA, col = col, ...)
}


#' @title Add Linear Model Trendline
#' @author Thomas Bryce Kelly
#' @description A helper function to plot the trendline determined from a base::lm model.
#' @keywords Statistics
#' @export
add.lm.trendline = function(model, col = 'black', lwd = 1, lty = 1, ...) {
  adline(model, col = col, lwd = lwd, lty = lty, ...)
}


#' @title Add Bootstrapped Confidence Intervals
#' @author Thomas Bryce Kelly
#' @description Add confidence bands to a figure based on results of a bootstrap.
#' @keywords Statistics Uncertainty
#' @param model A bootstrap object (i.e. a dataframe) containing bootstrapped estimates of m and b.
#' @param new.x The x values for which predictions are required.
#' @param col The desired color of the confidence band. We recomend colors with transparency for plotting.
#' @param conf The quantile ranges for the plotting (default is 95% two-tail).
#' @param border Do you want a border on the shaded confidence interval?
#' @param trendline Should the maximum liklihood values be plotted as a trendline?
#' @export
add.boot.conf = function(model, x = NULL, col = '#55555540', conf = c(0.025, 0.975),
                         border = FALSE, trendline = FALSE, n = 1e3, ...) {
  warning('Depreciated function, please use add.regess.conf instead.')
  if (is.null(x)) {
    x.new = seq(min(pretty(model$data$x)), max(pretty(model$data$x)), length.out = n)
  } else {
    x.new = x
  }
  y.upper = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = conf[2])})
  y.lower = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = conf[1])})
  y.mid = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = 0.5)})

  polygon(x = c(x.new, rev(x.new)), y = c(y.upper, rev(y.lower)), col = col, border = border, ...)

  if (trendline) {
    add.boot.trendline(model, ...)
  }
}

#' @title Add Bootstrapped Confidence Intervals
#' @author Thomas Bryce Kelly
#' @description Add confidence bands to a figure based on results of a bootstrap.
#' @keywords Statistics Uncertainty
#' @param model A bootstrap object (i.e. a dataframe) containing bootstrapped estimates of m and b.
#' @param new.x The x values for which predictions are required.
#' @param col The desired color of the confidence band. We recomend colors with transparency for plotting.
#' @param conf The quantile ranges for the plotting (default is 95% two-tail).
#' @param border Do you want a border on the shaded confidence interval?
#' @param trendline Should the maximum liklihood values be plotted as a trendline?
#' @export
add.regress.conf = function(model, x = NULL, col = '#55555540', conf = c(0.025, 0.975),
                            border = FALSE, trendline = FALSE, n = 1e3, ...) {
  if (is.null(x)) {
    x.new = seq(min(pretty(model$data$x)), max(pretty(model$data$x)), length.out = n)
  } else {
    x.new = x
  }
  y.upper = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = conf[2])})
  y.lower = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = conf[1])})
  y.mid = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = 0.5)})

  polygon(x = c(x.new, rev(x.new)), y = c(y.upper, rev(y.lower)), col = col, border = border, ...)

  if (trendline) {
    add.regress.trendline(model, ...)
  }
}
