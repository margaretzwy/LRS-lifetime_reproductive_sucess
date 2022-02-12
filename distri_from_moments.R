# Optimize a distribution which is best fit of the given raw moments using the Gram-Charlier method.

distri_from_moments <- function(raw_moments){
  #cumulants <- PDQutils::moment2cumulant(moments)
  xvals <- seq(0, 15, length.out = 5000) # pick up 15 for human data
  output <- PDQutils::dapx_gca(xvals, raw_moments)
  return(data.frame(mids = xvals, density = output))
}
