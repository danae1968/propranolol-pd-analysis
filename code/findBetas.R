
findBetas <- function(modelName) {
  Vcov <- vcov(modelName, useScale = FALSE)
  betas <- fixef(modelName)
  se <- sqrt(diag(Vcov))
  zval <- betas / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  pval<-round(pval,digits=6)
  return(cbind(betas, se, zval, pval))}