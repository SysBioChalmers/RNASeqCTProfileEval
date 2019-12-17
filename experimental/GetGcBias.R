GetScBias <- function(expr1, expr2, lb=-3, ub=12) {
  divvar = log2((expr1 + 0.05) / (expr2 + 0.05))
  #meanExprVar = log2((expr1 + expr2)/2 + 0.05)
  meanExprVar = log2(sqrt(expr1*expr2) + 0.05)
  #meanExprVar = log2(expr2 + 0.05)
  
  step = 0.5
  aimedxes = seq(lb, ub, by=step)
  uppers = aimedxes + step/2
  lowers = aimedxes - step/2
  
  #allocate of the correct length
  x = aimedxes
  y = aimedxes
  
  for (i in 1:length(aimedxes)) {
    ind = (meanExprVar >= lowers[i]) & (meanExprVar <= uppers[i])
    x[i] = mean(meanExprVar[ind])
    y[i] = mean(divvar[ind])
  }
  
  return (cbind(x, y))
}

