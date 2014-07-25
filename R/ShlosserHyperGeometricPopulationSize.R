 ShlosserHyperGeometricPopulationSize <- function(counts_hist, sampling_frac){
estimated_pop_size = sum(counts_hist)
numerator = counts_hist[1]*(1 - sampling_frac)*counts_hist[1]
denominator = sampling_frac*counts_hist[1]
for(i in 2:length(counts_hist)){
  numerator = numerator + counts_hist[1]*counts_hist[i]*(1 - sampling_frac)^i
  denominator = denominator + counts_hist[i]*i*sampling_frac*(1 - sampling_frac)^(i-1)
}
return(sum(counts_hist) + numerator/denominator)
}

ShlosserHyperGeometricPopulationSizeVariance <- function(counts_hist, sampling_prop){
p = sampling_prop
delta_f0 = mat.or.vec(nc = 1, nr = length(counts_hist) + 1)
denom = 0
numer = 0
numer_p_deriv = 0                
denom_p_deriv = 0
for(i in 1:length(counts_hist)){
  numer = numer + counts_hist[i]*(1 - p)^i
  denom = denom + counts_hist[i]*i*p*(1 - p)^(i - 1)
  numer_p_deriv = numer_p_deriv - i*counts_hist[i]*(1 - p)^(i - 1)
  denom_p_deriv = denom_p_deriv + i*counts_hist[i]*((1 - p)^(i - 1) - p*(i - 1)*(1 - p)^(i - 2))
}
delta_f0[1] = counts_hist[1]*(numer_p_deriv*denom - numer*denom_p_deriv)/(denom^2)
delta_f0[2] = numer/denom + counts_hist[1]*( ((1 - p)*denom - p*numer)/(denom^2) )
for(i in 2:length(counts_hist)){
  delta_f0[i + 1] = counts_hist[1]*( denom*(1 - p)^i - numer*i*p*(1 - p)^(i - 1))/(denom^2)
}
sigma = diag(c(p*(1 - p)/27900, counts_hist))
return(t(delta_f0) %*% sigma %*% delta_f0)
}
