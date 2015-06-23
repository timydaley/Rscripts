# counts_hist is a vector of count frequencies
# counts_hist[i,2] = # classes observed exactly counts_hist[i,1] times

# see Chiu, Wang, Walther, & Chao, Biometrics (2014)
ZeltermanLowerBound <- function(counts_hist){
  distinct = sum(counts_hist[,2])
  return( distinct/(1 - exp( -2*counts_hist[2,2]/counts_hist[1,2])) )
}

