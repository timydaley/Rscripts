# counts_hist is a vector of count frequencies
# counts_hist[i,2] = # classes observed exactly counts_hist[i,1] times

# see Chiu, Wang, Walther, & Chao, Biometrics (2014)
ChaoImprovedLowerBound <- function(counts_hist){
        total_obs = sum(counts_hist[,1]*counts_hist[,2])
        return((counts_hist[1,2]^2)/(2*counts_hist[2,2])
                + ((total_obs - 3)*counts_hist[3,2]/(4*total_obs*counts_hist[4,2]))
                   *max(0, counts_hist[1,2] - (total_obs - 3)*counts_hist[2,2]*counts_hist[3,2]
                                              /(2*(total_obs - 1)*counts_hist[4,2])) )
}

