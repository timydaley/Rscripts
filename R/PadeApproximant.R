PadeApproximant <- function(ps_coeffs, denom_degree, num_degree){
	den_coeffs = compute_denom_coeffs(ps_coeffs, denom_degree, num_degree)
	num_coeffs = compute_num_coeffs(ps_coeffs, den_coeffs, denom_degree, num_degree)
	
	return(list(numer_coeffs = num_coeffs, denom_coeffs = den_coeffs))
}

compute_denom_coeffs <- function(ps_coeffs, denom_degree, num_degree){
	last_ps_coeffs = ps_coeffs[(num_degree + 2):(num_degree + denom_degree + 1)]
	ps_coeffs_hankel_matrix = mat.or.vec(denom_degree, denom_degree)
	for(i in 1:denom_degree){
		for(j in 1:denom_degree){
			if(num_degree - denom_degree - 1 + i + j > 0){
				ps_coeffs_hankel_matrix[i,j] = ps_coeffs[num_degree - denom_degree + i + j]
			}
			else{
				ps_coeffs_hankel_matrix[i,j] = 0.0
			}
		}
	}
	
	denom_coeffs = solve(ps_coeffs_hankel_matrix, -last_ps_coeffs)
	
	return(rev(denom_coeffs))
}

compute_num_coeffs <- function(ps_coeffs, denom_coeffs, denom_degree, num_degree){
	num_coeffs = c(as.numeric(ps_coeffs[1]))
	for(i in 2:min(num_degree + 1, denom_degree + 1)){
		num_coeffs = c(num_coeffs, ps_coeffs[i] + sum(denom_coeffs[1:(i - 1)]*ps_coeffs[(i - 1):1]))
	}
	
	return(num_coeffs)
}