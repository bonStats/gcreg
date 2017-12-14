#' Transform beta to gamma polynomial coefficients
#' 
#' \code{b} or beta is standard polynomial coefficients w.r.t. raw data. \code{g} or
#' gamma is the orthogonalised polynomial coefficients w.r.t. scaled data (\eqn{-1 < x, y < 1}).
#' See \code{\link{untransform_gam}} for inverse transformation.
#' 
#' @param b beta, standard polynomial coefficients w.r.t. raw data.
#' @param g gamma, orthogonalised polynomial coefficients w.r.t. scaled data (\eqn{-1 < x, y < 1}).
#' @param sc_x_f scale functions for x data created by \code{\link{gen_scale_data_funs}}.
#' @param sc_y_f scale functions for y data created by \code{\link{gen_scale_data_funs}}.
#' @param cv orthogonal basis converter functions created by \code{\link{gen_poly_basis_converters}}.
#' @export

transform_beta <- function(b, sc_x_f, sc_y_f, cv){
  
  beta_par <- as.numeric(b)
  
  len <- length(beta_par)
  
  #scale x (range of x data was set to between -1 and 1)
  b_unsc <- attr(sc_x_f,"b") ^ (0:(length(beta_par)-1))
  sc_x_beta_par <- beta_par / b_unsc
  
  sc_x_beta <- as.numeric(
    change.origin(polynomial(sc_x_beta_par), o = - attr(sc_x_f,"c"))
  )
  
  sc_x_beta_par <- rep(0, times = len)
  sc_x_beta_par[1:length(sc_x_beta)] <- sc_x_beta
  
  # scale y
  sc_xy_beta_par <- as.numeric(sc_x_beta_par) * attr(sc_y_f,"b")
  sc_xy_beta_par <- replace(sc_xy_beta_par, 1, sc_xy_beta_par[1] + attr(sc_y_f,"c"))
  
  # un-orthonormalise (discrete polynomial orthogonalisation)
  gam_par <- as.numeric(cv$to_ortho(sc_xy_beta_par))
  
  return(gam_par)
  
}

#' Transform gamma to beta polynomial coefficients
#' 
#' \code{b} or beta is standard polynomial coefficients w.r.t. raw data. \code{g} or
#' gamma is the orthogonalised polynomial coefficients w.r.t. scaled data (\eqn{-1 < x, y < 1}).
#' See \code{\link{transform_beta}} for inverse transformation.
#' 
#' @inheritParams transform_beta
#' @export

untransform_gam <- function(g, sc_x_f, sc_y_f, cv){
  
  # un-orthonormalise (discrete polynomial orthogonalisation)
  sc_xy_beta_par <- as.numeric(cv$to_mono(g))
  
  len <- length(sc_xy_beta_par)
  
  # unscale y
  sc_x_beta_par <- replace(sc_xy_beta_par, 1, sc_xy_beta_par[1] - attr(sc_y_f,"c")) / attr(sc_y_f,"b")
  
  # un-scale (range of x data was set to between -1 and 1).
  b_unsc <- attr(sc_x_f,"b") ^ (0:(length(sc_xy_beta_par)-1))
  sc_x_origin_pl <- polynomial(sc_x_beta_par * b_unsc)
  
  beta_pl <- change.origin(sc_x_origin_pl, o = attr(sc_x_f,"c")/attr(sc_x_f,"b"))
  
  beta_num <- as.numeric(beta_pl)
  
  beta_par <- rep(0, times = len)
  beta_par[1:length(beta_num)] <- beta_num
  
  attr(beta_par, "beta_pl") <- beta_pl
  
  return(beta_par)
  
}
