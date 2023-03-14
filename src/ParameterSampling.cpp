#include <RcppArmadillo.h>
#include <cmath>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' @description Draws the transition probabilities for models that are governed by one regime process
//' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @param heteroscRegime boolean. If true, the variable names r and s are used instead of q and p
//' @return the filtered output
// [[Rcpp::export]]
rowvec TransProbs_fctn(vec regimeVec, List paramList, bool heteroscRegime = false)
{
  List Probs = paramList["Probs"];
  double Backup1;
  double Backup2;
  if (heteroscRegime == false)
  {
    Backup1 = Probs["q"];
    Backup2 = Probs["p"];
  }
  else
  {
    Backup1 = Probs["r"];
    Backup2 = Probs["s"];
  }

  // discard first 10 observations because of the high state vector variation due to the filter initialization
  regimeVec = regimeVec.subvec(10, regimeVec.n_elem - 1);
  // Count the transitions
  int S_lt = regimeVec[0];
  double n_00 = 0;
  double n_10 = 0;
  double n_01 = 0;
  double n_11 = 0;
  int S_t;
  double q;
  double p;
  rowvec oututVec;

  for (int i = 1; i < regimeVec.n_elem; i++)
  {
    S_t = regimeVec[i];
    if (S_lt == S_t)
    {
      if (S_t == 0)
      {
        n_00++;
      }
      else
      {
        n_11++;
      }
    }
    else
    {
      if (S_t == 1)
      {
        n_01++;
      }
      else
      {
        n_10++;
      }
    }
    S_lt = S_t;
  }

  // Set up the priors
  double u_00 = 0;
  double u_11 = 0;
  double u_01 = 0;
  double u_10 = 0;

  // Draw the sample
  q = R::rbeta(u_00 + n_00, u_01 + n_01);
  int threshold = 1e4;
  int counter;
  double lowerLimit = .9;
  double upperLimit = .997;

  counter = 0;
  if (heteroscRegime == false)
  {
    while ((q < lowerLimit || q > upperLimit) && counter <= threshold)
    {
      counter++;
      q = R::rbeta(u_00 + n_00, u_01 + n_01);
    }
    if (q < lowerLimit || q > upperLimit)
    {
      q = Backup1;
    }
  }
  else
  {
    while ((q < .01 || q > upperLimit) && counter <= threshold)
    {
      counter++;
      q = R::rbeta(u_00 + n_00, u_01 + n_01);
    }
    if (q < .01 || q > upperLimit)
    {
      q = Backup1;
    }
  }

  p = R::rbeta(u_11 + n_11, u_10 + n_10);
  if (heteroscRegime == false)
  {
    counter = 0;
    while ((p < lowerLimit || p > upperLimit) && counter <= threshold)
    {
      counter++;
      p = R::rbeta(u_11 + n_11, u_10 + n_10);
    }
    if (p < lowerLimit || p > upperLimit)
    {
      p = Backup2;
    }
  }
  else
  {
    counter = 0;
    while ((p < .01 || p > upperLimit) && counter <= threshold)
    {
      counter++;
      p = R::rbeta(u_11 + n_11, u_10 + n_10);
    }
    if (p < .01 || p > upperLimit)
    {
      p = Backup2;
    }
  }
  oututVec = {q, p};
  return oututVec;
}

//' @description Draws the added negative drift for regime 1 for models governed by two regime processes
//' @param stateVec array with draws of the state vector
//' @param regimeDriftVec vector with regime realizations governing the drift. Vector takes either values 1 or 0
//' @param regimeHeteroscVec vector with regime realizations governing the heteroscedasticity. Vector takes either values 1 or 0.
//' If NULL, the regimeDriftVec will be employed here as well
//' @param paramList list with the additional parameters
//' @return the filtered output
// [[Rcpp::export]]
double Nu_1_fctn(mat stateVec, vec regimeDriftVec, List paramList, Nullable<arma::vec> regimeHeteroscVec = R_NilValue)
{
  double nu_1_backup = paramList["nu_1"];

  // discard first 10 observations because of the high state vector variation due to the filter initialization
  int nPeriods = stateVec.n_cols;
  mat stateVecRed = stateVec.cols(10, nPeriods - 1);
  vec regimeDriftVecRed = regimeDriftVec.subvec(10, nPeriods - 1);
  vec regimeHeteroscVecRed;
  if (regimeHeteroscVec.isNull())
  {
    regimeHeteroscVecRed = regimeDriftVecRed;
  }
  else
  {
    regimeHeteroscVecRed = as<arma::vec>(regimeHeteroscVec).subvec(10, nPeriods - 1);
  }
  List Xi = paramList["xi"];
  double nu_1Backup = paramList["nu_1"];
  double SdXi_0 = Xi["xi_0"];
  double SdXi_1 = Xi["xi_1"];
  // Set up the priors
  double zeta_0 = 0;
  double delta_0 = 1e6;
  // Compute the SSR from the regression of of the trend on its lag and drift terms, adjusted for xi so that the error is standard normal
  rowvec mu = stateVecRed.row(0);
  int nPeriodsRed = mu.n_elem;
  rowvec mu_t = mu.tail(nPeriodsRed - 1);
  rowvec mu_lt = mu.head(nPeriodsRed - 1);
  rowvec nu_0 = stateVecRed.row(1).tail(nPeriodsRed - 1);
  vec SdXiVec = regimeHeteroscVecRed.tail(nPeriodsRed - 1) * SdXi_1;
  SdXiVec.elem(find(SdXiVec == 0)).fill(SdXi_0);
  mat Y = mat((mu_t - mu_lt - nu_0) / SdXiVec.t()).t();
  mat X = mat(regimeDriftVecRed.tail(nPeriodsRed - 1) / SdXiVec);
  mat transX = trans(X);
  double invXtX = as_scalar(inv(1 / delta_0 + transX * X));
  double zeta_1 = invXtX * as_scalar(transX * Y);
  // Constrain the parameter space
  double lowerLimit = -10;
  double upperLimit = -.01;
  if (zeta_1 >= upperLimit | zeta_1 <= lowerLimit)
  {
    zeta_1 = nu_1_backup;
  }
  double delta_1 = invXtX;
  double draw;
  draw = zeta_1 + as_scalar(randn(1, 1)) * delta_1;
  // if none out of 1000 draws produces a negative value, the function does not return an output
  int counter = 1;
  int threshold = 1e3;
  while ((draw >= upperLimit || draw <= lowerLimit) && counter <= threshold)
  {
    draw = zeta_1 + as_scalar(randn(1, 1)) * delta_1;
    counter++;
  }
  if (draw < upperLimit && draw > lowerLimit)
  {
    return draw;
  }
  else
  {
    return nu_1_backup;
  }
}

//' @description Draws the innovation standard deviation for the seasonal unit root process
//' @param stateVec array matrix with the draw of the state vector
//' @param paramList list with the additional parameters
//' @param seasPos integer with the position of the seasonal term within the state vector
//' @return the sampled standard deviation
// [[Rcpp::export]]
double SdOmega_fctn(mat stateVec, List paramList, int seasPos)
{
  List omega = paramList["omega"];
  // discard first 10 observations because of the high state vector variation due to the filter initialization
  mat stateVecRed = stateVec.cols(10, stateVec.n_cols - 1);
  int Dimens = stateVecRed.n_rows;
  int nPeriods = stateVecRed.n_cols;
  // Set the priors
  double alpha_0 = 0;
  double beta_0 = 0;
  // Constrain the parameter space
  double lowerLimit = 1e-5;
  double upperLimit = math::inf();
  // Compute the SSR based on the seasonal unit root (regression of the seasonal term on its determ dummies and the unit root)
  rowvec gamma = stateVecRed.row(seasPos - 1);
  rowvec gamma_t = gamma.subvec(6, nPeriods - 1);
  rowvec gamma_tl1 = gamma.subvec(5, nPeriods - 2);
  rowvec gamma_tl2 = gamma.subvec(4, nPeriods - 3);
  rowvec gamma_tl3 = gamma.subvec(3, nPeriods - 4);
  rowvec gamma_tl4 = gamma.subvec(2, nPeriods - 5);
  rowvec gamma_tl5 = gamma.subvec(1, nPeriods - 6);
  rowvec gamma_tl6 = gamma.subvec(0, nPeriods - 7);
  rowvec x_t = gamma_t + gamma_tl1 + gamma_tl2 + gamma_tl3 + gamma_tl4 + gamma_tl5 + gamma_tl6;
  rowvec x_tl7 = x_t.head(x_t.n_elem - 7);
  rowvec omega_t = x_t.tail(x_t.n_elem - 7) - x_tl7;
  double SSR = accu(square(omega_t.head(omega_t.n_elem - 1)));
  // Adjust the priors
  double alpha_1 = alpha_0 + nPeriods;
  double beta_1 = beta_0 + SSR;
  // Draw the new value
  double omega_var = 1 / randg(distr_param(alpha_1 / 2.0, 2.0 / beta_1));
  double omega_sd_candidate = sqrt(omega_var);
  double omega_sd_draw;
  if (omega_sd_candidate < lowerLimit || omega_sd_candidate > upperLimit)
  {
    omega_sd_draw = omega["omega_0"];
  }
  else
  {
    omega_sd_draw = omega_sd_candidate;
  }
  return omega_sd_draw;
}

//' @description Function to draw the innovation standard deviation of the trend component
//' @param stateVec matrix of the state vector draw
//' @param paramList list with the additional parameters
//' @param regimeVec vector with ones and zeros of the regime draw
//' @return draw of the standard deviation
// [[Rcpp::export]]
double SdXiConst_fctn(mat stateVec, List paramList, vec regimeVec)
{
  // discard first 10 observations because of the high state vector variation due to the filter initialization
  mat stateVecRed = stateVec.cols(10, stateVec.n_cols - 1);
  vec regimeVecRed = regimeVec.subvec(10, regimeVec.n_elem - 1);
  // Set the priors
  double alpha_0 = 0;
  double beta_0 = 0;
  // Constrain the parameter space
  double lowerLimit = 1e-5;
  double upperLimit = math::inf();
  // Compute the SSR from the regression of the trend on its lag and drift terms
  double nu_1 = paramList["nu_1"];
  List xi = paramList["xi"];
  rowvec mu = stateVecRed.row(0);
  int nPeriods = mu.n_elem;
  rowvec mu_t = mu.tail(nPeriods - 1);
  rowvec nu_0 = stateVecRed.row(1).tail(nPeriods - 1);
  rowvec mu_lt = mu.head(nPeriods - 1);
  vec S_t = regimeVecRed.tail(nPeriods - 1);
  rowvec v_t = mu_t - mu_lt - nu_0 - S_t.t() * nu_1;
  double SSR = accu(square(v_t));
  // Adjust the priors
  double alpha_1 = (alpha_0 + nPeriods);
  double beta_1 = beta_0 + SSR;
  // Draw the new value
  double xi_var = 1 / randg(distr_param(alpha_1 / 2.0, 2.0 / beta_1));
  double xi_sd_candidate = sqrt(xi_var);
  double xi_sd_draw;
  if (xi_sd_candidate < lowerLimit || xi_sd_candidate > upperLimit)
  {
    xi_sd_draw = xi["xi_0"];
  }
  else
  {
    xi_sd_draw = xi_sd_candidate;
  }
  return xi_sd_draw;
}

//' @description Function to draw the innovation standard deviation of the trend component in case of regime induced heterocedasticity
//' @param stateVec matrix of the state vector draw
//' @param regimeDriftVec vector with regime realizations governing the trend drift. Vector takes either values 1 or 0
//' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @return vector of the draw of the two standard deviations
// [[Rcpp::export]]
vec SdXiSwitch_fctn(mat stateVec, vec regimeDriftVec, vec regimeHeteroscVec, List paramList)
{
  // discard first 10 observations because of the high state vector variation due to the filter initialization
  mat stateVecRed = stateVec.cols(10, stateVec.n_cols - 1);
  vec regimeDriftVecRed = regimeDriftVec.subvec(10, regimeDriftVec.n_elem - 1);
  vec regimeHeteroscVecRed = regimeHeteroscVec.subvec(10, regimeHeteroscVec.n_elem - 1);

  // Set the priors
  double alpha_0 = 0;
  double beta_0 = 0;
  double alpha_2 = .1;
  double beta_2 = 0;
  // Define variables
  int nPeriods = stateVecRed.n_cols;
  double nu_1 = paramList["nu_1"];
  rowvec mu = stateVecRed.row(0);
  rowvec mu_t = mu.tail(nPeriods - 1);
  rowvec mu_lt = mu.head(nPeriods - 1);
  rowvec nu_0 = stateVecRed.row(1).tail(nPeriods - 1);
  vec Drift_S_t = regimeDriftVecRed.tail(nPeriods - 1);
  vec Heterosc_S_t = regimeHeteroscVecRed.tail(nPeriods - 1);
  // Get multiplicative specification of the two variances
  List xi = paramList["xi"];
  double xi_0_sq = pow(xi["xi_0"], 2);
  double xi_1_sq = pow(xi["xi_1"], 2);
  double x = xi_1_sq / xi_0_sq - xi_0_sq;
  vec correcFactor = 1 / sqrt(1 + x * Heterosc_S_t);
  // Compute the SSR from the regression of the trend on its lag and drift terms, corrected for S_t = 1
  rowvec resid = mu_t - mu_lt - nu_0 - Drift_S_t.t() * nu_1;
  rowvec v_t_1 = resid % correcFactor.t();
  double SSR_1 = accu(square(v_t_1));
  // Adjust the priors
  double alpha_1 = alpha_0 + nPeriods;
  double beta_1 = beta_0 + SSR_1;
  // Constrain the parameter space
  double lowerLimit = 1e-5;
  double upperLimit = math::inf();
  // Draw the new value xi_0
  double xi_0_var = 1 / randg(distr_param(alpha_1 / 2.0, 2.0 / beta_1));
  double xi_0_candidate = sqrt(xi_0_var);
  double xi_0_draw;
  if (xi_0_candidate < lowerLimit || xi_0_candidate > upperLimit)
  {
    xi_0_draw = xi["xi_0"];
  }
  else
  {
    xi_0_draw = sqrt(xi_0_var);
  }
  // Compute the SSR from the regression of the trend on its lag and drift terms, corrected for (heterosc.) S_t = 0
  vec v_t_2 = (Heterosc_S_t % resid.t()) / xi_0_draw;
  double SSR_2 = accu(square(v_t_2));
  // Adjust the priors
  double alpha_3 = alpha_2 + accu(Heterosc_S_t);
  double beta_3 = beta_2 + SSR_2;
  // Draw the new value xi_0
  double x_draw = 1 / randg(distr_param(alpha_3 / 2.0, 2.0 / beta_3));
  double xi_1_candidate = sqrt(xi_0_var * x_draw);
  double xi_1_draw;
  if (xi_1_candidate < lowerLimit || xi_1_candidate > upperLimit)
  {
    xi_1_draw = xi["xi_1"];
  }
  else
  {
    xi_1_draw = xi_1_candidate;
  }
  // Construct the output object
  vec xiVec = {xi_0_draw, xi_1_draw};
  return xiVec;
}
