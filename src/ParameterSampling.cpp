#include <RcppArmadillo.h>
#include <cmath>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' @description Draws the transition probabilities for models that are governed by one regime process
//' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @param secondRegime boolean. If true, the variable names r and s are used instead of q and p
//' @return the filtered output
// [[Rcpp::export]]
rowvec TransProbs_fctn(vec regimeVec, List paramList, bool secondRegime = false)
{
  List Probs = paramList["Probs"];
  double Backup1;
  double Backup2;
  if (secondRegime == false)
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
  // double u_00 = 300;
  // double u_11 = 100;
  double u_00 = 0;
  double u_11 = 0;
  // double u_01 = 0;
  // double u_10 = 0;
  double u_01 = 50;
  double u_10 = 50;
  // Draw the sample
  q = R::rbeta(u_00 + n_00, u_01 + n_01);
  int threshold = 1e04;
  int counter;
  double probLimit = .9;
  counter = 0;
  while (q < probLimit && counter <= threshold)
  {
    counter++;
    q = R::rbeta(u_00 + n_00, u_01 + n_01);
  };
  if (q < probLimit)
  {
    q = Backup2;
  }

  p = R::rbeta(u_11 + n_11, u_10 + n_10);
  counter = 0;
  while (p < probLimit && counter <= threshold)
  {
    counter++;
    p = R::rbeta(u_11 + n_11, u_10 + n_10);
  };
  if (p < probLimit)
  {
    p = Backup2;
  }
  oututVec = {q, p};
  return oututVec;
}

//' @description Draws the added negative drift for regime 1 for models governed by one regime processes
//' @param stateVec array with draws of the state vector
//' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @return the filtered output
// [[Rcpp::export]]
double Nu_1_fctn(mat stateVec, vec regimeVec, List paramList)
{
  double nu_1_backup = paramList["nu_1"];
  // discard first 10 observations because of the high state vector variation due to the filter initialization
  mat stateVecRed = stateVec.cols(10, stateVec.n_cols - 1);
  vec regimeVecRed = regimeVec.subvec(10, regimeVec.n_elem - 1);
  List Xi = paramList["xi"];
  double nu_1Backup = paramList["nu_1"];
  double SdXi = Xi["xi_1"];
  // Set up the priors
  double zeta_0 = 0;
  double delta_0 = 1e4;
  // Compute the SSR from the regression of of the trend on its lag and drift terms, adjusted for xi so that the error is standard normal
  rowvec mu = stateVecRed.row(0);
  rowvec mu_t = mu.subvec(1, mu.size() - 1);
  rowvec mu_lt = mu.subvec(0, mu.size() - 2);
  rowvec nu_0 = stateVecRed.row(1).subvec(1, mu.size() - 1);
  mat Y = mat((mu_t - mu_lt - nu_0) / SdXi).t();
  mat X = mat(regimeVecRed.subvec(1, mu.size() - 1) / SdXi);
  mat transX = trans(X);
  double invXtX = as_scalar(inv(1 / delta_0 + transX * X));
  double zeta_1 = invXtX * as_scalar(transX * Y);
  if (zeta_1 > 0)
  {
    zeta_1 = 0;
  }
  double delta_1 = invXtX;
  double draw;
  draw = zeta_1 + as_scalar(randn(1, 1)) * delta_1;
  // if none out of 1000 draws produces a negative value, the function does not return an output
  int counter = 1;
  int threshold = 1e3;
  while (draw >= 0 && counter <= threshold)
  {
    draw = zeta_1 + as_scalar(randn(1, 1)) * delta_1;
    counter++;
  }
  if (counter != threshold)
  {
    return draw;
  }
  else
  {
    return nu_1_backup;
  }
}

//' @description Draws the added negative drift for regime 1 for models governed by two regime processes
//' @param stateVec array with draws of the state vector
//' @param regimeDriftVec vector with regime realizations governing the drift. Vector takes either values 1 or 0
//' @param regimeHeteroscVec vector with regime realizations governing the heteroscedasticity. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @return the filtered output
// [[Rcpp::export]]
double Nu_1_TwoRegime_fctn(mat stateVec, vec regimeDriftVec, vec regimeHeteroscVec, List paramList)
{
  double nu_1_backup = paramList["nu_1"];
  // discard first 10 observations because of the high state vector variation due to the filter initialization
  mat stateVecRed = stateVec.cols(10, stateVec.n_cols - 1);
  vec regimeDriftVecRed = regimeDriftVec.subvec(10, regimeDriftVec.n_elem - 1);
  vec regimeHeteroscVecRed = regimeHeteroscVec.subvec(10, regimeHeteroscVec.n_elem - 1);
  List Xi = paramList["xi"];
  double nu_1Backup = paramList["nu_1"];
  double SdXi_0 = Xi["xi_0"];
  double SdXi_1 = Xi["xi_1"];
  // Set up the priors
  double zeta_0 = 0;
  double delta_0 = 1e4;
  // Compute the SSR from the regression of of the trend on its lag and drift terms, adjusted for xi so that the error is standard normal
  rowvec mu = stateVecRed.row(0);
  rowvec mu_t = mu.subvec(1, mu.size() - 1);
  rowvec mu_lt = mu.subvec(0, mu.size() - 2);
  rowvec nu_0 = stateVecRed.row(1).subvec(1, mu.size() - 1);
  rowvec SdXiVec = regimeHeteroscVecRed.subvec(1, mu.size() - 1) * SdXi_1;
  SdXiVec.elem(find(SdXiVec == 0)).fill(SdXi_0);
  mat Y = mat((mu_t - mu_lt - nu_0) / SdXiVec).t();
  mat X = mat(regimeDriftVecRed.subvec(1, mu.size() - 1) / SdXiVec);
  mat transX = trans(X);
  double invXtX = as_scalar(inv(1 / delta_0 + transX * X));
  double zeta_1 = invXtX * as_scalar(transX * Y);
  if (zeta_1 > 0)
  {
    zeta_1 = 0;
  }
  double delta_1 = invXtX;
  double draw;
  draw = zeta_1 + as_scalar(randn(1, 1)) * delta_1;
  // if none out of 1000 draws produces a negative value, the function does not return an output
  int counter = 1;
  int threshold = 1e3;
  while (draw >= 0 && counter <= threshold)
  {
    draw = zeta_1 + as_scalar(randn(1, 1)) * delta_1;
    counter++;
  }
  if (counter != threshold)
  {
    return draw;
  }
  else
  {
    return nu_1_backup;
  }
}