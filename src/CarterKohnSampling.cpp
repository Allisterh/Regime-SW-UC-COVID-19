#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' @description Function that generates random values from a multivariate normal distribution
//' @param mean vector of means
//' @param cov positive semi definite variance covariance matrix
//' @return vector with random values

vec rnormVec_fctn(int n, vec mu, mat sigma)
{
  int ncols = sigma.n_cols;
  mat iid = randn(n, ncols);
  mat draws = repmat(mu, 1, n).t() + iid * chol(sigma);
  return draws.row(0).t();
}

//' @description Function that returns the position of zeros on the diagonal of a diagonal matrix
//' @param A diagonal matrix
//' @return vector with the positions
uvec zeroDiagPos_fctn(mat A)
{
  uvec B;
  for (unsigned int i = 0; i < A.n_cols; i++)
  {
    if (A(i, i) != 0)
    {
      B = join_cols(B, uvec({i}));
    }
  }
  return B;
}

//' @description Normal pdf
//' @param double x a quantile
//' @returns double the pdf

double pnorm_fctn(double x)
{
  double dist = R::pnorm5(x, 0.0, 1.0, 1, 0);
  return dist;
}

//' @description Function to execute the Carter Kohn forward filtering backwards sampling algorithm for the state vector
//' @param filterOutput list with the updated state vector and var-cov matrix of the kalman filter
//' @param systemList list with the system matrices
//' @param paramList list with the additional parameters
//' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
//' @return the sampled state vectors
// [[Rcpp::export]]
mat BackwardsStateSampling_fctn(List filterOutput, List systemList, List paramList, vec regimeVec)
{
  // unpack the filter output
  mat a_t_ct_mat = filterOutput["a_ct"];
  cube P_t_ct_array = filterOutput["P_ct"];
  // System matrices
  mat Tt = systemList["Tt"];
  mat R = systemList["R"];
  // Pre transpose system matrices
  mat transTt = trans(Tt);
  mat transR = trans(R);

  // dimension of state vector and number of time nPeriods
  int dimens = Tt.n_cols;
  int nPeriods = a_t_ct_mat.n_cols;

  // Var-Cov Matrix of innovations in transition equation
  List xi = paramList["xi"];
  double xi_0 = xi["xi_0"];
  double xi_1 = xi["xi_1"];
  List omega = paramList["omega"];
  double omega_0 = omega["omega_0"];
  double omega_1 = omega["omega_1"];
  List eta = paramList["eta"];
  double eta_0 = eta["eta_0"];
  double eta_1 = eta["eta_1"];
  vec Q_0Vec = {pow(xi_0, 2), pow(eta_0, 2), pow(omega_0, 2)};
  vec Q_1Vec = {pow(xi_1, 2), pow(eta_1, 2), pow(omega_1, 2)};
  mat Q_0 = diagmat(Q_0Vec);
  mat Q_1 = diagmat(Q_1Vec);
  mat expandQ_0 = R * Q_0 * transR;
  mat expandQ_1 = R * Q_1 * transR;

  // Vector holding all the positions where R * Q * R is not singular
  uvec singularPos = zeroDiagPos_fctn(expandQ_0);
  uvec singularPos1 = zeroDiagPos_fctn(expandQ_1);

  // Extract the parts of the state vector where it is not deterministic given the last iteration
  mat expandQ_0_star = expandQ_0.submat(singularPos, singularPos);
  mat expandQ_1_star = expandQ_1.submat(singularPos, singularPos);
  mat Tt_star = Tt.rows(singularPos);
  mat transTt_star = trans(Tt_star);

  // Vector holding additional Drift in Regime 1
  double nu_1 = paramList["nu_1"];
  vec lambdaRaw(dimens);
  lambdaRaw.fill(0);
  vec lambdaNullVec = lambdaRaw.elem(singularPos);
  lambdaRaw(0) = nu_1;
  vec lambdaVec = lambdaRaw.elem(singularPos);

  // Initializes output objects
  mat a_t_draw_mat(dimens, nPeriods, fill::zeros);
  // and the sampler
  vec a_t1_cT = a_t_ct_mat.col(nPeriods - 1);
  mat P_t1_cT = P_t_ct_array.slice(nPeriods - 1);
  vec a_T_draw = rnormVec_fctn(1, a_t1_cT, P_t1_cT);
  a_t_draw_mat.col(nPeriods - 1) = a_T_draw;
  vec a_t_draw = a_T_draw.elem(singularPos);

  // Initialize further values
  vec a_t_ct;
  mat P_t_ct;
  vec lambda;
  mat Trans_P_t_ct_expand;
  mat Inv_Helper;
  mat P_t1_cT_sym;
  vec a_t_draw_output;
  vec a_t1_cT_star;
  mat P_t1_cT_star;

  // Execute the backwards recursions
  for (int i = nPeriods - 2; i >= 0; i--)
  {
    // Pick the updated values from the Kalman filter
    a_t_ct = a_t_ct_mat.col(i);
    P_t_ct = P_t_ct_array.slice(i);

    // Select the system matrices based on the regime
    Trans_P_t_ct_expand = P_t_ct * transTt_star;
    if (regimeVec[i] == 1)
    {
      Inv_Helper = inv(Tt_star * Trans_P_t_ct_expand + expandQ_1_star);
      lambda = lambdaVec;
    }
    else
    {
      Inv_Helper = inv(Tt_star * Trans_P_t_ct_expand + expandQ_0_star);
      lambda = lambdaNullVec;
    }
    // Compute the backwards updates
    a_t1_cT = a_t_ct + Trans_P_t_ct_expand * Inv_Helper * (a_t_draw - lambda - Tt_star * a_t_ct);
    P_t1_cT = P_t_ct - Trans_P_t_ct_expand * Inv_Helper * Tt_star * P_t_ct;
    // To prevent errors
    P_t1_cT = symmatl(P_t1_cT);
    // Draw the state vector from a multivariate normal
    a_t1_cT_star = a_t1_cT.elem(singularPos);
    P_t1_cT_star = P_t1_cT.submat(singularPos, singularPos);
    a_t_draw = rnormVec_fctn(1, a_t1_cT_star, P_t1_cT_star);

    // The output vector consists of the randomly drawn state vector and the updated state vector from the Kalman filter
    // for positions where the state vector is deterministic given the last iteration
    a_t_draw_output = a_t_ct;
    a_t_draw_output.elem(singularPos) = a_t_draw;
    // Store the sampled state vectors
    a_t_draw_mat.col(i) = a_t_draw_output;
  }
  return a_t_draw_mat;
}

//' @description Function to execute the Carter Kohn forward filtering backwards sampling algorithm for the state vector
//' and models that are governed by two regime processes
//' @param filterOutput list with the updated state vector and var-cov matrix of the kalman filter
//' @param systemList list with the system matrices
//' @param paramList list with the additional parameters
//' @param regimeDriftVec vector with regime realizations governing the trend drift. Vector takes either values 1 or 0
//' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
//' @return the sampled state vectors
// [[Rcpp::export]]
mat BackwardsTwoRegimeStateSampling_fctn(List filterOutput, List systemList, List paramList, vec regimeDriftVec, vec regimeHeteroscVec)
{
  // unpack the filter output
  mat a_t_ct_mat = filterOutput["a_ct"];
  cube P_t_ct_array = filterOutput["P_ct"];
  // System matrices
  mat Tt = systemList["Tt"];
  mat R = systemList["R"];
  // Pre transpose system matrices
  mat transTt = trans(Tt);
  mat transR = trans(R);

  // dimension of state vector and number of time nPeriods
  int dimens = Tt.n_cols;
  int nPeriods = a_t_ct_mat.n_cols;

  // Var-Cov Matrix of innovations in transition equation
  List xi = paramList["xi"];
  double xi_0 = xi["xi_0"];
  double xi_1 = xi["xi_1"];
  List omega = paramList["omega"];
  double omega_0 = omega["omega_0"];
  double omega_1 = omega["omega_1"];
  List eta = paramList["eta"];
  double eta_0 = eta["eta_0"];
  double eta_1 = eta["eta_1"];
  vec Q_0Vec = {pow(xi_0, 2), pow(eta_0, 2), pow(omega_0, 2)};
  vec Q_1Vec = {pow(xi_1, 2), pow(eta_1, 2), pow(omega_1, 2)};
  mat Q_0 = diagmat(Q_0Vec);
  mat Q_1 = diagmat(Q_1Vec);
  mat expandQ_0 = R * Q_0 * transR;
  mat expandQ_1 = R * Q_1 * transR;

  // Vector holding all the positions where R * Q * R is not singular
  uvec singularPos = zeroDiagPos_fctn(expandQ_0);
  uvec singularPos1 = zeroDiagPos_fctn(expandQ_1);

  // Extract the parts of the state vector where it is not deterministic given the last iteration
  mat expandQ_0_star = expandQ_0.submat(singularPos, singularPos);
  mat expandQ_1_star = expandQ_1.submat(singularPos, singularPos);
  mat Tt_star = Tt.rows(singularPos);
  mat transTt_star = trans(Tt_star);

  // Vector holding additional Drift in Regime 1
  double nu_1 = paramList["nu_1"];
  vec lambdaRaw(dimens);
  lambdaRaw.fill(0);
  vec lambdaNullVec = lambdaRaw.elem(singularPos);
  lambdaRaw(0) = nu_1;
  vec lambdaVec = lambdaRaw.elem(singularPos);

  // Initializes output objects
  mat a_t_draw_mat(dimens, nPeriods, fill::zeros);
  // and the sampler
  vec a_t1_cT = a_t_ct_mat.col(nPeriods - 1);
  mat P_t1_cT = P_t_ct_array.slice(nPeriods - 1);
  vec a_T_draw = rnormVec_fctn(1, a_t1_cT, P_t1_cT);
  a_t_draw_mat.col(nPeriods - 1) = a_T_draw;
  vec a_t_draw = a_T_draw.elem(singularPos);

  // Initialize further values
  vec a_t_ct;
  mat P_t_ct;
  vec lambda;
  mat Trans_P_t_ct_expand;
  mat Inv_Helper;
  mat P_t1_cT_sym;
  vec a_t_draw_output;
  vec a_t1_cT_star;
  mat P_t1_cT_star;

  // Execute the backwards recursions
  for (int i = nPeriods - 2; i >= 0; i--)
  {
    // Pick the updated values from the Kalman filter
    a_t_ct = a_t_ct_mat.col(i);
    P_t_ct = P_t_ct_array.slice(i);

    // Select the system matrices based on the regime
    Trans_P_t_ct_expand = P_t_ct * transTt_star;
    if (regimeDriftVec[i] == 1)
    {
      lambda = lambdaVec;
    }
    else
    {
      lambda = lambdaNullVec;
    }
    if (regimeHeteroscVec[i] == 1)
    {
      Inv_Helper = inv(Tt_star * Trans_P_t_ct_expand + expandQ_1_star);
    }
    else
    {
      Inv_Helper = inv(Tt_star * Trans_P_t_ct_expand + expandQ_0_star);
    }
    // Compute the backwards updates
    a_t1_cT = a_t_ct + Trans_P_t_ct_expand * Inv_Helper * (a_t_draw - lambda - Tt_star * a_t_ct);
    P_t1_cT = P_t_ct - Trans_P_t_ct_expand * Inv_Helper * Tt_star * P_t_ct;
    // To prevent errors
    P_t1_cT = symmatl(P_t1_cT);
    // Draw the state vector from a multivariate normal
    a_t1_cT_star = a_t1_cT.elem(singularPos);
    P_t1_cT_star = P_t1_cT.submat(singularPos, singularPos);
    a_t_draw = rnormVec_fctn(1, a_t1_cT_star, P_t1_cT_star);

    // The output vector consists of the randomly drawn state vector and the updated state vector from the Kalman filter
    // for positions where the state vector is deterministic given the last iteration
    a_t_draw_output = a_t_ct;
    a_t_draw_output.elem(singularPos) = a_t_draw;
    // Store the sampled state vectors
    a_t_draw_mat.col(i) = a_t_draw_output;
  }
  return a_t_draw_mat;
}

//' @description Executes the hamilton filter for a given state vector
//' @param stateVec array with draws of the state vector
//' @param paramList list with the additional parameters
//' @param endogen boolean. If TRUE, the function includes endogeneity
//' @return a matrix with the sampled regime realizations
// [[Rcpp::export]]
vec RegimeSampling_fctn(mat filterOutput, List paramList, bool endogen = false)
{
  double p;
  double q;
  double beta_0;
  double beta_1;
  double varrho;
  if (endogen == false)
  {
    // Setting up the transition probabilities
    List Probs = paramList["Probs"];
    p = Probs["p"];
    q = Probs["q"];
  }
  else
  {
    List beta = paramList["beta"];
    beta_0 = beta["beta_0"];
    beta_1 = beta["beta_1"];
    varrho = paramList["varrho"];
    // Transfer probit coefficients into probabilities
    p = 1 - pnorm_fctn(-beta_0 - beta_1);
    q = pnorm_fctn(-beta_0);
  }

  // Length of the observational period
  int nPeriods = filterOutput.n_rows;

  // Initialize output object
  vec S_t_vec(nPeriods, fill::zeros);
  // Initialize the routine
  double Pr_t1_cT_1 = filterOutput(nPeriods - 1, 1);
  double S_t_prob_draw = runif(1)(0);
  int S_t_draw;
  if (S_t_prob_draw <= Pr_t1_cT_1)
  {
    S_t_draw = 1;
  }
  else
  {
    S_t_draw = 0;
  }
  S_t_vec[nPeriods - 1] = S_t_draw;

  // Initialize further values
  double prob;

  for (int i = nPeriods - 2; i >= 0; i--)
  {
    // Pull the filtered probability
    double Pr_ct_1 = filterOutput(i, 1);
    // Pick the relevant transition probabilities
    vec P_vec;
    if (S_t_draw == 0)
    {
      P_vec = {q, 1 - q};
    }
    else
    {
      P_vec = {1 - p, p};
    }
    // Calculate the conditional probability for the next iteration
    Pr_t1_cT_1 = (P_vec[1] * Pr_ct_1) / (P_vec[1] * Pr_ct_1 + P_vec[0] * (1 - Pr_ct_1));

    // Draw a regime probability
    prob = runif(1)(0);
    if (prob <= Pr_t1_cT_1)
    {
      S_t_draw = 1;
    }
    else
    {
      S_t_draw = 0;
    }
    // Record the draw
    S_t_vec[i] = S_t_draw;
  }
  return S_t_vec;
}

// [[Rcpp::export]]
mat OldBackwardsStateSampling_fctn(List filterOutput, List systemList, List paramList, vec regimeVec)
{
  // unpack the filter output
  mat a_t_ct_mat = filterOutput["a_ct"];
  cube P_t_ct_array = filterOutput["P_ct"];
  // System matrices
  mat Tt = systemList["Tt"];
  mat R = systemList["R"];
  // Pre transpose system matrices
  mat transTt = trans(Tt);
  mat transR = trans(R);

  // dimension of state vector and number of time nPeriods
  int dimens = Tt.n_cols;
  int nPeriods = a_t_ct_mat.n_cols;

  // Var-Cov Matrix of innovations in transition equation
  List xi = paramList["xi"];
  double xi_0 = xi["xi_0"];
  double xi_1 = xi["xi_1"];
  List omega = paramList["omega"];
  double omega_0 = omega["omega_0"];
  double omega_1 = omega["omega_1"];
  List eta = paramList["eta"];
  double eta_0 = eta["eta_0"];
  double eta_1 = eta["eta_1"];
  vec Q_0Vec = {pow(xi_0, 2), pow(eta_0, 2), pow(omega_0, 2)};
  vec Q_1Vec = {pow(xi_1, 2), pow(eta_1, 2), pow(omega_1, 2)};
  mat Q_0 = diagmat(Q_0Vec);
  mat Q_1 = diagmat(Q_1Vec);
  mat expandQ_0 = R * Q_0 * transR;
  mat expandQ_1 = R * Q_1 * transR;

  // Vector holding all the positions where R * Q * R is not singular
  uvec singularPos = zeroDiagPos_fctn(expandQ_0);
  uvec singularPos1 = zeroDiagPos_fctn(expandQ_1);

  // Extract the parts of the state vector where it is not deterministic given the last iteration
  mat expandQ_0_star = expandQ_0.submat(singularPos, singularPos);
  mat expandQ_1_star = expandQ_1.submat(singularPos, singularPos);
  mat Tt_star = Tt.rows(singularPos);
  mat transTt_star = trans(Tt_star);

  // Vector holding additional Drift in Regime 1
  double nu_1 = paramList["nu_1"];
  vec lambdaRaw(dimens);
  lambdaRaw.fill(0);
  vec lambdaNullVec = lambdaRaw.elem(singularPos);
  lambdaRaw(0) = nu_1;
  vec lambdaVec = lambdaRaw.elem(singularPos);

  // Initializes output objects
  mat a_t_draw_mat(dimens, nPeriods, fill::zeros);
  // and the sampler
  vec a_t1_cT = a_t_ct_mat.col(nPeriods - 1);
  mat P_t1_cT = P_t_ct_array.slice(nPeriods - 1);
  vec a_T_draw = rnormVec_fctn(1, a_t1_cT, P_t1_cT);
  a_t_draw_mat.col(nPeriods - 1) = a_T_draw;
  vec a_t_draw = a_T_draw.elem(singularPos);

  // Initialize further values
  vec a_t_ct;
  mat P_t_ct;
  vec lambda;
  mat Trans_P_t_ct_expand;
  mat Inv_Helper;
  mat P_t1_cT_sym;
  vec a_t_draw_output;
  vec a_t1_cT_star;
  mat P_t1_cT_star;

  // Execute the backwards recursions
  for (int i = nPeriods - 2; i >= 0; i--)
  {
    // Pick the updated values from the Kalman filter
    a_t_ct = a_t_ct_mat.col(i);
    P_t_ct = P_t_ct_array.slice(i);

    // Select the system matrices based on the regime
    Trans_P_t_ct_expand = P_t_ct * transTt_star;
    if (regimeVec[i] == 1)
    {
      Inv_Helper = inv(Tt_star * Trans_P_t_ct_expand + expandQ_1_star);
      lambda = lambdaVec;
    }
    else
    {
      Inv_Helper = inv(Tt_star * Trans_P_t_ct_expand + expandQ_0_star);
      lambda = lambdaNullVec;
    }
    // Compute the backwards updates
    a_t1_cT = a_t_ct + Trans_P_t_ct_expand * Inv_Helper * (a_t_draw - lambda - Tt_star * a_t_ct);
    P_t1_cT = P_t_ct - Trans_P_t_ct_expand * Inv_Helper * Tt_star * P_t_ct;
    // To prevent errors
    P_t1_cT = symmatl(P_t1_cT);
    // Draw the state vector from a multivariate normal
    a_t1_cT_star = a_t1_cT.elem(singularPos);
    P_t1_cT_star = P_t1_cT.submat(singularPos, singularPos);
    a_t_draw = rnormVec_fctn(1, a_t1_cT_star, P_t1_cT_star);

    // The output vector consists of the randomly drawn state vector and the updated state vector from the Kalman filter
    // for positions where the state vector is deterministic given the last iteration
    a_t_draw_output = a_t_ct;
    a_t_draw_output.elem(singularPos) = a_t_draw;
    // Store the sampled state vectors
    a_t_draw_mat.col(i) = a_t_draw_output;
  }
  return a_t_draw_mat;
}

// [[Rcpp::export]]
vec OldRegimeSampling_fctn(mat filterOutput, List paramList, bool endogen)
{
  double p;
  double q;
  double beta_0;
  double beta_1;
  double varrho;
  if (endogen == false)
  {
    // Setting up the transition probabilities
    List Probs = paramList["Probs"];
    p = Probs["p"];
    q = Probs["q"];
  }
  else
  {
    List beta = paramList["beta"];
    beta_0 = beta["beta_0"];
    beta_1 = beta["beta_1"];
    varrho = paramList["varrho"];
    // Transfer probit coefficients into probabilities
    p = 1 - pnorm_fctn(-beta_0 - beta_1);
    q = pnorm_fctn(-beta_0);
  }

  // Length of the observational period
  int nPeriods = filterOutput.n_rows;

  // Initialize output object
  vec S_t_vec(nPeriods, fill::zeros);
  // Initialize the routine
  double Pr_t1_cT_1 = filterOutput(nPeriods - 1, 1);
  double S_t_prob_draw = runif(1)(0);
  int S_t_draw;
  if (S_t_prob_draw <= Pr_t1_cT_1)
  {
    S_t_draw = 1;
  }
  else
  {
    S_t_draw = 0;
  }
  S_t_vec[nPeriods - 1] = S_t_draw;

  // Initialize further values
  double prob;

  for (int i = nPeriods - 2; i >= 0; i--)
  {
    // Pull the filtered probability
    double Pr_ct_1 = filterOutput(i, 1);
    // Pick the relevant transition probabilities
    vec P_vec;
    if (S_t_draw == 0)
    {
      P_vec = {q, 1 - q};
    }
    else
    {
      P_vec = {1 - p, p};
    }
    // Calculate the conditional probability for the next iteration
    Pr_t1_cT_1 = (P_vec[1] * Pr_ct_1) / (P_vec[1] * Pr_ct_1 + P_vec[0] * (1 - Pr_ct_1));

    // Draw a regime probability
    prob = runif(1)(0);
    if (prob <= Pr_t1_cT_1)
    {
      S_t_draw = 1;
    }
    else
    {
      S_t_draw = 0;
    }
    // Record the draw
    S_t_vec[i] = S_t_draw;
  }
  return S_t_vec;
}