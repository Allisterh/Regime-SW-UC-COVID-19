#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' @description Normal pdf
//' @param double x a quantile
//' @returns double the pdf

double pnorm_fctn(double x)
{
   double dist = R::pnorm5(x, 0.0, 1.0, 1, 0);
   return dist;
}

//' @description Derives the normal density
//' @param x diagonal mat with the quantile values
//' @param var diagonal mat with the variances
//' @returns row vector with the densities

double dnorm_fctn(double x, double sd)
{
   double dens = R::dnorm4(x, 0, sd, FALSE);
   double densFloor = fmax(dens, 1e-10);
   return densFloor;
}

//' @description Function that computes the state vector via the Kalman filter given a vector of drawn regime realizations and additional parameters
//' @param data vector of the data
//' @param Ini initial value for the trend component
//' @param systemList list with the system matrices
//' @param paramList list with the additional parameters
//' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
//' @param endogen logical. If true, the function includes endogeneity
//' @return the filtered output
// [[Rcpp::export]]
List CondKalman_fctn(vec data, double Ini, List systemList, List paramList, vec regimeVec, bool endogen = false)
{
   // System matrices
   mat Tt = systemList["Tt"];
   mat Z = systemList["Z"];
   mat R = systemList["R"];
   // Pre transpose system matrices
   mat transTt = trans(Tt);
   mat transR = trans(R);
   mat transZ = trans(Z);

   // dimension of state vector and number of time nPeriods
   int dimens = Tt.n_cols;
   int nPeriods = data.n_elem;

   // Residual variance in the measurement eq. 1e-6 assures a positive definite var-cov variance of the state vector
   double epsilon = paramList["epsilon"];
   double epsilonSq = fmax(pow(epsilon, 2), 1e-6);

   // Vector holding additional Drift in Regime 1
   double nu_1 = paramList["nu_1"];
   vec lambda(dimens);
   lambda.fill(0);
   lambda(0) = nu_1;

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

   // Initial values for state vector
   vec a_t_clt = zeros(dimens);
   a_t_clt(0) = Ini;
   if (regimeVec[0] == 1)
   {
      a_t_clt(0) = Ini + nu_1;
   }
   // Initial variance matrix for state vector (diffuse initialization)
   vec P_t_vec(dimens, fill::value(1e3));
   mat P_t_clt = diagmat(P_t_vec);

   // Initializes output objects
   mat a_ct_mat(dimens, nPeriods, fill::zeros);
   cube P_ct_array(dimens, dimens, nPeriods, fill::zeros);
   vec v_t_vec(nPeriods, fill::zeros);
   vec F_t_vec(nPeriods, fill::zeros);

   // Initialize further values for the filter
   double v_t;
   double F_t;
   vec a_t_ct;
   mat P_t_ct;
   mat K_t;

   for (int i = 0; i < nPeriods; i++)
   {
      //---------------------//
      // Kalman part 1/2     //
      //---------------------//

      // One step ahead prediction error with state vector prediction from
      // (t|t-1)
      v_t = (data[i] - Z * a_t_clt)[0];
      v_t_vec[i] = v_t;
      // Variance of prediction error
      F_t = (Z * P_t_clt * transZ + epsilonSq)[0];
      F_t_vec[i] = F_t;

      // Updating step
      // (t|t, S_t, S_t-1)
      // Kalman gain
      K_t = transZ * (1 / F_t);
      // State vector
      a_t_ct = a_t_clt + P_t_clt * K_t * v_t;
      // State vector Var-Matrix
      P_t_ct = P_t_clt - P_t_clt * K_t * Z * P_t_clt;

      //---------------------//
      // Kalman part 2/2     //
      //---------------------//

      // Store updated values
      // (t|t)
      a_ct_mat.col(i) = a_t_ct;
      P_ct_array.slice(i) = P_t_ct;

      // pull the regime realization for the next period
      if (i < nPeriods - 1)
      {
         if (regimeVec[i + 1] == 0)
         {
            // Prediction step with approximated updates to complete loop
            // (t+1|t)
            // Regime 0 (high drift)
            a_t_clt = Tt * a_t_ct;
            P_t_clt = Tt * P_t_ct * transTt + expandQ_0;
         }
         else
         {
            // Regime 1 (low drift)
            a_t_clt = Tt * a_t_ct + lambda;
            P_t_clt = Tt * P_t_ct * transTt + expandQ_1;
         }
      }
   }

   List Output = List::create(
       Named("a_ct") = a_ct_mat, // updated state vector
       _["P_ct"] = P_ct_array,
       _["v_t_vec"] = v_t_vec,
       _["F_t_vec"] = F_t_vec // updated state vector var
   );
   return Output;
}

//' @description Function that computes the state vector via the Kalman filter given additional parameters and two vectors of drawn regime realizations
//' @param data vector of the data
//' @param Ini initial value for the trend component
//' @param systemList list with the system matrices
//' @param paramList list with the additional parameters
//' @param regimeDriftVec vector with regime realizations governing the trend drift. Vector takes either values 1 or 0
//' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
//' @param endogen logical. If true, the function includes endogeneity
//' @return the filtered output
// [[Rcpp::export]]
List CondTwoRegimeKalman_fctn(vec data, double Ini, List systemList, List paramList, vec regimeDriftVec, vec regimeHeteroscVec, bool endogen = false)
{
   // System matrices
   mat Tt = systemList["Tt"];
   mat Z = systemList["Z"];
   mat R = systemList["R"];
   // Pre transpose system matrices
   mat transTt = trans(Tt);
   mat transR = trans(R);
   mat transZ = trans(Z);

   // dimension of state vector and number of time nPeriods
   int dimens = Tt.n_cols;
   int nPeriods = data.n_elem;

   // Residual variance in the measurement eq
   double epsilon = fmax(paramList["epsilon"], 1e-3);
   double epsilonSq = pow(epsilon, 2);

   // Vector holding additional Drift in Regime 1
   double nu_1 = paramList["nu_1"];
   vec lambda(dimens);
   lambda.fill(0);
   lambda(0) = nu_1;

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

   // Initial values for state vector
   vec a_t_clt = zeros(dimens);
   a_t_clt(0) = Ini;
   if (regimeDriftVec[0] == 1)
   {
      a_t_clt(0) = Ini + nu_1;
   }
   // Initial variance matrix for state vector (Exact diffuse initialization)
   vec P_t_vec(dimens, fill::value(1000));
   mat P_t_clt = diagmat(P_t_vec);

   // Initializes output objects
   mat a_ct_mat(dimens, nPeriods, fill::zeros);
   cube P_ct_array(dimens, dimens, nPeriods, fill::zeros);
   vec v_t_vec(nPeriods, fill::zeros);
   vec F_t_vec(nPeriods, fill::zeros);

   // Initialize further values for the filter
   double v_t;
   double F_t;
   vec a_t_ct;
   mat P_t_ct;
   mat K_t;

   for (int i = 0; i < nPeriods; i++)
   {
      //---------------------//
      // Kalman part 1/2     //
      //---------------------//

      // One step ahead prediction error with state vector prediction from
      // (t|t-1)
      v_t = (data[i] - Z * a_t_clt)[0];
      v_t_vec[i] = v_t;
      // Variance of prediction error
      F_t = (Z * P_t_clt * transZ + epsilonSq)[0];
      F_t_vec[i] = F_t;

      // Updating step
      // (t|t, S_t, S_t-1)
      // Kalman gain
      K_t = transZ * (1 / F_t);
      // State vector
      a_t_ct = a_t_clt + P_t_clt * K_t * v_t;
      // State vector Var-Matrix
      P_t_ct = P_t_clt - P_t_clt * K_t * Z * P_t_clt;

      //---------------------//
      // Kalman part 2/2     //
      //---------------------//

      // Store updated values
      // (t|t)
      a_ct_mat.col(i) = a_t_ct;
      P_ct_array.slice(i) = P_t_ct;

      if (i < nPeriods - 1)
      {
         // pull the regime realization for the next period
         if (regimeDriftVec[i + 1] == 0)
         {
            // Prediction step with approximated updates to complete loop
            // (t+1|t)
            // Regime with high drift and xi_0
            if (regimeHeteroscVec[i + 1] == 0)
            {
               a_t_clt = Tt * a_t_ct;
               P_t_clt = Tt * P_t_ct * transTt + expandQ_0;
            }
            else
            // Regime with high drift and xi_1
            {
               a_t_clt = Tt * a_t_ct;
               P_t_clt = Tt * P_t_ct * transTt + expandQ_1;
            }
         }
         else
         // Regime 1 (low drift)
         {
            if (regimeHeteroscVec[i + 1] == 0)
            {
               // Regime with low drift and xi_0
               a_t_clt = Tt * a_t_ct + lambda;
               P_t_clt = Tt * P_t_ct * transTt + expandQ_0;
            }
            else
            {
               // Regime with low drift and xi_1
               a_t_clt = Tt * a_t_ct + lambda;
               P_t_clt = Tt * P_t_ct * transTt + expandQ_1;
            }
         }
      }
   }

   List Output = List::create(
       Named("a_ct") = a_ct_mat, // updated state vector
       _["P_ct"] = P_ct_array,
       _["v_t_vec"] = v_t_vec,
       _["F_t_vec"] = F_t_vec // updated state vector var
   );
   return Output;
}

//' @description Executes the hamilton filter for a given state vector
//' @param stateVec array with draws of the state vector
//' @param paramList list with the additional parameters
//' @param endogen boolean. If true, the function includes endogeneity
//' @return the filtered output
// [[Rcpp::export]]
mat CondHamilton_fctn(mat stateVec, List paramList, bool endogen = false)
{
   int nPeriods = stateVec.n_cols;

   // Setting up the transition probabilities
   double p;
   double q;
   double beta_0;
   double beta_1;
   double varrho;
   if (endogen == false)
   {
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

   // Additional drift in regime 1
   double nu_1 = paramList["nu_1"];

   // Sd of innovations to the drift process
   List xi = paramList["xi"];
   double xi_0 = xi["xi_0"];
   double xi_1 = xi["xi_1"];

   // Initializes output object
   mat Prob_t_ct(nPeriods, 2, fill::zeros);
   // and the sampler
   vec Pr_ct;
   double Pr_ct_0 = (1 - p) / (2 - p - q);
   double Pr_ct_1 = (1 - q) / (2 - p - q);
   Pr_ct = {Pr_ct_0, Pr_ct_1};
   Prob_t_ct.row(0) = Pr_ct.t();

   // Initialize further values for the filter
   double Pr_clt_00;
   double Pr_clt_10;
   double Pr_clt_01;
   double Pr_clt_11;
   double mu_t;
   double mu_lt;
   double nu_0;
   double v_t_0;
   double v_t_1;
   double y_dens_clt_0;
   double y_dens_clt_1;
   double y_dens_clt;
   double Pr_ct_00;
   double Pr_ct_10;
   double Pr_ct_01;
   double Pr_ct_11;
   double Pr_ct_0_floor;
   double Pr_ct_1_floor;
   double y_dens_clt_floor;

   // nu_0 = stateVec.row(1).col(nPeriods - 1)(0);

   for (int i = 1; i < nPeriods; i++)
   {
      // Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
      // (t|t-1)
      Pr_clt_00 = q * Pr_ct_0_floor;
      Pr_clt_10 = (1 - p) * Pr_ct_1_floor;
      Pr_clt_01 = (1 - q) * Pr_ct_0_floor;
      Pr_clt_11 = p * Pr_ct_1_floor;

      // Compute the error term of the RW trend for each regime
      mu_t = stateVec.row(0).col(i)(0);
      mu_lt = stateVec.row(0).col(i - 1)(0);
      nu_0 = stateVec.row(1).col(i)(0);
      v_t_0 = mu_t - mu_lt - nu_0;
      v_t_1 = mu_t - mu_lt - nu_0 - nu_1;

      // Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
      y_dens_clt_0 = dnorm_fctn(v_t_0, xi_0);
      y_dens_clt_1 = dnorm_fctn(v_t_1, xi_1);

      // Adjust the densities in case of endogenous switching

      // Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
      // (t|t-1))
      y_dens_clt = Pr_clt_00 * y_dens_clt_0 + Pr_clt_10 * y_dens_clt_0 + Pr_clt_01 * y_dens_clt_1 + Pr_clt_11 * y_dens_clt_1;
      y_dens_clt_floor = fmax(y_dens_clt, 1e-10);

      // Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
      // (t|t)
      Pr_ct_00 = Pr_clt_00 * y_dens_clt_0 / y_dens_clt_floor;
      Pr_ct_10 = Pr_clt_10 * y_dens_clt_0 / y_dens_clt_floor;
      Pr_ct_01 = Pr_clt_01 * y_dens_clt_1 / y_dens_clt_floor;
      Pr_ct_11 = Pr_clt_11 * y_dens_clt_1 / y_dens_clt_floor;

      // Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
      // Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
      Pr_ct_0 = Pr_ct_00 + Pr_ct_10;
      Pr_ct_1 = Pr_ct_01 + Pr_ct_11;
      Pr_ct_0_floor = fmax(Pr_ct_0, 1e-10);
      Pr_ct_1_floor = fmax(Pr_ct_1, 1e-10);
      // Pr_ct_0_floor = fmin(fmax(Pr_ct_0, .01), .99);
      // Pr_ct_1_floor = fmin(fmax(Pr_ct_1, .01), .99);

      // Records updated probability
      // (t|t)
      Pr_ct = {Pr_ct_0_floor, Pr_ct_1_floor};
      Prob_t_ct.row(i) = Pr_ct.t();
   }

   // Output of filtered values
   return Prob_t_ct;
}

//' @description Executes the hamilton filter for a regime process governing the drift process
//' @param stateVec array with draws of the state vector
//' @param regimeVec vector with regime realizations governing the heteroscedasticity. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @param endogen boolean. If true, the function includes endogeneity
//' @return the filtered output
// [[Rcpp::export]]
mat CondDriftHamilton_fctn(mat stateVec, vec regimeVec, List paramList, bool endogen = false)
{
   int nPeriods = stateVec.n_cols;

   // Setting up the transition probabilities
   double p;
   double q;
   double beta_0;
   double beta_1;
   double varrho;
   if (endogen == false)
   {
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

   // Additional drift in regime 1
   double nu_1 = paramList["nu_1"];

   // Sd of innovations to the drift process
   List xi = paramList["xi"];
   double xi_0 = xi["xi_0"];
   double xi_1 = xi["xi_1"];

   // Initializes output object
   mat Prob_t_ct(nPeriods, 2, fill::zeros);
   // and the sampler
   vec Pr_ct;
   double Pr_ct_0 = (1 - p) / (2 - p - q);
   double Pr_ct_1 = (1 - q) / (2 - p - q);
   Pr_ct = {Pr_ct_0, Pr_ct_1};
   Prob_t_ct.row(0) = Pr_ct.t();

   // Initialize further values for the filter
   double Pr_clt_00;
   double Pr_clt_10;
   double Pr_clt_01;
   double Pr_clt_11;
   double mu_t;
   double mu_lt;
   double nu_0;
   double v_t_0;
   double v_t_1;
   double y_dens_clt_0;
   double y_dens_clt_1;
   double y_dens_clt;
   double Pr_ct_00;
   double Pr_ct_10;
   double Pr_ct_01;
   double Pr_ct_11;
   double Pr_ct_0_floor;
   double Pr_ct_1_floor;
   double y_dens_clt_floor;
   int S_t;

   for (int i = 1; i < nPeriods; i++)
   {
      // Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
      // (t|t-1)
      Pr_clt_00 = q * Pr_ct_0;
      Pr_clt_10 = (1 - p) * Pr_ct_1;
      Pr_clt_01 = (1 - q) * Pr_ct_0;
      Pr_clt_11 = p * Pr_ct_1;

      // Compute the error term of the RW trend for each regime
      mu_t = stateVec.row(0).col(i)(0);
      mu_lt = stateVec.row(0).col(i - 1)(0);
      nu_0 = stateVec.row(1).col(i)(0);
      v_t_0 = mu_t - mu_lt - nu_0;
      v_t_1 = mu_t - mu_lt - nu_0 - nu_1;

      // Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
      S_t = regimeVec[i];
      // for periods where the regime governing the sd is 0
      if (S_t == 0)
      {
         y_dens_clt_0 = dnorm_fctn(v_t_0, xi_0);
         y_dens_clt_1 = dnorm_fctn(v_t_1, xi_0);
      }
      else
      {
         // for periods where the regime governing the sd is 1
         y_dens_clt_0 = dnorm_fctn(v_t_0, xi_1);
         y_dens_clt_1 = dnorm_fctn(v_t_1, xi_1);
      }

      // Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
      // (t|t-1))
      y_dens_clt = Pr_clt_00 * y_dens_clt_0 + Pr_clt_10 * y_dens_clt_0 + Pr_clt_01 * y_dens_clt_1 + Pr_clt_11 * y_dens_clt_1;
      y_dens_clt_floor = fmax(y_dens_clt, 1e-10);

      // Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
      // (t|t)
      Pr_ct_00 = Pr_clt_00 * y_dens_clt_0 / y_dens_clt_floor;
      Pr_ct_10 = Pr_clt_10 * y_dens_clt_0 / y_dens_clt_floor;
      Pr_ct_01 = Pr_clt_01 * y_dens_clt_1 / y_dens_clt_floor;
      Pr_ct_11 = Pr_clt_11 * y_dens_clt_1 / y_dens_clt_floor;

      // Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
      // Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
      Pr_ct_0 = Pr_ct_00 + Pr_ct_10;
      Pr_ct_1 = Pr_ct_01 + Pr_ct_11;
      Pr_ct_0_floor = fmax(Pr_ct_0, 1e-10);
      Pr_ct_1_floor = fmax(Pr_ct_1, 1e-10);

      // Records updated probability
      // (t|t)
      Pr_ct = {Pr_ct_0_floor, Pr_ct_1_floor};
      Prob_t_ct.row(i) = Pr_ct.t();
   }

   // Output of filtered values
   return Prob_t_ct;
}

//' @description Executes the hamilton filter for a regime process governing the heteroscedasticity
//' @param stateVec array with draws of the state vector
//' @param regimeVec vector with regime realizations governing the heteroscedasticity. Vector takes either values 1 or 0
//' @param paramList list with the additional parameters
//' @param endogen boolean. If true, the function includes endogeneity
//' @return the filtered output
// [[Rcpp::export]]
mat CondHeteroscHamilton_fctn(mat stateVec, vec regimeVec, List paramList, bool endogen = false)
{
   int nPeriods = stateVec.n_cols;

   // Setting up the transition probabilities
   double p;
   double q;
   double beta_0;
   double beta_1;
   double varrho;
   if (endogen == false)
   {
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

   // Additional drift in regime 1
   double nu_1 = paramList["nu_1"];

   // Sd of innovations to the drift process
   List xi = paramList["xi"];
   double xi_0 = xi["xi_0"];
   double xi_1 = xi["xi_1"];

   // Initializes output object
   mat Prob_t_ct(nPeriods, 2, fill::zeros);
   // and the sampler
   vec Pr_ct;
   double Pr_ct_0 = (1 - p) / (2 - p - q);
   double Pr_ct_1 = (1 - q) / (2 - p - q);
   Pr_ct = {Pr_ct_0, Pr_ct_1};
   Prob_t_ct.row(0) = Pr_ct.t();

   // Initialize further values for the filter
   double Pr_clt_00;
   double Pr_clt_10;
   double Pr_clt_01;
   double Pr_clt_11;
   double mu_t;
   double mu_lt;
   double nu_0;
   double v_t;
   double y_dens_clt_0;
   double y_dens_clt_1;
   double y_dens_clt;
   double Pr_ct_00;
   double Pr_ct_10;
   double Pr_ct_01;
   double Pr_ct_11;
   double Pr_ct_0_floor;
   double Pr_ct_1_floor;
   double y_dens_clt_floor;
   int S_t;

   for (int i = 1; i < nPeriods; i++)
   {
      // Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
      // (t|t-1)
      Pr_clt_00 = q * Pr_ct_0;
      Pr_clt_10 = (1 - p) * Pr_ct_1;
      Pr_clt_01 = (1 - q) * Pr_ct_0;
      Pr_clt_11 = p * Pr_ct_1;

      // Compute the error term of the RW trend for each regime
      mu_t = stateVec.row(0).col(i)(0);
      mu_lt = stateVec.row(0).col(i - 1)(0);
      nu_0 = stateVec.row(1).col(i)(0);
      S_t = regimeVec[i];
      // for periods where the regime governing the sd is 0
      if (S_t == 0)
      {
         v_t = mu_t - mu_lt - nu_0;
      }
      else
      {
         // for periods where the regime governing the sd is 1
         v_t = mu_t - mu_lt - nu_0 - nu_1;
      }

      // Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
      y_dens_clt_0 = dnorm_fctn(v_t, xi_0);
      y_dens_clt_1 = dnorm_fctn(v_t, xi_1);

      // Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
      // (t|t-1))
      y_dens_clt = Pr_clt_00 * y_dens_clt_0 + Pr_clt_10 * y_dens_clt_0 + Pr_clt_01 * y_dens_clt_1 + Pr_clt_11 * y_dens_clt_1;
      y_dens_clt_floor = fmax(y_dens_clt, 1e-10);

      // Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
      // (t|t)
      Pr_ct_00 = Pr_clt_00 * y_dens_clt_0 / y_dens_clt_floor;
      Pr_ct_10 = Pr_clt_10 * y_dens_clt_0 / y_dens_clt_floor;
      Pr_ct_01 = Pr_clt_01 * y_dens_clt_1 / y_dens_clt_floor;
      Pr_ct_11 = Pr_clt_11 * y_dens_clt_1 / y_dens_clt_floor;

      // Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
      // Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
      Pr_ct_0 = Pr_ct_00 + Pr_ct_10;
      Pr_ct_1 = Pr_ct_01 + Pr_ct_11;
      Pr_ct_0_floor = fmax(Pr_ct_0, 1e-10);
      Pr_ct_1_floor = fmax(Pr_ct_1, 1e-10);

      // Records updated probability
      // (t|t)
      Pr_ct = {Pr_ct_0_floor, Pr_ct_1_floor};
      Prob_t_ct.row(i) = Pr_ct.t();
   }

   // Output of filtered values
   return Prob_t_ct;
}
