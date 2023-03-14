#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' @description Derives the normal density
//' @param x diagonal mat with the quantile values
//' @param var diagonal mat with the variances
//' @returns row vector with the densities

double dnormVar_fctn(double x, double var)
{
   double dens = R::dnorm4(x, 0, sqrt(var), FALSE);
   return dens;
}

//' @description Normal pdf
//' @param double x a quantile
//' @returns double the pdf

double pnorm_fctn(double x)
{
   double dist = R::pnorm5(x, 0.0, 1.0, 1, 0);
   return dist;
}

//' @description Runs the Kim filter
//' @param data vector with the data
//' @param Ini double with the value initializing the trend in the state vector
//' @param systemList list with system matrices
//' @param paramList list with the additional parameters
//' @param outLogLik boolean. If FALSE, function returns the filtered output. If TRUE, function
//' retuns the log likelihood value of the prediction error decomposition
//' @param endogen boolean. If TRUE, the filter includes endogeneity between trend and regime
//' @returns list with the filtered output or the likelihood values
// [[Rcpp::export]]
List FilterRecursions_fctn(vec data, double Ini, List systemList, List paramList,
                           bool outLogLik = true, bool endogen = false)
{
   // System matrices
   mat Tt = systemList["Tt"];
   mat Z = systemList["Z"];
   mat R = systemList["R"];
   // Pre transpose system matrices
   mat transTt = trans(Tt);
   mat transZ = trans(Z);
   mat transR = trans(R);
   // Dimensions of state vector and number of time periods
   int Dimens = Tt.n_cols;
   int periods = data.n_elem;

   // Residual variance in the measurement eq
   double epsilon = paramList["epsilon"];
   double epsilonSq = fmax(pow(epsilon, 2), 1e-6);

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
   // Steady state probabilities
   double Pr_ct_0 = (1 - p) / (2 - p - q);
   double Pr_ct_1 = (1 - q) / (2 - p - q);

   // Vector holding additional Drift in Regime 1
   double nu_1 = paramList["nu_1"];
   vec lambda(Dimens);
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
   vec a_t_clt_00 = zeros(Dimens);
   a_t_clt_00(0) = Ini;
   vec a_t_clt_10 = a_t_clt_00;

   vec a_t_clt_01 = zeros(Dimens);
   a_t_clt_01(0) = Ini + nu_1;
   vec a_t_clt_11 = a_t_clt_01;

   // Initial variance matrix for state vector (Exact diffuse initialization)
   vec P_t_vec(Dimens, fill::value(1000));
   mat P_t_clt_00 = diagmat(P_t_vec);
   mat P_t_clt_10 = diagmat(P_t_vec);
   mat P_t_clt_01 = diagmat(P_t_vec);
   mat P_t_clt_11 = diagmat(P_t_vec);

   // Initializes log-likelihood to record values
   double logLik_T;
   // Initializes output objects
   mat a_t_clt;
   mat P_t_clt;
   mat a_t_ct;
   mat P_t_ct;
   logLik_T = 0;
   cube a_clt_array(Dimens, 4, periods, fill::zeros);
   cube P_clt_array(Dimens, Dimens * 4, periods, fill::zeros);
   cube a_ct_array(Dimens, 2, periods, fill::zeros);
   cube P_ct_array(Dimens, Dimens * 2, periods, fill::zeros);
   mat Prob_t_clt(periods, 2, fill::zeros);
   mat Prob_t_ct(periods, 2, fill::zeros);
   mat v_ct_mat(periods, 4, fill::zeros);
   mat F_ct_mat(periods, 4, fill::zeros);

   // Initialize further values for the filter
   double v_t_00;
   double v_t_10;
   double v_t_01;
   double v_t_11;
   double F_t_00;
   double F_t_10;
   double F_t_01;
   double F_t_11;
   vec a_t_ct_00;
   vec a_t_ct_10;
   vec a_t_ct_01;
   vec a_t_ct_11;
   mat P_t_ct_00;
   mat P_t_ct_10;
   mat P_t_ct_01;
   mat P_t_ct_11;
   double Pr_clt_00;
   double Pr_clt_10;
   double Pr_clt_01;
   double Pr_clt_11;
   double y_dens_clt_00;
   double y_dens_clt_10;
   double y_dens_clt_01;
   double y_dens_clt_11;
   double y_dens_clt;
   double Pr_ct_00;
   double Pr_ct_10;
   double Pr_ct_01;
   double Pr_ct_11;
   double Pr_ct_0_raw;
   double Pr_ct_1_raw;
   vec a_t_ct_0;
   vec a_t_ct_1;
   mat P_t_ct_0;
   mat P_t_ct_1;
   rowvec v_t;
   rowvec F_t;
   double Pr_clt_0;
   double Pr_clt_1;
   rowvec Pr_clt;
   rowvec Pr_ct;
   mat K_t_00;
   mat K_t_10;
   mat K_t_01;
   mat K_t_11;
   double logLik_t;

   // Begin the filtering routine
   for (int i = 0; i < periods; i++)
   {

      //---------------------//
      // Kalman part 1/2    //
      //--------------------//

      // One step ahead prediction error with state vector prediction from
      // (t|t-1)
      v_t_00 = (data[i] - Z * a_t_clt_00)[0];
      v_t_10 = (data[i] - Z * a_t_clt_10)[0];
      v_t_01 = (data[i] - Z * a_t_clt_01)[0];
      v_t_11 = (data[i] - Z * a_t_clt_11)[0];

      // Variance of prediction error
      F_t_00 = (Z * P_t_clt_00 * transZ + epsilonSq)[0];
      F_t_10 = (Z * P_t_clt_10 * transZ + epsilonSq)[0];
      F_t_01 = (Z * P_t_clt_01 * transZ + epsilonSq)[0];
      F_t_11 = (Z * P_t_clt_11 * transZ + epsilonSq)[0];

      // Updating step
      // (t|t, S_t, S_t-1)
      // Kalman gain
      K_t_00 = P_t_clt_00 * transZ * (1 / F_t_00);
      K_t_10 = P_t_clt_10 * transZ * (1 / F_t_10);
      K_t_01 = P_t_clt_01 * transZ * (1 / F_t_01);
      K_t_11 = P_t_clt_11 * transZ * (1 / F_t_11);

      // State vector
      a_t_ct_00 = a_t_clt_00 + K_t_00 * v_t_00;
      a_t_ct_10 = a_t_clt_10 + K_t_10 * v_t_10;
      a_t_ct_01 = a_t_clt_01 + K_t_01 * v_t_01;
      a_t_ct_11 = a_t_clt_11 + K_t_11 * v_t_11;

      // State vector Var-Matrix
      P_t_ct_00 = P_t_clt_00 - K_t_00 * Z * P_t_clt_00;
      P_t_ct_10 = P_t_clt_10 - K_t_10 * Z * P_t_clt_10;
      P_t_ct_01 = P_t_clt_01 - K_t_01 * Z * P_t_clt_01;
      P_t_ct_11 = P_t_clt_11 - K_t_11 * Z * P_t_clt_11;

      //-------------------//
      // Hamilton part     //
      //-------------------//

      // Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
      // (t|t-1)
      Pr_clt_00 = q * Pr_ct_0;
      Pr_clt_10 = (1 - p) * Pr_ct_1;
      Pr_clt_01 = (1 - q) * Pr_ct_0;
      Pr_clt_11 = p * Pr_ct_1;

      // Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
      y_dens_clt_00 = dnormVar_fctn(v_t_00, F_t_00);
      y_dens_clt_10 = dnormVar_fctn(v_t_10, F_t_10);
      y_dens_clt_01 = dnormVar_fctn(v_t_01, F_t_01);
      y_dens_clt_11 = dnormVar_fctn(v_t_11, F_t_11);
      // Adjust the densities in case of endogenous switching
      if (endogen == true)
      {
         y_dens_clt_00 = y_dens_clt_00 * pnorm_fctn((-beta_0 - varrho * (v_t_00 / sqrt(F_t_00))) / sqrt(1 - pow(varrho, 2))) / q;
         y_dens_clt_10 = y_dens_clt_10 * pnorm_fctn((-beta_0 - beta_1 - varrho * (v_t_10 / sqrt(F_t_10))) / sqrt(1 - pow(varrho, 2))) / (1 - p);
         y_dens_clt_01 = y_dens_clt_01 * pnorm_fctn((beta_0 - varrho * (v_t_01 / sqrt(F_t_01))) / sqrt(1 - pow(varrho, 2))) / (1 - q);
         y_dens_clt_11 = y_dens_clt_11 * pnorm_fctn((beta_0 + beta_1 - varrho * (v_t_11 / sqrt(F_t_11))) / sqrt(1 - pow(varrho, 2))) / p;
      }

      // Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
      // (t|t-1))
      y_dens_clt = y_dens_clt_00 * Pr_clt_00 + y_dens_clt_10 * Pr_clt_10 + y_dens_clt_01 * Pr_clt_01 + y_dens_clt_11 * Pr_clt_11;

      // Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
      // (t|t)
      Pr_ct_00 = Pr_clt_00 * y_dens_clt_00 / y_dens_clt;
      Pr_ct_10 = Pr_clt_10 * y_dens_clt_10 / y_dens_clt;
      Pr_ct_01 = Pr_clt_01 * y_dens_clt_01 / y_dens_clt;
      Pr_ct_11 = Pr_clt_11 * y_dens_clt_11 / y_dens_clt;

      // Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
      // Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
      Pr_ct_0_raw = Pr_ct_00 + Pr_ct_10;
      Pr_ct_1_raw = Pr_ct_01 + Pr_ct_11;
      Pr_ct_0 = fmax(Pr_ct_0_raw, 1e-10);
      Pr_ct_1 = fmax(Pr_ct_1_raw, 1e-10);

      //------------------------//
      // Approximation part     //
      //------------------------//

      // Approximate updated values to break exponential growth of required values
      // (t|t, S_t)
      a_t_ct_0 = (a_t_ct_00 * Pr_ct_00 + a_t_ct_10 * Pr_ct_10) / Pr_ct_0;
      a_t_ct_1 = (a_t_ct_01 * Pr_ct_01 + a_t_ct_11 * Pr_ct_11) / Pr_ct_1;

      P_t_ct_0 = ((Pr_ct_00 * (P_t_ct_00 + (a_t_ct_0 - a_t_ct_00) * trans(a_t_ct_0 - a_t_ct_00))) + (Pr_ct_10 * (P_t_ct_10 + (a_t_ct_0 - a_t_ct_10) * trans(a_t_ct_0 - a_t_ct_10)))) / Pr_ct_0;
      P_t_ct_1 = ((Pr_ct_01 * (P_t_ct_01 + (a_t_ct_1 - a_t_ct_01) * trans(a_t_ct_1 - a_t_ct_01))) + (Pr_ct_11 * (P_t_ct_11 + (a_t_ct_1 - a_t_ct_11) * trans(a_t_ct_1 - a_t_ct_11)))) / Pr_ct_1;

      //---------------------//
      // Kalman part 2/2     //
      //---------------------//

      if (outLogLik == false)
      {
         // Store prediction errors and variances
         v_t = {v_t_00, v_t_10, v_t_01, v_t_11};
         v_ct_mat.row(i) = v_t;
         F_t = {F_t_00, F_t_10, F_t_01, F_t_11};
         F_ct_mat.row(i) = F_t;
         // Pr of S_t = j conditional on t = t-1 is given by summarizing over all regimes at t = t-1
         // (t|t-1)
         Pr_clt_0 = Pr_clt_00 + Pr_clt_10;
         Pr_clt_1 = Pr_clt_01 + Pr_clt_11;
         // Records predicted probabilities
         // (t|t-1)
         Pr_clt = {Pr_clt_0, Pr_clt_1};
         Prob_t_clt.row(i) = Pr_clt;
         // Records updated probability
         // (t|t)
         Pr_ct = {Pr_ct_0, Pr_ct_1};
         Prob_t_ct.row(i) = Pr_ct;
         // Store predicted values
         // (t|t-1)
         a_t_clt = join_rows(a_t_clt_00, a_t_clt_10, a_t_clt_01, a_t_clt_11);
         a_clt_array.slice(i) = a_t_clt;
         P_t_clt = join_rows(P_t_clt_00, P_t_clt_10, P_t_clt_01, P_t_clt_11);
         P_clt_array.slice(i) = P_t_clt;
         // Store updated values
         // (t|t)
         a_t_ct = join_rows(a_t_ct_0, a_t_ct_1);
         a_ct_array.slice(i) = a_t_ct;
         P_t_ct = join_rows(P_t_ct_0, P_t_ct_1);
         P_ct_array.slice(i) = P_t_ct;
      }
      else
      {
         // Store approximate likelihood
         logLik_t = -log(y_dens_clt);
         // Sum up log-likelihood over all iterations
         logLik_T = logLik_T + logLik_t;
      }

      // Prediction step with approximated updates to complete loop
      // (t+1|t)
      // Regime 0 (high drift)
      a_t_clt_00 = Tt * a_t_ct_0;
      a_t_clt_10 = Tt * a_t_ct_1;
      P_t_clt_00 = Tt * P_t_ct_0 * transTt + expandQ_0;
      P_t_clt_10 = Tt * P_t_ct_1 * transTt + expandQ_0;

      // Regime 1 (low drift)
      a_t_clt_01 = Tt * a_t_ct_0 + lambda;
      a_t_clt_11 = Tt * a_t_ct_1 + lambda;
      P_t_clt_01 = Tt * P_t_ct_0 * transTt + expandQ_1;
      P_t_clt_11 = Tt * P_t_ct_1 * transTt + expandQ_1;
   }

   if (outLogLik == false)
   {
      List filterOutput = List::create(
          Named("a_clt") = a_clt_array, // predicted state vector
          _["P_clt"] = P_clt_array,     // predicted state vector var
          _["a_ct"] = a_ct_array,       // updated state vector
          _["P_ct"] = P_ct_array,       // updated state vector var
          _["Prob_clt"] = Prob_t_clt,   // predicted regime probs
          _["Prob_ct"] = Prob_t_ct,     // updated regime probs
          _["Pred_err"] = v_ct_mat,     // one-step-ahead prediction error
          _["Pred_err_Var"] = F_ct_mat  // Var of pred error
      );
      return filterOutput;
   }
   else
   {
      // Impose constraint that p,q > 0.9
      if (p <= 0.9)
      {
         logLik_T = 1e4 * (1 - p) + logLik_T;
      }
      else if (q <= 0.9)
      {
         logLik_T = 1e4 * (1 - q) + logLik_T;
      }
      List logLikList = List::create(logLik_T);
      return logLikList;
   }
}
