# Modeling the COVID-19 Infection Rates by Regime Switching Unobserved Components Models

This repository contains all necessary code to replicate the results, figures and tables of
**Haimerl, P. and Hartl T. 2023. Modeling the COVID-19-Infection Rates by Regime Switching Unobserved Components Models. ...**

If you use any parts of this code, please cite this paper. 

## Contents

1. **Replicate_paper.Rmd**: R-Notebook that guides you through all the necessary steps to replicate the paper
2. **R**: Folder containing R scripts
3. **src**: Folder holding Rcpp functions
4. **Output**: Model outputs

## Model specifications

To set up the regime-switching unobserved components (UC) model, let $i_t$ denote the daily reported COVID-19 cases and transform $y_t = log(i_t)$. The measurement equation of the state space form is then given by

$y_t = \mu_t + \gamma_t + c_t + \epsilon_t, \quad \epsilon_t \sim N(0, \sigma^2_{\epsilon})$.

The **trend** $\mu_t$ is defined as

$\mu_t = \mu_{t-1} + \nu_t + \xi_t, \quad \xi_t \sim i.i.d.N(0, \sigma^2_{\xi})$  with

$\nu_t = S_t \nu_1 + (1 - S_t) \nu_0$.

The **seasonal component** $\gamma_t$ admits to

$\gamma_t = -\sum^6_{j=1} \gamma_{t-j}$.

The **cyclical component** $c_t$ is given by

$(1 - L\phi_1 - L^2\phi_2)c_t = \eta_t, \quad \eta_t \sim i.i.d.N(0, \sigma^2_{\eta})$

where $L$ denotes the lag operator and all characteristic roots of the lag polynomial lie outside the unit circle.

The **regime indicator** $S_t \in \{ 0, 1 \}$ governs $\nu_t$. By imposing $\nu_1 < 0$ we declare $S_t = 1$ as the infection down-turning regime. The regimes follow a first order stationary Markov chain with the time homogeneous transition probabilities  

$\Pr(S_t = 0 | S_t = 0) = q, \ \Pr(S_t = 1 | S_t = 1) = p$.

To robustify our findings, we further introduce several extensions to this approach (see section 4). 

Among other, we include a specification that substitutes the deterministic seasonal component $\gamma_t$ with a **stochastic seasonal unit root process**

$\gamma^{UR}_t = -\sum^6_{j=1} \gamma^{UR}_{t-j} + x_t, \quad (1 - L^7) x_t = \omega_t, \quad \omega_t \sim i.i.d.N(0, \sigma^2_{\omega})$.

Furthermore, we consider a regime process that incorporates a **third state** $S_t = 2$, which is characterized by $\nu_t = 0$ during time periods where this state is active.

Lastly, even though not shown in the paper, we investigate specifications that allow for endogeneity between the trend and regime processes,  regime-dependent heteroscedasticity in either the innovations $\xi_t$ or $\eta_t$, as well as a model where a second independent Markov chain $S_t^*$ governs this regime-induced innovation heteroscedasticity.

Turning to the endogenous specification, innovations to the trend and regime processes may exhibit a non-zero correlation $\rho$ as proposed in [Kim et al. (2008)](#Kim); [Kang (2014)](#Kang). For a previous contribution that employs regime dependent heteroscedasticity, i.e. $\sigma^2_{\xi, S_t^*} = S_t^* \sigma^2_{\xi, S_t^*=1} + (1 - S_t^*) \sigma^2_{\xi, S_t^*=0}$ where $S_t^*$ may or may not equal $S_t$ (same applying to $\sigma^2_{\eta, S_t^*}$), we refer to [Engel and Kim (1999)](#Engel).

For further details, see sections 2 and 4.

The following table provides an overview of all implemented models as well as their respective specifications. A :white_check_mark: indicates a non-zero value or the respective extension being implemented. To use any of these specifications, define `modelSpec` in Chunk 2 of the notebook as one of the abbreviations. Note that the trailing dot must be included. Model specific functions can be found in R/Models/


|                        | $\gamma_t$         | $\gamma^{UR}_t$    | $\epsilon_t$       | $c_t$    | $S_t = 2$          | $\rho$             | 2 Regime processes ($S_t \neq S_t^*$) | $\sigma^2_{\xi, S_t^* = 0} \neq \sigma^2_{\xi, S_t^* = 1}$ | $\sigma^2_{\eta, S_t^* = 0} \neq \sigma^2_{\eta, S_t^* = 1}$ |
|:---------------------- |:------------------:|:------------------:|:------------------:|:------------------:|:------------------:|:------------------:|:------------------:|:------------------:|:-:|
| **D.Seas.**            | :white_check_mark: |                    | :white_check_mark: |                    |                    |                    |                    |                    |   |
| **D.Seas.C.**          | :white_check_mark: |                    |                    | :white_check_mark: |                    |                    |                    |                    |   |
| **D.Seas.C.3St.**      | :white_check_mark: |                    |                    | :white_check_mark: | :white_check_mark: |                    |                    |                    |   |
| **D.Seas.2P.TVarSw.**  | :white_check_mark: |                    | :white_check_mark: |                    |                    |                    | :white_check_mark: | :white_check_mark: |   |
| **UR.Seas.**           |                    | :white_check_mark: | :white_check_mark: |                    |                    |                    |                    |                    |   |
| **UR.Seas.C.**         |                    | :white_check_mark: |                    | :white_check_mark: |                    |                    |                    |                    |   |
| **UR.Seas.En.**        |                    | :white_check_mark: | :white_check_mark: |                    |                    | :white_check_mark: |                    |
| **UR.Seas.TVarSw.**    |                    | :white_check_mark: | :white_check_mark: |                    |                    |                    |                    | :white_check_mark: |   |
| **UR.Seas.CVarSw.**    |                    | :white_check_mark: |                    | :white_check_mark: |                    |                    |                    | | :white_check_mark:   |
| **UR.Seas.2P.TVarSw.** |                    | :white_check_mark: | :white_check_mark: |                    |                    |                    | :white_check_mark: | :white_check_mark: |   |
| **UR.Seas.2P.CVarSw.** |                    | :white_check_mark: |                    | :white_check_mark: |                    |                    | :white_check_mark: | | :white_check_mark:   |


## References

<a name="Engel"></a>Engel, C., & Kim, C.-J. (1999). The Long-Run U.S./U.K. Real Exchange Rate. Journal of Money, Credit and Banking, 31(3), 335–356. https://doi.org/10.2307/2601115

<a name="Kang"></a>Kang, K. H. (2014). Estimation of state-space models with endogenous Markov regime-switching parameters. The Econometrics Journal, 17(1), 56–82. http://www.jstor.org/stable/43697655

<a name="Kim"></a>Kim, C.-J., Piger, J., & Startz, R. (2008). Estimation of Markov regime-switching regression models with endogenous switching. Journal of Econometrics, 143(2), 263–273. https://doi:10.1016/j.jeconom.2007.10.002




