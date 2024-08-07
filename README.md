# Semicompeting Risks Model

This repository contains code for a semicompeting risks model with two nonterminal 
events and one terminal event. The model utilizes Bayesian inference using Markov chain 
Monte Carlo (MCMC) sampling, with options for direct sampling, adaptive rejection sampling, 
or the Metropolis-Hastings algorithm.

## Data Structure

The data structure used in this model includes the following variables:

- $t_0$: Time to death without developing any nonterminal event.
- $t_1$: Time to the development of the first nonterminal event.
- $t_2$: Time to death following the diagnosis of the first nonterminal event.
- $t_3$: Time to the development of a subsequent nonterminal event following the first nonterminal event.
- $t_4$: Time to death following the subsequent nonterminal event.
- $d$: Indicator for experiencing the first nonterminal event (1 if experienced, 0 otherwise).
- $u$: Indicator for experiencing the second nonterminal event (1 if experienced, 0 otherwise).
- $v$: Indicator for death occurrence (1 if death occurred, 0 otherwise).

### Cases

The six cases and corresponding survival time variables are as follows:

1. **Case 1**: $I(d=0, u=0, v=1)$; $t_0$
2. **Case 2**: $I(d=1, u=0, v=1)$; $t_1, t_2$
3. **Case 3**: $I(d=1, u=1, v=1)$; $t_1, t_3, t_4$
4. **Case 4**: $I(d=1, u=1, v=0)$; $t_1, t_3, t_4$
5. **Case 5**: $I(d=1, u=0, v=0)$; $t_1, t_2$
6. **Case 6**: $I(d=0, u=0, v=0)$; $t_0$

### Missing Values

For cases where $t_0$, $t_1$, $t_2$, $t_3$, and $t_4$ have missing values, these are indicated by 9999.
 These missing values or unavailable survival times are handled appropriately in the analysis.

## Covariates

The model uses two covariates $x_1$ and $x_2$. The same set of covariates is used in both the 
logistic regressions and the Cox regression.

## Initialization

It is suggested to run $t_0$, $t_1$, $t_2$, $t_3$, and $t_4$ separately and use those values as initial values 
for the MCMC algorithm. However, random initial values should also work.

## Model Parameters

- **n**: Sample size
- **nbin**: Number of burn-in iterations
- **nrep**: Number of iterations after burn-in
- **nthin**: Number of thinning iterations

## Baseline Hazard Function

The model considers different baseline hazard functions for each event time. 
For each $k$, $n_{g_k}$ indicates the number of pieces for the baseline hazard function, for $k=0, 1, 2, 3, 4$.

## Usage

- The main files
  + `main.f`: main functions
  + `gibbs.f`: functions for MCMC sampling
  + `exampledat.txt`: example data to run the model
- Additional files
  + `optim1.f`: for optimization
  + `gilks2.f`: function to compute adaptive rejection sampling
  + `hpd.f`: for computing highest posterior density
  
Keeping all the files in the same directory, run the `main.f` function. 
The IMSL library is required to run this code. Please ensure that the library is installed and properly 
configured before running the model.


