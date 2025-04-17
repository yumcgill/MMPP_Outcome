# The exact Gibbs sampler for the Markov-modulated Poisson process with an outcome process

This repository contains the R code supporting the paper by Yu Luo and Chris Sherlock (2025) ["Bayesian inference for the Markov-modulated Poisson process with an outcome process"](https://academic.oup.com/jrsssc/advance-article/doi/10.1093/jrsssc/qlaf021/8090125), published in <em><strong>Journal of the Royal Statistical Society Series C</strong></em>. 

The folder contains the R code for 
- <em><strong>data_Ex1</strong></em>: code to generate the simulated data in Example 1 in Section 5;
-	<em><strong>functions</strong></em>: Functions include forward backward algorithm, simulating the path given the start and end states, and the exact Gibbs sampler for one iteration.
-	<em><strong>Gibbs_run</strong></em>: The actual run for the Gibbs sampler.

