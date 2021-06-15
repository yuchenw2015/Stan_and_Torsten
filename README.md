# Stan_and_Torsten
This repository contains a few examples from 'Advanced Use of Stan, RStan and Torsten for Pharmacometric Applications' (https://www.metrumrg.com/course/advanced-use-stan-rstan-torsten-pharmacometric-applications/). 
For the Stan models that requires Torsten/Torsten built-in functions, it is common to fail when building those Stan models using 'rstan'. The possible reasons could be that 'rstan' package cannot use the latest Stan or a few built-in Torsten functions have been changed but the example scripts have not.

Here are the updates made in each example script:
1. In each .stan file, the Torsten built-in function names are updated and the output matrix dimensions are updated accordingly.
2. In each .R file, the model fitting codes using 'rstan' function 'stan()' are updated to using 'cmdstanr' functions, and the cmdstan outputs are converted to stan fit output.

Stan_and_Torsten.html is a step-by-step document to show how to modify scripts and compare the non-Torsten-used model to the Torsten-used model.
