# BEAST 2 xmls for use as templates or ready-to-go xmls with BEAST_pype.

This README contains an explanation (lustification of priors, etc) for the xmls in this folder. 

## coal_exp_fixed_clock_1_2e-3.xml & coal_exp_fixed_clock_1_9e-3.xml

If not explicitly stated all other priors are the default values for BEAST 2.7.7

### Site model

HKY gamma catagory 4, based on 

Dellicour, Simon, Guy Baele, Gytis Dudas, et al. 2018. “Phylodynamic Assessment of Intervention Strategies for the West African Ebola Virus Outbreak.” Nature Communications 9 (1): 2222. https://doi.org/10.1038/s41467-018-03763-2.

uncheck estimate

scale = 1

### Clock rate.

Strict clock rate. coal_exp_fixed_clock_1_2e-3.xml is fixed at 1.2e-3 & coal_exp_fixed_clock_1_9e-3.xml 1.9e-3. This is based on:

Abbott, Sam, Katharine Sherratt, Samuel Brand, and Sebastian Funk. 2026. “Estimating the Current Size of the 2026 DRC Bundibugyo Virus Outbreak: A Joint Bayesian Re-Analysis of the McCabe et al. Report.” Zenodo, May 27. https://epiforecasts.io/BVDOutbreakSize/dev/.

prior 
uniform with lower = 0.0011999, and upper = 0.0012 (as a cheat to force the fixed rate)
or      
lower = 0.0018999, and upper = 0.00190001 


### Pop size (ePopSize)

default:  1/X , initial 0.3, [0, inf], offset =0


### growth rate

default:  Laplace, initial = 0.0003,  [-inf,  inf]
mu = 0.001, scale = 0.5, offset = 0
--> in years
--> gives -2 to 2 / year




