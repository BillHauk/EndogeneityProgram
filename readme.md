# Simulation to Accompany "Endogeneity Bias and Growth Regressions"
## *Journal of Macroeconomics*, 2017

### General Info

This repository consists of 15 files -- 9 Stata do files and 6 Stata data files.  In order for the simulation to run, all files must be put in the same working directory.  When the simulation is run, an output file called montecarloresults.dta will be created.  Running the simulation will also create a few extra files in the directory, which can be deleted as needed once the relevant information has been extracted.

### Changing parameters

The "main" program file is MCprogram.do.  All of the parameters that vary in the published paper can be changed in the first 53 lines of code in this file.  The parameters that can be changed are:

#### "True" structural parameters

**alpha** -- the elasticity of output with respect to physical capital in the Solow model used in the data generating process, set by default to 1/3.

**beta** -- the elasticity of output with respect to human capital in the Solow model used in the data generating process, set by default to 1/3.

#### "True" reduced form parameters

**gamma1** -- the reduced form regression parameter on the physical capital in the Solow model, set by default to 0.0579475

**gamma2** -- the reduced form regression parameter on the human capital in the Solow model, set by default to 0.0407447

**gamma3** -- the reduced form regression parameter on population growth and depreciation in the Solow model, set by default to -0.2141621

**gamma4** -- the reduced form regression parameter on lagged GDP in the Solow model, set by default to 0.7962023

#### Variance and correlation scaling for the Fixed Effect term

**fevar** -- the ratio of the variance of the generated FE term to the observed FE term (scaling explained in paper), set by default to 0.5

**fecorr** -- the ratio of the correlations of the generated FE term regressors to the correlations of the observed FE term (scaling explained in paper), set by default to 0.5

#### Number of countries in Monte Carlo sample draw

**N** -- number of countries in a Monte Carlo sample draw

#### Controlled residual correlations

**control** -- one of 7 options in terms of the residual correlations that can be controlled in the parameters below

**skresidcorr** -- correlation between generated residual term to physical capital regressor

**shresidcorr** -- correlation between generated residual term to human capital regressor

**ndgresidcorr** -- correlation between generated residual term to population growth and depreciation regressor

#### Monte Carlo sample draws per simulation

**drawnum** -- number of Monte Carlo sample draws per simulation
