# PrecipitationFields

This repository contains files with R-code used to generate the plots and experiments for the paper 'Generating Calibrated Ensembles of Physically Realistic, High-Resolution Precipitation Forecast Fields Based on GEFS Model Output' submitted to the Journal of Hydrometeorology. The following gives a brief description of the individual files:


- AuxiliaryFunctions.r
A collection of small functions (calculating ensemble statistics, plotting reliability diagrams, etc.) used by the other R-scripts 


- CaseStudy-ECC.r
R-code for generating ECC-Q and ECC-mQ-SNP ensemble forecasts for a single date and for depicting these ensembles graphically


- CaseStudy-MDSS.r
R-code for generating MDSS-RO and MDSS-SDA ensemble forecasts for a single date and for depicting these ensembles graphically


- CaseStudy-StSS.r
R-code for generating StSS ensemble forecasts for a single date and for depicting these ensembles graphically


- CodeForGraphics.r
R-code for generating the plots used to illustrate and motivate the MDSS-SDA and ECC-mQ-SNP method


- GenerateFcstFields-ECC-StSS.r
R-code for generating StSS, ECC-Q and ECC-mQ-SNP ensemble forecasts for the entire verification period


- GenerateFcstFields-MDSS.r
R-code for generating MDSS-SDA and MDSS-RO ensemble forecasts for the entire verification period


- LookAtForecastFields.r
R-code for depicting and quickly browsing through ensemble forecasts generated via 'GenerateFcstFields-ECC-StSS.r' or 'GenerateFcstFields-MDSS.r'


- Verification.r
R-Code for calculating verification statistics and depicting verification results graphically

