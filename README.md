# glmm_floralvisitors

This repository contains the R script and data used to evaluate the effects of livestock management systems on floral visitor communities in dry Chaco forests of Argentina.

## data
nuevos_analisis_polinizadores.ods — abundance matrix of floral visitor species recorded across 12 sites

rasgosfuncionales_apiformes.ods — functional trait database for Apiformes species, including sociality, nesting location, and diet breadth

## glmm_negativebinomial

All statistical analyses were performed in R version 4.5.3. To evaluate the effects of livestock management on floral visitor communities, we fitted generalized linear mixed models (GLMMs) using the glmmTMB package.

We modelled total floral visitor abundance, Shannon diversity, species richness, the abundance of insect orders (Coleoptera, Hymenoptera, Diptera, and Lepidoptera), bee functional groups (stingless bees, Halictidae, other bees, and honeybee), and bee functional traits (sociality, nesting location, and diet breadth) as a function of livestock management treatment (two levels: peasant and silvopastoral). Total abundance was modelled with a negative binomial distribution, Shannon diversity with a Gaussian distribution, and species richness with a Conway-Maxwell-Poisson distribution to account for underdispersion. All remaining models were fitted with a negative binomial distribution to account for overdispersion in count data. For functional trait models, the interaction between treatment and each trait category was included to evaluate whether the effect of livestock management differed among trait groups.  Model diagnostics were performed using simulated residuals from the DHARMa package. 

## inext_analysis

Sampling completeness using iNEXT package. We calculated rarefaction and extrapolation curves to evaluate sampling completeness for each management system (peasant and silvopastoral) 
