# Random-Thin-Marginal
SCR with random thinning samplers marginalizing out latent individual IDs. Poisson observation model only. 
To speed up computation, I use the approach of Herliansyah et al. (2024, section 4.3) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

There are 4 types of models: 
1) Single session
2) Multisession
3) Single session with density covariates and habitat mask
4) Multisession with density covariates and habitat mask


Random thinning models that allow observation models other than Poisson and/or categorical partial IDs can be found here:

https://github.com/benaug/RandomThinIDCov

These are more limited (e.g., no habitat mask, density covariates), but can be modified.

Analogous repositories in the "Marginal Unmarked Trilogy" can be found here:
Unmarked SCR: https://github.com/benaug/Unmarked-SCR-Marginal
Spatial mark-resight: https://github.com/benaug/Spatial-Mark-Resight-Marginal
