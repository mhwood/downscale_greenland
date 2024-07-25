# Downscaled Ocean Models around Greenland using MITgcm

This repository was created to document downscaled ocean model configurations around Greenland using the MIT General Circulation Model. 

All models here are labeled L0, L1, and L2 to denote the extent to which they have been downscaled from existing models. The models here are derived from the Estimating the Circulation and Climate of the Ocean (ECCO) consoritum ECCOv5 Alpha ([LLC270](https://github.com/MITgcm-contrib/llc_hires/tree/master/llc_270)) state estimate. Global (L0) models are identical to the ECCO LLC270 solution with slight modifications. All regional L0 and L1 models are constructed on subsets of the "Lat-Lon-Cap" (LLC) grid and are forced with ECCO-derived data. The regional L2 models are built on manually-constructed grids which are separate from the LLC grid. 


### Manuscripts
1. [Decadal Evolution of Ice-Ocean Interactions at a Large East Greenland Glacier Resolved at Fjord Scale With Downscaled Ocean Models and Observations](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL107983)

   Wood et al 2024

   The models associated with this manuscript are the L0 [global](https://github.com/mhwood/downscale_greenland/tree/main/L0/global) model, the [L1_CE_Greenland](https://github.com/mhwood/downscale_greenland/tree/main/L1/L1_CE_Greenland) model, and the [L2_Scoresby_Sund](https://github.com/mhwood/downscale_greenland/tree/main/L2/L2_Scoresby_Sund) model.

2. Feedbacks between fjord circulation, m´elange melt, and the subglacial discharge plume at Kangerlussuaq Glacier, East Greenland

   Wood et al 2024 b, Under Review

   The models associated with this manuscript are listed under [L2_Kanger](https://github.com/mhwood/downscale_greenland/tree/main/L2/L2_Kanger) model.


