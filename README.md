# Designs

This repository holds all the used designs as GitHub releases.

## Latest versions

- SB Night
  - [v3.0.0](https://github.com/JaakkoAhola/DESIGN/releases/tag/v3.0.0)
- SB Day  
  - [v3.1.0](https://github.com/JaakkoAhola/DESIGN/releases/tag/v3.1.0)
- SALSA Night
  - [v3.2.4](https://github.com/JaakkoAhola/DESIGN/releases/tag/v3.2.4)
- SALSA Day
  - [v3.3.3](https://github.com/JaakkoAhola/DESIGN/releases/tag/v3.3.3)


## Variables

### Meteorological variables
- q_inv
  - inversion strength of total water mass mixing ratio at the boundary layer top
  - (g kg^-1)
- tpot_pbl
  - liquid water potential temperature in the boundary layer
  - (K)
- tpot_inv
  - inversion strength of liquid water potential temperature (K)
- lwp
  - liquid water path for the cloud
  - (g m^-2)
- pblh
  - planetary boundary layer height described as a pressure difference from the surface
  - (hPa)

### Aerosol related variables
- cdnc
  - number of particles with diameter greater than 50 nm in the PBL (SB only)
  - (kg^-1)
- rdry_AS_eff
  - effective dry radius of accumulation mode (SALSA only)
  - (m)
- ks
  - aerosol number concentration in Aitken mode (SALSA only)
  - (kg^-1)
- as
  - aerosol number concentration in accumulation mode (SALSA only)
  - (kg^-1)
- cs
  - aerosol number concentration in coarse mode (SALSA only)
  - (kg^-1)
### Solar radiation related variable
- cos_mu
  - cosine of solar zenith angle (day time only)

## Authors
Designs made by: Muzzaffer Ege Alper
Readme written by: Jaakko Ahola


## git tags
i.j.k

i
- major update i.e. change of design purpose, first one would be proof-of-concept

j
- minor update i.e. change of variables or constraints

k
- patching i.e. bug fixes but NO changes in variables or constraints

