# elasticmapper
Python tool to analyze and visualize elastic maps.

## Installation

`> git clone https://github.com/uafgeotools/elasticmapper.git`

`> cd elasticmapper`

`> conda env create -f environment.yml`

## Notes

Matrices are represented by capitalized first letters in functions and 
variables.\
C --> Stiffness matrix in the Voigt representation\
T --> stiffness matrix represented in the TapeTape basis\
tsp --> theta, sigma and phi angles representing 3D rotation

## References

If you use this package in your own research, please cite the following paper:

- Gupta, A., Tape, C., Brown, M.J., & Becker, T.W., (2025): “Navigating the 
  space of seismic anisotropy for crystal and whole-Earth scales”.
  [Geophysical Journal International](https://doi.org/10.1093/gji/ggaf134), 242(1), 1-18.

The following papers can also be cited in relation to this package:

- Tape, W., & Tape, C., (2024): "A Reformulation of the Browaeys and Chevrot 
  Decomposition of Elastic Maps". [Journal of Elasticity](https://doi.org/10.1007/s10659-024-10056-x), 156, 415-454.
- Tape, W., & Tape, C., (2022): "Two complementary methods of inferring 
  elastic symmetry". [Journal of Elasticity](https://doi.org/10.1007/s10659-022-09898-0), 150, 91-118.
- Tape, W., & Tape, C., (2021): "Elastic symmetry with beachball pictures". 
  [Geophysical Journal International](https://doi.org/10.1093/gji/ggab183), 
  227(2), 970-1003.