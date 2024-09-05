# Realistic Microstructure Simulator (RMS): Monte Carlo simulations of diffusion in 3D cells of permeable membrane (CUDA C++)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10651814.svg)](https://doi.org/10.5281/zenodo.10651814)

The code implements 3d Monte Carlo simulations originally developed in [Lee, et al., Journal of Neuroscience Methods, 2021](https://doi.org/10.1016/j.jneumeth.2020.109018), demonstrating the effect of water exchange and axonal beading on the time-dependent kurtosis in [Wu et al., Science Advances 2024](https://doi.org/10.1126/sciadv.adk1817).

* **Demo 0, one sphere:** We simulate water diffusion initialized at the center of a permeable sphere, and validate the permeaion probility equation based on the calculated particle density, as detailed in [Fieremans & Lee, NeuroImage, 2018](https://doi.org/10.1016/j.neuroimage.2018.06.046) and [Lee, et al., Journal of Neuroscience Methods, 2021](https://doi.org/10.1016/j.jneumeth.2020.109018).
* **Demo 1, random stick:** We perform simulations of diffusion in a medium composed of randomly positioned, randomly oriented cylinders with or without caliber variations (beadings) with a fixed membrane permeability.
* **Demo 2, random stick:** We analyze the simulation result in Demo 1.
* **Demo 3, random stick, fixed exchange time:** We perform simulations of diffusion in a medium composed of randomly positioned, randomly oriented cylinders with or without caliber variations (beadings) with a fixed exchange time.
* **Demo 4, random stick, fixed exchange time:** We analyze the simulation result in Demo 3.
* **Demo 5, random stick, fixed exchange time in two values:** We perform simulations of diffusion in a medium composed of randomly positioned, randomly oriented cylinders with or without caliber variations (beadings) with a fixed exchange time in two different values.
* **Demo 6, random stick, fixed exchange time in two values:** We analyze the simulation result in Demo 5.
* **Demo 7, plot geometry:** We plot the 3-dimensional geometries.

## References
* **Monte Carlo simulation**
  - [Fieremans, et al., NMR Biomed, 2010](https://doi.org/10.1002/nbm.1577)
  - [Novikov, et al., Nature Physics, 2011](https://doi.org/10.1038/nphys1936)
  - [Fieremans and Lee, NeuroImage 2018](https://doi.org/10.1016/j.neuroimage.2018.06.046)
  - [Lee, et al., Communications Biology 2020](https://doi.org/10.1038/s42003-020-1050-x)
  - [Lee, et al., NeuroImage 2020](https://doi.org/10.1016/j.neuroimage.2020.117228)
  - [Lee, et al., Journal of Neuroscience Methods 2021](https://doi.org/10.1016/j.jneumeth.2020.109018)
  - [Lee, et al., NMR in Biomedicine 2024](https://doi.org/10.1002/nbm.5087)

## Authors
* [Hong-Hsi Lee](http://www.diffusion-mri.com/people/hong-hsi-lee)
* [Dmitry S Novikov](http://www.diffusion-mri.com/people/dmitry-novikov)
* [Els Fieremans](http://www.diffusion-mri.com/people/els-fieremans)

## Acknowledgement
We would like to thank [Ricardo Coronado-Leija](https://scholar.google.com/citations?user=V5hykxgAAAAJ&hl=en) for the fruitful discussion of simulation implementation.

## License
This project is licensed under the [LICENSE](https://github.com/leehhtw/monte-carlo-simulation-3D-RMS-exchange/blob/main/LICENSE).
