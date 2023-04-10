<!--
Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL)

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# Convergent Inverse Scattering using Optimization and Regularization (CISOR)

This software package implements the CISOR reconstruction algorithm along with other benchmark algorithms that attempt to recover the distribution of refractive indices of an object in a multiple scattering regime. The problem of reconstructing an object from the measurements of the light it scatters is common in numerous imaging applications. While the most popular formulations of the problem are based on linearizing the object-light relationship, there is an increased interest in considering nonlinear formulations that can account for multiple light scattering. Our proposed algorithm for nonlinear diffractive imaging, called Convergent Inverse Scattering using Optimization and Regularization (CISOR), is based on our new variant of fast iterative shrinkage/thresholding algorithm (FISTA) and total variation (TV) regularization. The code allows for a systematic comparison between CISOR, SEAGLE, and other state-of-the-art methods on simulated data. We also provide preprocessing scripts that allow for running the algorithm on experimentally measured data that is publicly available from the Fresnel Institute.

This software package reproduces the results from the manuscript:

[1] Y. Ma, H. Mansour, D. Liu, P. T. Boufounos and U. S. Kamilov, "Accelerated Image
 Reconstruction for Nonlinear Diffractive Imaging," 2018 IEEE International Conference
 on Acoustics, Speech and Signal Processing (ICASSP), 2018, pp. 6473-6477.

 and the accompanying extended manuscript published on arxiv:

 [2] https://arxiv.org/pdf/1708.01663v2.pdf

The code includes scripts that implement the CISOR algorithm for imaging from inverse scattering as well as benchmark algorithms for performance comparison in the 2D setting.


## Requirements
Matlab 8 and higher with Image Processing Toolbox.


## Installation and Usage

Select the algorithm that you wish to run by modifying the setting for "algo" in the main_*.m files that are available in each
directory. These will perform the measurement and reconstruction from simulated data.

To run the code on real measurements, you may wish to download into the data/ directory the public dataset provided by
the Fresnel Institute:

- J.-M. Geffrin, P. Sabouroux, and C. Eyraud, "Free space experimental scattering database continuation:
experimental set-up and measurement precision," Inv. Probl., vol. 21, no. 6, pp. S117?S130, 2005.

and

- J.-M. Geffrin and P. Sabouroux, "Continuing with the Fresnel database:
experimental setup and improvements in 3D scattering measurements," Inv. Probl., vol. 25, no. 2, p. 024001, 2009.

We have included a script 'PreProcessFresnelData.m' to preprocess the Fresnel Institute dataset.


## Citation

If you use this code please cite our publications as follows:

```Bibtex
@inproceedings{Ma2018apr,
  author = {Ma, Yanting and Mansour, Hassan and Liu, Dehong and Boufounos, Petros T. and Kamilov, Ulugbek},
  title = {Accelerated Image Reconstruction for Nonlinear Diffractive Imaging},
  booktitle = {IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP)},
  year = 2018,
  pages = {6473--6477},
  month = apr,
  doi = {10.1109/ICASSP.2018.8462400},
  url = {https://www.merl.com/publications/TR2018-008}
}

@article{Ma2018cisor,
  title={Accelerated Image Reconstruction for Nonlinear Diffractive Imaging},
  author={Ma, Yanting and Mansour, Hassan and Liu, Dehong and Boufounos, Petros T. and Kamilov, Ulugbek},
  journal={arXiv preprint},
  year={2017}
}
```


## Contact

Yanting Ma: (yma@merl.com); Hassan Mansour: (mansour@merl.com)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for our policy on contributions.

## License

Released under `AGPL-3.0-or-later` license, as found in the [LICENSE.md](LICENSE.md) file.

All files:

```
Copyright (C) 2022-2023 Mitsubishi Electric Research Laboratories (MERL).

SPDX-License-Identifier: AGPL-3.0-or-later
```
