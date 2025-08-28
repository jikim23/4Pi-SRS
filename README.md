# 4Pi-SRS Deconvolution & PSF Simulation

This repository contains MATLAB implementations for 4Pi-stimulated Raman scattering (4Pi-SRS) image restoration and point-spread function (PSF) simulation. 

## 4Pi-SRS deconvolution
The deconvolution script implements a hybrid L₀–TV optimization framework designed to restore sidelobe-removed 4Pi-SRS image stacks.
Reconstructions are first initialized with Richardson–Lucy (RL) deconvolution, providing a physically consistent starting point and reducing the risk of noise-driven artifacts. 
An experimentally measured PSF is incorporated, capturing the sidelobes characteristic of 4Pi-SRS optics. 
The hybrid deconvolution algorithm then iteratively refines the reconstruction, with convergence monitored through relative residuals between predicted and observed data. This approach ensures stable, artifact-free results while preventing over-iteration.


The script has been adapted from: Dong, W.; Tao, S.; Xu, G.; Chen, Y. Blind Deconvolution for Poissonian Blurred Image With Total Variation and L0 -Norm Gradient Regularizations. IEEE Trans. Image Process. 2021, 30, 1030–1043. https://doi.org/10.1109/TIP.2020.3038518.

## 4Pi-SRS PSF Simulation

The PSF simulation module models the three-dimensional point spread function of the 4Pi-SRS microscope using a vectorial diffraction framework based on Richards–Wolf theory. 
The simulation accounts for key optical parameters, including numerical aperture, refractive index, excitation wavelength, and polarization state, to accurately represent field propagation through the dual-objective 4Pi geometry. 
By incorporating interference between counter-propagating beams, the model reproduces the axial sidelobes characteristic of 4Pi microscopy and enables systematic evaluation of their dependence on imaging conditions. 
The simulated PSFs can be directly compared with bead-based experimental measurements or used as input for deconvolution, providing a flexible tool for both microscope design and image reconstruction.

A detailed code description can be found of the supplemental information.
