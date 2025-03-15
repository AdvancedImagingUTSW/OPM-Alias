# OPM-Alias

Accompanying code for [Increasing the acquisition speed in oblique plane microscopy via aliasing](https://www.biorxiv.org/content/10.1101/2024.12.25.630337v1.full). 



#### Requirements:
- Installed version of [PetaKit5D](https://github.com/abcucberkeley/PetaKit5D)
- MATLAB 2023a
- Python 3.11
- Conda env with `scikit-image`, `numpy` and `matplotlib`
<br><br>
#### [OPM_aliasing_diagCut.m](https://github.com/AdvancedImagingUTSW/OPM-Alias/blob/main/OPM_aliasing_diagCut.m)

Main MATLAB script used to perform aliasing reconstruction on a single stack. Uses the diagonal top-hat filter on the deskewed stack before rotation. `fullySampled = true` allows user to either simulate downsampling from a ground truth critically sampled stack, or `fullySampled = false` allows loading of real undersampled data, as long as the downsample factor `dsFactor` is known.

Used to generate main Figures 5-6 as well as Supplementary Figures S7, S11.
<br><br>
#### [OPM_aliasing_diagCut_timeSeries.m](https://github.com/AdvancedImagingUTSW/OPM-Alias/blob/main/OPM_aliasing_diagCut_timeSeries.m)

Essentially performs **OPM_aliasing_diagCut** on a full time series. Include the path to the folder as `dataPath` and all Tif/Tiff files will be loaded and processed automatically, with optional saving (`saveData`). Only the first and last frames will be plotted as Figures to save on processing time. Outputs two datasets: aliasing reconstruction vs basic interpolation (PetaKit5D).

Used to generate main Figure 7 and Supplementary Figure S12.
<br><br>
#### [OPM_aliasing_timingCompare.m](https://github.com/AdvancedImagingUTSW/OPM-Alias/blob/main/OPM_aliasing_timingCompare.m)

MATLAB script comparing computational performance for OPM data processing via 4 methods:
- PetaKit5D fast deskew-rotation (DSR) with interpolation
- PetaKit5D interpolation before rotation
- Aliasing reconstruction with top-hat filter before rotation
- Aliasing reconstruction with top-hat filter after rortation (uses fast DSR)

Used to generate Supplementary Figures S9-10.
<br><br>
#### [opm_aliasing_sim_psf.ipynb](https://github.com/AdvancedImagingUTSW/OPM-Alias/blob/main/opm_aliasing_sim_psf.ipynb)

Python simulation of an undersampled, tilted PSF in 2D being recovered using our comb upsampling + top-hat filtering method. Adjustable parameters include the image size `(m,n)`, `tiltAngle`, undersampling factor `dsFactor` and the axial/lateral size of the psf `pSize`.

Used to generate Supplementary Figures S1-2.
