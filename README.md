# fk3c - Three-component FK analysis in MATLAB

Author: Nima Riahi
Date: 2014-JUL-14


## Description

Small MATLAB package to perfrom three-component array processing. The code was written with seismic data in mind, but should work for any kind of 2D array where each location has three orthogonal motion sensors.

Feedback, suggestions, and bug reports are most welcome: nimariahizrh@gmail.com


## Contents

Important files:

- SynthDat.m: A trace gather data structure with synthetic data for testing.
- FK3C_Fourier.m: MATLAB script that computes short-time spectrograms from a trace gather data structure.
- FK3C_FK.m: MATLAB script that computes 3C wavenumber spectra and makes some simple visualization.
	
FK3C directory: the package functions

  - compFK3C: computes 3C wavenumber spectra from a time-series of SDM matrices.
  - compSDM: computes SDM matrices from a batch of spectrograms.
  - polpar2cmplx: converts 'human-readable' parameterization of polarization ellipse to complex 3-vector.
  - extrema/extrema2: find extrema in vectors and matrices (downloaded from the web, author Carlos Adrin Vargas Aguilera, Uni de Guadalajara, 2005)
		

	
## Getting started

You can run the two scripts `FK3C_Fourier` and `FK3C_FK` sequentially to test drive the package. Those scripts are also a useful starting point for your own script development and modifications. To get started:

Copy the package content into a working directory of your choice (refered to as ./). Add the package functions to your path:

``` {matlab}
addpath('./FK3C');
````


Open the script `FK3C_Fourier` and make sure that the input directory (top of script) points to where the MAT file `SynthDat.mat` is located. Run the script.

Now open the script `FK3C_FK` and make sure that the input directory is the same as the output directory from `FK3C_Fourier`. If everything worked out you should see two plots showing you the results for the synthetic test dataset stored in SynthDat.mat.

Once this sanity test has passed, you can move on to real data or your own synthetic data (nice exercise). I recommend to prepare a data loading routine that prepares your data in the format of the data structure that the synthetic test data is stored in (see next section).

Once you've familiarized yourself with the package we should schedule a 1-2 hr skype session to discuss it.




## Test dataset

SynthDat.mat contains a trace gather from a recangular 3C array with 3x3 km aperture and 500m receiver spacing. Two synthetic polarized plane waves are contained in the data with a SNR of 0.5 (defined as std(noise)/std(signal) in time domain).

The two plane waves are:

- A Rayleigh wave traveling in positive x-direction at 1900 m/s, with the radial amplitude being 0.4*vertical amplitude.
- A Love wave travling Northwest at 2500 m/s.

Here's the description of the MATLAB structured variable that contains the trace gather and all relevant meta data (copied from FK3C_Fourier). 


'DAT' is a structure that contains a trace gather from all 3*n receivers
of the n element 3C array. It has the following contents:

```
> DAT.h             :  Header 
> DAT.h.coords      : (n, 2) orthogonal projection (if possible) of coordinates of receivers
> DAT.h.stations    : (n, p) custom table containing metadata about the receivers. This field is ignored in processing but may be useful for analysis and troubleshooting
> DAT.h.t0          : (1, 1) Start time of trace gather in MATLAB days
> DAT.h.dt          : (1, 1) Sampling period of trace gather

> DAT.procpars      : (struct) Processing parameters that went into the data, e.g. bandpass filter, etc. (same applies as with the field 'h.stations')


> DAT.data          : (L, 3*n) Actual trace gather. Each column contains the time recording of a receivers. There are 3*n receivers for an n element 3C array. THE ORDERING IS KEY!
>                       - First n MUST be the East components,
>                       - Next n MUST be the North components with the
>                       same order as the East
>                       - Last n MUST be the vertical components with same
>                       ordering as East
>                       
>                       It is CRUCIAL that the field DAT.h.coords ALSO has
>                       the identical ordering as the columns in DAT.data

```


## References

**N. Riahi** and E. H. Saenger (2013); "Rayleigh and Love wave anisotropy in Southern California using seismic noise"; *Geophysical Research Letters* 41 363--369; *doi: 10.1002/2013GL058518*   
[Link](http://dx.doi.org/10.1002/2013GL058518)

**N. Riahi** and G. Bokelmann and P. Sala and E. H.
Saenger (2013); "Time-lapse analysis of ambient surface wave anisotropy: A three-component array study above an underground gas storage";
*Journal of Geophysical Research: Solid Earth* 118 5339--5351; *doi: 10.1002/jgrb.50375*   
[Link](http://onlinelibrary.wiley.com/doi/10.1002/jgrb.50375/abstract)