# Implementation of the image source method
This repository contains the implementation of the image source method [AB79,BEW18] as used in the PhD thesis "Wave Field Synthesis in a listening room" by Vera Erbes. You can use the script ```Test_calculate_ISM.m``` to try it out.

## Software requirements
* **Matlab**:\
This Matlab code has been tested with Matlab R2015a.
* **Sound Field Synthesis Toolbox**:\
The scripts use the Matlab version of the <a rel="license" href="http://www.sfstoolbox.org">Sound Field Synthesis Toolbox</a>, <a rel="license" href="http://doi.org/10.5281/zenodo.2597212">release 2.5.0</a>.
* **Spatially Oriented Format for Acoustics (SOFA)**:\
When simulating with a binaural receiver, you need HRTFs given in the format <a rel="license" href="http://www.sofaconventions.org">SOFA</a>. To read this data, you need the <a rel="license" href="http://sourceforge.net/projects/sofacoustics/files/">SOFA API for Matlab/Octave</a>. The code has been tested with release 1.0.4.

## Licenses
The Matlab code is licensed under the MIT license.

This Material is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Literature
[AB79] J. B. Allen and D. A. Berkley (1979): Image source method for efficiently simulating small-room acoustics. J. Acoust. Soc. Am. 65(4):943-950
         
[BEW18] F. Brinkmann, V. Erbes and S. Weinzierl (2018): Extending the closed form image source model for source directivity. Proc. of the 44th German Annual Conf. on Acoustics (DAGA)
