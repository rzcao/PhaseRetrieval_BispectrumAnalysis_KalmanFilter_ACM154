# ACM154
Course project, ACM 154
Phase estimation algorithm is developed based on the work by Wu et. al [1]. 
Tengfei Wu, Ori Katz, Xiaopeng Shao, and Sylvain Gigan, Opt. Lett. 41, 5003-5006 (2016)

To run the Kalman filter based phase estimation algorithm, please download the files provided here:
https://dx.doi.org/10.6084/m9.figshare.3593370.v1
(Required files: Digit 4 experimental data.tif & PolarToCartesian.m)

In speckle correlation imaging [2-3], the unkown object is required to be smaller the memory effect range [5]. This is a prior we can use to obtain a better phase estimate. Here, we improved the phase estimate of Bispecturm Analysis by cooperating with such prior on the object size with the help of Kalman filter.

[1] Wu,T., Katz, O., Shao, X. & Gigan, S. Single-shot diffraction-limited imaging through scattering layers via bispectrum analysis. Opt. Lett. 41, 5003-5006 (2016)
[2] Bertolotti, J. et al. Non-invasive imaging through opaque scattering layers. Nature 491, 232–234 (2012).
[3] Katz, O., Heidmann, P., Fink, M. et al. Non-invasive single-shot imaging through scattering layers and around corners via speckle correlations. Nature Photon 8, 784–790 (2014).
[4] Freund, I., Rosenbluh, M. & Feng, S. Memory effects in propagation of optical waves through disordered media. Phys. Rev. Lett. 61, 2328–2331 (1988).
