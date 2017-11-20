This package performs the Bayesian model selection on different types of FCS data, 
including multiple TACF curves, single intensity traces, and imaging FCS data. Usage of the package 
for each data type is demonstrated in the scripts 'run_FCS_bayes_multi_curves.m'
, 'run_FCS_bayes_single_trace.m', and 'run_imagingFCS_bayes.m'. The FCS data sets for 
the demonstrations can be downloaded from:

http://fcs-bayes.org/examples/examples.zip

To run the demonstration, unpack the folder 'example' in the same folder of the package and run the MATLAB scripts. 

This package was tested on MATLAB R2011a with Statistics Toolbox and Image Processing Toolbox installed.

###################################################################
Example 1 Analysis of multiple TACF curves without raw intensity traces

Proper calculation of model probabilities requires incorporating the noise covariance matrix. 
When only TACF curves are available instead of raw intensity traces, the noise covariance matrix must be 
estimated from multiple independent TACF curves measured from the same underlying physical process.

Example file 'FCS_64cur_vary_Dr.mat' contains a cell 'corrFCS_cell' of TACF curves data sets, 
a lag time vector 't', and the aspect ratio of the focal volume 'ks'. 'corrFCS_cell' contains 15 sets
of simulated two-component diffusion TACF curves data with varying D2/D1 ratios. Each data set has 64 
individual curves of 2-component diffusion simulated with the same parameters. 
Please refer to the paper for details about the simulation. The script 'run_FCS_bayes_multi_curves.m' 
runs the Bayesian model selection on the provided example data.

Typically > 4 TACF curves with 128 data points are needed for the regularized covariance matrix to be 
non-singular, but the minimum required number of TACFs increases as the number of data points increases.

When the raw photon-count trace or photon arrival time is available, the noise in a single TACF can be 
calculated from its underlying photon-count products (see the paper for details).

######################################################################################
Example 2 Analysis of raw intensity traces

Example files 'FCS_8int_vary_Dr_D2=#.mat' contains a series simulated photon count traces of two-component 
diffusion with varying D2/D1 ratios. Variable 'im' is an array of eight simulated 8-bit photon count traces with D2=# and the sampling time 'dt'. Note that the photon count trace contains more than 100 million counts, so a machine equipped with at least 8GB ram is required to run this analysis. The run-script 'run_FCS_bayes_single_trace.m' performs the analysis through the following steps:

(1) Running the blocking analysis to find the minimum block-time (dot) in the blocking curve to generate 
independent samples required for noise covariance estimation. When the minimum block-time of the blocking 
curve cannot be found by the algorithm (gray curves), the maximal block-time available in the blocking curve 
is used to estimate the noise.

(2) Computing the TACF and its covariance matrix using the multi-tau algorithm. If dual detection scheme is 
used (two photon count traces recorded simultaneously using two detectors), use the function 'compute_CCF2' 
instead of 'compute_ACF2' for computing TACF.

(3) Running the Bayesian model selection on the computed TACF curves, same as the analysis of multiple curves.

##############################################################################################
Example 3 Analysis of fluorescence movies

For imaging FCS in which each pixel acts as a single detector, the intensity trace at each pixel may be 
analyzed individually in the same manner as in Example 2, yielding maps of model probabilities and diffusion 
coefficients.

The example file is a simulated 16-bit fluorescence movie of a partitioned microdomain in the membrane 
with lower diffusivity than the outer, non-domain region (D1 = 4 um2/s and D2 = 0.04 um2/s). Particles 
undergo free diffusion and switch diffusivity instantaneously upon crossing the microdomain boundary. 
Two-component diffusion is detected near the domain boundary because the PSF averages intensity from 
particles inside and outside the domains. The run-script 'run_imagingFCS_bayes.m' performs the analysis 
through the following steps:

(1) Running the blocking analysis at each pixel to find the minimum block-time (plateau) in the blocking 
curve to generate independent samples required for noise covariance estimation. When the minimum block-time 
cannot be found by the algorithm, the pixel is labeled and the maximal block-time available in the blocking 
curve is used to estimate the noise.

(2) Computing the TACF and its covariance matrix using the multi-tau algorithm and the specified block-time.

(3) Running the Bayesian model selection on the computed TACF curves and generating the map of model 
probabilities, the map and histogram of diffusion coefficients.

For more information, visit http://fcs-bayes.org