
# Minimax NMF Version 1.0 #
Copyright (c) <2020> <University of Mons>  
COLORAMAP Group, Univ. of Mons, Belgium  
(See individual files for license)  


# Description #

This matlab package contains files implementing an approximate subgradient 
method for solving a minimax nonnegative matrix factorization problem,

       min               max   || X_i - (W)(H_i) ||_F
    W in R_{>0}^{b x r}   i 
    H in R_{>0}^{r x p}

for some data, {X_i} in R_{>0}^{b x p} for i = 1, ... , n, as described 
in the paper: 
"Hyperspectral Unmixing with Rare Endmembers via Minimax Nonnegative Matrix 
Factorization." by TIM MARRINAN and NICOLAS GILLIS.


# Contents #

The matlab package includes the following directories and files:  
\  
01. README.txt  
02. demo.m

\src\  
03. main.m  
04. minimaxNMF_DataGen.m  
05. minimaxNMF.m  


(The code in this directory accompanies [3] and was downloaded from https://sites.google.com/site/nicolasgillis/code on 22/10/2019 )

\src\minvolNMF\  
06. compareWs.m  
07. FGMfcnls.m  
08. FGMqpnonneg.m  
09. minvolNMF.m  
10. munkres.m  
11. nnlsHALSupdt.m  
12. ReadMe.m  
13. sample_dirichlet.m  
14. SimplexProj.m  
15. SNPA.m  
16. synthetic_data.m  

\examples\  
17. minimaxNMF_Fig1.m  
18. minimaxNMF_Fig2.m   
19. minimaxNMF_Fig3.m  
20. minimaxNMF_ScenarioSpecification.m  
21. USGS_Library.mat  


# Abstract #

Hyperspectral images are used for ground-cover classification because many 
materials can be identified by their spectral signature, even in images 
with low spatial resolution. Pixels in such an image are often modeled as 
a convex combination of vectors, called endmembers, that correspond to the 
reflectance of a material to different wavelengths of light. This is the 
so-called linear mixing model. Since reflectance is inherently nonnegative, 
the task of unmixing hyperspectral pixels can be posed as a low-rank 
nonnegative matrix factorization (NMF) problem, where the data matrix is 
decomposed into the product of the estimated endmembers and their 
abundances in the scene. The standard NMF problem then minimizes the 
residual of the decomposition. Thus, using NMF works well when materials 
are present in similar amounts, but if some materials are 
under-represented, they may be missed with this formulation. Alternatively, 
we propose a novel hyperspectral unmixing model using a collection of NMF 
subproblems solved for patches of the original image. The endmembers are 
estimated jointly, such that the the maximum residual across all patches 
is minimized. In this paper we estimate the solution to the patch-based 
minimax NMF model, and show that it can estimate rare endmembers with 
superior accuracy.


# File Usage #

The file demo.m can be run with no modifications to add the folders to the 
path and run a synthetic example using endmembers taken from the USGS 
spectral library [4]. It computes the minimax NMF of the data and the 
traditional minimum-volume NMF from [3], displays the estimated endmembers, 
and compares the accuracy of the factorizations.  

Alternatively, the file main.m can be run with no modifications to 
recreate the figures from the associated paper.  This will call each of the 
functions in the 'examples' directory, whose parameters are set by 
minimaxNMF_ScenarioSpecification.m. The variables 'queued' and 'nRuns' can 
be modified to control which experiments are run and for how many 
iterations.

To re-run the experiments and generate the plots from the paper:
01. Run main.m

In the paper the experiments are each run for nRuns = 50 Monte Carlo 
trials, however this takes numerous hours to complete. Alternatively, the 
behavior of the methods can be assessed from a much smaller number of 
independent trials or with less noise variances tested to avoid the long 
computation time. 

On a laptop with a 2.6 GHz Intel Core i7-8850H processor and 16 GB of RAM, 
the experiments from the paper run for nRuns = 5 take approximately 80 
minutes. To run the complete experiment set with nRuns = 50 takes 
approximately 13 hours.

To run your own experiments:
01. In the file minimaxNMF_ScenarioSpecification.m, modify the parameters 
of the 'custom' scenario as desired.
02. Depending on your desired outcome, look to the files demo.m, 
minimaxNMF_Fig1.m, or  minimaxNMF_Fig2.mfor examples of how to use the data 
generation function and the minimax NMF algorithm.


# Release Notes #

Version 1.0 of the code does not currently provide a standalone function 
to cut data into patches.



# References #

 [1]   Marrinan, Timothy and Nicolas Gillis. 
       "Hyperspectral Unmixing with Rare Endmembers via Minimax
       Nonnegative Matrix Factorization." In 28th European Signal 
       Processing Conference (EUSIPCO) , pp. 69-78. IEEE, 2020.

 [2]   Gillis, N., 2014. Successive nonnegative projection algorithm for 
       robust nonnegative blind source separation. SIAM Journal on Imaging
       Sciences, 7(2), pp.1420-1450.

 [3]    Leplat, Valentin, Andersen MS Ang, and Nicolas Gillis. 
       "Minimum-volume rank-deficient nonnegative matrix factorizations." 
       In ICASSP 2019-2019 IEEE International Conference on Acoustics, 
       Speech and Signal Processing (ICASSP), pp. 3402-3406. IEEE, 2019.

 [4]   Clark, Roger N., Gregg A. Swayze, A. Gallagher, T. V. V. King, 
       and W. M. Calvin. "The us geological survey digital spectral 
       library." US Geological Survey Open File Report (1993): 93-592.



# Citation: #

If you find this code useful in your research, please cite:  
       
       Hyperspectral unmixing with rare endmembers via minimax nonnegative matrix factorization.
       T. Marrinan and N. Gillis.
       Proc. IEEE 28th European Signal Processing Conference (EUSIPCO), (2021): 1015-1019.



# Contact #

In case of questions, suggestions, problems etc. please send an email.

Tim Marrinan:  
marrinat@oregonstate.edu

This matlab package is hosted at:  
http://www.tmarrinan.com/research-interests-code/  
https://sites.google.com/site/nicolasgillis/code
