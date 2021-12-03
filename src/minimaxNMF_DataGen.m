%   File: minimaxNMF_DataGen.m
%   Copyright (c) <2020> <University of Mons>
%   COLORAMAP Group, Univ. of Mons, Belgium
%   https://sites.google.com/site/nicolasgillis/projects/overview
%
%   Permission is hereby granted, free of charge, to any person
%   obtaining a copy of this software and associated documentation
%   files (the "Software"), to deal in the Software without restriction,
%   including without limitation the rights to use, copy, modify and
%   merge the Software, subject to the following conditions:
%
%   1.) The Software is used for non-commercial research and
%       education purposes.
%
%   2.) The above copyright notice and this permission notice shall be
%       included in all copies or substantial portions of the Software.
%
%   3.) Publication, Distribution, Sublicensing, and/or Selling of
%       copies or parts of the Software requires special agreements
%       with the University of Mons and is in general not permitted.
%
%   4.) Modifications or contributions to the software must be
%       published under this license. The University of Mons is granted 
%       the non-exclusive right to publish modifications or contributions 
%       in future versions of the Software free of charge.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%   OTHER DEALINGS IN THE SOFTWARE.
%
%   Persons using the Software are encouraged to notify the
%   COLORAMAP Group, Univ. of Mons, Belgium about bugs. Please reference 
%   the Software in your publications if it was used for them.
% ------------------------------------------------------------------------
% SYNTAX:  
% [ data ] = minimaxNMF_DataGen( param );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% Generates samples as a convex combination of nonnegative endmembers plus
% noise.
%
% ------------------------------------------------------------------------
% INPUTS:
% param     MATLAB structure that contains the parameters for data
%           generation.
%
%   Fields:
%   'scenario'      string - specificies the scenario used                      
%   'nTotalEMs'     scalar - # of endmembers. total needs to be equal to 
%                   (nDominantEMs + nRareEMs)
%   'nBands'        scalar - # of dimensions/spectral bands 
%   'nSamples'      scalar - # samples/pixels to generate. total needs to
%                   be equal to (nRows x nColumns)
%   'nRows'         scalar - # of rows in the "image"
%   'nColumns'      scalar - # of columns in the "image"
%   'W'             nTotalEMs x nBands matrix - groundtruth endmembers
%   'rarePercent'   scalar - percent of pixels containing nonzero Dirichlet 
%                   parameter for the rare endmembers
%   'nDominantEMs'  scalar - # of dominant endmembers, always have nonzero
%                   Dirichlet parameters
%   'nRareEMs'      scalar - # of rare endmembers, have nonzero Dirichlet 
%                   parameter in (nSamples x rarePercent) samples
%   'maxEMsPerPix'  scalar - pixel sparsity
%   'alpha'         nTotalEMs x 1 vector - parameter of the Dirichlet 
%                   distribution
%   'maxPurity'     scalar - maximum abundance of a single endmember in any
%                   sample/pixel
%   'sigmaN'        scalar - noise variance
%
%   Addtional Fields for scenarios 'custom'/'Fig.2_Exp.IVb'/'Fig.3':
%   'usedEMs'       1 x nTotalEMs vector - indices of library spectra   
%   'EM_names'      nTotalEMs x 29 character array - names of used spectra
%
% ------------------------------------------------------------------------
% OUTPUTS: 
% data      MATLAB structure with fields for noisy data matrix and ground 
%           truth.
%
%   Fields:
%   'groundTruth_mat'   nSamples x nBands matrix - indicates  which 
%                       endmembers have nonzero Dirichlet parameter in 
%                       which sample
%   'groundTruth_cube'  nBands x nRows x nColumns tensor - data tensor 
%                       (with noise)  
%   'abundance_mat'     nSamples x nBands matrix - indicates the abundance 
%                       of each endmember in each sample
%   'abundance_cube'    nBands x nRows x nColumns tensor - data tensor 
%                       (with noise)  
%   'synthetic_mat'     nSamples x nBands matrix - data matrix (with noise)
%   'synthetic_cube'    nBands x nRows x nColumns tensor - data tensor 
%                       (with noise)   
%
% ------------------------------------------------------------------------
% DEPENDENCIES:
% [a] sample_dirichlet.m
%       Not included within. 
%       Written by Valentin Leplat, Andersen MS Ang, and Nicolas Gillis. 
%       See file for documentation and reference [2] for details. 
%
% ------------------------------------------------------------------------
% DETAILS:
% The purpose of this data generation function is to create samples for 
% which the ground truth is known
%
% ------------------------------------------------------------------------
% CITATION:
% If you find this code useful in your research, please cite
%
% [1]   Marrinan, Timothy and Nicolas Gillis. 
%       "Hyperspectral Unmixing with Rare Endmembers via Minimax
%       Nonnegative Matrix Factorization." In 28th European Signal 
%       Processing Conference (EUSIPCO) , pp. 69-78. IEEE, 2020.
%
% ------------------------------------------------------------------------
% REFERENCES:
% [1]   Marrinan, Timothy and Nicolas Gillis. 
%       "Hyperspectral Unmixing with Rare Endmembers via Minimax
%       Nonnegative Matrix Factorization." In 28th European Signal 
%       Processing Conference (EUSIPCO) , pp. 69-78. IEEE, 2020.
%
% [2]   Leplat, Valentin, Andersen MS Ang, and Nicolas Gillis. 
%       "Minimum-volume rank-deficient nonnegative matrix factorizations." 
%       In ICASSP 2019-2019 IEEE International Conference on Acoustics, 
%       Speech and Signal Processing (ICASSP), pp. 3402-3406. IEEE, 2019.
% 
% ------------------------------------------------------------------------
% CREATED:      20/07/2020 by Tim Marrinan
%
% LAST EDITED:  28/07/2020 by Tim Marrinan
%
% NOTES: 
%
% ------------------------------------------------------------------------
function [ data ] = minimaxNMF_DataGen( param )

    %% Set parameters
    % More condensed variable names because I am lazy
    Wt              = param.W'; % W is now transposed
    nb              = param.nBands;
    N               = param.nSamples;
    nr              = param.nRows;
    nc              = param.nColumns;
    nEMs            = param.nTotalEMs;
    nDom            = param.nDominantEMs;
    nRare           = param.nRareEMs;
    rarePercent     = param.rarePercent;
    maxEMsPerPix    = param.maxEMsPerPix;
    alpha           = param.alpha;
    purity          = param.maxPurity;
    sigmaN          = param.sigmaN;
    nRarePix        = floor((rarePercent)*N);    

    %% Randomly select pixels to contaminate with rare endmember
    groundTruth = false(nEMs,nr,nc);
    for i = 1 : nr
        for j = 1 : nc
            groundTruth(randperm(nDom,maxEMsPerPix),i,j) = true;
        end
    end

    % Pixels with nonzero Dirichlet parameters for the rare endmembers are 
    % chosen as a rectangle with random dimensions whose area is nRarePix
    rRare = zeros(nRare,1);
    cRare = zeros(nRare,1);
    rStart= zeros(nRare,1);
    cStart = zeros(nRare,1);
    for i = 1 : nRare
        ops = divisors(nRarePix);
        if length(ops)>2
            rRare(i) = ops(randi(length(ops)-2)+1);
        else
            rRare(i) = ops(randi(length(ops)));
        end
        cRare(i) = nRarePix/rRare(i);
        rStart(i) = randi(nr-rRare(i));
        cStart(i) = randi(nc-cRare(i));
        groundTruth(nDom+i,rStart(i):rStart(i)+rRare(i)-1,...
            cStart(i):cStart(i)+cRare(i)-1) = true(1,rRare(i),cRare(i));
    end
    groundTruth = reshape(groundTruth,[nEMs,nr*nc]); 


    %%
    abundances = zeros(nEMs,N);    
    for j = 1 : N
       abundances(groundTruth(:,j),j) =  sample_dirichlet(alpha(groundTruth(:,j)),1)'; 
       while max( abundances(:,j) > purity )
           abundances(groundTruth(:,j),j) =  sample_dirichlet(alpha(groundTruth(:,j)),1)'; 
       end
    end

    abundances = abundances';
    u = reshape(abundances,[nr,nc,nEMs]);
    u = permute(u,[3,1,2]);
    groundTruth = abundances>10^-6;
    v = reshape(abundances,[nr,nc,nEMs]);
    v = permute(v,[3,1,2]);
    vec = zeros(N,nb);
    for i = 1 : N
        if sigmaN > 0
            vec(i,:) = max(0,abundances(i,:)*Wt + sqrt(sigmaN)*(1/sqrt(nb))*randn(1,nb)); %
        else
            vec(i,:) = max(0,abundances(i,:)*Wt); %
        end
    end
    x = reshape(vec,[nr,nc,nb]);
    x = permute(x,[3,1,2]);
    
    
    data.abundance_mat = abundances';
    data.abundance_cube = u;
    data.groundTruth_mat = groundTruth';
    data.groundTruth_cube = v;
    data.synthetic_mat = vec';
    data.synthetic_cube = x;
