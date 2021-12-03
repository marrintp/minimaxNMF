%   File: main.m
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
%       published under this license. The University of Mons
%       is granted the non-exclusive right to publish modifications
%       or contributions in future versions of the Software free of charge.
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
%
% ------------------------------------------------------------------------
% SYNTAX:  
% This is a script not a function.  Parameters are set within. 
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This script runs all of the experiments and reproduces the plots from the
% associated paper.
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
% [2]    Leplat, Valentin, Andersen MS Ang, and Nicolas Gillis. 
%       "Minimum-volume rank-deficient nonnegative matrix factorizations." 
%       In ICASSP 2019-2019 IEEE International Conference on Acoustics, 
%       Speech and Signal Processing (ICASSP), pp. 3402-3406. IEEE, 2019.
%
% ------------------------------------------------------------------------
% CREATED:      04/11/2019 by Tim Marrinan
%
% LAST EDITED:  28/07/2020 by Tim Marrinan
%
% NOTES: 
%
% 16/07/2020:
% The Experiments IV-a & IV-b in the paper are each run for 50 Monte Carlo 
% trials, but to reproduce this will take time, like hours, for all plots. 
% Just be warned.  You can see that the behavior is consistent from a 
% smaller number of trials.
%
% Back of the napkin estimate:
% 'Fig.2_ExpIVa' takes about 4 seconds per run with min-vol NMF and 5 
% seconds with minimaxNMF
% 'Fig.2_ExpIVb' takes about 20 seconds per run with min-vol NMF and 30 
% seconds with minimaxNMF
%
% So,
% 'Fig.2_ExpIVa' takes 16*50*(4+5) = ~2hrs 
% 'Fig.2_ExpIVa' takes 16*50*(20+30) = ~11hrs 
% ------------------------------------------------------------------------
path1 = genpath('src');
path2 = genpath('examples');
addpath(path1);
addpath(path2);

scen{1}     = 'Fig.1';
scen{2}     = 'Fig.2_Exp.IVa';
scen{3}     = 'Fig.2_Exp.IVb';
scen{4}     = 'Fig.3';
nRuns       = 10;
queued      = 1:4;


%% This will be updated when the code is finished.
for i = queued
    scenario = scen{i};
    if i == 1
        [ param{i}, mvOptions{i}, mmOptions{i} ] = minimaxNMF_Fig1( scenario );
    elseif i == 2 || i == 3
        [ mvResults{i}, mmResults{i}, nV{i}, param{i}, mvOptions{i}, mmOptions{i} ] = minimaxNMF_Fig2( scenario, nRuns );
    elseif i == 4
        minimaxNMF_Fig3( scenario );
    end
end
fprintf('\tFinished.\n');


