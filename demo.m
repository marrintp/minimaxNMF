%   File: demo.m
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
% This script demonstrates how to use the minimaxNMF package by generating
% data from convex combinations of spectral signatures + noise, and trying
% to recover the original spectral endmembers.  The results are compared to
% minimum volume NMF [3] computed in the traditional way, for the entire 
% data set, versus the proposed method that minimizes the maximum minimum
% volume error over a collection of subsample sets.
%
% ------------------------------------------------------------------------
% INPUTS:
% n/a
% ------------------------------------------------------------------------
% OUTPUTS: 
% n/a
%
% ------------------------------------------------------------------------
% DEPENDENCIES:
% [a] minimaxNMF_ScenarioSpecification.m
%       Not included.
%       Sets scenario parameters.
%       Written by T. Marrinan; see file for documentation.
%
% [b] minimaxNMF_DataGen.m
%       Not included.
%       Written by T. Marrinan, see file for documentation.
%
% [c] minimaxNMF.m
%       Not included.
%       Written by T. Marrinan, see file for documentation.
%
% [d] minvolNMF.m
%       Not included within. 
%       Written by Valentin Leplat, Andersen MS Ang, and Nicolas Gillis. 
%       See file for documentation and reference [2] for details.
%
% [e] CompareWs.m
%       Not included within. 
%       Written by Valentin Leplat, Andersen MS Ang, and Nicolas Gillis. 
%       See file for documentation and reference [3] for details.
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
% [2]   Clark, Roger N., Gregg A. Swayze, A. Gallagher, T. V. V. King, 
%       and W. M. Calvin. "The us geological survey digital spectral 
%       library." US Geological Survey Open File Report (1993): 93-592.
%
% [2]    Leplat, Valentin, Andersen MS Ang, and Nicolas Gillis. 
%       "Minimum-volume rank-deficient nonnegative matrix factorizations." 
%       In ICASSP 2019-2019 IEEE International Conference on Acoustics, 
%       Speech and Signal Processing (ICASSP), pp. 3402-3406. IEEE, 2019.
%
% ------------------------------------------------------------------------
% CREATED:      21/03/2020 by Tim Marrinan
%
% LAST EDITED:  08/07/2020 by Tim Marrinan
%
% NOTES: 
%
% ------------------------------------------------------------------------
path1 = genpath('src');
path2 = genpath('examples');
addpath(path1);
addpath(path2);



%% Data setup
scenario                        = 'custom'; % Choose scenario
[ param, mmOptions, mvOptions ] = minimaxNMF_ScenarioSpecification( scenario );
[ data ]                        = minimaxNMF_DataGen( param );

% Some parameters are needed in the experiment
r           = param.nTotalEMs;
nb          = param.nBands;
nr          = param.nRows;
nc          = param.nColumns;
W           = param.W;
ww          = param.patch_width;
wh          = param.patch_height;
nEMs        = param.nTotalEMs;
EM_names    = param.EM_names;
usedb       = find(param.usedBands);
printout    = param.printout;  

switch lower(printout)
    case 'verbose'
        clc        
        fprintf('\n\tminimax NMF demo\n');
        fprintf('\t-------------------------------------------------------------\n');
        fprintf('\tScenario:\t\t\t\t\t\t%s\n',scenario);
        fprintf('\tNumber of samples:\t\t\t\t%g\n',param.nSamples);
        fprintf('\tNumber of endmembers:\t\t\t%g\n',param.nTotalEMs);
        fprintf('\tRare endmembers:\t\t\t\t%g\n',param.nRareEMs);
        fprintf('\tSamples with rare endmembers:\t%g%%\n',param.rarePercent*100);
        fprintf('\tNoise variance:\t\t\t\t\t%g\n',param.sigmaN);
        fprintf('\t-------------------------------------------------------------\n');
    otherwise
end

%% Initial solution via [3]
normalize = true;
[J,Hfull] = SNPA(data.synthetic_mat,r,normalize);
if length(J) < r
    warning('SNPA recovered less than r basis vectors.');
    warning('The data poins have less than r vertices.');
    r = length(J);
    fprintf('The new value of r is %2.0d.\n',r);
end

%% Create patches for minimax NMF
% Some samples will not be used if the patches don't tessellate the image.
%
% In reality, we have not investigated the relationship between the spatial
% location of pixels in the same "patch" so this whole tessellation thing
% is more or less unnecessary.  But it it was what do for the examples in
% the paper, so it is what we do here.
%
% For your own experiments, you can just make patches of an appropriate
% size so that: # pixels per patch x # patches = nSamples.  The pixels in
% each patch can be randomly selected without replacement or whatever.  For
% HSI applications we might expect the rare pixels to be collocated, but it
% remains to be seen whether endmember extraction improves with all rare
% pixels in one patch or one rare pixel in many patches.  I suspect the
% former is better.
nhor    = floor(nc/ww); 
nvert   = floor(nr/wh);
nwin    = nhor*nvert;

% Chop image into patches
windowed_data = cell(nwin,1);
for i = 1 : nhor
    for j = 1 : nvert
        temp = data.synthetic_cube(:,(j-1)*wh+1:j*wh,(i-1)*ww+1:i*ww);
        windowed_data{((i-1)*nvert)+j} = reshape(temp,[nb,wh*ww]);
    end 
end

% Find patches of initial H corresponding to patches of data
holder = reshape(Hfull,[r,nr,nc]);
windowed_H = cell(nwin,1);
for i = 1 : nhor
    for j = 1 : nvert
        temp = holder(:,(j-1)*wh+1:j*wh,(i-1)*ww+1:i*ww);
        windowed_H{((i-1)*nvert)+j} = reshape(temp,[r,wh*ww]);
    end 
end

%% Minvol NMF
% Initial solution for minvolNMF
mvOptions.W = data.synthetic_mat(:,J);
mvOptions.H = Hfull;

[minvolNMF_W,minvolNMF_H,~,~,~,~] = minvolNMF(data.synthetic_mat,r,mvOptions);


%% Minimax NMF
% Initial solution for minimaxNMF
mmOptions.W = data.synthetic_mat(:,J);
mmOptions.H = windowed_H;

[minimaxNMF_W, windowed_minimaxNMF_H, iterates, options] = minimaxNMF(windowed_data, r, mmOptions);

%% Match estimated endmembers with best fit groundtruth endmembers
% minvolNMF
[minvolEM_error,minvol_assignment] = compareWs( W,minvolNMF_W );

% minimaxNMF
[minimaxEM_error,minimax_assignment] = compareWs( W, minimaxNMF_W );


%% Compute reconstruction error
% minvolNMF
minvolX_error = norm(data.synthetic_mat-minvolNMF_W*minvolNMF_H,'fro')/norm(data.synthetic_mat,'fro');

% minimaxNMF (first reconstruct H)
minimaxNMF_H = zeros(r,nr,nc);
for j = 1 : nvert
    for i = 1 : nhor
        % If samples were lost due to a mismatch in window size, they are
        % still gone here and reconstruction error will be off
        minimaxNMF_H(:,(j-1)*wh+1:j*wh,(i-1)*ww+1:i*ww) = reshape(windowed_minimaxNMF_H{((i-1)*nvert)+j},[r,wh,ww]);
    end 
end
minimaxNMF_H = reshape(minimaxNMF_H,[r,nr*nc]);
minimaxX_error = norm(data.synthetic_mat-minimaxNMF_W*minimaxNMF_H,'fro')/norm(data.synthetic_mat,'fro');

%% Plotting
fig_width = 18;
fig_height = 12;

% I want the plots to show the missing bands
Wnan                    = NaN(224,nEMs);
Wnan(usedb,:)           = W;

mmW                     = NaN(224,nEMs);
mmW(usedb,:)            = minimaxNMF_W(:,minimax_assignment);

mvW                     = NaN(224,nEMs);
mvW(usedb,:)            = minvolNMF_W(:,minvol_assignment);

temp                    = NaN(224,1);
temp(usedb)             = usedb;
usedb                   = temp;


%%

switch lower(scenario)
    % -------------------------------------------------------------------------
    case 'custom'
        % Plot groundtruth spectra
        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        plot(usedb,Wnan(:,1),'DisplayName',strrep(EM_names(1,:), '_', '\_'));
        hold on;
        for i = 2 : nEMs
            plot(usedb,Wnan(:,i),'DisplayName',strrep(EM_names(i,:), '_', '\_'));
        end
        lgd = legend({},'Location','southoutside','FontSize',10);
        lgd.NumColumns = 2;
        title('Groundtruth endmember spectra','FontSize',12 );
        
        % Plot estimated endmembers from Algorithm 1
        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        plot(usedb,mmW(:,1),'DisplayName',strrep(EM_names(1,:), '_', '\_'));
        hold on;
        for i = 2 : nEMs
            plot(usedb,mmW(:,i),'DisplayName',strrep(EM_names(i,:), '_', '\_'));
        end
        lgd = legend({},'Location','southoutside','FontSize',10);
        lgd.NumColumns = 2;
        title('Algorithm 1: Estimated endmember spectra','FontSize',12 );
        
        % Plot estimated endmembers from min-vol NMF [3]
        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        plot(usedb,mvW(:,1),'DisplayName',strrep(EM_names(1,:), '_', '\_'));
        hold on;
        for i = 2 : nEMs
            plot(usedb,mvW(:,i),'DisplayName',strrep(EM_names(i,:), '_', '\_'));
        end
        lgd = legend({},'Location','southoutside','FontSize',10);
        lgd.NumColumns = 2;
        title('min-vol NMF [12]: estimated endmember spectra','FontSize',12 );
        
        % Compare rare endmember estimates
        for i = 1 : param.nRareEMs
            figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
            set(0, 'DefaultLineLineWidth', 2);
            plot(usedb,Wnan(:,end-i+1),'DisplayName','Groundtruth');
            hold on;
            
            mmSA = acos(W(:,end-i+1)'*minimaxNMF_W(:,minimax_assignment(end-i+1))...
                /(norm(W(:,end-i+1))*norm(minimaxNMF_W(:,minimax_assignment(end-i+1)))));
            mmSTR=sprintf('%.2f',mmSA);
            plot(usedb,mmW(:,end-i+1),'DisplayName',strcat(['Algorithm 1: $\theta$ = ',mmSTR]));
            
            mvSA = acos(W(:,end-i+1)'*minvolNMF_W(:,minvol_assignment(end-i+1))...
                /(norm(W(:,end-i+1))*norm(minvolNMF_W(:,minvol_assignment(end-i+1)))));
            mvSTR=sprintf('%.2f',mvSA);
            plot(usedb,mvW(:,end-i+1),'DisplayName',strcat(['min-vol NMF [12]: $\theta$ = ',mvSTR]));
        
            lgd = legend({},'Location','southoutside','FontSize',10,'Interpreter','latex');
            lgd.NumColumns = 3;
            title('Rare endmember comparison','FontSize',12 );
        end
end        

%% Print results
switch lower(printout)
    case 'verbose'       
        fprintf('\n\tResults\n');
        fprintf('\t-------------------------------------------------------------\n');        
        fprintf('\tEndmember matrix error:\n');
        fprintf('\tAlgorithm 1:\t\t\t%1.2f%%\n',100*minimaxEM_error);
        fprintf('\tmin-vol NMF [12]:\t\t%1.2f%%\n\n',100*minvolEM_error);
        
        fprintf('\tReconstruction error:\n');
        fprintf('\tAlgorithm 1:\t\t\t%1.2f%%\n',100*minimaxX_error);
        fprintf('\tmin-vol NMF [12]:\t\t%1.2f%%\n',100*minvolX_error);
        fprintf('\t-------------------------------------------------------------\n');
    otherwise
end
fprintf('\n\tdone.\n');






