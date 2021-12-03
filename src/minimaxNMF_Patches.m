%   File: minimaxNMF_Patches.m
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
% ------------------------------------------------------------------------
% SYNTAX:
% minimaxNMF_Fig1( scenario );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This experiment compares the accuracy the endmembers extracted using 
% minimum volume NMF [2] computed in the traditional way, for the entire 
% data set, versus the proposed method that minimizes the maximum minimum
% volume error over a collection of subsample sets.
%
% Displays illustrative example corresponding to Fig. 1 in the paper.
%
% ------------------------------------------------------------------------
% INPUTS:
% scenario  String that specificies the figure to reproduce.
%           'Fig.1'
%
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
%       Not included.
%       Written by V. Leplat, see file for documentation.
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
% ------------------------------------------------------------------------
% CREATED:      01/01/2020 by Tim Marrinan
%
% LAST EDITED:  30/07/2020 by Tim Marrinan
%
% NOTES: 
% ------------------------------------------------------------------------
%clear all; clc;
function minimaxNMF_Fig1( scenario )


%% Illustrative Example
% Simultaneous hyperspectral unmixing and rare target detection via NMF minimum enclosing ball
clear variables
% scen                            = 'Fig.1';
printout                        = 'verbose';

%% Data setup
[ param, mmOptions, mvOptions ] = minimaxNMF_ScenarioSpecification( scen );
[ data, param ]                 = minimaxNMF_DataGen( param );

% Some parameters are needed in the experiment
r       = param.nTotalEMs;
nb      = param.nBands;
nr      = param.nRows;
nc      = param.nColumns;
W       = param.W;
ww      = param.patch_width;
wh      = param.patch_height;

switch lower(printout)
    case 'verbose'
        clc        
        fprintf('\n\tminimax NMF demo\n');
        fprintf('\t--------------------------------------------------------\n');
        fprintf('\tScenario:\t\t\t\t\t\t%s\n',scen);
        fprintf('\tNumber of samples:\t\t\t\t%g\n',param.nSamples);
        fprintf('\tNumber of endmembers:\t\t\t%g\n',param.nTotalEMs);
        fprintf('\tRare endmembers:\t\t\t\t%g\n',param.nRareEMs);
        fprintf('\tSamples with rare endmembers:\t%g%%\n',param.rarePercent*100);
        fprintf('\tNoise variance:\t\t\t\t\t%g\n',param.sigmaN);
        fprintf('\t--------------------------------------------------------\n');
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

% The output for this function is gonna change soon
[iterates, ~] = NMF_MEB_dual_subgrad_mod(windowed_data, r, mmOptions);


minimaxNMF_W = iterates.W{iterates.lastIter};
windowed_minimaxNMF_H = iterates.H(iterates.lastIter,:);

%% Match estimated endmembers with best fit groundtruth endmembers
% minvolNMF
[minvolEM_error,minvol_assignment] = compareWs( W,minvolNMF_W );

% minimaxNMF
[minimaxEM_error,minimax_assignment] = compareWs( W, minimaxNMF_W );


%% Compute reconstruction error
% minvolNMF
minvolX_error = norm(data.synthetic_mat-minvolNMF_W*minvolNMF_H,'fro')/norm(data.synthetic_mat,'fro');

% minimaxNMF 
minimaxNMF_H = zeros(r,nr,nc);
for j = 1 : nvert
    for i = 1 : nhor
        % Need to put the windows back in the right order
        % If samples were lost due to a mismatch in window size, they are
        % still gone here and reconstruction error will be off
        minimaxNMF_H(:,(j-1)*wh+1:j*wh,(i-1)*ww+1:i*ww) = reshape(windowed_minimaxNMF_H{((i-1)*nvert)+j},[r,wh,ww]);
    end 
end
minimaxNMF_H = reshape(minimaxNMF_H,[r,nr*nc]);
minimaxX_error = norm(data.synthetic_mat-minimaxNMF_W*minimaxNMF_H,'fro')/norm(data.synthetic_mat,'fro');

%% Display data, groundtruth, and estimated endmembers
dColor = [0.98431,0.50196,0.44706];
gtColor = [0.59608,0.30588,0.63922];
mmColor = [0,0,0];
mvColor = [0.55294,0.82745,0.78039];

fig_width = 12;
fig_height = 12;
figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
set(0, 'DefaultLineLineWidth', 2);
minvolW     = [ minvolNMF_W(:,minvol_assignment), minvolNMF_W(:,minvol_assignment(1)) ];
minimaxW    = [ minimaxNMF_W(:,minimax_assignment),minimaxNMF_W(:,minimax_assignment(1)) ];
d = plot(data.synthetic_mat(2,:), data.synthetic_mat(3,:),'o',...   
        'MarkerEdgeColor',dColor,...
        'DisplayName','Data'); 
hold on
gt = plot(W(2,:), W(3,:), 'o--','MarkerSize',10,...
        'Color',gtColor,'MarkerEdgeColor',gtColor,...
        'DisplayName','Groundtruth');
mm = plot(minimaxW(2,:),minimaxW(3,:),'d-','MarkerSize',10,...   
        'Color',mmColor,'MarkerEdgeColor',mmColor,...
        'DisplayName','Algorithm 1');
mv = plot(minvolW(2,:),minvolW(3,:),'x-','MarkerSize',10,...
        'Color',mvColor,'MarkerEdgeColor',mvColor,...
        'DisplayName','min-vol NMF [12]'); 
grid on; 
title(scen)
l = legend('Location','best');
set(l, 'Interpreter', 'latex')
hold off


%% Print results
switch lower(printout)
    case 'verbose'       
        fprintf('\n\tResults\n');
        fprintf('\t--------------------------------------------------------\n');
        fprintf('\tEndmember matrices:\n');
        fprintf('\tGroundtruth:\t\t Algorithm 1:\t\t   minvol NMF [12]:\n');
        fprintf('\t|');
        fprintf( num2str(param.W(1,:)) );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minimaxW(1,1:4),'% 1.2f') );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minvolW(1,1:4),'% 1.2f') );
        fprintf('|\n');
        
        fprintf('\t|');
        fprintf( num2str(param.W(2,:)) );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minimaxW(2,1:4),'% 1.2f') );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minvolW(2,1:4),'% 1.2f') );
        fprintf('|\n');
        
        fprintf('\t|');
        fprintf( num2str(param.W(3,:)) );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minimaxW(3,1:4),'% 1.2f') );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minvolW(3,1:4),'% 1.2f') );
        fprintf('|\n');
        
        fprintf('\t|');
        fprintf( num2str(param.W(4,:)) );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minimaxW(4,1:4),'% 1.2f') );
        fprintf('|');
        fprintf('\t|');
        fprintf( num2str(minvolW(4,1:4),'% 1.2f') );
        fprintf('|\n\n');
        
        fprintf('\tEndmember matrix error:\n');
        fprintf('\tAlgorithm 1:\t\t\t%1.2f%%\n',100*minimaxEM_error);
        fprintf('\tMinvol NMF [12]:\t\t%1.2f%%\n\n',100*minvolEM_error);
        
        fprintf('\tReconstruction error:\n');
        fprintf('\tAlgorithm 1:\t\t\t%1.2f%%\n',100*minimaxX_error);
        fprintf('\tMinvol NMF [12]:\t\t%1.2f%%\n',100*minvolX_error);
        fprintf('\t--------------------------------------------------------\n');
    otherwise
end
