%   File: minimaxNMF_Fig2.m
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
% [ mvResults, mmResults, nV ] = minimaxNMF_Fig2( scenario, nRuns );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This experiment compares the accuracy the endmembers extracted using 
% minimum volume NMF [12] computed in the traditional way, for the entire 
% data set, versus Algorithm 1 [1] that minimizes the maximum error over a 
% collection of subsample sets.
%
% Displays illustrative example corresponding to Fig. 2 in the paper.
%
% ------------------------------------------------------------------------
% INPUTS:
% scenario  String that specificies the figure to reproduce.
%           'Fig.2_Exp.IVa'/'Fig.2_Exp.IVb'
%
% nRuns     Scalar that specficies the number of Monte Carlo trials to run.
%           (Paper uses 50 - but it will take hours)
%
% ------------------------------------------------------------------------
% OUTPUTS: 
% mvResults, mmResults, nV
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
% [2]    Leplat, Valentin, Andersen MS Ang, and Nicolas Gillis. 
%       "Minimum-volume rank-deficient nonnegative matrix factorizations." 
%       In ICASSP 2019-2019 IEEE International Conference on Acoustics, 
%       Speech and Signal Processing (ICASSP), pp. 3402-3406. IEEE, 2019.
%
% ------------------------------------------------------------------------
% CREATED:      03/01/2020 by Tim Marrinan
%
% LAST EDITED:  31/07/2020 by Tim Marrinan
%
% NOTES: 
% ------------------------------------------------------------------------
%clear all; clc;
function [ mvResults, mmResults, nV, param, mvOptions, mmOptions ] = minimaxNMF_Fig2( scenario, nRuns )

%% Data setup
[ param, mmOptions, mvOptions ] = minimaxNMF_ScenarioSpecification( scenario );

% These fields just need to exist so I can remove them
mvOptions.W = 'n/a';
mvOptions.H = 'n/a';
mmOptions.W = 'n/a';
mmOptions.H = 'n/a';

% Some parameters are needed in the experiment
nb          = param.nBands;
nr          = param.nRows;
nc          = param.nColumns;
W           = param.W;
ww          = param.patch_width;
wh          = param.patch_height;
printout    = param.printout;  

% Noise variance is the independent variable in the experiments
% nV      =  [  1, .5, .1, 0.05, 0.01, 0.005, 0.001 ]; % Noise variance 
nV      =  [  1, .8, .6, .4, .2, .1, 0.08, 0.06 0.04,0.02, 0.01, 0.008,...
           0.006, 0.004, 0.002, 0.001 ]; % Noise variance 
nVar    = size(nV,2);

% Initialize variables for results
mvResults.EM_error = zeros(nRuns,nVar);
mvResults.X_error= zeros(nRuns,nVar);

mmResults.EM_error = zeros(nRuns,nVar);
mmResults.X_error = zeros(nRuns,nVar);

switch lower(printout)
    case 'verbose'       
        fprintf('\n\tminimax NMF experiment\n');
        fprintf('\t-------------------------------------------------------------\n');
        fprintf('\tScenario:\t\t\t\t\t\t%s\n',scenario);
        fprintf('\tNumber of samples:\t\t\t\t%g\n',param.nSamples);
        fprintf('\tNumber of endmembers:\t\t\t%g\n',param.nTotalEMs);
        fprintf('\tRare endmembers:\t\t\t\t%g\n',param.nRareEMs);
        fprintf('\tSamples with rare endmembers:\t%g%%\n',param.rarePercent*100);
        fprintf('\tNumber of variances to test:\t%1.0d\n',nVar);
        fprintf('\tNumber of runs per variance:\t%1.0d\n',nRuns);
        fprintf('\t-------------------------------------------------------------\n');
   case 'terse'
       clc        
        fprintf('\n\tminimax NMF experiment\n');
        fprintf('\t-------------------------------------------------------------\n');
        fprintf('\tScenario:\t\t\t\t\t\t%s\n',scenario);
        fprintf('\t-------------------------------------------------------------\n');
   case 'none'
    otherwise
end

for iVar = 1 : nVar
    switch lower(printout)
        case 'verbose'       
            fprintf('\n\tIndependent variable:\t\t\t%1.0d/%1.0d\n',iVar,nVar);
            fprintf('\tTrial:\t');
       case 'terse'
       case 'none'
        otherwise
    end
    for iter = 1 : nRuns
        switch lower(printout)
            case 'verbose'       
                fprintf('%1.0d\t',iter);
           case 'terse'
           case 'none'
            otherwise
        end

        
        %% Generate new data for each run
        [ param, mmOptions, mvOptions ] = minimaxNMF_ScenarioSpecification( scenario );
        % Update noise variance
        param.sigmaN = nV(iVar);
        % These fields just need to exist so I can remove them
        mvOptions.W = 'n/a';
        mvOptions.H = 'n/a';
        mmOptions.W = 'n/a';
        mmOptions.H = 'n/a';
        [ data ]    = minimaxNMF_DataGen( param );
        
        % Reset initial iterates
        r           = param.nTotalEMs;
        mvOptions   = rmfield(mvOptions,{'W','H'});
        mmOptions   = rmfield(mmOptions,{'W','H'});
        
        %% Initial solution via [3]
        normalize   = true;
        [J,Hfull]   = SNPA(data.synthetic_mat,r,normalize);
        if length(J) < r
            warning('SNPA recovered less than r basis vectors.');
            warning('The data poins have less than r vertices.');
            r       = length(J);
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
        
        tstart= tic;
        [minvolNMF_W,minvolNMF_H,~,~,~,~] = minvolNMF(data.synthetic_mat,r,mvOptions);
        mvResults.time(iter,iVar) = toc(tstart);

        %% Minimax NMF
        % Initial solution for minimaxNMF
        mmOptions.W = data.synthetic_mat(:,J);
        mmOptions.H = windowed_H;
        
        tstart= tic;
        [minimaxNMF_W, windowed_minimaxNMF_H, ~, mmOptions] = minimaxNMF(windowed_data, r, mmOptions);
        mmResults.time(iter,iVar) = toc(tstart);
        %% Match estimated endmembers with best fit groundtruth endmembers
        % minvolNMF
        [mvResults.EM_error(iter,iVar),~] = compareWs( W, minvolNMF_W );

        % minimaxNMF
        [mmResults.EM_error(iter,iVar),~] = compareWs( W, minimaxNMF_W );
        
        %% Compute reconstruction error
        % minvolNMF
        mvResults.X_error(iter,iVar) = norm(data.synthetic_mat-minvolNMF_W*minvolNMF_H,'fro')/norm(data.synthetic_mat,'fro');

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
        mmResults.X_error(iter,iVar) = norm(data.synthetic_mat-minimaxNMF_W*minimaxNMF_H,'fro')/norm(data.synthetic_mat,'fro');
    end
    fprintf('\n\tTime for per ind. var.:\t\t\t%1.0d\n',sum(mmResults.time(:,iVar)+mvResults.time(:,iVar)));
end

%% Plotting
dColor = [0.98431,0.50196,0.44706];
gtColor = [0.59608,0.30588,0.63922];
mmColor = [0,0,0];
mvColor = [0.55294,0.82745,0.78039];

fig_width = 18;
fig_height = 12;

switch lower(scenario)
    % -------------------------------------------------------------------------
    case 'fig.2_exp.iva'
        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        loglog(nV,mean(mmResults.EM_error),'d-','MarkerSize',10,...   
        'Color',mmColor,'MarkerEdgeColor',mmColor,...
        'DisplayName','Algorithm 1');
        hold on;
        loglog(nV,mean(mvResults.EM_error),'x-','MarkerSize',10,...
        'Color',mvColor,'MarkerEdgeColor',mvColor,...
        'DisplayName','min-vol NMF [12]')
        xlabel('$\sigma_N^2$','Interpreter','latex')
        ylabel('$\frac{\|W-W_{est}\|_F}{\|W\|_F}$','Interpreter','latex')
        title({['Exp. IV-a: Mean relative endmember matrix error']});
        legend('Location','northwest')


        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        semilogx(nV,mean(mmResults.X_error),'d-','MarkerSize',10,...   
        'Color',mmColor,'MarkerEdgeColor',mmColor,...
        'DisplayName','Algorithm 1');
        hold on;
        semilogx(nV,mean(mvResults.X_error),'x-','MarkerSize',10,...
        'Color',mvColor,'MarkerEdgeColor',mvColor,...
        'DisplayName','min-vol NMF [12]')
        xlabel('$\sigma_N^2$','Interpreter','latex')
        ylabel('$\frac{\|X-W_{est}H_{est}\|_F}{\|X\|_F}$','Interpreter','latex')
        title({'Exp. IV-a: Mean relative reconstruction error'});
        legend('Location','northwest')
        
    % -------------------------------------------------------------------------
    case 'fig.2_exp.ivb'
        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        semilogx(nV,mean(mmResults.EM_error),'d--','MarkerSize',10,...   
        'Color',gtColor,'MarkerEdgeColor',gtColor,...
        'DisplayName','Algorithm 1');
        hold on;
        semilogx(nV,mean(mvResults.EM_error),'x--','MarkerSize',10,...
        'Color',dColor,'MarkerEdgeColor',dColor,...
        'DisplayName','min-vol NMF [12]')
        xlabel('$\sigma_N^2$','Interpreter','latex')
        ylabel('$\frac{\|W-W_{est}\|_F}{\|W\|_F}$','Interpreter','latex')
        title({'Exp. IV-b: Mean relative endmember matrix error'});
        legend('Location','northwest')

        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        loglog(nV,mean(mmResults.X_error),'d-','MarkerSize',10,...   
        'Color',gtColor,'MarkerEdgeColor',gtColor,...
        'DisplayName','Algorithm 1');
        hold on;
        loglog(nV,mean(mvResults.X_error),'x-','MarkerSize',10,...
        'Color',dColor,'MarkerEdgeColor',dColor,...
        'DisplayName','min-vol NMF [12]')
        xlabel('$\sigma_N^2$','Interpreter','latex')
        ylabel('$\frac{\|X-W_{est}H_{est}\|_F}{\|X\|_F}$','Interpreter','latex')
        title({'Exp. IV-b: Mean relative reconstruction error'});
        legend('Location','northwest')
        
end


tp = mvResults.EM_error - mmResults.EM_error;
mmWinner = sum(sum(tp>0))./(nVar*nRuns);

tq = mmResults.X_error - mvResults.X_error;
mvWinner = sum(sum(tq>0))./(nVar*nRuns);

switch lower(printout)
    case 'verbose'       
        fprintf('\n\tResults\n');
        fprintf('\t-------------------------------------------------------------\n');   
        fprintf('\tPercent of trials in which...\n');
        fprintf('\tAlgorithm 1 has lower endmember matrix error:\t\t%1.0d/%1.0d=%1.1f%%\n',sum(sum(tp>0)), nVar*nRuns, 100*mmWinner);
        
        fprintf('\tmin-vol NMF [12] has lower reconstruction error:\t%1.0d/%1.0d=%1.1f%%\n',sum(sum(tq>0)), nVar*nRuns, 100*mvWinner);
        fprintf('\t-------------------------------------------------------------\n');
        fprintf('\ndone.\n');
    case 'terse'
        fprintf('\ndone.\n');
    case 'none'
    otherwise
end
fprintf('\n\tdone.\n');
