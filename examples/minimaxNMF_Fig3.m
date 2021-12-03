%   File: minimaxNMF_Fig3.m
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
% minimaxNMF_Fig3( scenario );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This function just displays the groundtruth endmembers used for 
% Experiment IV-b in the paper. The spectra come from the USGS spectral
% library [2].
%
% Displays illustrative example corresponding to Fig. 3 in the paper.
%
% ------------------------------------------------------------------------
% INPUTS:
% scenario  String that specificies the figure to reproduce.
%           'Fig.3'
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
% ------------------------------------------------------------------------
% CREATED:      02/03/2020 by Tim Marrinan
%
% LAST EDITED:  04/08/2020 by Tim Marrinan
%
% NOTES: 
% ------------------------------------------------------------------------
%clear all; clc;
function minimaxNMF_Fig3( scenario )

%% Experimental setup:

%% Data setup
[ param, ~, ~ ] = minimaxNMF_ScenarioSpecification( scenario );

% Some parameters are needed for the plot
W           = param.W;
nEMs        = param.nTotalEMs;
EM_names    = param.EM_names;
usedb       = find(param.usedBands);
printout    = param.printout;




switch lower(printout)
    case 'verbose'       
        fprintf('\n\tminimax NMF experiment\n');
        fprintf('\t-------------------------------------------------------------\n');
        fprintf('\tScenario:\t\t\t\t\t\t%s\n',scenario);
        fprintf('\t-------------------------------------------------------------\n');
    otherwise
end

%% Plotting
mycolor{1} = [0.00000,0.44700,0.74100];
mycolor{2} = [0.85000,0.32500,0.09800];
mycolor{3} = [0.92900,0.69400,0.12500];
mycolor{4} = [0.49400,0.18400,0.55600];
mycolor{5} = [0.46600,0.67400,0.18800];
mycolor{6} = [0.30100,0.74500,0.93300];
mycolor{7} = [0.63500,0.07800,0.18400];
mycolor{8} = [0.99216,0.70588,0.38431];

fig_width = 18;
fig_height = 12;

% I want the plots to show the missing bands
Wnan                    = NaN(224,nEMs);
Wnan(usedb,:)           = W;

temp                    = NaN(224,1);
temp(usedb)             = usedb;
usedb                   = temp;


switch lower(scenario)
    % -------------------------------------------------------------------------
    case 'fig.3'
        figure('units','centimeters','Position',[4,4,fig_width,fig_height]);
        set(0, 'DefaultLineLineWidth', 2);
        
        % Dominant endmembers
        plot(usedb,Wnan(:,1),'-','Color',mycolor{1},'DisplayName',strrep(EM_names(1,:), '_', '\_'));
        hold on;
        for i = [3 5 7]
            plot(usedb,Wnan(:,i),'-','Color',mycolor{i},'DisplayName',strrep(EM_names(i,:), '_', '\_'));
        end
        
        % Rare endmember #1
        plot(usedb,Wnan(:,9),'--','Color',mycolor{2},'DisplayName',EM_names(9,:));
        
        % Dominant endmembers
        for i = [2 4 6 8]
            plot(usedb,Wnan(:,i),'-','Color',mycolor{i},'DisplayName',strrep(EM_names(i,:), '_', '\_'));
        end
        
        % Rare endmember #2
        plot(usedb,Wnan(:,10),'--','Color',mycolor{3},'DisplayName',EM_names(10,:));
        lgd = legend({},'Location','southoutside','FontSize',10);
        lgd.NumColumns = 2;
        title('Exp. IV-b: Groundtruth endmember spectra','FontSize',12 );
end        
        
        fprintf('\n\tdone.\n');