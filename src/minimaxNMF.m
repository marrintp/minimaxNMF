%   File: minimax_NMF.m
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
%       with the University of Paderborn and is in general not permitted.
%
%   4.) Modifications or contributions to the software must be
%       published under this license. The University of Paderborn
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
% [ W, H, iterates, options ] = minimax_NMF( data, r, options )
%
% ------------------------------------------------------------------------
% OVERVIEW: 
% We wish to solve the low-rank nonnegative factorization problem,
%
%       min               max   || X_i - (W)(H_i) ||_F
%    W in R_{>0}^{b x r}   i 
%    H in R_{>0}^{r x p}
%
% for some data, {X_i} in R_{>0}^{b x p} for i = 1, ... , n.  
%
% This code solves the Lagrange dual of this problem, which can be written
% as
%
%       max         g(lambda)
%   lambda in R^m
%
%       s.t.        g(lambda)    =   min                sum ( lambda_i || X_i - (W)(H_i) ||_F )
%                                W in R_{>0}^{b x r}    i
%                                H_i in R_{>0}^{r x m}
%                   lambda_i    >=   0 for i = 1,...,n
%                || lambda ||_1  =   1.
% 
% This dual can be solved via an 'approximate subgradient' algorithm, where 
% the scare quotes are there because we can't actually compute the
% subgradient exactly because of the NP-hard NMF subproblems.  For a fixed 
% lambda, let
%
%   f_i(W, H_i) = lambda_i || X_i - (W)(H_i) ||_F. 
%
% A vector included in the subdifferential is then the vector of function
% values at the minimum of the sum.  That is 
%
%   [f_1(W*,H*_1), f_2(W*,H*_2), ..., f_n(W*,H*_n)]^T is in g(lambda).
%
% We approximate this vector by alternating between fixing H_i and 
% minimizing f_i with respect to W, and fixing W and minimizing f_i with
% respect to H_i. After a few iterations of each, we use the values of W
% and H_i for i = 1, ..., n to compute this approximate subgradient, update
% lambda with the approximate subgradient, and repeat the process.
%
% The subgradient algorithm in general terms is then:
% (1) Intialize the weights
%       lambda^(0) = unit vector with equal elements
%       t = 0
%
% Repeat (2) - (6) until convergence:
%
%   At step t repeat (2) & (3) for a fixed number of iterations:
%   (2) W* =  argmin           lambda_i || X_i - (W)(H_i) ||_F
%           W in R_{>0}^{b x r}  
%   (3) H*_i =  argmin           lambda_i || X_i - (W)(H_i) ||_F
%           H_i in R_{>0}^{r x m}  
%           for i = 1, ..., n
%
% (4) Pick a subgradient
%       del g^(t) = [f_1^(t)(W*,H*_1), f_2^(t)(W*,H*_2), ..., f_n^(t)(W*,H*_n)]^T
%
% (5) Take a step in the direction of the subgradient with stepsize alpha
%       \hat{lambda}^(t+1) = lambda^(t) - alpha del g^(t)
%
% (6) Project onto the feasible set using [a]
%       lambda^(t+1) = P( \hat{lambda^(t+1)} )
%
% ------------------------------------------------------------------------
% INPUTS:
% data      nPatches x 1 cell array of matrices - Cells contain 
%           nonnegative matrices of size b x m where (nPatches x m = p, the
%           total number of samples)
%
% r         scalar - number of endmembers to extract
%
% options   MATLAB structure - contains the parameters for algorithm. 
%           If not included, default parameters will be used.
%   Fields and [default values]:
%   'tol'               scalar - convergence threshold for the duality gap
%                       and for lack of movement of the dual variable 
%                       [10^(-4)]
%   'lambda'            nPatches x 1 vector - initial value of the dual variable
%                       [1/nPatches, 1/nPatches, ... , 1/nPatches]
%   'beta'              scalar - volume penalty parameter
%                       [0.001]
%   'delta'             scalar - eigenvalue weight
%                       [0.1]
%   'maxiter'           scalar - maximum number of subgradient steps
%                       [100]
%   'inneriter'         scalar - maximum iterations of fast gradient 
%                       method per call
%                       [25]
%   'W'                 b x r matrix - initial iterate of the endmember 
%                       matrix
%                       [default value computed with SNPA using [2]]
%   'H'                 nPatches x 1 cell array of matrices - initial
%                       value of abundance matrices.  each cell is a matrix
%                       of size r x m
%                       [default value computed with SNPA using [2]]
%   'a'                 scalar - step length coefficient
%                       [ 1 / min_i || X_i - W H_i||_F^2 ]
%   'display'           boolean - flag that determines whether to show the
%                       progress of the fast projected gradient method of
%                       not
%                       [false]
%
% ------------------------------------------------------------------------
% OUTPUTS:
% 'W'       b x r matrix - final iterate of the endmember matrix
%
% 'H'       nPatches x 1 cell array of matrices - final value of abundance 
%           matrices.  each cell is a matrix of size r x m
%
% iterates  MATLAB structure - contains the value the variables in the dual
%           subgradient algorithm at each iteration.
%   Fields:
%   'W'                 maxiter x 1 cell array - cell t contains the 
%                       endmember matrix at iteration t, a b x r matrix.
%   'H'                 maxiter x nPatches cell array - cell (t,i) contains 
%                       the ith abundance matrix at iteration t, an r x m
%                       matrix.
%   'lambda'            maxiter x nPatches matrix - row t of this matrix
%                       contains the dual variable at iteration t
%   'g'                 maxiter x nPatches matrix - row t of this matrix
%                       contains the subgradient at iteration t.
%   'dual_cost'         maxiter x 1 vector - row t of this vector
%                       contains the dual cost at iteration t
%   'primal_cost'       maxiter x 1 vector - row t of this vector
%                       contains the primal cost at iteration t
%   'e'                 maxiter x 1 vector - row t of this vector contains
%                       the objective value of the min-volume problem,
%                       e = ||X - WH||_F^2 + beta*logdet( W^T W + delta I),
%                       at iteration t
%   'err1'              maxiter x 1 vector - row t of this vector contains
%                       the reconstruction error ||X - WH||_F^2 at
%                       iteration t
%   'err2'              maxiter x 1 vector - row t of this vector contains 
%                       the volume of the convex hull of the columns of W
%                       and the origin, logdet( W^T W + delta I), at
%                       iteration t
%
% options   MATLAB structure - contains the parameters for algorithm. 
%           Only differs from input structure if fields were left blank.
% ------------------------------------------------------------------------
% DEPENDENCIES:
%
% [a] SimplexProj.m
%       Not included within. 
%       Written by Nicolas Gillis. 
%       See file for documentation and reference [2] for details.
%
% [b] FGMfcnls.m
%       Not included within. 
%       Written by Valentin Leplat, Andersen MS Ang, and Nicolas Gillis. 
%       See file for documentation and reference [3] for details.
%
% [c] FGMgpnonneg.m
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
% [2]   Gillis, N., 2014. Successive nonnegative projection algorithm for 
%       robust nonnegative blind source separation. SIAM Journal on Imaging
%       Sciences, 7(2), pp.1420-1450.
%
% [3]    Leplat, Valentin, Andersen MS Ang, and Nicolas Gillis. 
%       "Minimum-volume rank-deficient nonnegative matrix factorizations." 
%       In ICASSP 2019-2019 IEEE International Conference on Acoustics, 
%       Speech and Signal Processing (ICASSP), pp. 3402-3406. IEEE, 2019.
%
% ------------------------------------------------------------------------
% CREATED:      01/01/2020 by Tim Marrinan
%
% LAST EDITED:  07/08/2020 by Tim Marrinan
%
% NOTES: 
%
% ------------------------------------------------------------------------

function [W, H, iterates, options] = minimaxNMF(data, r, options)
% assumes data is a cell array of size nDataSets x 1 and each cell
% is an orthonormal matrix of size nSamples x nDims

% path to called functions
myFolder = 'C:\Users\Timothy\Documents\MATLAB\minimaxNMF\';
p = genpath(myFolder);
addpath(p);

%% Optional inputs and default parameter values
    [nPatches,~] = size(data);
    [nr,nc] = size(data{1});
    cols = zeros(nPatches,1);
    for i = 1 : nPatches
        cols(i) = size(data{i},2);
    end
    options.r = r;
    if nargin <= 2
        options = [];
    end
    if ~isfield(options,'tol') % Duality gap threshold
        options.tol = 10^(-4);
        tol = options.tol;
    else
        tol = options.tol;
    end
    if ~isfield(options,'lambda') % initial value of the dual variable
        options.lambda = (1/nPatches).*(ones(1,nPatches));
        lambda = options.lambda;
    else
        lambda = options.lambda;
    end
    if ~isfield(options,'beta') % volume penalty parameter
        options.beta=0.001;
        beta = options.beta;
    else
        beta = options.beta;
    end
    if ~isfield(options,'delta') % eigenvalue weight
        options.delta = 0.1;
        delta=options.delta;
    else
        delta=options.delta;
    end
    if ~isfield(options,'maxiter') % maximum number of subgradient steps
        options.maxiter=100;
        maxiter = options.maxiter;
    else
        maxiter = options.maxiter;
    end
    if ~isfield(options,'inneriter') % maximum iterations of fast gradient method per call
        options.inneriter=25;
        inneriter = options.inneriter;
    else
        inneriter = options.inneriter;
    end
    if isfield(options,'W') && isfield(options,'H') % initial primal variables
        W = options.W;
        H = options.H;
    else
        % Use SNPA [2] if no intial solution is provided
        M = [data{:}];
        normalize = true;
        [J,Hfull] = SNPA(M,r,normalize);
        options.W = M(:,J);
        options.H = cell(nPatches,1);
        for i = 1 : nPatches
            options.H{i} = Hfull(:,sum(cols(1:i-1))+1:sum(cols(1:i)));
        end
        W = options.W;
        H = options.H;
        if length(J) < r
            warning('SNPA recovered less than r basis vectors.');
            warning('The data poins have less than r vertices.');
            r = length(J);
            fprintf('The new value of r is %2.0d.\n',r);
        end
    end
    if ~isfield(options,'a')
        initial_steps = zeros(nPatches,1);
        for i = 1 : nPatches
           initial_steps(i) = norm(data{i} - options.W*options.H{i},'fro')^2; 
        end
        options.a = 2/min(initial_steps);
        a = options.a;
        specifiedA = false;
    else
        specifiedA = true;
        a = options.a;
    end
    if ~isfield(options,'display')
        options.display=1;
    end
    

%% Initialize variables
    iterates.W = cell(maxiter,1);
    iterates.H = cell(maxiter,nPatches);
    iterates.lambda = zeros(maxiter,nPatches);
    iterates.g = zeros(maxiter,nPatches);
    iterates.primal_cost = zeros(maxiter,1);
    iterates.dual_cost = zeros(maxiter,1);
    iterates.e = zeros(maxiter,1);
    iterates.err1 = zeros(maxiter,1);
    iterates.err2 = zeros(maxiter,1);
    
    iterates.W{1} = W;
    iterates.H(1,:) = H(:);
    iterates.lambda(1,:) = lambda; 

%% Compute initial costs
    g = zeros(1,nPatches);
    for i = 1 : nPatches
        g(i) = norm(data{i}-W*H{i},'fro')^2;
    end
    iterates.g(1,:) = g;
    iterates.primal_cost(1) = max(g);
    iterates.dual_cost(1) = lambda*g';
    
    
    W = iterates.W{1};
    H = iterates.H(1,:)';    
    % First round needs the unweighted data because SNPA was computed with
    % equal weights -TM 10/08/2020
    unweightedX= [data{:}];
    unweightedH = [H{:}];
    normX2 = norm(unweightedX,'fro')^2;
    WtW = W'*W;
    WtX = W'*unweightedX;
    err1 = max(0,normX2-2*sum(sum(WtX.*unweightedH))+sum(sum( WtW.*(unweightedH*unweightedH'))));
    
    
    err2 = log( det (WtW  + delta*eye(r) ) );
    beta = beta * max(1e-6,err1) / abs( err2 );
    e =  err1 + beta * err2 ;
    iterates.e(1) = e;
    iterates.err1(1) = err1;
    iterates.err2(1) = err2;
    
    % Main loop 
    for iter = 2 : maxiter        
        % ***** Take a subgradient step *****
        %alpha = a/sqrt(iter)
        alpha = a/iter;
        y = lambda + alpha*g;
        % Project onto the nonnegative orthant.  Code from [2]
        lambda = SimplexProj(y');
        lambda = lambda';
        support = lambda>0;
        iterates.lambda(iter,:) = lambda;
        
        W = iterates.W{iter-1};
        H = iterates.H(iter-1,:)';
        result = cellfun(@times,data(support),num2cell(sqrt(lambda(support)')),'uni',false); % weight the data and remove the zero elements
        weightedX = [result{:}];
        
        % ***** Update W, H_i's a few times *****
        % Most of the time this makes the approximate subgradient an
        % actual subgradient.  If each term is only updated once we can end
        % up moving in the wrong direction.  If there is only one patch in
        % the support of the minimum enclosing ball, we necessarily have
        % not found a subgradient so we do extra inner iterations
        if sum(support)<2
            twice = 1;
        else
            twice = 0;
        end
        for t = 1 : inneriter + (inneriter*twice)            
            % The H_i's need to be weighted before updating W each time
            % H_i's with dual variables (weights) equal to zero are removed
            result = cellfun(@times,H(support),num2cell(sqrt(lambda(support)')),'uni',false);
            weightedH = [result{:}];
            
            % *** Update W ***
            % W is updated using the weighted data and the weighted H_i's
            % Uses the fast gradient method formulated as a quadratic
            % program from [3]
            % Precompute the Hessian, A, for each row of W
            XHt = weightedX*weightedH';
            HHt = weightedH*weightedH';
            Y = inv( ( W'*W + delta*eye(r) ) );
            A = beta*Y + HHt;
            W = FGMqpnonneg(A,XHt,W,inneriter);

            % ***** Update the H_i's in the support *****
            % The H_i's are updated without the weights
            % Uses the fast gradient method code from [3]
            for i = 1 : nPatches
                if support(i)
                    [H{i},~,~,~] = FGMfcnls(data{i},W,H{i},inneriter);
                end
            end
            
        end

        % ****** Update all H_i's not in the support *****
        % Compute the error
        % Uses the fast gradient method code from [3]
        for i = 1 : nPatches
            if ~support(i)
                [H{i},~,~,~] = FGMfcnls(data{i},W,H{i},inneriter);
            end
            g(i) = norm(data{i}-W*H{i},'fro')^2;
        end  
       
        % Usually the initial step size is not a great estimate of the
        % minimum error, so we recompute after the first outer iteration
        if ~specifiedA && iter == 2
            options.a = 2/min(g);
            a = options.a; 
        end
        
      
        iterates.W{iter} = W; 
        iterates.H(iter,:) = H(:);
        iterates.g(iter,:) = g;
        iterates.primal_cost(iter) = max(g);
        iterates.dual_cost(iter) = lambda*g';
     
        W = iterates.W{iter};
        H = iterates.H(iter,:)';
        weightedX = [data{:}];      
        weightedH = [H{:}];
        normX2 = norm(weightedX,'fro')^2;
        WtW = W'*W;
        WtX = W'*weightedX;
        err1 = max(0,normX2-2*sum(sum(WtX.*weightedH))+sum(sum( WtW.*(weightedH*weightedH'))));
        err2 = log( det (WtW  + delta*eye(r) ) );
        e =  err1 + beta * err2 ;
        iterates.e(iter) = e;
        iterates.err1(iter) = err1;
        iterates.err2(iter) = err2;

        if options.display == 1
            fprintf('%1.0d ...', iter);
            if mod(iter,10) == 0
                fprintf('\n');
            end
        end
        if iterates.primal_cost(iter) - iterates.dual_cost(iter) < tol
            W = iterates.W{iter};
            H = iterates.H(iter,:)';
            iterates.lastIter = iter;        
           return 
        end
        W = iterates.W{iter};
        H = iterates.H(iter,:)';
        iterates.lastIter = iter;    
    end
end % of minimaxNMF   