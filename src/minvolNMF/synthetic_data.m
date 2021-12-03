% Synthetic data set from 
% Minimum-Volume Rank-Deficient Nonnegative Matrix Factorizations, 
% Valentin Leplat, Andersen M.S. Ang, Nicolas Gillis, 2018. 
clear all; clc; close all; 
% True generating W, with rank_+(W) = 4 > rank(W) = 3
disp('True W:') 
Wt = [1 0 0 1; 1 0 1 0; 0 1 1 0; 0 1 0 1]' 
r = 4; % nonnegative rank 
% Generate H with n columns
n = 500; 
purity = 0.8; % ~min distance between the columns of W and the data points
alpha = 0.05*ones(r,1);  % Parameter of the Dirichlet distribution
Ht = [sample_dirichlet(alpha,n)']; 
for j = 1 : n
    while max( Ht(:,j) ) > purity
        Ht(:,j) = sample_dirichlet(alpha,1)';
    end
end
% Generate X = Wt*Ht + noise
epsilon = 0.01; % standard deviation 
Xt = Wt*Ht; 
X = max( 0 , Xt + epsilon*randn(size(Xt)) ); 
% Run min-vol NMF (default parameters: delta=0.1, maxiter=100)
options.lambda=0.01;
[W,H,e,er1,er2,lambda] = minvolNMF(X,r,options); 
% Display results 
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
figure; 
plot3(X(1,:), X(2,:), X(3,:),'bo'); hold on; 
plot3(Wt(1,:), Wt(2,:), Wt(3,:), 'ro', 'MarkerSize', 12); 
plot3(W(1,:), W(2,:), W(3,:),'mx', 'MarkerSize', 12); 
grid on; 
legend('Data points', 'True W', 'Estimated W'); 
figure; 
plot(er1); hold on; 
plot(er2); 
legend('||X-WH||_F^2','logdet(W^TW+ \delta I)'); 
xlabel('Iterations'); 
disp('Computed W:')
W 
fprintf('Error ||W-Wt||/||Wt|| = %2.2f%%.\n', 100*compareWs( Wt, W ) ); 