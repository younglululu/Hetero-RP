ddpath('cvx')
addpath('cvx/structures')
addpath('cvx/lib')
addpath('cvx/functions')
addpath('cvx/commands')
addpath('cvx/builtins')
%run('cvx/cvx_startup')

% XURL1 = '../data/rbp/X1.csv'; X1 = csvread(XURL1);
% XURL2 = '../data/rbp/X2.csv'; X2 = csvread(XURL2);
% XURL3 = '../data/rbp/X3.csv'; X3 = csvread(XURL3);
% XURL4 = '../data/rbp/X4.csv'; X4 = csvread(XURL4);
% XURL5 = '../data/rbp/X5.csv'; X5 = csvread(XURL5);
% 
% supAdjMatURL = '../data/rbp/matrix_truth.csv';
% tmpArr = load(supAdjMatURL); tmpArr = [tmpArr ones(size(tmpArr,1),1)];
% weightMat = spconvert(tmpArr); 
% if size(weightMat,1) ~= n || size(weightMat,2) ~= n, weightMat(n,n) = 0; end
% A_attr = weightMat + weightMat'; clear tmpArr weightMat;
% A_repul = ones(n)-A_attr;
% 
% save('../data/rbp.mat', 'X1','X2','X3','X4','X5','A_attr','A_repul');

load('../data/rbp.mat');
X = [X1; X2; X3; X4; X5]; [m,n]=size(X); size(X)

[X_reweigh, weightVec] = heteroRP(X, A_attr, A_repul);