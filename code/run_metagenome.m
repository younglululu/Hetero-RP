addpath('cvx')
addpath('cvx/structures')
addpath('cvx/lib')
addpath('cvx/functions')
addpath('cvx/commands')
addpath('cvx/builtins')
%run('cvx/cvx_startup')

XURL1 = '../data/metagenome/X1.csv';
XURL2 = '../data/metagenome/X2.csv';
alignAdjMatURL = '../data/metagenome/alignMat.csv';
originAdjMatURL = '../data/metagenome/originMat.csv';
dataURL = '../data/metagenome.mat';

% X1 = csvread(XURL1); X2 = csvread(XURL2); X = [X1; X2]; [m,n]=size(X);
% A_origin = buildKnnGraph(X, originAdjMatURL, top);
% 
% tmp = load(originAdjMatURL);
% weightMat = spconvert([tmp(:,1:2) ones(size(tmp,1),1)]);
% if size(weightMat,1) ~= n || size(weightMat,2) ~= n, weightMat(n,n) = 0; end
% A_sup_bin = (weightMat > 0); A_sup_bin = A_sup_bin + A_sup_bin';
% A_align = weightMat + weightMat';  
% A_align(A_align>0) = A_align(A_align>0) ./ A_sup_bin(A_sup_bin>0); clear weightMat A_sup_bin tmp;
% save(dataURL, 'X1', 'X2', 'A_origin', 'A_align');

load(dataURL);
X = [X1; X2]; [m,n]=size(X); 

gamma = trace(A_origin'*A_align)/trace(A_origin'*A_origin);
A_attr = A_align + gamma*A_origin;

[X_reweigh, weightVec] = heteroRP(X, A_attr, []);





