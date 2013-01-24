function model = kpca(X,K,options)
% KPCA Kernel Principal Component Analysis.
%  
% Synopsis:
%  model = kpca(X,K)
%  model = kpca(X,K,options)
%
% Description:
%  This function is implementation of Kernel Principal Component 
%  Analysis (KPCA) [Schol98b]. The input data X are non-linearly
%  mapped to a new high dimensional space induced by prescribed
%  kernel function. The PCA is applied on the non-linearly mapped 
%  data. The result is a model describing non-linear data projection.
%  See 'help kernelproj' for info how to project data.
%
% Input:
%  X [dim x num_data] Training data.
%  
%  options [struct] Decribes kernel and output dimension:
%   .ker [string] Kernel identifier (see 'help kernel'); 
%     (default 'linear').
%   .arg [1 x narg] kernel argument; (default 1).
%   .new_dim [1x1] Output dimension (number of used principal 
%     components); (default dim).
%
% Output:
%  model [struct] Kernel projection:
%   .Alpha [num_data x new_dim] Multipliers.
%   .b [new_dim x 1] Bias.
%   .sv.X [dim x num_data] Training vectors.
%  
%   .nsv [1x1] Number of training data.
%   .eigval [1 x num_data] Eigenvalues of centered kernel matrix.
%   .mse [1x1] Mean square representation error of maped data.
%   .MsErr [dim x 1] MSE with respect to used basis vectors;
%      mse=MsErr(new_dim).
%   .kercnt [1x1] Number of used kernel evaluations.
%   .options [struct] Copy of used options.
%   .cputime [1x1] CPU time used for training.
%
% Example:
%  X = gencircledata([1;1],5,250,1);
%  model = kpca( X, struct('ker','rbf','arg',4,'new_dim',2));
%  XR = kpcarec( X, model );
%  figure; 
%  ppatterns( X ); ppatterns( XR, '+r' );
%  
% See also 
%  KERNELPROJ, PCA, GDA.
%   

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-may-2004, VF
% 10-july-2003, VF, computation of kercnt added
% 22-jan-2003, VF
% 11-july-2002, VF, mistake "Jt=zeros(N,L)/N" repared 
%              (reported by SH_Srinivasan@Satyam.com).
% 5-July-2001, V.Franc, comments changed
% 20-dec-2000, V.Franc, algorithm was implemented

% timer
start_time = cputime;

% gets dimensions
[dim,num_data] = size(X);  

% process input arguments
%-----------------------------------
if nargin < 3, options = []; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'new_dim'), options.new_dim = dim; end


% Centering kernel matrix (non-linearly mapped data).
J = ones(num_data,num_data)/num_data;
JK = J*K;
JKJ = bsxfun(@minus,JK,K)*J;
Kc = K - JK + JKJ;

% eigen decomposition of the kernel marix
% [U,D] = eig(Kc);
% Faster way is given here, in case where f
[U,D] = symeig(Kc, options.new_dim);
% Lambda=real(diag(D));

Lambda = D;

% normalization of eigenvectors to be orthonormal 
Lambda(Lambda==0) = eps;
% U = U ./ sqrtLambda(ones(size(U,1), 1),:);
U = bsxfun(@rdivide,U,sqrt(Lambda)');

% Sort the eigenvalues and the eigenvectors in descending order.
% [Lambda,ordered]=sort(-Lambda);    
% Lambda=-Lambda;
% U=U(:,ordered);                    

% use first new_dim principal components
% A=U(:,1:options.new_dim);              

% compute Alpha and compute bias (implicite centering)
% of kernel projection
% Mehrdad have changed Matrix A to U for his faster code
model.Alpha = (eye(num_data,num_data)-J)*U;
Jt=ones(num_data,1)/num_data;
KJt=K*Jt;
model.b = U'*(J'*KJt-KJt);

% fill output structure
model.sv.X = X;
model.nsv = num_data;
model.options = options;
model.eigval = D;   % change Lambda to D by Mehrdad
model.kercnt = num_data*(num_data+1)/2;
% model.MsErr = triu(ones(num_data,num_data),1)*model.eigval/num_data;
% model.mse = model.MsErr(options.new_dim);
model.cputime = cputime - start_time;
model.fun = 'kernelproj';

return;
