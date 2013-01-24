function [evecs, evals] = symeig(A, k)
%SLSYMEIG Compute the eigenvalues and eigenvectors for symmetric matrix
%
% $ Syntax $
%   - evals = slsymeig(A)
%   - [evals, evecs] = slsymeig(A)
%   - ... = slsymeig(A, k)
%   - ... = slsymeig(A, k, ord)
%
% $ Arguments $
%   - A:        the target symmetrix matrix to be decomposed
%   - k:        the number of eigenvalues to be solved
%   - ord:      the order of the sorting ('ascend' | 'descend')
%               default = 'descend';
%   - evals:    the column vector of eigenvalues of A (in descending order)
%   - evecs:    the matrix of eigenvectors of A
%
% $ Description $
%   - evals = slsymeig(A) computes only the eigenvalues of symmetric 
%     matrix A.
%
%   - [evals, evecs] = slsymeig(A) computes both the eigenvalues and
%     eigenvectors of matrix A. If A is a d x d matrix, then evals
%     is a d x 1 column vector with the eigenvalues in descending order.
%     evecs is a d x d matrix with each column vector corresponding to
%     an eigenvalue in evals.
%
%   - ... = slsymeig(A, k) the user can specify the number of eigenvalues
%     to preserve, if not specified, all eigenvalues will be preserved.
%
%   - ... = slsymeig(A, k, ord) the user can specify how to sort the
%     eigenvalues. By default, they are sorted in descending order.
%
% $ History $
%   - Created by Dahua Lin on Apr 21, 2006
%   - Modified by Dahua Lin on Sep 9, 2006
%       - The user now can specify the number of eigenvalues to preserve
%         if necessary.  
%       - An auto-selection scheme, which will use eigs when k is small
%         in order to achieve higher efficiency.
%         
%

%% parse and verify input arguments

[d, d2] = size(A);


%% compute

% Matrix should be symetric

% decide scheme and compute
if k > d / 3 && ~issparse(A)
    [evecs, evals] = eig(A);
else            
    opts = struct('disp', 0, 'issym', true);
    [evecs, evals] = eigs(A, k, 'LM', opts);
end

evals = diag(evals);

% enforce real
if ~isreal(evecs)
    evecs = real(evecs);
end

if ~isreal(evals)
    evals = real(evals);
end

%% sort: make the eigenvalues in descending order

[evals, si] = sort(evals, 1, 'descend');
evecs = evecs(:, si);

k0 = length(evals);
if k < k0
    evals = evals(1:k);
    evecs = evecs(:, 1:k);
end



