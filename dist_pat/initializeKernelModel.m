function [clusterModel, Z] = initializeKernelModel(Y, K, sigma, par)

Pat = par.Pat;
Patz= par.Patz;
MDS = par.MDS;

fprintf('Initialize Kernel PCA    ......................');

% Initialization of the model
clusterModel = kpca(Y', K, struct('ker','rbf','arg',sigma,'new_dim',MDS));
% clusterModel = kpca(Y', K, struct('ker','rbf','arg',sigma,'new_dim',Pat^2*Patz));
% clusterModel = greedykpca(Y',struct('ker','rbf','arg',sigma,'new_dim',Pat^2*Patz)); 
Z = kernelproj(Y', clusterModel)';

fprintf(' Done!\n');
end

