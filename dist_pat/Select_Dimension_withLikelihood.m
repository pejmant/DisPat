function dim = Select_Dimension_withLikelihood(evals)

% This function uses the concept of scree plot to find the elbow or a big
% gap in the plot. This is done by maximizing a simple profile likelighood
% function. 
%   evals : descending eignenvalues of the covariance matrix.
%
%                                           Written by: Mehrdad Honarkhah

global l;

% number of eigenvalues
p = length(evals);

if p < 3
    dim = p;
    return;
end

% vector to hold the profile log-Likelihood of different dimensions
l = zeros(p-2,1);

for i = 2:p-1
    
    % assign two sets of samples
    Set1 = evals(1:i);
    Set2 = evals(i+1:end);
    
    % Calculate their means
    mu1  = sum(Set1)/i;
    mu2  = sum(Set2)/(p-i);
    
    % Calculate sample variance of sets
    var1 = var(Set1);
    var2 = var(Set2);
    
    % Calculate variance used for both samples
    varT = ((i-1)*var1 + (p-i-1)*var2)/(p-2);
    
    % find the log-Likelihood function for this i
    lSet1= log(1/sqrt(2*pi*var1)) - (Set1 - mu1)./(2*var1);
    lSet2= log(1/sqrt(2*pi*var2)) - (Set2 - mu2)./(2*var2);
    
    l(i-1) = sum(lSet1) + sum(lSet2);
end

[dummy, index] = max(l);

dim = index+1;
end