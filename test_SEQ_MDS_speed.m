function test_SEQ_MDS_speed(d1, d2)
% This function tests the speed of classical MDS and SEQ-MDS for different
% data sets. The dimensionality of the data in the original space is given
% by d1. and we test the results and plot them for 2000-10000 data points.
% The results are mapped to MDS space with the dimension given by d2.


%number of tests
n=30;


time_CMDS   = zeros(1,n-2);
time_SEQMDS = zeros(1,n-2);


for i=3:n
    
    % doing repetition for averaging the time over some attemps
    for j=1:4
        
        
        % generate a set of samples
        num = i*100;
        X = rand(num,d1);

        
        %CMDS test
        tic;
        disMat = pdist(X);
        [Y1,e] = cmdscale(double(disMat));
        Y = Y1(:,1:d2);
        time_CMDS(i-2) = time_CMDS(i-2) + toc;

        
        %SEQ-MDS test
        tic;
        [Y,totaltime] = scmdscale_withDistAndEigs(X,d2,200,floor(3*d2/2),1);
        time_SEQMDS(i-2) = time_SEQMDS(i-2) + toc;
        
        
    end
    
    % averaging
    time_CMDS(i-2)   = time_CMDS(i-2)/4;
    time_SEQMDS(i-2) = time_SEQMDS(i-2)/4;
    
end

% plotting the results
plot(100*[3:n],time_CMDS,'k-');
hold on;
plot(100*[3:n],time_SEQMDS,'r--');

end