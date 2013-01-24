% Some codes



%% making the training image to be blurred (using Log(distance transform)

% out = training image

outbw = bwdist(out);                % Do distance Transform
outbw = log(outbw+1);               % +1 is to make sure we have positive values
outbw = 1-outbw./max(max(outbw));   % inversely rescale it to [0,1]
drawPattern(outbw);                 % show resulting TI


% Now for different multi-grids we find a correspodning TI by resizing
outbwScaled_3 = imresize(outbw, 1/3);
drawPattern(outbwScaled_3);           % show the TI for multi-grid = 3
outbwScaled_2 = imresize(outbw, 1/2);
drawPattern(outbwScaled_2);           % show the TI for multi-grid = 2


% Show the reconstruction of the generated multi-grid over full scale
outbw_3 = imresize(outbwScaled_3, 3);
drawPattern(outbw_3);                 % show the reconstructed original TI
outbw_2 = imresize(outbwScaled_2, 2);
drawPattern(outbw_2);                 % show the reconstructed original TI




%% How to automatically select template size using entropy

% maximum size to be tested (only odd values)
MAX_SIZE = 43;

ent = zeros(1,(MAX_SIZE-1)/2);      % only odd values will be tested
var1 = zeros(1,(MAX_SIZE-1)/2);      % only odd values will be tested
for i = 3:2:MAX_SIZE
%     f = @(x) entropy(x);
%     I2 = nlfilter(out,[i i],f);
    I2 = entropyfilt(out,ones(i));  % faster implementation than two previous lines
    
    figure; imagesc(I2);
    ent((i-1)/2) = mean(mean(I2));
    var1((i-1)/2) = var(I2(:));
end
global l;

% ent = diff(ent,2);

Tsize = 2*Select_Dimension_withLikelihood(ent) + 1; % considering odd ordering
fprintf('Template Size = %d\n',Tsize);

%visualize the results 
figure;
subplot(1,2,1);plot(3:2:length(ent)*2+1,ent,':ro','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);
title('Entropy Plot');xlabel('Template size');ylabel('Mean Entropy');axis square;set(gca,'XTick',1:2:MAX_SIZE);
subplot(1,2,2);plot(5:2:length(ent)*2+1-2,l,'-ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);
title('Profile Log-likelihood');xlabel('Template size');ylabel('log-likelihood'); axis square;set(gca,'XTick',1:2:MAX_SIZE);
hold on; line([Tsize,Tsize],get(gca,'Ylim'),'LineStyle',':','Color','k')
plot(Tsize,l((Tsize-1)/2-1),'ks','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
figure;
plot(3:2:length(ent)*2+1,var1,':ro','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);
title('Variance of Entropy Plot');xlabel('Template size');ylabel('Variance Entropy');axis square;set(gca,'XTick',1:2:MAX_SIZE);




%% Kmeans speed analysis

% fastet
tic;[centroid, labels] = slkmeans(Y','K',100,'clsfunc','ann');toc

% matlab
opts = statset('MaxIter',1000);
tic;kmeans(Y, par.clus,'Options',opts);toc

% kernel kmeans
tic;K = slkernel(Y', 'gauss', sigma);idx = dualkmeansFast(K,par.clus);toc

% fast kmeans by triangular
tic;[centers,mincenter,mindist,q2,quality] =Fkmeans(Y, par.clus);toc

% mpi-kmeans (triangular coded in c++) needs to be mexed
 tic;[CX, sse, assignment]=mpi_kmeans(Y', par.clus);toc

 
 
 %% PCA analysis
 
% statistical learning toolbox
% tic; S   = slpca(X','preserve',{'number',3}); toc
tic; S   = slpca(X','preserve',{'Everrit',0}); toc
tic;Ynew = slapplypca(S,Xnew');toc
tic;Xe = slpcarecon(S, Ynew);toc
subplot(1,2,2);drawPattern(reshape(Xe(:,1)'>0.5,9,9));
subplot(1,2,1);drawPattern(reshape(Xnew,9,9));


% dimensionality reduction toolbox
tic;[Y mapping]= compute_mapping(X,'PCA',12);toc
tic;Ynew = out_of_sample(Xnew,mapping);toc



%% out of sample with kmedoid

% get represantative samples
c = min(size(X,1), 5*par.MDS);
repSample = kmedoid(X, c);

% calculate the distance matrix of kmedoids
repSample.dist = zeros(5*par.MDS + 1);
repSample.dist(1:end-1,1:end-1) = squareform(pdist(repSample.v));

% get the out-of-sample point for Xsample
Ysample   = addOutOfSample(Xsample, repSample, Y, par);

% or to test the patterns in the database themselves, do this code:
nbp = 140;      % index of pattern to be tested
Ysample   = addOutOfSample(X(nbp,:), repSample, Y, par);
parallelcoords(Y,'color',[0.8,0.8,0.8]); hold on
parallelcoords(Y(nbp,:),'color','black','Linewidth',2)
parallelcoords(Ysample,'color','red','Linewidth',2)



%% Kmeans ++    
% to find initial centroids for kmeans (faster)

centroid = zeros(par.clus, par.MDS);

centroid(1,:) = Y(ceil(size(Y,1)*rand(1,1)),:);
for i = 2:par.clus
   dist = slmetric_pw(Y', centroid(1:i-1,:)', 'sqdist');
   dummy = min(dist,[],2);
   [dummy, inx]  = max(dummy);
   centroid(i,:) = Y(inx,:);
end
tic;[centroid, labels] = slkmeans(Y','K',100,'init_means',centroid','clsfunc','ann');toc



%% another k means method (hierarchical)
tic;
centroid = [];
nn= size(Y,1);
rr = [1,2];
for i = rr
    if i ~= 2
        Y1 = Y(ceil(nn.*rand(ceil(nn*i/2),1)),:);
    else
        Y1 = Y;
    end
[centroid, labels] = slkmeans(Y1','K',clus,'init_means',centroid','clsfunc','ann','verbose',true);
centroid = centroid';
end
toc

%% Hierarchical K-means test on the proportion of data points for first sampling

Y = rand(10000,15);
nn= size(Y,1);
clus = 80;

% number of prportion test (2:MAX)
MAX = 20;
time  = zeros(1,MAX-1); %for hierarchical
time1 = zeros(1,MAX-1); %without hierarchical

for prop=2:MAX
    for test_count = 1:10
        tic;
        [centroid, labels] = slkmeans(Y(ceil(nn.*rand(ceil(nn/prop),1)),:)','K',clus,'init_means',[],'clsfunc','ann','verbose',false);
        [centroid, labels] = slkmeans(Y','K',clus,'init_means',centroid,'clsfunc','ann','verbose',false);
        time(prop-1) = time(prop-1) + toc/30;
        tic;
        [centroid, labels] = slkmeans(Y','K',clus,'clsfunc','ann','verbose',false);
        time1(prop-1) = time1(prop-1) + toc/30;
    end
end

figure;
plot(2:MAX,time ,'-ko','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);
hold on;
plot(2:MAX,time1,':ko','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3);
    



% another test on changing proportation a lot like 3,2,1.5,1
figure;hold on;
for clus=20:10:100
    time2=0; time3=0;
    for i=1:30
        tic;
        prop=10;[centroid, labels] = slkmeans(Y(ceil(nn.*rand(ceil(nn/prop),1)),:)','K',clus,'init_means',[],'clsfunc','ann','verbose',false);
        prop=1;[centroid, labels] = slkmeans(Y','K',clus,'init_means',centroid,'clsfunc','ann','maxiter',300,'verbose',false);
        time2 = time2 + toc/30;
        tic;prop=1;[centroid, labels] = slkmeans(Y','K',clus,'init_means',[],'clsfunc','ann','maxiter',300,'verbose',false);
        time3 = time3 + toc/30;
    end
    fprintf('\ntime for hierarchical = %g\n',time2);
    fprintf('time for normal = %g\n',time3);
    plot(clus,time2 ,'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3);
    plot(clus,time3 ,'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);
end



%% procedure to test the speed of matlab kmeans and Ann kmeans

clus = 100;
nn = size(Y,1);

opts = statset('MaxIter',1000);
time_k = []; time_a = [];
for i = 500:200:2209
    i
    tt1 = 0;
    tt2 = 0;
    for j=1:30
        Yttt = Y(1:i,:);
        tic;kmeans(Yttt, clus,'Options',opts);tt1 = tt1 + toc/30;
        tic;[centroid, labels] = slkmeans(Yttt','K',clus,'clsfunc','ann','verbose',false,'maxiter',500);tt2 = tt2 + toc/30;
    end
    time_k = horzcat(time_k, tt1);
    time_a = horzcat(time_a, tt2);
end

plot( 500:200:2209,time_k ,':ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3);
hold on;
plot( 500:200:2209,time_a ,':ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);




%% comparison of speed of dualkmeansFast with (dualkmeansFast + Ann-k-means)

clus = 80;
time_k = []; time_a = [];

for i = 500:100:2209
    i
    tt1 = 0;
    tt2 = 0;
    for j=1:40
        rr= randperm(i);
        II = K(rr,rr);
        tic;idx = dualkmeansFast(II,clus);tt1 = tt1 + toc/40;
        
        tic;[centroid, labels] = slkmeans(Y(rr,:)','K',clus,'clsfunc','ann','verbose',false,'maxiter',300);
        idx = dualkmeansFast(II,clus,labels');tt2 = tt2 + toc/40;
    end
    time_k = horzcat(time_k, tt1);
    time_a = horzcat(time_a, tt2);
end

plot( 500:100:2209,time_k ,':ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3);
hold on;
plot( 500:100:2209,time_a ,':ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',4);






%% procedure to find the number of clusters in kernel kmeans

% the maximum number that we want to analyze upto
MAX = 200;

% number of data points is NUM
NUM = size(K,1);

% fast eigenvalue decomposition 
tic;[EigV,EigD] = eigs(K,MAX);

int = zeros(1,MAX);
for i=1:MAX
    int(i)=EigD(i,i)*(1/NUM*ones(1,NUM)*EigV(:,i))^2;
end

toc
% showing the result
figure;plot(1:MAX,int,'-ko','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3);
title('Eigenvalue of centered Kernel');

% smoothing the result steadily
figure1 = figure;
axes('Parent',figure1,'YScale','log','YMinorTick','on','XScale','log','XMinorTick','on');
title('Eigenvalue of centered Kernel - SMOOTHED');
for i = 0.001:0.05:0.2
    yy2 = smooth(1:MAX,int,i,'rloess');
    loglog(1:MAX,yy2,'-ko','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3);
    title(['Eigenvalue of centered Kernel - SMOOTHED with i = ',num2str(i)]);
    ylim([min(int) max(int)]);
    drawnow expose;
    pause(0.3);
end

% hold on;hh = fill([1,1:80,80],[3*10^(-8),yy2(1:80)',3*10^(-8)],[0.8 0.8 0.8]);



%% visualizations

% whole MDS points
plot3(Y(:,1),Y(:,2),Y(:,3),'.');grid on;hold on; 

% single MDS point
c=25;
plot3(Y(c,1),Y(c,2),Y(c,3),'ro');grid on;hold on;

% some reconstructed point
point = Ynew';
plot3(point(1,1),point(1,2),point(1,3),'ko');grid on;hold on;
plot3(point(1,1),point(1,2),point(1,3),'k*');grid on;hold on;


%% greedy kernel PCA on patterns and reconstructing one sample pattern

tic;model = greedykpca(X',struct('ker','rbf','arg',sigma,'new_dim',25,'m',100));toc


% pattern index for testing
testP = 150;

% one simple test on reconstruction of original pattern
% ---------------------------------------------------------
XR = kpcarec(X(testP,:)',model);
% visualizing;
figure;
subplot(1,2,1); drawPattern(X(testP,:)); title('Original');
subplot(1,2,2); drawPattern(XR'); title('Reconstructed');


% another simple test on reconstruction of partial(=0.5) original pattern
% ----------------------------------------------------------------------
testP= 199;
Xtest = X(testP,:);
Xtest = reshape(Xtest,par.Pat,par.Pat);
Xtest(1:9,3:7)=0.5;
Xtest = reshape(Xtest,1,81);
XR = kpcarec(Xtest',model);
% visualization

drawPattern(X(testP,:)); title('Original');
drawPattern(Xtest); title('partial');
drawPattern(XR'); title('Reconstructed');




%% Plotting the realization beautifully

if sum(realization(:)) ~= sum(round(realization(:)))
    figure1   = figure('Color',[0.8 0.8 0.8]);  %categorical
else
    figure1 = figure('Colormap',[0.23 0.23 0.23;0.2422 0.2422 0.2422;0.2544 0.2544 0.2544;0.2667 0.2667 0.2667;0.2789 0.2789 0.2789;0.2911 0.2911 0.2911;0.3033 0.3033 0.3033;0.3156 0.3156 0.3156;0.3278 0.3278 0.3278;0.34 0.34 0.34;0.3522 0.3522 0.3522;0.3644 0.3644 0.3644;0.3767 0.3767 0.3767;0.3889 0.3889 0.3889;0.4011 0.4011 0.4011;0.4133 0.4133 0.4133;0.4256 0.4256 0.4256;0.4378 0.4378 0.4378;0.45 0.45 0.45;0.4622 0.4622 0.4622;0.4744 0.4744 0.4744;0.4867 0.4867 0.4867;0.4989 0.4989 0.4989;0.5111 0.5111 0.5111;0.5233 0.5233 0.5233;0.5356 0.5356 0.5356;0.5478 0.5478 0.5478;0.56 0.56 0.56;0.5722 0.5722 0.5722;0.5844 0.5844 0.5844;0.5967 0.5967 0.5967;0.6089 0.6089 0.6089;0.6211 0.6211 0.6211;0.6333 0.6333 0.6333;0.6456 0.6456 0.6456;0.6578 0.6578 0.6578;0.67 0.67 0.67;0.6822 0.6822 0.6822;0.6944 0.6944 0.6944;0.7067 0.7067 0.7067;0.7189 0.7189 0.7189;0.7311 0.7311 0.7311;0.7433 0.7433 0.7433;0.7556 0.7556 0.7556;0.7678 0.7678 0.7678;0.78 0.78 0.78;0.7922 0.7922 0.7922;0.8044 0.8044 0.8044;0.8167 0.8167 0.8167;0.8289 0.8289 0.8289;0.8411 0.8411 0.8411;0.8533 0.8533 0.8533;0.8656 0.8656 0.8656;0.8778 0.8778 0.8778;0.89 0.89 0.89;0.9022 0.9022 0.9022;0.9144 0.9144 0.9144;0.9267 0.9267 0.9267;0.9389 0.9389 0.9389;0.9511 0.9511 0.9511;0.9633 0.9633 0.9633;0.9756 0.9756 0.9756;0.9878 0.9878 0.9878;1 1 1], 'Color',[0.8 0.8 0.8]);
end
axes1 = axes('Visible','off','Parent',figure1,'DataAspectRatio',[1 1 1],'CLim',[0 1]);
% annotation(figure1,'rectangle',[0.2135 0.1119 0.5151 0.6857],'LineWidth',2,'FaceColor','flat');
box('on');
hold('all');
imagesc(realization);
axis square off;
axis ij;
% saveas(gcf,'old_mg_3.fig')



%% plotting the realization from SGEMS beautifully

Dim = 101;      % Set this value manually, also the location of the file
Dimz = 1;

file_location = 'C:\Documents and Settings\mehrdadh\Desktop\eee.dat';
[datain, colnames, line1]=loadgeoeas(file_location);
out = geoeas2matlab(datain,[Dim Dim Dimz]);
if max(max(out)) > 1
    figure1   = figure('Color',[0.8 0.8 0.8]);
    out = out./max(max(out));
else
    figure1 = figure('Colormap',[0.3137 0.3176 0.3137;0.314 0.3177 0.3194;0.3142 0.3178 0.3251;0.3145 0.3179 0.3308;0.3147 0.318 0.3365;0.315 0.318 0.3421;0.3152 0.3181 0.3478;0.3155 0.3182 0.3535;0.3157 0.3183 0.3592;0.316 0.3184 0.3649;0.3162 0.3184 0.3706;0.3165 0.3185 0.3762;0.3167 0.3186 0.3819;0.317 0.3187 0.3876;0.3172 0.3187 0.3933;0.3175 0.3188 0.399;0.3177 0.3189 0.4047;0.318 0.319 0.4103;0.3182 0.3191 0.416;0.3184 0.3191 0.4217;0.3187 0.3192 0.4274;0.3189 0.3193 0.4331;0.3192 0.3194 0.4388;0.3194 0.3194 0.4444;0.3333 0.3385 0.4583;0.3472 0.3576 0.4722;0.3611 0.3767 0.4861;0.375 0.3958 0.5;0.3889 0.4149 0.5139;0.4028 0.434 0.5278;0.4167 0.4531 0.5417;0.4306 0.4722 0.5556;0.4444 0.4913 0.5694;0.4583 0.5104 0.5833;0.4722 0.5295 0.5972;0.4861 0.5486 0.6111;0.5 0.5677 0.625;0.5139 0.5868 0.6389;0.5278 0.6059 0.6528;0.5417 0.625 0.6667;0.5556 0.6441 0.6806;0.5694 0.6632 0.6944;0.5833 0.6823 0.7083;0.5972 0.7014 0.7222;0.6111 0.7205 0.7361;0.625 0.7396 0.75;0.6389 0.7587 0.7639;0.6528 0.7778 0.7778;0.6745 0.7917 0.7917;0.6962 0.8056 0.8056;0.7179 0.8194 0.8194;0.7396 0.8333 0.8333;0.7613 0.8472 0.8472;0.783 0.8611 0.8611;0.8047 0.875 0.875;0.8264 0.8889 0.8889;0.8481 0.9028 0.9028;0.8698 0.9167 0.9167;0.8915 0.9306 0.9306;0.9132 0.9444 0.9444;0.9349 0.9583 0.9583;0.9566 0.9722 0.9722;0.9783 0.9861 0.9861;1 1 1],'Color',[0.8 0.8 0.8]);
end
axes1 = axes('Visible','off','Parent',figure1,'DataAspectRatio',[1 1 1],'CLim',[0 1]);
% annotation(figure1,'rectangle',[0.2135 0.1119 0.5151 0.6857],'LineWidth',2,'FaceColor','flat');
box('on');
hold('all');
imagesc(out);
set(gcf,'Colormap',[0 0 0.005208;0.01998 0.01998 0.01779;0.03995 0.03995 0.03037;0.05993 0.05993 0.04296;0.0799 0.0799 0.05554;0.09988 0.09988 0.06812;0.1199 0.1199 0.0807;0.1398 0.1398 0.09328;0.1598 0.1598 0.1059;0.1798 0.1798 0.1184;0.1998 0.1998 0.131;0.2197 0.2197 0.1436;0.2397 0.2397 0.1562;0.2597 0.2597 0.1688;0.2797 0.2797 0.1814;0.2996 0.2996 0.1939;0.3196 0.3196 0.2065;0.3396 0.3396 0.2191;0.3596 0.3596 0.2317;0.3795 0.3795 0.2443;0.3995 0.3995 0.2569;0.4195 0.4195 0.2694;0.4395 0.4395 0.282;0.4594 0.4594 0.2946;0.4794 0.4794 0.3072;0.4994 0.4994 0.3198;0.5194 0.5194 0.3323;0.5393 0.5393 0.3449;0.5593 0.5593 0.3575;0.5793 0.5793 0.3701;0.5993 0.5993 0.3827;0.6192 0.6192 0.3953;0.6392 0.6392 0.4078;0.6512 0.6512 0.4142;0.6633 0.6633 0.4207;0.6753 0.6753 0.4271;0.6873 0.6873 0.4335;0.6993 0.6993 0.4399;0.7114 0.7114 0.4463;0.7234 0.7234 0.4527;0.7354 0.7354 0.4591;0.7475 0.7475 0.4655;0.7595 0.7595 0.4719;0.7715 0.7715 0.4783;0.7835 0.7835 0.4847;0.7956 0.7956 0.4911;0.8076 0.8076 0.4975;0.8196 0.8196 0.5039;0.8316 0.8316 0.5103;0.8437 0.8437 0.5167;0.8557 0.8557 0.5231;0.8677 0.8677 0.5295;0.8797 0.8797 0.5359;0.8918 0.8918 0.5424;0.9038 0.9038 0.5488;0.9158 0.9158 0.5552;0.9278 0.9278 0.5616;0.9399 0.9399 0.568;0.9519 0.9519 0.5744;0.9639 0.9639 0.5808;0.9759 0.9759 0.5872;0.988 0.988 0.5936;1 1 0.6;1 1 0.6]);
axis square;
axis off;
% saveas(gcf,'crack_sgems_15_11_3.fig')


%% Matlab realization/TI    ------>     SGems *.dat file


fid = fopen([dirName,'\Training Images\meander.dat'],'w+');
fileData = matlab2geoeas(out);
[Iout,Jout,Kout] = size(out);
fprintf(fid,'Grid (%dx%dx%d)\n1\nMatlab\n',Iout,Jout,Kout);
fprintf(fid,'%g\n',fileData);
fclose(fid);
pause(1);fclose all;
fprintf('Done.\n');



%% make some 3D slices beautifully


v=out1;
clrlmt=[]; varargin=[];
[ny,nx,nz] = size(v); x=1:nx; y=1:ny; z=1:nz;

hax3d=figure('Colormap',[0.3137 0.3176 0.3137;0.314 0.3177 0.3194;0.3142 0.3178 0.3251;0.3145 0.3179 0.3308;0.3147 0.318 0.3365;0.315 0.318 0.3421;0.3152 0.3181 0.3478;0.3155 0.3182 0.3535;0.3157 0.3183 0.3592;0.316 0.3184 0.3649;0.3162 0.3184 0.3706;0.3165 0.3185 0.3762;0.3167 0.3186 0.3819;0.317 0.3187 0.3876;0.3172 0.3187 0.3933;0.3175 0.3188 0.399;0.3177 0.3189 0.4047;0.318 0.319 0.4103;0.3182 0.3191 0.416;0.3184 0.3191 0.4217;0.3187 0.3192 0.4274;0.3189 0.3193 0.4331;0.3192 0.3194 0.4388;0.3194 0.3194 0.4444;0.3333 0.3385 0.4583;0.3472 0.3576 0.4722;0.3611 0.3767 0.4861;0.375 0.3958 0.5;0.3889 0.4149 0.5139;0.4028 0.434 0.5278;0.4167 0.4531 0.5417;0.4306 0.4722 0.5556;0.4444 0.4913 0.5694;0.4583 0.5104 0.5833;0.4722 0.5295 0.5972;0.4861 0.5486 0.6111;0.5 0.5677 0.625;0.5139 0.5868 0.6389;0.5278 0.6059 0.6528;0.5417 0.625 0.6667;0.5556 0.6441 0.6806;0.5694 0.6632 0.6944;0.5833 0.6823 0.7083;0.5972 0.7014 0.7222;0.6111 0.7205 0.7361;0.625 0.7396 0.75;0.6389 0.7587 0.7639;0.6528 0.7778 0.7778;0.6745 0.7917 0.7917;0.6962 0.8056 0.8056;0.7179 0.8194 0.8194;0.7396 0.8333 0.8333;0.7613 0.8472 0.8472;0.783 0.8611 0.8611;0.8047 0.875 0.875;0.8264 0.8889 0.8889;0.8481 0.9028 0.9028;0.8698 0.9167 0.9167;0.8915 0.9306 0.9306;0.9132 0.9444 0.9444;0.9349 0.9583 0.9583;0.9566 0.9722 0.9722;0.9783 0.9861 0.9861;1 1 0.6],'Color',[0.8 0.8 0.8]);
slice(x,y,z,v,[],[],[7,20,39]);
shading flat; 
rotate3d on; box on; hold all;
if ~isempty(clrlmt), set(hax3d,'clim',clrlmt); end;
if isempty(clrlmt), clrlmt=get(hax3d,'clim'); end;
set(gca,'zdir','reverse','Projection','perspective','LineWidth',1,...
    'color','black','Xcolor','white','Ycolor','white','Zcolor','white'); axis tight;

set(hax3d,'tag','ax3d','userdata',{x,y,z,v});
view([-34.5 18]);
drawnow;
c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off','zcolor', 'k', 'zgrid', 'off');
xlabel('X'); ylabel('Y'); zlabel('Z'); 




%% To test 100 realization and calculating the MPH of each one of them (mine)

% number of realizations
NUM = 300;

% mph = zeros(NUM,2^16);
Real100 = zeros(101,101,NUM);
for re = 1:NUM
    fprintf('%d\n',re);
    distPat;
    Real100(:,:,re) = realization;
%     mph(re,:) = calculateMPH(Real100(:,:,re));
end

mphTI = calculateMPH(out);
error_mine = sum(abs(bsxfun(@minus,mph,mphTI)),2);

[mdsMPH, eMPH] = cmdscale(pdist(mph));
plot3(mdsMPH(:,1),mdsMPH(:,2),mdsMPH(:,3),'k.');
grid on;

% hausdorf distance
dismatrix = zeros(100,100);
for i = 1:100
    for j = i:100
        dismatrix(i,j) = haussdorff(Real100(:,:,i),Real100(:,:,j));
        dismatrix(j,i) = dismatrix(i,j);
    end
end

%% To test 100 realization and calculating the MPH of each one of them (SGeMS)

% number of realizations
NUM  = 100;
Dim  = 111;
Dimz = 1;

mphRe = zeros(NUM,2^16);
real100 = zeros(Dim,Dim,NUM);
for re = 0:NUM-1
    fprintf('%d\n',re);
    file_location = ['C:\Documents and Settings\mehrdadh\Desktop\', num2str(re),'.dat'];
    [datain, colnames, line1]=loadgeoeas(file_location);
    out = geoeas2matlab(datain,[Dim Dim Dimz]);
    real100(:,:,re+1) = out;
    mphRe(re+1,:) = calculateMPH(out);
end

[mdsMPHre, eMPHre] = cmdscale(pdist(mphRe));
plot3(mdsMPHre(:,1),mdsMPHre(:,2),mdsMPHre(:,3),'r.');
grid on; hold on
    

dismatriX = zeros(100,100);
for i = 1:100
    for j = i:100
        dismatriX(i,j) = haussdorff(real100(:,:,i),real100(:,:,j));
        dismatriX(j,i) = dismatriX(i,j);
    end
end

[rrr,ttt]=cmdscale(dismatrix);
[rrr1,ttt1]=cmdscale(dismatriX);

plot3(rrr(:,1),rrr(:,2),rrr(:,3),'ko')
grid on; hold on;
plot3(rrr1(:,1),rrr1(:,2),rrr1(:,3),'r*')

figure;parallelcoords(rrr1(:,1:12),'color',[0.6,0.6,0.6]); hold on;parallelcoords(rrr(:,1:12),'color','r')
figure;parallelcoords(rrr1(:,1:12),'color',[0.3,0.3,0.3],'Quantile',0.1); hold on;parallelcoords(rrr(:,1:12),'color','r','Quantile',0.1)


% Doing the distance function of JS divergence
jsmatrix1 = zeros(100,100);
jsmatrix2 = zeros(100,100);
for i = 1:100
    for j =i:100
        jsmatrix1(i,j) = kldiv(1:size(mph,2),(mph(i,:)+10^-15)./sum(mph(i,:)+10^-15),(mph(j,:)+10^-15)./sum(mph(j,:)+10^-15),'js');
        jsmatrix1(j,i) = jsmatrix1(i,j) ;
        jsmatrix2(i,j) = kldiv(1:size(mphRe,2),(mphRe(i,:)+10^-15)./sum(mphRe(i,:)+10^-15),(mphRe(j,:)+10^-15)./sum(mphRe(j,:)+10^-15),'js');
        jsmatrix2(j,i) = jsmatrix2(i,j) ;
        fprintf('%d\n',i*j);
    end
end
figure;
[rrrjs,tttjs]=cmdscale(jsmatrix1);
[rrr1js,ttt1js]=cmdscale(jsmatrix2);
plot3(rrrjs(:,1),rrrjs(:,2),rrrjs(:,3),'ko')
grid on; hold on;
plot3(rrr1js(:,1),rrr1js(:,2),rrr1js(:,3),'r*')

figure;parallelcoords(rrr1js(:,1:20),'color',[0.6,0.6,0.6]); hold on;parallelcoords(rrrjs(:,1:20),'color','r')
figure;parallelcoords(rrr1js(:,1:20),'color',[0.3,0.3,0.3],'Quantile',0.1); hold on;parallelcoords(rrrjs(:,1:20),'color','r','Quantile',0.1)


%% making the axis to be copied

c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off','zcolor', 'k', 'zgrid', 'off');
