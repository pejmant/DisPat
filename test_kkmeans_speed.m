% test time required wrt. number of data points
% the results are O(N^3)
clus = 80;
figure;hold on;
for i=200:50:1800
   tt=0;
   for j=1:20
       rr= randperm(i);
       II = K(rr,rr);
       tic;
       idx = dualkmeansFast(II,clus);
       tt= tt+toc;
   end
   plot(i,tt/20,'ko');
end



% test on time required wr.t. number of clusters
figure;hold on;
k=1;
for clus=3:3:300
   tt=0;
   for j=1:30
       tic;
       idx = dualkmeansFast(K,clus);
       tt= tt+toc;
   end
   time(k)=tt/30;
   k=k+1;
   
end
plot(3:3:300,time,':ko','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3);