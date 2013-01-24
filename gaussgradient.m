function [gx,gy]=gaussgradient(IM,sigma)
%GAUSSGRADIENT Gradient using first order derivative of Gaussian.

%determine the appropriate size of kernel. The smaller epsilon, the larger
%size.
epsilon=1e-2;
halfsize=ceil(sigma*sqrt(-2*log(sqrt(2*pi)*sigma*epsilon)));
size=2*halfsize+1;
%generate a 2-D Gaussian kernel along x direction
for i=1:size
    for j=1:size
        u=[i-halfsize-1 j-halfsize-1];
        hx(i,j)=gauss(u(1),sigma)*dgauss(u(2),sigma);
    end
end
hx=hx/sqrt(sum(sum(abs(hx).*abs(hx))));
%generate a 2-D Gaussian kernel along y direction
hy=hx';
%2-D filtering
gx=imfilter(IM,hx,'replicate','conv');
gy=imfilter(IM,hy,'replicate','conv');

%Plot
figure('name','Gradient');
imshow(IM(1:20,1:20),'InitialMagnification','fit');
hold on;
quiver(gx(1:20,1:20),gy(1:20,1:20));
title(['sigma=',sigma]);


%calculate orientation
hg=im2col(gx(1:20,1:20),[1,1],'sliding')';
hv=im2col(gy(1:20,1:20),[1,1],'sliding')';
[u,s,v]=svd([hg,hv]);
figure;
plot(hg,hv,'.');
R = (s(1,1)-s(2,2))/(s(1,1)+s(2,2))
theta = atand(-v(2,1)/v(1,1)) + 90


function y = gauss(x,sigma)
%Gaussian
y = exp(-x^2/(2*sigma^2)) / (sigma*sqrt(2*pi));

function y = dgauss(x,sigma)
%first order derivative of Gaussian
y = -x * gauss(x,sigma) / sigma^2;