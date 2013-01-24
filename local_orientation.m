function local_orientation(I)


seh = strel([1 1 1]);
sev = strel([1;1;1]);

horizontal_gradient = imdilate(I,seh) - imerode(I,seh);
vertical_gradient = imdilate(I,sev) - imerode(I,sev);

subplot(1,2,1)
imshow(horizontal_gradient, []), title('Horizontal gradient')

subplot(1,2,2)
imshow(vertical_gradient, []), title('Vertical gradient')

hg=im2col(horizontal_gradient(1:20,1:20),[1,1],'sliding')';
hv=im2col(vertical_gradient(1:20,1:20),[1,1],'sliding')';
[u,s,v]=svd([hg,hv]);
figure;
plot(hg,hv,'.');
R = (s(1,1)-s(2,2))/(s(1,1)+s(2,2))
theta = atand(v(2,1)/v(1,1)) + 90

end
