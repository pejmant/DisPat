function Y = Yrec(Z,model)

% same function as kpcarec, except that Z is already the projections no projection is done


[dim, num_data] = size(Z);

% allocate memory
Y = zeros(size(model.sv.X,1),num_data);
img = model;

fprintf('Computing preimages');
for i=1:num_data,
   fprintf('.'); 
   img.Alpha = model.Alpha*(Z(:,i) - model.b);
   Y(:,i) = rbfpreimg(img);       % Schoelkopf's fix-point algorithm
   %   Y(:,i) = rbfpreimg2(img);  % Gradient method
    %    Y(:,i) = rbfpreimg3(img,7);     % Kwok & Tsang
    
end
fprintf('done\n');

return;
