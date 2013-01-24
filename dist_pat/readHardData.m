function hardData = readHardData(file_location, hardData, par)
% reads the hard data from a SGeMS dat file


temp = hardData((par.szRealization-par.Dim)/2+1:(par.szRealization-par.Dim)/2+par.Dim,(par.szRealization-par.Dim)/2+1:(par.szRealization-par.Dim)/2+par.Dim);


[datain, colnames, line1]=loadgeoeas(file_location);
if sum(datain(:,3)) == 0
    for i = 1:size(datain,1)
        temp(datain(i,1),datain(i,2),1) = datain(i,4);
    end
else
    for i = 1:size(datain,1)
        temp(datain(i,1),datain(i,2),datain(i,3)) = datain(i,4);
    end
end

temp = permute(temp, [2 1 3]);
temp = temp(par.Dim:-1:1,:,:);

hardData((par.szRealization-par.Dim)/2+1:(par.szRealization-par.Dim)/2+par.Dim,(par.szRealization-par.Dim)/2+1:(par.szRealization-par.Dim)/2+par.Dim) = temp;
