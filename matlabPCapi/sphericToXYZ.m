temp = lidarDataSpheric;

n = length(temp);
out = zeros(n,3);

%convert to spheric coordinates
for i=1:n
    out(i,1) = temp(i,1)*sin(temp(i,2))*cos(temp(i,3)); %x
    out(i,2) = temp(i,1)*sin(temp(i,2))*sin(temp(i,3)); %y
    out(i,3) = temp(i,1)*cos(temp(i,2)); %z
end

lidarData{1} = out;