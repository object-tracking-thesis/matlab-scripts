function cuboid = fourCorners2Cube(p1,p2,p3,p4,height)
%create a cuboid description from 
%in:
%p1,p2,p3,p4 - 1x3 vectors, corner points of the lower xy-plane face of the 
%              cuboid, recorded in clockwise direction
%out:
%origin - 1x3 vector, center point of the cube
%length - scalar, x
%width  - scalar, y
%height - scalar, z
%rotation - 1x3 vector, [roll pitch yaw]

p = p1;
q = p2;
r = p3-p1;
s = p4-p2;

%finding the intersection point according to
%https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
t = cross((q-p),s)/cross(r,s);

center = p+t.*r;
length = sqrt(sum([p4-p1].^2));
width = sqrt(sum([p2-p1].^2));
anglevec = p4-p1;
rotation = [0 0 atan(anglevec(2)/anglevec(1))];  

cuboid = [center, length, width, height, rotation];

end