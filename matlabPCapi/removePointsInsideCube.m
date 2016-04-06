function points = removePointsInsideCube(cube, points)
%return whether a point is inside a cube or not
%in:
%cube - 9x1 vector, x,y,z(center), length,width,height, roll,pitch,yaw 
%frame - nx3 vector, x,y,z
%out:
%remaining - nx3 vector, x,y,z the remaining points not inside the cube

%TODO why -4 in the z direction?
%cube = cube'
[px, py, pz] = cubePlanes(cube(1:3)',cube(4),cube(5),cube(6),cube(7:9)',-2,cube(3)-4);
nFaces = size(px,2);
origin = repmat(cube(1:3)',1,nFaces);

j=1:nFaces;
%3 points for each of the faces of the cube
p11 = [px(1,j); py(1,j); pz(1,j)];
p12 = [px(2,j); py(2,j); pz(2,j)];
p13 = [px(3,j); py(3,j); pz(3,j)];

%normal vectors for each of the faces of the cube
normalVectors = cross(p12-p11,p13-p11);

%calculate distance for the center point of the cube from all 6 faces,
%normalize it and build a repeat matrix over the number of all points
%(6x6 matrix for each point, row repeated)
od = (origin-p11)'*normalVectors;
od = od./abs(od);
od = repmat(od,size(points,1),1);

%calculate distance for a repeat matrix of all points and normalize it
%(6x6 matrix for each point, row repeated)
pointMat = repelem(points',1,nFaces);
pd = (pointMat-repmat(p11,1,size(points,1)))'*normalVectors;
pd = pd./abs(pd);

%if pd and od are equal then both points have to lie on the same side of 
%each of the 6 faces of the cube, otherwise the point lies outside the cube
%subtract od from pd and sum up the 6x6 submatrices, t is finally a vector
%of length(points) that is equal to zero if the point lies within the cube
t = pd-od;
t = abs(t);
t = sum(t,2);
t = sum(vec2mat(t,6),2);

%remove all points that are within the cube
points(t==0,:) = [];