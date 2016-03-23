function remainingPoints = removePointsInsideCube(cube, points)
%return whether a point is inside a cube or not
%in:
%cube - 9x1 vector, x,y,z(center), length,width,height, roll,pitch,yaw 
%frame - nx3 vector, x,y,z
%out:
%remaining - nx3 vector, x,y,z the remaining points not inside the cube

remainingPoints = zeros(60000,3);

%TODO why -2 in the z direction?
[px, py, pz] = cubePlanes(cube(1:3)',cube(4),cube(5),cube(6),cube(7:9)',-1.5,cube(3)-2.5);
nFaces = size(px,2);
origin = repmat(cube(1:3)',1,nFaces);

j=1:nFaces;
%3 points for each of the faces of the cube
p11 = [px(1,j); py(1,j); pz(1,j)];
p12 = [px(2,j); py(2,j); pz(2,j)];
p13 = [px(3,j); py(3,j); pz(3,j)];

%normal vectors for each of the faces of the cube
normalVectors = cross(p12-p11,p13-p11);

%calculate distance for the center point of the cube from all 6 faces and
%normalize it
od = (origin-p11)'*normalVectors;
od = od./abs(od);

for i=1:size(points,1)
    point = repmat(points(i,:)',1,nFaces);
    %calculate distance for the test point from all 6 faces and normalize it
    pd = (point-p11)'*normalVectors;
    pd = pd./abs(pd);
    
    %both points have to lie on the same side of each of the 6 faces of the
    %cube, otherwise the point lies outside the cube
    t = pd-od;
    
    %..if the test matrix is not all zeros, then the point and the cube origin
    %are on different sides of at least one face and therefore the point
    %lies outside of the cube, therefore keep...
    if ~(all(t(:) == 0))
        remainingPoints(i,:) = points(i,:);
    end
end

remainingPoints(all(remainingPoints == 0, 2),:) = [];

%remainingPoints = points;