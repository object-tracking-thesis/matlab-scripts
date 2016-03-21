function [] = plotCubes(origin, length, width, height, rot_angle, tuneCubeSize, groundPlaneZ)
%plot the static map cubes (mainly walls)
%in:
%origin - 1x3 vector, center point of the cube
%length - scalar, x
%width - scalar, y
%height - scalar, z
%rot_angle - 1x3 vector, [roll pitch yaw]
%tuneCubeSize - 1x3 vector, tune cube to be bigger/smaller
%groundPlaneZ - scalar, offset in z direction

default_x=[0 1 1 0 0 0;
    1 1 0 0 1 1;
    1 1 0 0 1 1;
    0 1 1 0 0 0]-0.5;

default_y=[0 0 1 1 0 0;
    0 1 1 0 0 0;
    0 1 1 0 1 1;
    0 0 1 1 1 1]-0.5;

default_z=[0 0 0 0 0 1;
    0 0 0 0 0 1;
    1 1 1 1 0 1;
    1 1 1 1 0 1]-0.5;

cubeSize = [length,width,height]-tuneCubeSize;

default_x = default_x.*cubeSize(1);
default_y = default_y.*cubeSize(2);
default_z = default_z.*cubeSize(3);


cube_rot_x = [1, 0, 0;...
    0, cos(rot_angle(1)), -sin(rot_angle(1));...
    0, sin(rot_angle(1)), cos(rot_angle(1))];

cube_rot_y = [cos(rot_angle(2)), 0, sin(rot_angle(2));...
    0, 1, 0;...
    -sin(rot_angle(2)), 0, cos(rot_angle(2))];

cube_rot_z = [cos(rot_angle(3)), -sin(rot_angle(3)), 0;...
    sin(rot_angle(3)), cos(rot_angle(3)), 0;...
    0, 0, 1];

for iRot = 1:6
    for jRot = 1:4
        tmpPoint = [default_x(jRot,iRot);...
            default_y(jRot,iRot);default_z(jRot,iRot)];
        tmpPoint = cube_rot_x*tmpPoint;
        tmpPoint = cube_rot_y*tmpPoint;
        tmpPoint = cube_rot_z*tmpPoint;
        default_x(jRot,iRot) = tmpPoint(1);
        default_y(jRot,iRot) = tmpPoint(2);
        default_z(jRot,iRot) = tmpPoint(3);
    end
end

cubeSize = cube_rot_x*cubeSize';
cubeSize = cube_rot_y*cubeSize;
cubeSize = cube_rot_z*cubeSize;

origin(3) = groundPlaneZ + abs(cubeSize(3)/2);

default_x = default_x + origin(1);
default_y = default_y + origin(2);
default_z = default_z + origin(3);

for iPatch=1:6
    h=fill3(default_x(:,iPatch),default_y(:,iPatch),default_z(:,iPatch),'y','facecolor', [1,0.7,0],'facealpha',0.7);
    set(h,'edgecolor','k');
end

end
