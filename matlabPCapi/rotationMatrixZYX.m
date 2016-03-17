function rotationMatrixZYX = rotationMatrixZYX(angleZ, angleY, angleX)
% forms the rotation matrix in order Z-Y-X from the given angles

% Rotation functions
rotMatZ = @(theta) [cos(theta) -sin(theta) 0;
                    sin(theta)  cos(theta) 0;
                    0           0          1];
                
rotMatY = @(theta) [cos(theta) 0 sin(theta)  ;
                    0          1            0;
                    -sin(theta) 0 cos(theta)];

rotMatX = @(theta) [1          0            0;
                    0  cos(theta) -sin(theta);
                    0  sin(theta) cos(theta)];                

%first around Z, then around Y, then around X
%http://www.oxts.com/Downloads/Products/Inertial2/inertialPlusman.pdf
rotationMatrixZYX = rotMatZ(angleZ)*rotMatY(angleY)*rotMatX(angleX);

end