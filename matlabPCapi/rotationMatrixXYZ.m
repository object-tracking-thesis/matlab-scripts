function rotationMatrixXYZ = rotationMatrixXYZ(angleX, angleY, angleZ)
% forms the rotation matrix in order X-Y-Z from the given angles

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

rotationMatrixXYZ = rotMatX(angleX)*rotMatY(angleY)*rotMatZ(angleZ);

end