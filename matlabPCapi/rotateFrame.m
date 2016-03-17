function rotatedFrame = rotateFrame(lidarFrame, angleZ, angleY, angleX, varargin)
% Takes a lidarFrame and rotates all coordinates around Z-axis 
% in the frame with 'angle' (input in radians). An additional argument may
% be passed, specifying any angle offset to be used when rotating.
% lidarframe is of Nx3 size, angle is scalar, offset is also scalar. 

if length(varargin) == 1
    angleOffset = varargin{1};
elseif length(varargin) > 1
   error('Too many function inputs');
else
   angleOffset = 0;
end

% Rotation function 
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
rotatedFrame = (rotMatZ((angleZ+angleOffset))*lidarFrame')';
rotatedFrame = (rotMatY((angleY+angleOffset))*rotatedFrame')';
rotatedFrame = (rotMatX((angleX+angleOffset))*rotatedFrame')'; 

end