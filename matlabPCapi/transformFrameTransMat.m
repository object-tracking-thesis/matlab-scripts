function transformedFrame = transformFrameTransMat(lidarFrame, transMat)
% Takes a lidarFrame and translates+rotates all coordinates according to the
% 4x4 affine matrix 'transMat'

%transform the padded positions of each point via affine transformation
lidarFrame = [lidarFrame ones(size(lidarFrame,1),1)];
%transMat(1:3,1:3) = eye(3); %remove rotation
transformedFrame = (transMat*lidarFrame')';                

%rotate, then translate!
%transformedFrame = lidarFrame;
%transformedFrame = (transMat(1:3,1:3)*transformedFrame')';
%transformedFrame = transformedFrame + repmat(transMat(1:3,4)',size(transformedFrame,1),1);

%remove the padding
transformedFrame = transformedFrame(:,1:3);
end