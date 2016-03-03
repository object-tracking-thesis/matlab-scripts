function transformedFrame = transformFrameTransMat(lidarFrame, transMat, staticPos)
% Takes a lidarFrame and translates+rotates all coordinates according to the
% 4x4 affine matrix 'transMat'

%position difference between static and live
transMat(1:3,4) = transMat(1:3,4) - staticPos';

%transform the padded positions of each point via affine transformation
lidarFrame = [lidarFrame ones(length(lidarFrame),1)];
transformedFrame = (transMat*lidarFrame')';                

%rotate, then translate!
%transformedFrame = lidarFrame;
%transformedFrame = (transMat(1:3,1:3)*transformedFrame')';
%transformedFrame = transformedFrame + repmat(transMat(1:3,4)',length(transformedFrame),1);

%remove the padding
transformedFrame = transformedFrame(:,1:3);
end