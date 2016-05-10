% A function which calls upon the powers of Xzibit to draw a plot of your
% car, based on the statevector 'state'. Returns plothandle.
% 
% drawMyRide(state)
% drawMyRide(state, colSpec)
% drawMyRide(state, colSpec, LineSpec)
%
% plotHandles = drawMyRide(state,..)
%
%
function varargout = drawMyRide(state, varargin)
    if length(varargin) == 1
        colSpec = varargin{1};
        lineSpec = '-';
    elseif length(varargin) == 2
        colSpec = varargin{1};
        lineSpec = varargin{2};
    else
        colSpec = [0 0 0];
        lineSpec = '-';
    end
    
    if length(colSpec) == 1
        switch colSpec
            case 'b'
                colSpec = [0 0 1];
            case 'r'
                colSpec = [1 0 0];
            case 'g'
                colSpec = [0 1 0];
            case 'c'
                colSpec = [0 1 1];
            case 'y'
                colSpec = [1 1 0];
            case 'm'
                colSpec = [1 0 1];
            case 'k'
                colSpec = [0 0 0];
        end
    end

    % Corner i, width and length MGPs definitions
    corner_i = @(st) [ (st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4))),...
                       (st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)))];

    % Corner ii, width and length MGPs definitions
    corner_ii = @(st) [ (st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2)),...
                        (st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2))];

    % Corner iii, width and length MGPs definitions
    corner_iii = @(st) [ (st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + st(7)*cos(st(4))),...
                         (st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + st(7)*sin(st(4)))];

    % Corner iv, width and length MGPs definitions

    corner_iv = @(st) [ (st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(7)*cos(st(4))),...
                        (st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(7)*sin(st(4)))];


    % Define each vector     
    Allcorners = [  corner_i(state);
                   corner_ii(state);
                  corner_iii(state);
                   corner_iv(state);
                    corner_i(state)];
    
    % Define heading & carMiddle 
    
    carMiddle = [state(1), state(2)];
    
    directionX = [state(1) state(1)+cos(state(4))];
    directionY = [state(2) state(2)+sin(state(4))];
    
    plotHandle = plot(Allcorners(:,1), Allcorners(:,2), carMiddle(1), carMiddle(2),'x',directionX, directionY,'Color',colSpec, 'LineStyle',lineSpec);
        
    if nargout > 0    
        varargout{1} = plotHandle;
    end
                
end