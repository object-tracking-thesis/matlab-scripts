



function corner = getCorner2(massVec1, massVec2, heading)
    
    % Make sure they are unit vectors
    massVec1 = massVec1./norm(massVec1,2);
    massVec2 = massVec2./norm(massVec2,2);

    vHeading = [cos(heading), sin(heading)];
    % Find which vector is along length & width
    if abs(dot(massVec1, vHeading)) < abs(dot(massVec2,vHeading))
        vLength =  massVec2;
        vWidth  =  massVec1;
    else
        vLength =  massVec1;
        vWidth  =  massVec2;
    end
    
    % Generate composite vectors     
    compV = vLength + vWidth;
            
    phi = heading;
    
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    
    % Negate heading 
    
    compV = compV*R;
    
    angle = mod(atan2(compV(2), compV(1)),2*pi);
    
    
    if (0 < angle) && (angle < pi/2)
        corner = 1;
    elseif (pi/2 < angle) &&  (angle < pi)
        corner = 4;
    elseif (pi < angle) && (angle < 3*pi/2)
        corner = 3;
    elseif (3*pi/2 < angle) && (angle < 2*pi)
        corner = 2;
    end        

end