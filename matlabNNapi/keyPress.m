function keyPress(src, e)
    switch e.Key
        case '1'
            evalin('base','isCarMat(i,j) = 1;');
        case '2'
            evalin('base','isCarMat(i,j) = 2;');
        case '3'
            evalin('base','isCarMat(i,j) = 3;');
        case '4'
            evalin('base','isCarMat(i,j) = 4;');
        case '5'
            evalin('base','isCarMat(i,j) = 5;');    
        case 'leftarrow'
            evalin('base','j=j-2;');
        case 'downarrow'
            evalin('base','i=i-1;');
    end
end