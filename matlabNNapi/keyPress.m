function keyPress(src, e)
    switch e.Key
        case '1'
            evalin('base','isCarMat(i,j) = 1;');
        case '2'
            evalin('base','isCarMat(i,j) = 2;');
    end
end