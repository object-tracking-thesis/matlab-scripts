function [] = waitforkey()
    now = 1;
    while now
        keydown = waitforbuttonpress;
        if keydown == 1
            now = 0;
        end
    end
end