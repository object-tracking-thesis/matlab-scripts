function btn = mydialog
% Makes a dialog box where you can choose to step the twoTracks.m
% simulation forward or backward. If the box i closed, the simulation will
% progress automatically.
%
% Code based on template, don't really know how it works, so it might be
% smart to not change stuff.
    d = dialog('Position',[1300 600 300 100],'Name','Time Stepper');

    uicontrol('Parent',d,...
        'Style','text',...
        'Position',[50 40 210 40],...
        'String','Click the close button when you''re done.');

    uicontrol('Parent',d,...
        'Position',[50 30 70 25],...
        'String','<<',...
        'Callback',@prev);

    uicontrol('Parent',d,...
        'Position',[190 30 70 25],...
        'String','>>',...
        'Callback',@next);
    
    btn = 10;
    
    function prev(~, ~)
        btn = -1; 
            delete(d)
    end

    function next(~, ~)
        btn = 1;
            delete(d)
    end
    uiwait(d);


end