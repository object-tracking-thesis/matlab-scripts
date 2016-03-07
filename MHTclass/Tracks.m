% Tracks class for storing the posterior of different tracks. Essentially a
% container, a list of tracks and their Id. Instantiated empty, with new
% tracks added by using the member function 'addTrack'. Used as building
% block in Hypotheses.
%
% Tracks              - creates a Tracks object
% Tracks.addTrack(p1) - add the track p1
% Tracks.track        - gets the list of all tracks stored
% Tracks.trackId      - the corresponding trackId for each track
%
classdef Tracks < handle
    properties
        track;
        trackId;
    end
    
    methods (Access = public)
        function addTrack(this,posteriors)
            if ~isempty(posteriors)
                this.track = [this.track, posteriors];
                this.trackId = 1:length(this.track);
            else
               error('No Posteriors'); 
            end
        end
    end
end