% Scan class for storing measurements gathered at a time step k. Stores the
% measurements in no specific order in a list, and assigns a Id to each
% measurement for easy access. Can be instantiated as an empty object or by
% directly inserting measurements.
%
% Scan(meas)                 - instantiates the Scan object with the measurements meas
% Scan                       - instantiates the Scan object with no measurements
% Scan.addMeasurements(meas) - adds 'meas' to the current list of measurements. 
%                              Dimension (rows) must be the same accross all stored
%                              measurements.
% Scan.measurements          - get all measurements.
% Scan.measId                - get id for stored measurements.
%
%

classdef Scan < handle
    properties
        measurements;
        measId;
    end
    
    methods (Access = public)
        function this = Scan(meas)
            if nargin  == 1
                addMeasurements(this,meas);
            elseif nargin > 1
                error('Wrong number of input arguments (1 expected)')
            end
        end
        
        function addMeasurements(this,meas)
            this.measurements = [this.measurements meas];
            [~, c] = size(this.measurements);
            this.measId = 1:c;
        end
    end
end