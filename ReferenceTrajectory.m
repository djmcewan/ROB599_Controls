classdef (HandleCompatible = true) ReferenceTrajectory < matlab.mixin.SetGet
    properties
        uRef;
        xRef;
        nRefStates;
        sRef;    
    end
    
    properties(Access = private)
        filename;
    end
    
    methods
		%% Constructor
		function [obj] = ReferenceTrajectory(trackObj)
            % Define the reference filename
            obj.filename = 'ROB599_ControlsProject_Team23.mat';
            
            % Part 1. ROB599_ControlsProject_part1_input (for the vector)
            % Part 2. ROB599_ControlsProject_part2_function (for the function)
            refData = load(obj.filename);
            obj.uRef = refData.ROB599_ControlsProject_part1_input;
            dataFieldnames = fieldnames(refData);
            if(any(strcmp(dataFieldnames,'xRef')))
                obj.xRef = refData.xRef;
            else
                xRef = forwardIntegrateControlInput(obj.uRef);
                save(obj.filename,'xRef','-append');
                obj.xRef = xRef;
            end
            obj.nRefStates = size(obj.xRef,1);
            
            % Parameterize the reference trajectory as a functions of s
            if(any(strcmp(dataFieldnames,'sRef')))
                obj.sRef = refData.sRef;
            else
                sRef = NaN(obj.nRefStates,1);
                hWaitbar = waitbar(0,'Parameterizing reference trajectory');
                for iState = 1:obj.nRefStates
                    waitbar(iState/obj.nRefStates);
                    sRef(iState) = trackObj.cartesian2Track(obj.xRef(iState,:));
                end
                close(hWaitbar);
                save(obj.filename,'sRef','-append');
                obj.sRef = sRef;
            end
        end
    end
end