classdef (HandleCompatible = true) NonLinearConstraintsClass < matlab.mixin.SetGet   
    properties(Constant)
        % Vehicle Params
        W = 13720;
        Nw = 2;
        f = 0.01;
        Iz = 2667;
        a = 1.35;
        b = 1.45;
        By = 0.27;
        Cy = 1.2;
        Dy = 2921;
        Ey = -1.6;
        Shy = 0;
        Svy = 0;
        m = 1400; 
    end
    
    properties
        track;  % Reference to the track object
        nStates;
        nInputs;
        nDecPerStep;
    end
    
	methods
		%% Constructor
		function [obj] = NonLinearConstraintsClass(track,nStates,nInputs,nDecPerStep)
            % Store the input data
            obj.track = track;
            obj.nStates = nStates;
            obj.nInputs = nInputs;
            obj.nDecPerStep = nDecPerStep;
            
        end
		
		function [c, ceq, gradC, gradCeq] = constraintFcn(obj,D)   
            % Extract the states / inputs from the decision vector
            U = [D((obj.nStates + 1):obj.nDecPerStep:end),D((obj.nStates + 2):obj.nDecPerStep:end)]';
            X = [D(1:obj.nDecPerStep:end),D(2:obj.nDecPerStep:end),D(3:obj.nDecPerStep:end),D(4:obj.nDecPerStep:end),D(5:obj.nDecPerStep:end),D(6:obj.nDecPerStep:end)]';
            nTotalStates = size(X,2);

            % Setup the inequality constraints
            nInEqualityConstraints = 4*nTotalStates;
            c = zeros(nInEqualityConstraints,1);
            gradC = zeros(size(D,1),nInEqualityConstraints);
            posMargin = 0.05;
            for iStep = 1:nTotalStates
                cIdxOffset = (iStep-1)*4;
                stepS = obj.track.cartesian2Track(X(:,iStep));
                stepSeg = obj.track.getTrackSegment(stepS);
                                
                % X Limits
                stepXLimCoeffs = [stepSeg.leftXCoeffs;stepSeg.rightXCoeffs];
                xLimVals = [polyval(stepXLimCoeffs(1,:),stepS);polyval(stepXLimCoeffs(2,:),stepS)];
                [~,sortXIdx] = sort(xLimVals);
                c(1+cIdxOffset) =  (xLimVals(sortXIdx(1)) + posMargin) - X(1,iStep); % Min X limit
                c(2+cIdxOffset) =  X(1,iStep) - (xLimVals(sortXIdx(2)) - posMargin); % Max X limit
                
                % Compute the X gradients
                
                % Y Limits
                stepYLimCoeffs = [stepSeg.leftYCoeffs;stepSeg.rightYCoeffs];
                yLimVals = [polyval(stepYLimCoeffs(1,:),stepS);polyval(stepYLimCoeffs(2,:),stepS)];
                [~,sortYIdx] = sort(yLimVals);
                c(3+cIdxOffset) =  (yLimVals(sortYIdx(1)) + posMargin) - X(3,iStep); % Min Y limit
                c(4+cIdxOffset) =  X(3,iStep) - (yLimVals(sortYIdx(2)) - posMargin); % Max Y limit    
                
                % Compute the Y gradients

            end
            
            % Setup the equality constraints
            nEqualityConstraints = 6*nTotalStates;
            ceq = zeros(nEqualityConstraints,1);
            gradCeq = zeros(size(D,1),nEqualityConstraints);
            
        end

    end
end