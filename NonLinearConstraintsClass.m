classdef (HandleCompatible = true) NonLinearConstraintsClass < matlab.mixin.SetGet   
    properties(Constant, Hidden)
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
        R = @(rotAngle) [cos(rotAngle),-sin(rotAngle);sin(rotAngle),cos(rotAngle)];
    end
    
    properties
        track;  % Reference to the track object
        nStates;
        nInputs;
        nDecPerStep;
    end
    
    properties(Hidden)
        dT;
        dX;
        
        % Partial derivatives
        dXdU;
        dSdX;
    end
    
	methods
		%% Constructor
		function [obj] = NonLinearConstraintsClass(track,nStates,nInputs,nDecPerStep,dT)
            % Store the input data
            obj.track = track;
            obj.nStates = nStates;
            obj.nInputs = nInputs;
            obj.nDecPerStep = nDecPerStep;
            obj.dT = dT;
            
            % Define the symbolic variables
            U = sym('U', [nInputs, 1], 'real'); % Control vector u_1 = d_f / u_2 = F_x
            X = sym('X', [nStates, 1], 'real'); % State vector
            theta = sym('theta','real'); % Angle of the track
            
            % Slip angle functions in degrees
            radians2Degrees = 180/pi;
            a_f = radians2Degrees*(U(1)-atan2(X(4) + obj.a*X(6),X(2)));
            a_r = radians2Degrees*(-atan2((X(4)-obj.b*X(6)),X(2)));

            % Nonlinear tire dynamics
            phi_yf = (1-obj.Ey)*(a_f+obj.Shy) + (obj.Ey/obj.By)*atan(obj.By*(a_f+obj.Shy));
            phi_yr = (1-obj.Ey)*(a_r+obj.Shy) + (obj.Ey/obj.By)*atan(obj.By*(a_r+obj.Shy));
            F_yf = obj.Dy*sin(obj.Cy*atan(obj.By*phi_yf)) + obj.Svy;
            F_yr = obj.Dy*sin(obj.Cy*atan(obj.By*phi_yr)) + obj.Svy;

            % Vehicle dynamics
            obj.dX = [X(2)*cos(X(5))-X(4)*sin(X(5));...
                      (-obj.f*obj.W+obj.Nw*U(2)-F_yf*sin(U(1)))/obj.m+X(4)*X(6);...
                      X(2)*sin(X(5))+X(4)*cos(X(5));...
                      (F_yf*cos(U(1))+F_yr)/obj.m-X(2)*X(6);...
                      X(6);...
                     (F_yf*obj.a*cos(U(1))-F_yr*obj.b)/obj.Iz];
               
            % Compute the dXdU jacobian
            evalFcn_dXdU = subs(jacobian(obj.dX,U),{'X1','X2','X3','X4','X5','X6','U1','U2'},{'X(1)','X(2)','X(3)','X(4)','X(5)','X(6)','U(1)','U(2)'});
            obj.dXdU = eval(sprintf('@(X,U) [[%s;%s;%s;%s;%s;%s],[%s;%s;%s;%s;%s;%s]]',...
                        char(evalFcn_dXdU(1,1)),char(evalFcn_dXdU(2,1)),...
                        char(evalFcn_dXdU(3,1)),char(evalFcn_dXdU(4,1)),...
                        char(evalFcn_dXdU(5,1)),char(evalFcn_dXdU(6,1)),...
                        char(evalFcn_dXdU(1,2)),char(evalFcn_dXdU(2,2)),...
                        char(evalFcn_dXdU(3,2)),char(evalFcn_dXdU(4,2)),...
                        char(evalFcn_dXdU(5,2)),char(evalFcn_dXdU(6,2))));
                    
            % Compute the dSdX jacobian
            Vxy = obj.R(X(5))*[X(2);X(4)];
            trackNorm = [cos(theta);sin(theta)];
            evalFcn_dSdX = subs(jacobian(Vxy'*trackNorm,X),{'X1','X2','X3','X4','X5','X6'},{'X(1)','X(2)','X(3)','X(4)','X(5)','X(6)'});
            obj.dSdX = eval(sprintf('@(X,theta) [%s,%s,%s,%s,%s,%s]',char(evalFcn_dSdX(1)),char(evalFcn_dSdX(2)),char(evalFcn_dSdX(3)),char(evalFcn_dSdX(4)),char(evalFcn_dSdX(5)),char(evalFcn_dSdX(6))));
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
                if(iStep < nTotalStates)
                    useUPartials = true;
                else
                    useUPartials = false;
                end
                
                cIdxOffset = (iStep-1)*4;
                dIdxOffset = (iStep-1)*obj.nDecPerStep;
                stepS = obj.track.cartesian2Track(X(:,iStep));
                stepSeg = obj.track.getTrackSegment(stepS);
                stepPathAngle = stepSeg.getPathAngle(stepS);
                stepPartialSdX = obj.partialSdX(X(:,iStep),stepPathAngle);
                if(useUPartials)
                    stepPartialSdU = obj.partialSdU(X(:,iStep),U(:,iStep),stepPathAngle);
                end
                                
                % X Limits
                stepXLimCoeffs = [stepSeg.leftXCoeffs;stepSeg.rightXCoeffs];
                stepdXLimCoeffs = [stepSeg.leftdXCoeffs;stepSeg.rightdXCoeffs];
                xLimVals = [polyval(stepXLimCoeffs(1,:),stepS);polyval(stepXLimCoeffs(2,:),stepS)];
                [~,sortXIdx] = sort(xLimVals);
                c(1+cIdxOffset) =  (xLimVals(sortXIdx(1)) + posMargin) - X(1,iStep); % Min X limit
                c(2+cIdxOffset) =  X(1,iStep) - (xLimVals(sortXIdx(2)) - posMargin); % Max X limit
                
                % Compute the X gradients         
                mindXdS = polyval(stepdXLimCoeffs(sortXIdx(1),:),stepS);
                maxdXdS = polyval(stepdXLimCoeffs(sortXIdx(2),:),stepS);
                gradC((1:6)+dIdxOffset,1+cIdxOffset) = mindXdS.*stepPartialSdX + [-1,0,0,0,0,0];
                gradC((1:6)+dIdxOffset,2+cIdxOffset) = [1,0,0,0,0,0] - maxdXdS.*stepPartialSdX;
                
                if(useUPartials)
                    gradC((7:8)+dIdxOffset,1+cIdxOffset) = mindXdS.*stepPartialSdU;
                    gradC((7:8)+dIdxOffset,2+cIdxOffset) = maxdXdS.*stepPartialSdU;
                end
                
                % Y Limits
                stepYLimCoeffs = [stepSeg.leftYCoeffs;stepSeg.rightYCoeffs];
                stepdYLimCoeffs = [stepSeg.leftdYCoeffs;stepSeg.rightdYCoeffs];
                yLimVals = [polyval(stepYLimCoeffs(1,:),stepS);polyval(stepYLimCoeffs(2,:),stepS)];
                [~,sortYIdx] = sort(yLimVals);
                c(3+cIdxOffset) =  (yLimVals(sortYIdx(1)) + posMargin) - X(3,iStep); % Min Y limit
                c(4+cIdxOffset) =  X(3,iStep) - (yLimVals(sortYIdx(2)) - posMargin); % Max Y limit    
                
                % Compute the Y gradients
                mindYdS = polyval(stepdYLimCoeffs(sortYIdx(1),:),stepS);
                maxdYdS = polyval(stepdYLimCoeffs(sortYIdx(2),:),stepS);
                gradC((1:6)+dIdxOffset,3+cIdxOffset) = mindYdS.*stepPartialSdX + [0,0,-1,0,0,0];
                gradC((1:6)+dIdxOffset,4+cIdxOffset) = [0,0,1,0,0,0] - maxdYdS.*stepPartialSdX;
                
                if(useUPartials)
                    gradC((7:8)+dIdxOffset,1+cIdxOffset) = mindYdS.*stepPartialSdU;
                    gradC((7:8)+dIdxOffset,2+cIdxOffset) = maxdYdS.*stepPartialSdU;
                end
            end
            
            % Setup the equality constraints
            nEqualityConstraints = obj.nStates*(nTotalStates-1);
            ceq = zeros(nEqualityConstraints,1);
            gradCeq = zeros(size(D,1),nEqualityConstraints);
            predValues = forwardIntegrateControlInput(U',X(:,1)')';
            for iConst = 1:(nTotalStates-1)
                eqIdxOffset = obj.nStates*(iConst-1);
                prevIdx = iConst;
                nextIdx = prevIdx + 1;
                ceq((1:6) + eqIdxOffset,1) = X(:,nextIdx) - predValues(:,prevIdx);

                % Set the gradients
                decIdxOffset = obj.nDecPerStep*(iConst-1);
%                 gradCeq((1:6) + decIdxOffset, (1:6) + eqIdxOffset) = -(eye(6) + obj.partialX(X(:,prevIdx),U(:,prevIdx))');  % dG/dX(k)
%                 gradCeq((7:8) + decIdxOffset, (1:6) + eqIdxOffset) = -obj.partialU(X(:,prevIdx),U(:,prevIdx))';             % dG/dU(k)
%                 gradCeq((1:6) + decIdxOffset + obj.nDecPerStep, (1:6) + eqIdxOffset) = eye(6);                              % dG/dX(k+1)
            end        
        end

        function [J] = partialXdU(obj,xVec,uVec)
           J = obj.dXdU(xVec(:),uVec(:));
        end
        
        function [J] = partialSdX(obj,xVec,trackAngle)
            J = obj.dSdX(xVec(:),trackAngle);
        end
        
        function [J] = partialSdU(obj,xVec,uVec,trackAngle)
            J = obj.partialSdX(xVec,trackAngle)*obj.partialXdU(xVec,uVec);
        end
        
        function [stop] = outputFcn(obj,x,optimValues,state)
            obj.track.hWorkingLine.XData = x(1:obj.nDecPerStep:end);
            obj.track.hWorkingLine.YData = x(3:obj.nDecPerStep:end);
            drawnow;
            stop = false;
            disp('Output function');
        end
        
        function [stop] = plotFcn(obj,x,optimValues,state)
            stop = false;
%             disp('Plot function');
        end
    end
end