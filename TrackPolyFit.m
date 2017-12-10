classdef (HandleCompatible = true) TrackPolyFit < matlab.mixin.SetGet
	properties(Constant, Hidden)
        % Fit parameters
        order = 2;
    end
    
    properties
        % Line coefficients
        leftXCoeffs;
        leftYCoeffs;
        centerXCoeffs;
        centerYCoeffs;
        rightXCoeffs;
        rightYCoeffs;
        
        % dLine coefficients
        leftdXCoeffs;
        leftdYCoeffs;
        centerdXCoeffs;
        centerdYCoeffs;
        rightdXCoeffs;
        rightdYCoeffs;
    end
		
    methods(Static)
        function [output] = interpPoly(Cx,Cy,sPts)
            xPts = polyval(Cx,sPts);
            yPts = polyval(Cy,sPts);
            output = [xPts(:)';yPts(:)'];
        end
        
        function [poly4Vals] = square2OrderPoly(c)
            assert(length(c) == 3);
            poly4Vals = [c(1)^2, 2*c(1)*c(2), 2*c(1)*c(3) + c(2)^2, 2*c(2)*c(3), c(3)^2];
        end
    end
    
	methods
		%% Constructor
		function [obj] = TrackPolyFit(track,fitS,trackPolySpacing)
            nPtsExtra = 1;
            firstIdx = max([find(track.arc_s <= (fitS - trackPolySpacing),1,'Last') - nPtsExtra,1]);
            lastIdx  = min([find(track.arc_s >= (fitS + trackPolySpacing),1,'First')+ nPtsExtra,length(track.arc_s)]);
            fitIdx = firstIdx:lastIdx;
            
            % Fit the polynomials for each component
            obj.leftXCoeffs     = polyfit(track.arc_s(fitIdx),track.bl(1,fitIdx),obj.order);
            obj.leftYCoeffs     = polyfit(track.arc_s(fitIdx),track.bl(2,fitIdx),obj.order);
            obj.centerXCoeffs   = polyfit(track.arc_s(fitIdx),track.cline(1,fitIdx),obj.order);
            obj.centerYCoeffs   = polyfit(track.arc_s(fitIdx),track.cline(2,fitIdx),obj.order);
            obj.rightXCoeffs    = polyfit(track.arc_s(fitIdx),track.br(1,fitIdx),obj.order);
            obj.rightYCoeffs    = polyfit(track.arc_s(fitIdx),track.br(2,fitIdx),obj.order);
            
            % Compute the derivatives
            obj.leftdXCoeffs = polyder(obj.leftXCoeffs);
            obj.leftdYCoeffs = polyder(obj.leftYCoeffs);
            obj.centerdXCoeffs = polyder(obj.centerXCoeffs);
            obj.centerdYCoeffs = polyder(obj.centerYCoeffs);
            obj.rightdXCoeffs = polyder(obj.rightXCoeffs);
            obj.rightdYCoeffs = polyder(obj.rightYCoeffs);
        end
		
		function [pos] = interpLeft(obj,sPts)
            pos = obj.interpPoly(obj.leftXCoeffs,obj.leftYCoeffs,sPts);
        end
        
        function [pos] = interpCenter(obj,sPts)
            pos = obj.interpPoly(obj.centerXCoeffs,obj.centerYCoeffs,sPts);
        end
        
        function [pos] = interpRight(obj,sPts)
            pos = obj.interpPoly(obj.rightXCoeffs,obj.rightYCoeffs,sPts);
        end
        
        function [angle] = getPathAngle(obj,sPts)
            angle = atan2(polyval(obj.centerdYCoeffs,sPts),polyval(obj.centerdXCoeffs,sPts));
        end
        
    end
end