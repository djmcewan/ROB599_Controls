classdef (HandleCompatible = true) SplineClass < matlab.mixin.SetGet
	properties(Hidden)
		XRaw;
        YRaw;
        SRaw;
        dSize = 2;
    end
		
	methods
		%% Constructor
		function [obj] = SplineClass(X,Y,S)
			obj.XRaw = X(:);
			obj.YRaw = Y(:);
            obj.SRaw = S(:);
        end
		
		function [pos] = interp(obj,pts)
            pos = [obj.interpX(pts),...
                   obj.interpY(pts)];
        end
        
        function [pos] = interpX(obj,pts)
            %pos = interp1(obj.SRaw,obj.XRaw,pts,'spline','extrap');
            pos = interp1(obj.SRaw,smooth(obj.SRaw,obj.XRaw),pts,'spline','extrap');
        end
        
        function [pos] = interpY(obj,pts)
            %pos = interp1(obj.SRaw,obj.YRaw,pts,'spline','extrap');
            pos = interp1(obj.SRaw,smooth(obj.SRaw,obj.YRaw),pts,'spline','extrap');
        end
        
        function [T] = theta(obj,pts)
            pts = pts(:)';
            allPts = bsxfun(@plus,obj.dSize.*[-1;0;1],pts);
            dX = mean(diff(obj.interpX(allPts))./obj.dSize);
            dY = mean(diff(obj.interpY(allPts))./obj.dSize);
            T = atan2(dY,dX);
        end
        
        function [kappa] = curvature(obj,pts)
            pts = pts(:)';
            allPts = bsxfun(@plus,obj.dSize.*[-1;0;1],pts);
            kappa = mean(wrapToPi(diff([obj.theta(allPts(1,:));obj.theta(allPts(2,:));obj.theta(allPts(3,:))]))./obj.dSize);
        end
        
        function [R] = radiusOfCurvature(obj,pts)
            R = abs(1./obj.curvature(pts));
        end
        
        function [minS] = findMinRadiusLocation(obj)
            pts = linspace(obj.SRaw(2),obj.SRaw(end-1),100);
            [~,minIdx] = min(obj.radiusOfCurvature(pts));
            minS = pts(minIdx);
        end
        
        function [direction] = curveDirection(obj)
            direction = mode(obj.curvature(obj.SRaw) < 0);
        end
    end
end
