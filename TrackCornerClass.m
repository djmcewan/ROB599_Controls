classdef (HandleCompatible = true) TrackCornerClass	< matlab.mixin.SetGet
	properties(Constant, Hidden)
        LEFT = 1;
        RIGHT = -1;
    end
    
    properties(Hidden)
        innerSpline;
        centerSpline;
        outerSpline;
        entryTheta;
        exitTheta;
        raceCornerTheta;
        cornerSVector;
        innerCurPoly;
        trackdTheta;
        track;
        rInner;
        rOuter;
        rRace;
    end
		
    properties
        cornerNum = [];
        apexS;
        startS;
        stopS;
        direction;
        innerEntryXY;
        innerExitXY;
        centerEntryXY;
        centerExitXY;
        outerEntryXY;
        outerExitXY;
    end
    
	methods
		%% Constructor
		function [obj] = TrackCornerClass(track,startS,apexS,stopS,direction)
            % Set object properties
            obj.track = track;
            obj.startS = startS;
            obj.apexS = apexS;
            obj.stopS = stopS;
            obj.direction = direction;
                        
            % Determine which points should be used to define the arc
            prevIdx = find(track.arc_s < obj.startS,1,'Last');
            nextIdx = find(track.arc_s > obj.stopS,1,'First');
            fitIdx = prevIdx:nextIdx;    
            fitLeft = track.bl(1:2,fitIdx);
            fitRight = track.br(1:2,fitIdx);
            fitCenter = track.cline(1:2,fitIdx);
            fitS = track.arc_s(fitIdx);
            obj.cornerSVector = linspace(obj.startS,obj.stopS,300);
            
            % Fit the centerline
            obj.centerSpline = SplineClass(fitCenter(1,:),fitCenter(2,:),fitS);
    
            % Fit the inner and outer lines
            if(obj.direction == obj.LEFT)
                obj.innerSpline = SplineClass(fitLeft(1,:),fitLeft(2,:),fitS);
                obj.outerSpline = SplineClass(fitRight(1,:),fitRight(2,:),fitS);
            else
                obj.innerSpline = SplineClass(fitRight(1,:),fitRight(2,:),fitS);
                obj.outerSpline = SplineClass(fitLeft(1,:),fitLeft(2,:),fitS);
            end
            
            % Look at the curvature and find the critical points
            nPoly = 2;
            obj.innerCurPoly = polyfit(obj.cornerSVector,obj.innerSpline.curvature(obj.cornerSVector),nPoly);           
            %obj.innerCurPoly = polyfit(fitS,track.fdtheta(fitS),nPoly);           
            fitCurvature = obj.direction.*obj.getCurvature(obj.cornerSVector);
            [~,apexIdx] = max(fitCurvature);
            obj.apexS = obj.cornerSVector(apexIdx);
            zeroS = real(roots(obj.innerCurPoly));
            obj.startS = max([min(zeroS),obj.startS]);
            obj.stopS = min([max(zeroS),obj.stopS ]);          
                        
            % Define corner properties
            obj.cornerSVector = obj.startS:0.3:obj.stopS;
            obj.entryTheta = obj.centerSpline.theta(obj.startS);
            obj.exitTheta = obj.centerSpline.theta(obj.stopS);
            obj.trackdTheta = track.fdtheta(obj.cornerSVector);
            
            % Determine the radius
            extraSpace = 3.5;
            obj.innerSpline.radiusOfCurvature(obj.apexS);
            obj.rInner = obj.direction.*obj.innerSpline.radiusOfCurvature(obj.apexS);
            obj.rOuter = obj.rInner + obj.direction.*(obj.getTrackWidth(obj.apexS) - extraSpace);
            
            % Determine the racing parameters            
            obj.raceCornerTheta = wrapToPi(obj.exitTheta - obj.entryTheta);
            obj.rRace = obj.rInner + (obj.rOuter - obj.rInner)/(1-cos(obj.raceCornerTheta/2));
            
            % Entry/Exit points
            obj.innerEntryXY = obj.innerSpline.interp(obj.startS);
            obj.innerExitXY = obj.innerSpline.interp(obj.stopS);
            obj.centerEntryXY = obj.centerSpline.interp(obj.startS);
            obj.centerExitXY = obj.centerSpline.interp(obj.stopS);
            obj.outerEntryXY = obj.outerSpline.interp(obj.startS);
            obj.outerExitXY = obj.outerSpline.interp(obj.stopS);
        end
        
        function [] = setCornerNumber(obj,num)
            obj.cornerNum = num;
        end
		
		function [] = plotCornerAnalysis(obj)
            % Record the previous working figure
            prevFigFocus = gcf;
            
            % Plot the corner being analyzed
            hFit = figure(200+obj.cornerNum); clf; hold on;
            set(hFit,'NumberTitle','Off','Name',['Corner Analysis - ',num2str(obj.cornerNum)]);
                        
            % Plot the fit curves
            subplot(2,1,1); hold on;
            plot(obj.innerSpline.interpX(obj.cornerSVector),obj.innerSpline.interpY(obj.cornerSVector),'r','LineWidth',3);
            plot(obj.centerSpline.interpX(obj.cornerSVector),obj.centerSpline.interpY(obj.cornerSVector),'c','LineWidth',3);
            plot(obj.outerSpline.interpX(obj.cornerSVector),obj.outerSpline.interpY(obj.cornerSVector),'g','LineWidth',3);
            
            % Plot the raw points
            plot(obj.innerSpline.XRaw,obj.innerSpline.YRaw,'.-k','MarkerSize',20);
            plot(obj.centerSpline.XRaw,obj.centerSpline.YRaw,'.-k','MarkerSize',20);
            plot(obj.outerSpline.XRaw,obj.outerSpline.YRaw,'.-k','MarkerSize',20);

            % Plot location of corner from manual defintion
            plot([obj.innerSpline.interpX(obj.apexS),obj.centerSpline.interpX(obj.apexS),obj.outerSpline.interpX(obj.apexS)],[obj.innerSpline.interpY(obj.apexS),obj.centerSpline.interpY(obj.apexS),obj.outerSpline.interpY(obj.apexS)],'-k','LineWidth',2,'MarkerSize',30);

            % Plot location of maximum curviture
            minS = obj.centerSpline.findMinRadiusLocation();
            plot([obj.innerSpline.interpX(minS),obj.outerSpline.interpX(minS)],[obj.innerSpline.interpY(minS),obj.outerSpline.interpY(minS)],'-m','LineWidth',2);
            plot(obj.innerSpline.interpX(minS),obj.innerSpline.interpY(minS),'.g','MarkerSize',25);
            plot(obj.centerSpline.interpX(minS),obj.centerSpline.interpY(minS),'.c','MarkerSize',25);
            plot(obj.outerSpline.interpX(minS),obj.outerSpline.interpY(minS),'.r','MarkerSize',25);
            title(sprintf('Inner Radius = %f',obj.rInner));
                        
            % Axes properties
            axis equal;  
            xlabel('X Position');
            ylabel('Y Position');
            
            % Curvature plot
            subplot(2,1,2); hold on;
            plot(obj.cornerSVector,zeros(size(obj.cornerSVector)),'-k');
            h(1) = plot(obj.cornerSVector,obj.innerSpline.curvature(obj.cornerSVector),'.-r','LineWidth',2,'DisplayName','Inner Curvature');
            h(2) = plot(obj.cornerSVector,obj.centerSpline.curvature(obj.cornerSVector),'.-c','DisplayName','Centerline Curvature');
            h(3) = plot(obj.cornerSVector,obj.outerSpline.curvature(obj.cornerSVector),'.-g','DisplayName','Outer Curvature');
            h(4) = plot(obj.cornerSVector,obj.getCurvature(obj.cornerSVector),'--r','LineWidth',2,'DisplayName','Fitted 2nd Order Poly.');
            h(5) = plot(obj.cornerSVector,obj.trackdTheta,'-m','DisplayName','Raw Track Data');
            legend(h);
            
            % Return focus to the previous figure
            figure(prevFigFocus);
        end
        
        function [] = addPlotToFigure(obj,hFig)
           % Record the previous working figure
            prevFigFocus = gcf;
            
            % Plot the corner being analyzed
            figure(hFig);
            
            % Plot the fit curves
            plot(obj.innerSpline.interpX(obj.cornerSVector),obj.innerSpline.interpY(obj.cornerSVector),'r','LineWidth',3);
            plot(obj.centerSpline.interpX(obj.cornerSVector),obj.centerSpline.interpY(obj.cornerSVector),'c','LineWidth',3);
            plot(obj.outerSpline.interpX(obj.cornerSVector),obj.outerSpline.interpY(obj.cornerSVector),'g','LineWidth',3);
            
            % Plot the raw points
%             plot(obj.innerSpline.XRaw,obj.innerSpline.YRaw,'.-k','MarkerSize',20);
%             plot(obj.centerSpline.XRaw,obj.centerSpline.YRaw,'.-k','MarkerSize',20);
%             plot(obj.outerSpline.XRaw,obj.outerSpline.YRaw,'.-k','MarkerSize',20);

            % Plot location of corner from manual defintion
            plot([obj.innerSpline.interpX(obj.apexS),obj.centerSpline.interpX(obj.apexS),obj.outerSpline.interpX(obj.apexS)],[obj.innerSpline.interpY(obj.apexS),obj.centerSpline.interpY(obj.apexS),obj.outerSpline.interpY(obj.apexS)],'-k','LineWidth',2,'MarkerSize',30);
               
            % Plot the radii            
            cornerMargin = 0.3;
            radiS = obj.apexS;
            trackWidth = obj.getTrackWidth(radiS);
            perpAngle = obj.innerSpline.theta(radiS) + (pi/2);
%             plotPhiArc = (linspace(-obj.raceCornerTheta,obj.raceCornerTheta,100)./2) - obj.direction.*perpAngle;
            plotPhiArc = linspace(0,2.*pi,100);
            startInnerX = interp1([0,trackWidth],[obj.innerSpline.interpX(radiS),obj.outerSpline.interpX(radiS)],cornerMargin);
            startInnerY = interp1([0,trackWidth],[obj.innerSpline.interpY(radiS),obj.outerSpline.interpY(radiS)],cornerMargin);
             
            % Inner arc
            centerOfInner = [startInnerX + obj.rInner.*cos(perpAngle); startInnerY + obj.rInner.*sin(perpAngle)];
            innerCircPlot = [centerOfInner(1) + obj.rInner.*cos(plotPhiArc); centerOfInner(2) + obj.rInner.*sin(plotPhiArc)];
            plot([startInnerX,centerOfInner(1)],[startInnerY,centerOfInner(2)],'-r','MarkerSize',30);
            plot(innerCircPlot(1,:),innerCircPlot(2,:),'--r');
%             text(centerOfInner(1),centerOfInner(2),num2str(obj.cornerNum),'FontWeight','bold','HorizontalAlignment','center','Color','r','FontSize',14);
                
            % Race arc
%             centerOfRace = [startInnerX + obj.rRace.*cos(perpAngle); startInnerY + obj.rRace.*sin(perpAngle)];
%             raceCircPlot = [centerOfRace(1) + obj.rRace.*cos(plotPhiArc); centerOfRace(2) + obj.rRace.*sin(plotPhiArc)];
% %             plot(raceCircPlot(1,:),raceCircPlot(2,:),'-m','LineWidth',2);
        
            % Outer arc
            outerCircPlot = [centerOfInner(1) + obj.rOuter.*cos(plotPhiArc); centerOfInner(2) + obj.rOuter.*sin(plotPhiArc)];
            plot(outerCircPlot(1,:),outerCircPlot(2,:),'--g'); 
            
            % Label the corner
            text(obj.centerSpline.interpX(obj.apexS),obj.centerSpline.interpY(obj.apexS),num2str(obj.cornerNum),'FontWeight','bold','HorizontalAlignment','center','Color','k','FontSize',14);

            % Return focus to the previous figure
            figure(prevFigFocus); 
        end
        
        function [width] = getTrackWidth(obj,pts)
            width = norm(obj.outerSpline.interp(pts) - obj.innerSpline.interp(pts));
        end
        
        function [vals] = getCurvature(obj,pts)
            vals = polyval(obj.innerCurPoly,pts); 
        end
    end
end
