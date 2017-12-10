classdef (HandleCompatible = true) TrackClass < matlab.mixin.SetGet    
    properties(Constant, Hidden)
        MAIN_FIGURE_NUM = 100;
    end
    
    properties
        % Provided track properties
        bl;
        br;
        cline;
        theta;
        
        % Derived track properties
        arc_s;
        dtheta;

        % Inline functions
        center;
        ftheta;
        fdtheta;
        fbl;
        fbr;

        % Random obstacle properties
        nRandObs = 40;
        randObs;
       
        % Track corner properties
        cornerObjs;
        nCorners;
        
        % Initial state
        X0;
        profileSPts;
        profileXYPts;
        profileThetaPts;
        
        % Plotting properties
        hMainFig;
        hMainAxes;
        hWithoutObsLine;
        hWithObsLine;
        hWorkingLine;
    end
    
    properties(Access = private)
        % Polynomial fitting properties
        trackPolySpacing = 7;
        trackPolyLocations;
        trackPolyLocationLims;
        trackPolyObjs;
    end
        
	methods
		%% Constructor
		function [obj] = TrackClass(trackFilename,X0)
            % Set the initial state
            obj.X0 = X0;
            
            % Load track Information
            RawTrackData = load(trackFilename);
            
            % Set the provided values
            obj.bl = RawTrackData.TestTrack.bl;
            obj.br = RawTrackData.TestTrack.br;
            obj.cline = RawTrackData.TestTrack.cline;
            obj.theta = RawTrackData.TestTrack.theta;

            % Compute the arclength
            obj.arc_s = cumsum([0,vecnorm(diff(obj.cline,[],2))]);
            obj.dtheta = [0,wrapToPi(diff(obj.theta))./diff(obj.arc_s)];

            % Define the inline functions
            obj.center =    @(s)[interp1(obj.arc_s,obj.cline(1,:),s);interp1(obj.arc_s,obj.cline(2,:),s)];
            obj.ftheta =    @(s)interp1(obj.arc_s,obj.theta,s);
            obj.fdtheta =   @(s)interp1(obj.arc_s,obj.dtheta,s);
            obj.fbl =       @(s)[interp1(obj.arc_s,obj.bl(1,:),s);interp1(obj.arc_s,obj.bl(2,:),s)];
            obj.fbr =       @(s)[interp1(obj.arc_s,obj.br(1,:),s);interp1(obj.arc_s,obj.br(2,:),s)];
        
            % Generate the random obstacles
            obj.randObs = generateRandomObstacles(obj.nRandObs,obj);
            
            % Decompose the track into segments
            obj.trackPolyLocations = (obj.arc_s(1):obj.trackPolySpacing:obj.arc_s(end))';
            obj.trackPolyLocationLims = [1,length(obj.trackPolyLocations)];
            warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
            obj.trackPolyObjs = arrayfun(@(x) TrackPolyFit(obj,x,obj.trackPolySpacing),obj.trackPolyLocations);
            warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
            
            % Define the track corners
            obj.defineTrackCorners();
        end
        
        function [seg] = getTrackSegment(obj,s)
            seg = obj.trackPolyObjs(obj.trackPolyIdx(s));
        end
        
        function [s] = cartesian2Track(obj,state)
            % Extract components from the state vector
            x = state(1);
            y = state(3);

            % Do a rough cut
            roughDistLimit = 30;
            roughDist = sqrt((obj.cline(1,:)-x).^2+(obj.cline(2,:)-y).^2); 
            firstIdx = find(roughDist <= roughDistLimit,1,'First');
            lastIdx = find(roughDist <= roughDistLimit,1,'Last');
            sLims = obj.arc_s([firstIdx,lastIdx]);

            % Get the polynomial the describes the track
            trackPoly = obj.getTrackSegment(mean(obj.arc_s(firstIdx:lastIdx)));

            % Compute the distance polynomials
            [minLeftDist,   minLeftS]   = findMinimumLineDistance(trackPoly.leftXCoeffs,    trackPoly.leftYCoeffs,  sLims);    
            [minCenterDist, minCenterS] = findMinimumLineDistance(trackPoly.centerXCoeffs,  trackPoly.centerYCoeffs,sLims);    
            [minRightDist,  minRightS]  = findMinimumLineDistance(trackPoly.rightXCoeffs,   trackPoly.rightYCoeffs, sLims);    

            % Determine the minimum s value
            minSVals = [minLeftS,minCenterS,minRightS];
            [~, minDistIdx] = min([minLeftDist,minCenterDist,minRightDist]);
            s = minSVals(minDistIdx);

            function [minDist, minS] = findMinimumLineDistance(xPolyCoeffs,yPolyCoeffs,sLims)
                distXPoly = xPolyCoeffs + [0,0,-x];
                distYPoly = yPolyCoeffs + [0,0,-y];
                sqDistPoly = TrackPolyFit.square2OrderPoly(distXPoly) + TrackPolyFit.square2OrderPoly(distYPoly);
                sqDistDerPoly = polyder(sqDistPoly);
                distRoots =  roots(sqDistDerPoly);
                [~,realRootIdx] = min(abs(imag(distRoots)));

                % Check the endpoints to ensure minimum
                possibleSVals = [sLims,real(distRoots(realRootIdx))];
                [~,sIdx] = min(polyval(sqDistPoly,possibleSVals));
                minS = possibleSVals(sIdx);
                minDist = sqrt(polyval(sqDistPoly,minS));
            end
        end
        
        function [] = plotTrackBase(obj)
            obj.hMainFig = figure(obj.MAIN_FIGURE_NUM); clf; hold on; box on; axis equal;
            set(obj.hMainFig,'NumberTitle','Off','Name','COTA Track');
            plot([obj.bl(1,1),obj.br(1,1)],[obj.bl(2,1),obj.br(2,1)],'g-','LineWidth',1);
            plot([obj.bl(1,end),obj.br(1,end)],[obj.bl(2,end),obj.br(2,end)],'r-','LineWidth',1);
            plot(obj.bl(1,:),obj.bl(2,:),'k-','LineWidth',1);
            plot(obj.br(1,:),obj.br(2,:),'k-','LineWidth',1);
            plot(obj.cline(1,:),obj.cline(2,:),'k--','LineWidth',1);
            obj.hWithoutObsLine = plot(obj.X0(1).*[1,1],obj.X0(3).*[1,1],'-b','LineWidth',2);
            obj.hWithObsLine    = plot(obj.X0(1).*[1,1],obj.X0(3).*[1,1],'-g','LineWidth',2);
            obj.hWorkingLine    = plot(obj.X0(1).*[1,1],obj.X0(3).*[1,1],'-b','LineWidth',2);
            obj.hMainAxes = gca;
        end
        
        function [] = plotObstacles(obj)
            % Draw the obstacles
            for iObs = 1:obj.nRandObs
                fill([obj.randObs{iObs}(:,1);obj.randObs{iObs}(1,1)],[obj.randObs{iObs}(:,2);obj.randObs{iObs}(1,2)],'r');
            end
        end
        
        function [] = plotCorners(obj,plotRadii)            
            % Add to the figures
            for iCorner = 1:obj.nCorners
                obj.cornerObjs(iCorner).addPlotToFigure(obj.hMainFig,plotRadii);
            end
        end
        
        function [] = plotPolynomials(obj)
            % Draw the track fits
            for iPoly = 1:length(obj.trackPolyLocations)
                plotSRange = (obj.trackPolySpacing/2).*[-1,1] + obj.trackPolyLocations(iPoly);
                plotSPts = plotSRange(1):0.025:plotSRange(2);
                trackSeg = obj.getTrackSegment(obj.trackPolyLocations(iPoly));
                plot(polyval(trackSeg.leftXCoeffs,plotSPts),polyval(trackSeg.leftYCoeffs,plotSPts),'c','LineWidth',2);
                plot(polyval(trackSeg.rightXCoeffs,plotSPts),polyval(trackSeg.rightYCoeffs,plotSPts),'c','LineWidth',2);
            end
        end
        
        function [] = plotCenterline(obj)
            plot(obj.profileXYPts(:,1),obj.profileXYPts(:,2),'m','LineWidth',3);
        end
        
        function [] = zoomInOnCorner(obj,cornerNumber)
            limits = [obj.cornerObjs(cornerNumber).outerEntryXY; obj.cornerObjs(cornerNumber).outerExitXY];
            displayLims = (max(diff(limits)+30)/2).*[-1,1];
            axis([displayLims + mean(limits(:,1)),displayLims + mean(limits(:,2))]); 
        end
        
        function [dist,gradDist] = trackNegDistanceTraveled(obj,D,nStates,nlcObj)
            % Get the first and last state
            beginState = D(1:nStates);
            finalState = D((end-nStates+1):end);
            finalS = obj.cartesian2Track(finalState);
            finalSeg = obj.getTrackSegment(finalS);

            % Compute the negative distance because we want to minimize
            dist = obj.cartesian2Track(beginState) - finalS;

            % Compute the gradients
            gradDist = zeros(size(D));
            gradDist((end-nStates+1):end) = -nlcObj.partialSdX(finalState,finalSeg.getPathAngle(finalS));
        end
    end
    
    methods(Hidden)
        function [polyIdx] = trackPolyIdx(obj,s)
            polyIdx = min([max([round((s/obj.trackPolySpacing) + 1),obj.trackPolyLocationLims(1)]),obj.trackPolyLocationLims(2)]);
        end
        
        function [] = defineTrackCorners(obj)
            % Define the corners manually and sort them before numbering
            obj.cornerObjs = [TrackCornerClass(obj,140,220,300,TrackCornerClass.RIGHT);...    % 1
                              TrackCornerClass(obj,450,480,510,TrackCornerClass.LEFT);...     % 2
                              TrackCornerClass(obj,545,560,605,TrackCornerClass.RIGHT);...    % 3
                              TrackCornerClass(obj,625,680,705,TrackCornerClass.LEFT);...     % 4
                              TrackCornerClass(obj,735,761,828,TrackCornerClass.RIGHT);...    % 5
                              TrackCornerClass(obj,845,849,925,TrackCornerClass.RIGHT);...    % 6
                              TrackCornerClass(obj,1015,1017,1055,TrackCornerClass.LEFT);...  % 7
                              TrackCornerClass(obj,1070,1100,1120,TrackCornerClass.LEFT);...  % 8
                              TrackCornerClass(obj,1150,1206,1267,TrackCornerClass.RIGHT);... % 9
                              TrackCornerClass(obj,1268,1274,1330,TrackCornerClass.LEFT);...  % 10
                              TrackCornerClass(obj,1440,1491,1530,TrackCornerClass.LEFT);...  % 11
                             ];
            % Sort the corners based on their position
            obj.nCorners = length(obj.cornerObjs);
            arrayfun(@(cornerObj,num) cornerObj.setCornerNumber(num),obj.cornerObjs,(1:obj.nCorners)');
            
            % Get the smooth centerline of the track
            [obj.profileSPts,obj.profileXYPts] = obj.defineCenterline();
            obj.profileThetaPts = [atan2(diff(obj.profileXYPts(:,2)),diff(obj.profileXYPts(:,1)));atan2(diff(obj.profileXYPts((end-1:end),2)),diff(obj.profileXYPts((end-1:end),1)))];            
        end
        
        function [profileSPts,profileXYPts] = defineCenterline(obj)
            % Connect all the turns together
            profileSSpacing = 0.5;
            profileSPts = 0;
            profileXYPts = obj.X0([1,3]);
            for iCorner = 1:obj.nCorners
                % Add the line segment portion
                prevS = profileSPts(end);
                nextS = obj.cornerObjs(iCorner).startS;
                prevX = profileXYPts(end,1);
                prevY = profileXYPts(end,2);
                nextX = obj.cornerObjs(iCorner).centerEntryXY(1);
                nextY = obj.cornerObjs(iCorner).centerEntryXY(2);

                % Add the line segment if necessary
                if((nextS - prevS) > 1)
                    lineSegPts = (prevS + profileSSpacing):profileSSpacing:(nextS - profileSSpacing);
                    profileSPts = [profileSPts;lineSegPts']; %#ok<AGROW>
                    lineSegXVals = interp1([prevS,nextS],[prevX,nextX],lineSegPts)';
                    lineSegYVals = interp1([prevS,nextS],[prevY,nextY],lineSegPts)';
                    profileXYPts = [profileXYPts;[lineSegXVals,lineSegYVals]]; %#ok<AGROW>
                end

                % Add the curve portion
                cornerS = obj.cornerObjs(iCorner).cornerSVector';
                profileSPts = [profileSPts;cornerS]; %#ok<AGROW>
                profileXYPts = [profileXYPts;obj.cornerObjs(iCorner).centerSpline.interp(cornerS)]; %#ok<AGROW>    
            end

            % Add the end point
            profileSPts = [profileSPts;obj.arc_s(end)];
            profileXYPts = [profileXYPts;obj.cline(:,end)'];   
        end
    end
end