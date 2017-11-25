%% Analysis of track
% Define the corners manually and sort them before numbering
offsetS = -400;
cornerObjs = [TrackCornerClass(Track,100,220,320,TrackCornerClass.RIGHT);...    % 1
              TrackCornerClass(Track,430,480,520,TrackCornerClass.LEFT);...     % 2
              TrackCornerClass(Track,520,560,610,TrackCornerClass.RIGHT);...    % 3
              TrackCornerClass(Track,610,680,710,TrackCornerClass.LEFT);...     % 4
              TrackCornerClass(Track,710,761,828,TrackCornerClass.RIGHT);...    % 5
              TrackCornerClass(Track,831,849,950,TrackCornerClass.RIGHT);...    % 6 
              TrackCornerClass(Track,1010,1017,1055,TrackCornerClass.LEFT);...  % 7
              TrackCornerClass(Track,1070,1100,1120,TrackCornerClass.LEFT);...  % 8
              TrackCornerClass(Track,1140,1206,1265,TrackCornerClass.RIGHT);... % 9
              TrackCornerClass(Track,1265,1274,1330,TrackCornerClass.LEFT);...  % 10
              TrackCornerClass(Track,1440,1491,1530,TrackCornerClass.LEFT);...  % 11
];
[cornerSs,sortIdx] = sort([cornerObjs(:).apexS]);
cornerObjs = cornerObjs(sortIdx);
nCorners = length(cornerObjs);
arrayfun(@(cornerObj,num) cornerObj.setCornerNumber(num),cornerObjs,(1:nCorners)');

% Add to the figures
pointSpread = 5;
cornerMargin = 1;
cornerPlotColors = {'g','r'};
cornerPlotsOn = true;
for iCorner = 1:nCorners
    if(cornerPlotsOn)
        cornerObjs(iCorner).addPlotToFigure(FIG_TRACK_WITHOUT_OBS);
        cornerObjs(iCorner).addPlotToFigure(FIG_TRACK_WITH_OBS);
    end
end

% Connect all the turns together
profileSSpacing = 0.5;
profileSPts = 0;
profileXYPts = X0([1,3]);
for iCorner = 1:length(cornerObjs)
    % Add the line segment portion
    prevS = profileSPts(end);
    nextS = cornerObjs(iCorner).startS;
    prevX = profileXYPts(end,1);
    prevY = profileXYPts(end,2);
    nextX = cornerObjs(iCorner).centerEntryXY(1);
    nextY = cornerObjs(iCorner).centerEntryXY(2);
        
    % Add the line segment if necessary
    if((nextS - prevS) > 1)
        lineSegPts = (prevS + profileSSpacing):profileSSpacing:(nextS - profileSSpacing);
        profileSPts = [profileSPts;lineSegPts']; %#ok<AGROW>
        lineSegXVals = interp1([prevS,nextS],[prevX,nextX],lineSegPts)';
        lineSegYVals = interp1([prevS,nextS],[prevY,nextY],lineSegPts)';
        profileXYPts = [profileXYPts;[lineSegXVals,lineSegYVals]]; %#ok<AGROW>
    end
    
    % Add the curve portion
    cornerS = cornerObjs(iCorner).cornerSVector';
    profileSPts = [profileSPts;cornerS]; %#ok<AGROW>
    profileXYPts = [profileXYPts;cornerObjs(iCorner).centerSpline.interp(cornerS)]; %#ok<AGROW>    
end

% Add the end point
profileSPts = [profileSPts;Track.arc_s(end)];
profileXYPts = [profileXYPts;Track.cline(:,end)'];   
    
% Dummy velocity profile
maxVel = 80;
minVel = 5;
velWayPoints = [...
    0,                          minVel;...
    cornerObjs(1).stopS+0.1,    10;...
    cornerObjs(2).startS,       10;...
    cornerObjs(2).stopS,        15;...
    cornerObjs(3).stopS,        15;...
    cornerObjs(8).startS-30,    minVel;...
    cornerObjs(8).stopS,        minVel;...
    cornerObjs(11).stopS,       15;...
    profileSPts(end),           maxVel
];
profileVPts = interp1(velWayPoints(:,1),velWayPoints(:,2),profileSPts);

%% Velocity profile with turn information
figure(200); clf; hold on;
set(figure(200),'Name','Velocity Profile','NumberTitle','Off');
plot(profileSPts,profileVPts,'b','LineWidth',2);
for iCorner = 1:length(cornerObjs)
    hC(1) = plot(cornerObjs(iCorner).startS.*[1,1],80.*[0,1],'g','DisplayName','Start of Corner');     
    hC(2) = plot(cornerObjs(iCorner).stopS.*[1,1],80.*[0,1],'r','DisplayName','End of Corner');
    text(cornerObjs(iCorner).apexS,40,sprintf('%d',cornerObjs(iCorner).cornerNum),'HorizontalAlignment','center','VerticalAlignment','middle');
    zoom xon;    
end
legend(hC);
xlabel('Centerline Distance [m]');
ylabel('Velocity [m/s]');

% Ensure unique grid vector
[~,uniqueIdx] = unique(profileSPts);
profileSPts = profileSPts(uniqueIdx);
profileXYPts = profileXYPts(uniqueIdx,:);
profileVPts = profileVPts(uniqueIdx);

% Make interpolation functions
fProfileXYPts = @(s) interp1(profileSPts,profileXYPts,s,'linear','extrap');
fProfileVPts = @(s) interp1(profileSPts,profileVPts,s,'linear','extrap');
