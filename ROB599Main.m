%% ROB 599 Controls Project: Formula 1 car racing main script
% Clear the workspace
clear; clc; close all;

% Load track Information
load('TestTrack.mat');
Track = TestTrack; clear TestTrack;

% Compute the arclength
Track.arc_s = cumsum([0,vecnorm(diff(Track.cline,[],2))]);
Track.dtheta = [0,wrapToPi(diff(Track.theta))./diff(Track.arc_s)];

% Define the inline functions
Track.center = @(s)[interp1(Track.arc_s,Track.cline(1,:),s);interp1(Track.arc_s,Track.cline(2,:),s)];
Track.ftheta = @(s)interp1(Track.arc_s,Track.theta,s);
Track.fdtheta = @(s)interp1(Track.arc_s,Track.dtheta,s);

% Generate the random obstacles
nRandObs = 40;
randObs = generateRandomObstacles(nRandObs,Track);

%% Set the initial conditions
X0 = [287,5,-176,0,2,0];

%% Draw the track with and without obstacles
FIG_TRACK_WITHOUT_OBS = 100;
FIG_TRACK_WITH_OBS = 101;

for iFig = [FIG_TRACK_WITHOUT_OBS,FIG_TRACK_WITH_OBS]
    figure(iFig); clf; hold on; box on; axis equal;
    plot(Track.bl(1,:),Track.bl(2,:),'k-','LineWidth',1);
    plot(Track.br(1,:),Track.br(2,:),'k-','LineWidth',1);
    plot(Track.cline(1,:),Track.cline(2,:),'k--','LineWidth',1);
    hAnimatedLine = animatedline(X0(1),X0(3));

    switch(iFig)
        case FIG_TRACK_WITHOUT_OBS
            set(iFig,'NumberTitle','Off','Name','Track w/o Obstacles');
            hWithObsAnimatedLine = hAnimatedLine;

        case FIG_TRACK_WITH_OBS
            set(iFig,'NumberTitle','Off','Name','Track w/Obstacles');
            hWithoutObsAnimatedLine = hAnimatedLine;

            % Draw the obstacles
            for iObs = 1:nRandObs
                fill([randObs{iObs}(:,1);randObs{iObs}(1,1)],[randObs{iObs}(:,2);randObs{iObs}(1,2)],'r');
            end
    end
    clear hAnimatedLine
end

% Propose a path and velocity profile
planRacePath();    

% Draw the proposed path on the track figures
for iFig = [FIG_TRACK_WITHOUT_OBS,FIG_TRACK_WITH_OBS]
    figure(iFig);
    plot(profileXYPts(:,1),profileXYPts(:,2),'m','LineWidth',3);
end

%% Generate the necessary control inputs to follow proposed path
% forwardIntegrateControlInput([[0,0];[0,0]],X0)
% checkTrajectory(profileXYPts,zeros(size(profileXYPts)))
% checkTrajectory(profileXYPts,zeros(size(profileXYPts)),randObs)
