%% ROB 599 Controls Project: Formula 1 car racing main script
% Clear the workspace
clear; clc; close all;

% Set the update time of used in the forwardIntegrateControlInput function
dT = 0.01;      

% Set the initial conditions
X0 = [287,5,-176,0,2,0];

% Load track Information
track = TrackClass('TestTrack.mat',X0);

% Load the reference trajectory to get started
ref = ReferenceTrajectory(track);

%% Draw the track with additional contect
track.plotTrackBase();
track.plotObstacles();
track.plotCorners(false);
track.plotPolynomials();
track.plotCenterline();

% Look at individual corners if necessary
% track.zoomInOnCorner(4);

%% Perform the optimization
% Define the spacing parameters                
lookAheadTime = 20; % [steps]
incrementTime = 10;  % [steps]
nStates = size(X0,2);
nInputs = 2;
nDecPerStep = (nStates+nInputs);
nDecTotal = (nStates + lookAheadTime*nDecPerStep);

% Objective is to maximize the distance travelled down the track, because fmincon minimizes we add a negative sign
objectiveFun = @(d) track.trackNegDistanceTraveled(d,nStates);

% Define the lower bounds of the decision vector
luPosMargin = 5;    % [m]
maxYawRate = 100;   % [rad/s]
lb = -Inf(nDecTotal,1);
lb(1:nDecPerStep:end) = min([track.bl(1,:),track.br(1,:)]) - luPosMargin;
lb(3:nDecPerStep:end) = min([track.bl(2,:),track.br(2,:)]) - luPosMargin;
lb(6:nDecPerStep:end) = -maxYawRate;
lb(7:nDecPerStep:end) = -0.5;
lb(8:nDecPerStep:end) = -10000;

% Define the upper bounds of the decision vector
ub = Inf(nDecTotal,1);
ub(1:nDecPerStep:end) = max([track.bl(1,:),track.br(1,:)]) + luPosMargin;
ub(3:nDecPerStep:end) = max([track.bl(2,:),track.br(2,:)]) + luPosMargin;
ub(6:nDecPerStep:end) = maxYawRate;
ub(7:nDecPerStep:end) = 0.5;
ub(8:nDecPerStep:end) = 5000;

% Define the non-linear constraint function
nlcObj = NonLinearConstraintsClass(track,nStates,nInputs,nDecPerStep);

% Set optimization options
opts = optimoptions('fmincon',...
                    'SpecifyConstraintGradient',true,...
                    'SpecifyObjectiveGradient',false,...
                    'CheckGradients',true,...
                    'UseParallel',false);

% Seed the initial vector
D = NaN(nDecTotal,1);
D(1:nStates,:) = ref.xRef(1,:);
for iStep = 1:lookAheadTime
    D((1:nDecPerStep) + nDecPerStep*(iStep-1) + nStates,1) = [ref.uRef(iStep,:),ref.xRef(iStep+1,:)]';
end           
                
% Create a waitbar
hWaitbar = waitbar(0,'Test');

currentS = track.cartesian2Track(X0);
% while()
for i = 1
    waitbar(currentS/track.arc_s(end));  
    optvals = fmincon(objectiveFun,D,[],[],[],[],lb,ub,@nlcObj.constraintFcn,opts);    
end

figure(1); clf; hold on;
plot([D,optvals]);

% Close the waitbar
close(hWaitbar);

%%
% - Use gradient checker that is part of fmincon
% - Use small discretation
% - Be very careful
% - Do trajectory design over a section of the track and then shift the horizon
