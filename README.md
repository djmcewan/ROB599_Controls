## ROB599_Controls
ROB 599 Controls Project

## File description
#### Provided via Canvas
*checkTrajectory.p* - The function checkTrajectory will tell you if your car leaves the track, crashes into any obstacles, or exceeds any input limits. It takes in your trajectory (as an Nx2 vector, where the first column is the x coordinates of your trajectory and the second column is the y coordinates), your input (as an Nx2 vector where the first column is SWA and the second column is Fx), and a cell array of obstacles as an optional argument. The cell array of obstacles should be in the same format as the output of generateRandomObstacles. This function returns false if your car did not crash and did not exceed any input limits.

*forwardIntegrateControlInput.m* - The function forwardIntegrateControlInput is exactly the same forward integration method that we will use to check the validity of your control input and resulting trajectory. It is based on ode45, and assumes a zero order hold on your control inputs. It also assumes that your control inputs are defined every 0.01s. This function takes in a control input and initial condition and returns a trajectory in state space.

*generateRandomObstacles.m* - The function generateRandomObstacles is the same function that we will use to generate obstacles when checking your code for part 6.2. This function takes in an integer Nobs, the number of obstacles along the track, and returns Xobs, a cell array of size 1xNobs, where each cell contains a 4x2 matrix describing each obstacle (the first column is each obstacle’s x-coordinates and the second column is its y-coordinates). This function spaces the obstacles roughly evenly along the track. When we test your code, we will use 10-25 obstacles.

#### Team provided code
*ROB599Main.m* - Main script that calls all the sub-functions

*planRacePath.m* - Function to determine the optimal race path and velocity. Currently this is just the center of the track at a slow pace to get started.

*TrackCornerClass.m* - Used to define the corners or turns of the track

*SplineClass.m* - Dependency of TrackCornerClass.m
