%% CALCULATE SPEED OF SOUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = speedsound(P,rho)
% Inputs: Pressure and Density

% speed of sound = sqrt(gamma*pressure/rho)
c = sqrt(1.4*P/rho);