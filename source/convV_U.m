%% CALCULATE CONSERVATIVE STATE VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function U = convV_U(V)
% Inputs: Primitive state vector from a cell

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% U = [rho rho*u_vel rho*v_vel rho*energy]
U(1) = rho;
U(2) = rho*u;
U(3) = rho*v;
U(4) = rho*e_0(V);