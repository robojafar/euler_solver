%% CALCULATE TOTAL ENTHALPY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = h_0(V)
% Inputs: Primitive state vector from a cell

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% Total enthalpy = total energy + pressure/rho
h = e_0(V)+P/rho;