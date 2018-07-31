%% CALCULATE PRIMITIVE STATE VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = convU_V(U)
% Inputs: Conservative state vector from a cell

% Parse out variables and apply more names
rho   = U(1);   % Density
rho_u = U(2);   % Density * u vel.
rho_v = U(3);   % Density * v vel.
rho_e = U(4);   % Density * total energy

% V = [rho u_vel v_vel pressure]
V(1) = rho;
V(2) = rho_u/rho;
V(3) = rho_v/rho;
V(4) = (rho_e - rho/2*(V(2)^2+V(3)^2))*(1.4-1);