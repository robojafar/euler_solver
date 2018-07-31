%% CALCULATE POSITIVE FLUX VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fpos(V,n)
% Inputs: Primitive state vector, V, from a cell
%         Normal vector, n, from a cell face

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

gamma = 1.4;                    % Gamma
a     = speedsound(P,rho);      % Speed of sound
%conV  = dot([u v],n);           % Contravariant velocity
conV  = u*n(1) + v*n(2);
M     = conV/a;                 % Contravarient Mach number

% Assemble flux vector
f(1)  = rho*a/4*(M+1)^2;
f(2) = f(1) * (u + n(1)*(-conV+2*a)/gamma);
f(3) = f(1) * (v + n(2)*(-conV+2*a)/gamma);
f(4) = f(1) * (h_0(V) - a^2*(M-1)^2/(gamma + 1));