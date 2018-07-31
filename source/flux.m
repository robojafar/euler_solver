%% CALCULATE FULL FLUX VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = flux(V,n)
% Inputs: Primitive state vector, V, from a cell
%         Normal vector, n, from a cell face

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

%conV  = dot([u v],n);           % Contravariant velocity
conV  = u*n(1) + v*n(2);

% Assemble flux vector
f(1) = rho*conV;
f(2) = rho*u*conV + P*n(1);
f(3) = rho*v*conV + P*n(2);
f(4) = rho*h_0(V)*conV;