%% CALCULATE RESIDUAL OF GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resid = calcResid(V, V_free, normals, areas, nci, ncj)
% Inputs: V       - Cell array containing the primitive state vector for
%                   each cell in the grid
%         V_free  - Primitive state vector containing freestream conditions
%         normals - Cell array containing all of the face normals for each
%                   cell
%         areas   - Cell array containing all of the face areas for each
%                   cell
%         nci,ncj - Number of cell in i and j

% Parse out normal vectors in smaller, specialized cell arrays
nE = normals{1};
nN = normals{2};
nW = normals{3};
nS = normals{4};

% Parse out areas vectors in smaller, specialized cell arrays
sE = areas{1};
sN = areas{2};
sW = areas{3};
sS = areas{4};

% Initialize cell array containing the flux for each face
e_flux = cell(nci,ncj);
w_flux = cell(nci,ncj);
n_flux = cell(nci,ncj);
s_flux = cell(nci,ncj);

% Initialize residual cell array
resid = cell(nci,ncj);

% Assemble freestream velocity vector
free_vel = [V_free(2) V_free(3)];

% Calculate freestream Mach number
free_a   = speedsound(V_free(4), V_free(1));

% Loop through all cells
for i = 1:nci
    for j = 1:ncj
        % Calculate speed of sound for cell
        cell_a = speedsound(V{i,j}(4),V{i,j}(1));
        
        % Assemble velocity vector of cell
        vel    = [V{i,j}(2) V{i,j}(3)];
        
        % Contravariant Mach for face
        %conM = dot(vel,nE{i,j})/cell_a;
        conM = (vel(1)*nE{i,j}(1) + vel(2)*nE{i,j}(2))/cell_a;
        
        % If east face is at the outflow boundary
        if i == nci
            % Boundary normal equal to face normal
            nB = nE{i,j};
            
            % Check to see is contravariant M >= 1
            if conM >= 1
                % Assign boundary flux to nci cell flux
                Vb = V{i,j};
                
                % Calculate east face flux
                e_flux{i,j} = flux(Vb,nB);
            else
                % Calculate positive & negative Riemann invariants
                rpos = (vel.*nB) + 2/(1.4-1)*cell_a;
                rneg = (free_vel.*nB) - 2/(1.4-1)*free_a;

                % Average invariants to get normal velocities
                un = (rpos+rneg)/2;

                % Obtain velocities at boundaries
                ub = vel + (un - vel).*nB;                 

                % Vb is primitive state vector at outflow boundary
                % Density is extrapolated
                % u,v are calculated
                % Pressure is extrapolated
                Vb = [V{i,j}(1) ub(1) ub(2) V{i,j}(4)];

                % Calculate East face flux
                e_flux{i,j} = flux(Vb,nB);  
            end
        % If east face is internal
        else
            % Check to see is contravariant M >= 1
            if conM >= 1
                % Calculate east face flux using full flux
                e_flux{i,j} = flux(V{i,j},nE{i,j});
            else
                % Calculate east face flux using f+ and f-
                e_flux{i,j} = fpos(V{i,j},nE{i,j}) + fneg(V{i+1,j},nE{i,j});
            end
        end
        
        
        % Contravariant Mach for face
        %conM = dot(vel,nN{i,j})/cell_a;
        conM = (vel(1)*nN{i,j}(1) + vel(2)*nN{i,j}(2))/cell_a;
        
        % If north face is invicid wall
        if j == ncj
            % Boundary normal equal to face normal
            nB = nN{i,j};
            
            % Contravariant velocity:
            % Dot product of the u,v velocity of cell(i,ncj)
            % and the wall normal vector
            conV = dot(vel, nB);
            
            % Wall velocity:
            % Velocity of cell(i,ncj) - contravariant vel. * wall normal
            % vector
            velB = vel - conV*nB;
            
            % Vb is primitive state vector at the wall
            % Density is extrapolated
            % u,v are calculated
            % Pressure is extrapolated
            Vb = [V{i,j}(1) velB(1) velB(2) V{i,j}(4)];
            
            % Calculate north face flux
            n_flux{i,j} = flux(Vb,nB);
        else
            % Check to see is contravariant M >= 1
            if conM >=1
                % Calculate north face flux using full flux
                n_flux{i,j} = flux(V{i,j},nN{i,j});
            else
                % Calculate north face flux using f+ and f-
                n_flux{i,j} = fpos(V{i,j},nN{i,j}) + fneg(V{i,j+1},nN{i,j});
            end
        end
        
        
        % Contravariant Mach for face
        %conM = dot(vel,nW{i,j})/cell_a;
        conM = (vel(1)*nW{i,j}(1) + vel(2)*nW{i,j}(2))/cell_a;
        
        % If west face is inflow boundary
        if i == 1
            % Boundary normal equal to face normal
            nB = nW{i,j}; 
            
            if conM >= 1
                % Assign boundary state vec to freestream vals
                Vb = V_free;

                % Calc flux for west face
                w_flux{i,j} = flux(Vb,nB);
            else
                % Calculate positive & negative Riemann invariant
                rpos = (vel.*nB) + 2/(1.4-1)*cell_a;
                rneg = (free_vel.*nB) - 2/(1.4-1)*free_a;

                % Average invariants to get normal velocities
                un = (rpos+rneg)/2;

                % Obtain velocities at boundaries
                ub = free_vel + (un + free_vel).*nB;

                % Vb is primitive state vector at inflow boundary
                % Density is freestream
                % u,v are calculated
                % Pressure is freesteam
                Vb = [V_free(1) ub(1) ub(2) V_free(4)];

                % Calculate west face flux
                w_flux{i,j} = flux(Vb,nB);
            end
        else
            % Check to see is contravariant M >= 1
            if conM >=1
                % Calculate west face flux using full flux
                w_flux{i,j} = flux(V{i,j},nW{i,j});
            else
                % Calculate west face flux using the negative
                % of the adjacent cell east face flux
                w_flux{i,j} = -1*e_flux{i-1,j};
            end
        end
        
        
        % Contravariant Mach for face
        %conM = dot(vel,nS{i,j})/cell_a;
        conM = (vel(1)*nS{i,j}(1) + vel(2)*nS{i,j}(2))/cell_a;
        
        % If south face is invicid wall
        if j == 1
            % Boundary normal equal to face normal
            nB = nS{i,j};
            
            % Contravariant velocity:
            % Dot product of the u,v velocity of cell(i,1)
            % and the wall normal vector
            conV = dot(vel, nB);
            
            % Wall velocity:
            % Velocity of cell(i,1) - contravariant vel. * wall normal
            % vector
            velB = vel - conV*nB;
            
            % Vb is primitive state vector at the wall
            % Density is extrapolated
            % u,v are calculated
            % Pressure is extrapolated
            Vb = [V{i,j}(1) velB(1) velB(2) V{i,j}(4)];
            
            % Calculate south face flux
            s_flux{i,j} = flux(Vb,nB);
        else
            % Check to see is contravariant M >= 1
            if conM >= 1
                % Calculate south face flux using full flux
                s_flux{i,j} = flux(V{i,j},nS{i,j});
            else
                % Calculate south face flux using the negative
                % of the adjacent cell north face flux
                s_flux{i,j} = -1*n_flux{i,j-1};                            
            end
        end
        
        % Assemble residual:
        % SUM(face flux * face area)
        resid{i,j} = e_flux{i,j}*sE(i,j) + ...
                     n_flux{i,j}*sN(i,j) + ...
                     w_flux{i,j}*sW(i,j) + ...
                     s_flux{i,j}*sS(i,j);
    end
end