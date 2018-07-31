% Jafar Mohammed
% AE516 Final Project

% Euler Solver
% First Order Van Leer Scheme
% M-stage timestepping

clear                               % Clear variables from memory
clc                                 % Clear the command window

%% LOAD/GENERATE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadGrid = 1;            %0 - create cartesian grid
                         %1 - load grid from file
tic
if loadGrid == 1
    % Open the file and assign it an ID
    %fid = fopen('BUMP/bump10perc.grd', 'r');
    fid = fopen('BUMP/bump04perc.grd', 'r');
    %fid = fopen('bump.grd', 'r');
    
    if fid >= 1
        % Read in file headers
        zones = fscanf(fid, '%d', 1);
    
        % Code only handles 1 zone
        % Therefore, check for number of zones
        if (zones == 1)
            % Read in number of i,j,k points
            npi = fscanf(fid, '%d', 1);
            npj = fscanf(fid, '%d', 1);
            npk = fscanf(fid, '%d', 1);
    
            % Retrieve i,j,k coordinates
            x = fscanf(fid, '%f', [npi,npj]);
            y = fscanf(fid, '%f', [npi,npj]);
            z = fscanf(fid, '%f', [npi,npj]);
            disp('Grid read successfully');
        end
        fclose(fid);
    end
else
    % Generate cartesian grid using the meshgrid function
    npi   = 21;                          % Num. of pts. in i-dir
    npj   = 11;                          % Num. of pts. in j-dir
    xc    = linspace(0,20,npi);          % Distribution in i-dir
    yc    = linspace(0,10,npj);          % Distribution in j-dir
    [y,x] = meshgrid(yc,xc);             % Creates 2D grid
end

% Visualizes the mesh in a plot window
%mesh(x,y,zeros(npi,npj))
%axis image;view(2);drawnow;

disp('Grid generated');
toc;disp(' ');

%% GRID METRICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
z_width = [0,0,1];                   % Unit vector in z-dir
nci = npi-1;                         % Number of cells (pts-1) in i dir
ncj = npj-1;                         % Number of cells (pts-1) in j dir

for i = 1:nci
    for j = 1:ncj
        % Assemble the face lengths
        e_xlen = x(i+1,j+1)-x(i+1,j);
        e_ylen = y(i+1,j+1)-y(i+1,j);
        
        n_xlen = x(i,j+1)-x(i+1,j+1);
        n_ylen = y(i,j+1)-y(i+1,j+1);

        w_xlen = x(i,j)-x(i,j+1);
        w_ylen = y(i,j)-y(i,j+1);

        s_xlen = x(i+1,j)-x(i,j);
        s_ylen = y(i+1,j)-y(i,j);
        
        % Compute midpoint of cell (for plotting)
        xmid(i,j) = (x(i,j) + x(i+1,j))/2;
        ymid(i,j) = (y(i,j) + y(i,j+1))/2;
        
        % Compute volume of cell using A.BxC (volume of parallelepiped)
        volume(i,j) = abs(dot(z_width,cross(-1*[s_xlen,s_ylen,0],...
                      [e_xlen,e_ylen,0])));
        
        % Compute area of cell
        sE(i,j) = sqrt((e_xlen)^2 + (e_ylen)^2);
        sN(i,j) = sqrt((n_xlen)^2 + (n_ylen)^2);
        sW(i,j) = sqrt((w_xlen)^2 + (w_ylen)^2);
        sS(i,j) = sqrt((s_xlen)^2 + (s_ylen)^2);
        
        % Compute outward normal of faces (return 3 component vector)
        temp_nE = cross([e_xlen,e_ylen, 0]/sE(i,j), z_width);
        temp_nN = cross([n_xlen,n_ylen, 0]/sN(i,j), z_width);
        temp_nW = cross([w_xlen,w_ylen, 0]/sW(i,j), z_width);
        temp_nS = cross([s_xlen,s_ylen, 0]/sS(i,j), z_width);
        
        % Truncate normal vector to 2 components
        nE{i,j} = [temp_nE(1) temp_nE(2)];
        nN{i,j} = [temp_nN(1) temp_nN(2)];
        nW{i,j} = [temp_nW(1) temp_nW(2)];
        nS{i,j} = [temp_nS(1) temp_nS(2)];
        
        % Clear unecessary variables
        clear temp_nE temp_nN temp_nW temp_nS
        clear e_xlen n_xlen w_xlen s_xlen
        clear e_ylen n_ylen w_ylen s_ylen
    end 
end
disp('Grid Metrics calculated');
toc;disp(' ');

%% SET FREESTREAM CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
free_rho      = 1.2;                            % kg/m^3
free_T        = 300.00;                         % Kelvin
free_P        = 100000;                         % Pa
free_M        = 1.5;                            % Mach
free_aoa      = 0;                              % deg

free_a        = speedsound(free_P,free_rho);    % m/s
free_u        = free_M*free_a*cosd(free_aoa);   % m/s
free_v        = free_M*free_a*sind(free_aoa);   % m/s
free_vel      = [free_u free_v];

% Freestream Primitive State Vector (PSV)
V_free        = [free_rho free_u free_v free_P];

diary('output.txt')
fprintf('Freestream Conditions:\n');
fprintf('Mach:         %10.2f\n',free_M);
fprintf('Flow AOA:     %10.2f deg\n',free_aoa);
fprintf('u Velocity:   %10.2f m/s\n',free_u);
fprintf('v Velocity:   %10.2f m/s\n',free_v);
fprintf('Pressure:     %10.2e Pa\n',free_P);
fprintf('Temperature:  %10.2f K\n',free_T);
fprintf('Density:      %10.2f kg/m^3\n\n',free_rho);

%% SET ITERATION VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration variables
iterations   = 100;             % Number of iterations to run

gbl_timestep = 0;               % 0 = local timestepping
                                % 1 = global timestepping
                                
timestep     = 1e-5;            % Timestep for global timestepping
CFL          = 0.5;             % Courant number

m_stage      = 4;               % m-stage time stepping
                                % e.g. 1 for Euler step
                                %      4 for 4th-order RK

% Output variables
freq         = 10;              % Reporting frequency for output

plotcontours = 0;               % Plot contours during iterations
                                % 0 - off
                                % 1 - Density
                                % 2 - U Velocity
                                % 3 - V Velocity
                                % 4 - Pressure

plots_on     = 1 ;              % Plot PSV contour plots after 
                                % completion
                                
fprintf('Iteration Variables:\n');
fprintf('Iterations:   %5d\n',iterations);
fprintf('CFL:          %5.2f\n',CFL);
fprintf('M-Stage:      %5d\n',m_stage);
if gbl_timestep == 0
    fprintf('Timestepping  %5s\n\n','local');
else
    fprintf('Timestepping %5s\n\n','global');
end

%% INITIALIZE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
resid_i      = 0;               % Iterative residual
resid_0      = 0;               % Step 0 residual
start_iter   = 0;               % Used for multiple runs
end_iter     = 0;               % Used for multiple runs
residReduced = 0;               % If divergence detected
fm = 1;                         % Used for capturing movies

% Combine normals and areas into big cell array which will be passed
% to the function which computes the residual
normals = {nE nN nW nS};
areas   = {sE sN sW sS};

% Initalize variables which will allow for visualization
% i.e. Plot Contours
con_density = zeros([nci,ncj]);
con_uvel    = zeros([nci,ncj]);
con_vvel    = zeros([nci,ncj]);
con_pres    = zeros([nci,ncj]);

% Loop through all cells and init PSV to freestream conditions
% Convert PSV to conservative state vector
% Init residual to 0
for i = 1:nci
    for j = 1:ncj
        V{i,j}     = V_free;
        U{i,j}     = convV_U(V{i,j});
        resid{i,j} = [0 0 0 0];
    end
end

% Add a disturbance at cell (10,5) if cartesian grid created
if loadGrid == 0
    V{10,5}(4) = 1.1*free_P;
    U{10,5}    = convV_U(V{10,5});
end

diary off
disp('Solution initialized');
toc

%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_iter = start_iter + 1;        % Start at iteration 1
diary on

% Create plot window
if plotcontours == 1
    hc = figure('name','Density Contour'); 
elseif plotcontours == 2
    hc = figure('name','U Velocity Contour'); 
elseif plotcontours == 3
    hc = figure('name','V Velocity Contour'); 
elseif plotcontours == 4
    hc = figure('name','Pressure Contour'); 
end

tic

% Main loop in time
for iter = start_iter:(end_iter + iterations)
    % Time variable used to measure time/iteration
    ti1 = cputime;
    
    % Initialize iteration residual to 0
    resid_i = 0;
    
    % Save CSV from this timestep (to be used in m-stage)
    U0 = U;
    
    % M-stage timestepping scheme
    for m = 1:m_stage
        % Calculate residual using function calcResid
        % Passes PSV, normals, areas cell array,
        %        freestream PSV, and nci and ncj
        resid = calcResid(V, V_free, normals, areas, nci, ncj);
        
        % Loop through all cells to update solution
        for i = 1:nci
            for j = 1:ncj
                
                % If local timestepping, calculate timestep of the cell
                if gbl_timestep == 0
                    vel = [V{i,j}(2) V{i,j}(3)];
                    cell_a = speedsound(V{i,j}(4),V{i,j}(1));
                    dt(1) = CFL * sE(i,j)/(abs(vel(1)*nE{i,j}(1) +...
                            vel(2)*nE{i,j}(2))+cell_a);
                    dt(2) = CFL * sN(i,j)/(abs(vel(1)*nN{i,j}(1) +...
                            vel(2)*nN{i,j}(2))+cell_a);
                    dt(3) = CFL * sW(i,j)/(abs(vel(1)*nW{i,j}(1) +...
                            vel(2)*nW{i,j}(2))+cell_a);
                    dt(4) = CFL * sS(i,j)/(abs(vel(1)*nS{i,j}(1) +...
                            vel(2)*nS{i,j}(2))+cell_a);
                    timestep = min(dt);
                end
                
                % Update solution using the saved CSV
                % Multiply by 'alpha' constant 
                U{i,j} = U0{i,j} -...
                         1/(m_stage-(m-1))*timestep/volume(i,j)*resid{i,j};
                
                % Update cell PSV
                V{i,j} = convU_V(U{i,j});
            
                % Update contour arrays used for plotting
                con_density(i,j) = V{i,j}(1);
                con_uvel(i,j)    = V{i,j}(2);
                con_vvel(i,j)    = V{i,j}(3);
                con_pres(i,j)    = V{i,j}(4);

                % Assemble first part of L2 norm for residual
                resid_i = resid_i + resid{i,j}.^2;
            end
        end
    end   
    
    % Assemble second part of L2 norm for residual
    resid_i = (resid_i).^.5/(nci*ncj);
    
    % Assign normalization value in first interation
    if iter == 1
        resid_0 = resid_i;
        
        % More detailed iurput
        fprintf('\nIter  cont. resid  x-mom resid  y-mom resid  energy resid  Time left\n');
    end
    
    % Detects extreme divergence (at the point of no return)
    % and shuts down simulation
    if isnan(resid_i/resid_0)
        break;
        disp('Solution corrupt.');
    end
    
    % Detects divergence happening in x-mom resid and cuts CFL in half
	if ((resid_i(2)/resid_0(2)) >= (1e1))
        if residReduced == 0
            CFL = CFL/2; 
            notice = sprintf('Divergence detected.  CFL reduced to %5.2f',CFL);
            disp(notice);
            residReduced = residReduced + 1;
        %elseif ((residReduced > 0) &&((resid_i(2)/resid_0(2)) >= (2e1)))
        %    CFL = CFL/2; 
        %    notice = sprintf('Divergence detected.  CFL reduced to %5.2f',CFL);
        %    disp(notice);
        end
    end
    
    % Computes time/iteration
    ti2 = cputime-ti1;
    
    % Displays output and updates plots at user specified time
    % interval
    if mod(iter,freq) == 0
        % Plot contours if wanted
        if plotcontours == 1
            % Assembles contour plot
            [C,h] = contourf(xmid',ymid',con_density');
            
            % Removes lines from contour plot
            set(h, 'LineStyle','none');
            
            % Adds a title to the plot
            s = sprintf('Density Contour at %4d iterations',iter);
            title(s)
            
            % Adds gridlines, corrects aspect ratio
            grid on;axis image;drawnow;
            
            % Assembles array for movie viewing
            mov(fm) = getframe(gca);
            fm=fm+1;
        elseif plotcontours == 2
            [C,h] = contourf(xmid',ymid',con_uvel');
            set(h, 'LineStyle','none');
            s = sprintf('U-Vel. Contour at %4d iterations',iter);
            title(s)
            colorbar('peer',gca,'SouthOutside'); 
            grid on;axis image;drawnow;
            mov(fm) = getframe(gca);
            fm=fm+1;
        elseif plotcontours == 3
            [C,h] = contourf(xmid',ymid',con_vvel');
            set(h, 'LineStyle','none');
            s = sprintf('V-Vel. Contour at %4d iterations',iter);
            title(s)
            colorbar('peer',gca,'SouthOutside'); 
            grid on;axis image;drawnow;
            mov(fm) = getframe(gca);
            fm=fm+1;
        elseif plotcontours == 4
            [C,h] = contourf(xmid',ymid',con_pres');
            set(h, 'LineStyle','none');
            s = sprintf('Pressure Contour at %4d iterations',iter);
            title(s)
            colorbar('peer',gca,'SouthOutside'); 
            grid on;axis image;drawnow;
            mov(fm) = getframe(gcf);
            fm=fm+1;
        end
            
        % More detailed output
        fprintf('%4d  %11.2e  %11.2e  %11.2e   %11.2e  %7s\n',iter,...
                resid_i(1)/resid_0(1),...
                resid_i(2)/resid_0(2),...
                resid_i(3)/resid_0(3),...
                resid_i(4)/resid_0(4),...
                fixTime(ti2*(end_iter+iterations-iter)));
    end
end

start_iter = iter;
end_iter = iter;
toc
diary off

%% PLOT PRIMITIVE STATE VECTOR CONTOURS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a figure with 4 contour subplots
% Plot 1: Density (kg/m^3) contour
% Plot 2: U Velocity (m/s) contour
% Plot 3: V Velocity (m/s) contour
% Plot 4: Pressure (Pa)

if plots_on == 1
    figure('name','Primitive State Variables');
    
    subplot(221);
    contourf(xmid',ymid',con_density');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('Density (kg/m^3)')

    subplot(222);
    contourf(xmid',ymid',con_uvel');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('U Velocity (m/s)')

    subplot(223);
    contourf(xmid',ymid',con_vvel');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('V Velocity (m/s)')

    subplot(224);
    contourf(xmid',ymid',con_pres');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('Pressure (Pa)')
end

% Plots movie in new figure window
%figure;
%movie(mov);
%movie2avi(mov,'movie.avi','fps',5)