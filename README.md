# Finite Volume Euler Solver
Read the [technical brief](https://github.com/robojafar/euler_solver/blob/master/euler_solver.pdf) for a full explantion of the mechanics behind this code.

## Overview
- Cell centered
- Van Leer flux splitting method
- m-stage time stepping
- Written in Matlab



## Euler Code
The code was divided into five major sections:
- Loading/Creating a grid
- Calculating grid metrics
- Setting freestream condition and loop variables
- Initializing the solution to freestream conditions
- Executing the main loop in time
  - Calculating the residual
  
### Loading the Grid
The code provides the option to create a Cartesian grid or load a grid in Plot3D format. When a grid is loaded, the number of zones is checked. If the number of zones is greater than one, then the grid will not be loaded because the code cannot handle multiple zones. If the grid passes this check, then the file is read and the x, y and z coordinates are stored in arrays. The code which completes this task is shown below:
```Matlab
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
```

### Calculating Grid Metrics
Once the grid is loaded, then the grid metrics can be calculated. Although it is feasible to calculate quantities during the iteration, it may be inefficient. Therefore, all necessary quantities will be calculated and stored in arrays before the iterations begin. Since this code is a cell-centered, finite volume solver, certain cell quantities will be needed such as the areas of each face (north, south, east, west), the outward normal vector of each face, and the volume of each cell.

The area of each face is calculated by taking the square root of the square of the x and y lengths of the face.

The volume of each cell is calculated using the formula for the volume of a parallelepiped (or the absolute value of the scalar triple product), which can be found in any calculus book.

The normal vector for each face was found by crossing the unit normal vector of the face by the unit vector in the z-direction. The result is an outward normal of each face.

The code which accomplishes these three tasks is shown below:

```Matlab
% Compute volume of cell using A.BxC (volume of parallelepiped)
volume(i,j) = abs(dot(z_width,cross(-1*[s_xlen,s_ylen,0],[e_xlen,e_ylen,0]))); 

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
```

### Setting Freestream Conditions and Loop Variables
Before the solution can be initialized, the freestream conditions need to be specified. The following freestream quantities are given:
- Density: 1.2 kg/m2
- Temperature: 300 K
- Pressure: 100,000 Pa
- Mach: 0.3, 0.7, 1.5 (depending on test case)

After the freestream conditions are defined, the speed of sound is calculated and used to find the u and v velocities. The primitive state vector, V, contains the primitive variables of the solution: 
- ρ is the density
- u is the velocity in the x-dir.
- v is the velocity in the y-dir.
- P is the pressure

Also, several looping variable are set. The number of iterations and the frequency of text output can are defined. The user can choose between local and global timestepping. If global timestepping is chosen, then a global timestep must be set. If local timestepping is chosen, then a CFL number must be set. The user also has the ability to plot contours of the primitive variables, which are updated during each iteration (or at a specified frequency).

This code employs the m-stage timestepping scheme, as outlined by Jameson. This allows the user to choose the order of the timestepping scheme. For example, if the order of the scheme is 1, then the forward Euler step is used. If the order of the scheme is 4, then a modified fourth-order Runge Kutta timestepping scheme is used. However, for non-linear equations, the modified fourth-order Runge Kutta scheme is only second-order accurate. The fourth-order scheme enables the user to use a higher CFL number and less iterations to achieve an accurate result.

### Initializing the Solution
Each cell in the grid is assigned the freestream primitive state vector, which is then converted to the conservative state vector via a function. Also, the residual of the cell is initialized to 0. The V, U, and resid variables are cell arrays. This is accomplished using the following code:

```Matlab
% Loop through all cells and init PSV to freestream conditions
% Convert PSV to conservative state vector
% Init residual to 0
for i = 1:nci
  for j = 1:ncj
    V{i,j} = V_free;
    U{i,j} = convV_U(V{i,j});
    resid{i,j} = [0 0 0 0];
  end
end
```

Please read the [technical brief](https://github.com/robojafar/euler_solver/blob/master/euler_solver.pdf) for further explanation of the code. GitHub markdown can only do so much.
