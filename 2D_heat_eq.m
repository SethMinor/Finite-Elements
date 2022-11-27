% COOK THE 2D TURKEY
clear, clc, clf;

% 2D heat equation, with forcing
% PDE: u_t - (u_xx + u_yy) = f(x,y) (with Dirichlet BC's)

% WEAK FORM (assuming no flux)
% int(u_t*v)dV + int(grad(u)*grad(v))dV = int(f*v)dV
% Then discretize u_t


% Nice examples?:
% https://www.mathworks.com/matlabcentral/fileexchange/77545-2d-heat-transfer-fem



% Define the rectangular domain
%               Lx
%      ----------------------
%      |                    |
%   Ly |                    |
%      |                    |
%      ----------------------

Nx = 10; % number of x nodes
Ny = 20; % number of y nodes

Lx = 1; % length of x side
Ly = 2; % length of y side

dx = Lx/(Nx-1); % x-axis discretization
dy = Ly/(Ny-1); % y-axis discretization



% Discretize 2D rectangular domain into triangular elements
%      ----------------------       NODE NUMBERING
%      |/\/\/\/\/\/\/\/\/\/\|       1    -> 2    -> ... -> Nx
%      |/\/\/\/\/\/\/\/\/\/\|       Nx+1 -> Nx+2 -> ... -> 2*Nx
%      |/\/\/\/\/\/\/\/\/\/\|                  ...
%      ----------------------       ... ->  ...  -> ... -> Nx*Ny

% Triangular elements are defined by three boundary nodes, e.g.:
%     1 -- 2                 2
%     |  /                  /|
%     | /         or       / |      etc.
%     |/                  /  |
%   (Nx+1)           (Nx+1)--(Nx+2)

nodes = zeros(Ny,Nx); % nodes for discretized domain

% Define the Cartesian coords of the nodes, e.g.:
% (x,y)=(0dx,0dy)       (x,y)=(1dx,0dy)             (x,y)=(Lx,0dy)
%     (1) ----------------- (2) ----- ...  ------------ (Nx)
% (x,y)=(0dx,1dy)       (x,y)=(1dx,1dy)             (x,y)=(Lx,1dy)
%    (Nx+1) ------------- (Nx+2) ----- ...  ---------- (2*Nx)
%                           ...
% (x,y)=(0dx,Ly)        (x,y)=(1dx,Ly)              (x,y)=(Lx,Ly)
%    ([Nx-1]*Ny) ------ ([Nx-1]*Ny + 1) --- ...  ----- (Nx*Ny)

x = 0:dx:Lx;
y = 0:dy:Ly;
[X,Y] = meshgrid(x,y); % Cartesian coordinates at 'nodes'

% Define the Cartesian coords of the vertices for each element
% (Listed in element order, 2*(Nx-1)*(Ny-1) total elements)
% vertices = [[X(1,1),Y(1,1)]; [X(1,2),Y(1,1)]; [X(1,1),Y(2,1)]];
% figure (1)
% plot([vertices(:,1);vertices(1,1)],[vertices(:,2);vertices(1,2)])
% grid on

% Initialize (3 vertices) x (2*(Nx-1)*(Ny-1) element) array
num_elements = 2*(Nx-1)*(Ny-1);
vertices = zeros(num_elements,3);
for i = 1:Nx
    for j = 1:Ny
        % Even-numbered elements face one way
        if mod(j,2) == 0
            vertices(j,:) = [[X(1,1),Y(1,1)]; [X(1,2),Y(1,1)]; [X(1,1),Y(2,1)]];
    
        % Odd-numbered elements face another way
        else
            vertices(j,:) = [[0, 0]; [0, 0]; [0, 0]];
        end
    end
end




%DT = delaunay(x,y);
%triplot(DT,x,y);



% Dirichlet boundary conditions
%               B2
%      ----------------------
%      |                    |
%   B1 |                    | B4
%      |                    |
%      ----------------------
%               B3

B1 = 0;
B2 = 3;
B3 = 0;
B4 = 1;

% Initial condition (t = 0)

% Time stepping (RK-4?)
t0 = 0;  % initial time
tf = 10; % final time

dt = 0.0015; % time discretization
Nt = ceil((tf-t0)/dt);
% Stability conidition is apparently dt <= 1/(2*((1/dx)^2+(1/dy)^2))

%Turk = questdlg('Cook the turkey?','God calling:','Oh yeah','Nay','Nay');
