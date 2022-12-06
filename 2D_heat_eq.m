%% COOK THE 2D TURKEY
% 2D heat equation, with forcing
clear, clc, clf;

% Font size, for plotting
fs = 14;


%% Define the rectangular domain
% Discretize 2D rectangular domain into triangular elements
% and choose a basis of tent functions for interior nodes
%               Lx
%      ----------------------
%      |                    |
%   Ly |                    |
%      |                    |
%      ----------------------

Nx = 12; % number of x nodes
Ny = 12; % number of y nodes

% Inline function for creating linear elements
myelement =@(x,y,a,b) max(1 - abs(x-a) - abs(y-b),0);

[x,y] = meshgrid(1:Nx,1:Ny);
T = delaunay(x,y);

%Plotting some example elements
z1 = myelement(x,y,6,5);
z2 = myelement(x,y,10,5);

figure (1)
trisurf(T, x, y, z1 + z2,'facecolor',[0.8500 0.3250 0.0980])
alpha 0.75; hold on;
trisurf(T,x,y,0.01+0*z2,'facecolor','w')
alpha 0.75;
axis equal; colormap jet;

title('Triangulation of Domain','Interpreter','latex','FontSize',fs)
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)


%% Dirichlet boundary conditions
%               B2
%      ----------------------
%      |                    |
%   B1 |                    | B4
%      |                    |
%      ----------------------
%               B3

B1 = 0;
B2 = 0;
B3 = 0;
B4 = 0;


%% Set up initial condition u0 = u(x,y,0)
% Speficy the values of u0 on the grid (xi,yj) = (i,j)

% Initial Gaussian
Amplitude = 10;
Damping = 0.01;
G =@(x,y,xo,yo) Amplitude*exp(-Damping*(x-xo).^2-Damping*(y-yo).^2);

% Center the Gaussian somewhere
xo = Nx/2;
yo = Ny/2;

% Fill a matrix with the discretized initial conditions
u0 = zeros(Ny,Nx);
for j = 1:Ny
    for i = 1:Nx
        u0(j,i) = G(i,j,xo,yo);
    end
end

% Plot the continuous version of the IC alongside u0
figure (2)

% Continuous
subplot(2,1,1)
[X,Y] = meshgrid(linspace(1,Nx,5*Nx), linspace(1,Ny,5*Ny));
Z = G(X,Y,xo,yo);
surf(X,Y,Z)
colorbar
shading interp;
colormap hot;
axis equal
title('$u_0 = u(x,y,0)$ (Continuous)','Interpreter','latex','FontSize',fs)
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)

% Piece-wise linear approximation
subplot(2,1,2)
s = mesh(u0);
s.FaceColor = 'interp';
colorbar
shading interp;
colormap hot;
axis equal
title('$u_0 = u(x,y,0)$ (Piecewise linear)','Interpreter','latex','FontSize',fs)
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)


%% Time stepping algorithm
% FORWARD EULER: Stability conidition, dt <= 1/(2*((1/dx)^2+(1/dy)^2))
t0 = 0;  % initial time
tf = 10; % final time
dt = 0.001; % time discretization
Nt = ceil((tf-t0)/dt);

%Turk = questdlg('Cook the turkey?','God calling:','Oh yeah','Nay','Nay');

% Fill up the M matrix
% (Elements M_ij: phi(i)*phi(j) is zero unless adjacent/identical)
Mii = @(x,y) max(1-abs(x)-abs(y), 0).*max(1-abs(x)-abs(y), 0);
Mij = @(x,y) max(1-abs(x-1)-abs(y), 0).*max(1-abs(x)-abs(y), 0);

% Integral for identical tent functions. int(phi_i^2)
identical = integral2(Mii, -2,2, -2,2);

% Integral for adjacent tent functions, int(phi_i*phi_j)
offset = integral2(Mij, -2,2, -2,2);
