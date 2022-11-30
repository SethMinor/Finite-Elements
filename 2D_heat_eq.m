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

Nx = 10; % number of x nodes
Ny = 15; % number of y nodes

% Inline function for creating linear elements
myelement =@(x,y,a,b) max(1 - abs(x-a) - abs(y-b),0);

[x,y] = meshgrid(1:Nx,1:Ny);
T = delaunay(x,y);

% Plotting some example elements
z1 = myelement(x,y,6,5);
z2 = myelement(x,y,5,5);
z3 = myelement(x,y,1,4);
z4 = myelement(x,y,1,3);
z5 = myelement(x,y,Nx,Ny);

figure (1)
trisurf(T,x,y,z1+1.3*z2+z3+0.2*z4+z5,'facecolor',[0.8500 0.3250 0.0980])
alpha 0.75; hold on;
trisurf(T,x,y,0.001+0*z2,'facecolor','w')
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
G =@(x,y,xo,yo) 2*exp(-(0.15)*(x-xo).^2-(0.15)*(y-yo).^2);

% Center the Gaussian somewhere
xo = 6;
yo = 7;

% Fill a matrix with the discretized initial conditions
u0 = zeros(Ny,Nx);
for j = 1:Ny
    for i = 1:Nx
        u0(j,i) = G(i,j,xo,yo);
    end
end

% Plot the continuous version of the IC alongside u0
figure (2)
subplot(2,1,1)
[xg,yg] = meshgrid(1:0.1:Nx,1:0.1:Ny);
zg = 2*exp(-(0.15)*(xg-xo).^2-(0.15)*(yg-yo).^2);
surf(xg,yg,zg)
colorbar
shading interp;
colormap hot;
axis equal
title('$u_0 = u(x,y,0)$ (Continuous)','Interpreter','latex','FontSize',fs)
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)

subplot(2,1,2)
s = mesh(u0);
s.FaceColor = 'interp';
colorbar
shading interp;
colormap hot;
axis equal
title('$u_0 = u(x,y,0)$ (Discrete)','Interpreter','latex','FontSize',fs)
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
% (Elements M_ij: phi(i,j) * phi(I,J) is zero unless i=I and/or j=J)
fun1 = @(x,y) max(1 - abs(x-2) - abs(y-2),0).*max(1 - abs(x-2) - abs(y-2),0);
fun2 = @(x,y) max(1 - abs(x-1) - abs(y-2),0).*max(1 - abs(x-2) - abs(y-2),0);

ontop = integral2(fun1, 0,4, 0,4);
offset = integral2(fun2, 0,4, 0,4);
