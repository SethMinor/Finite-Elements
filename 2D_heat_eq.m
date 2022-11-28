%% COOK THE 2D TURKEY
clear, clc, clf;

% Font size, for plotting
fs = 14;

% 2D heat equation, with forcing
% PDE: u_t - (u_xx + u_yy) = f(x,y) (with Dirichlet BC's)

% WEAK FORM
% int(u_t*v)dV + int(grad(u)*grad(v))dV + int(BC's)dV = int(f*v)dV
% Insert u = sum(u_i(t)v_i(x,y)) and f = sum(f_i(t)v_i(x,y))
% Then discretize u_t

% STEADY STATE: KU=F
% DYNAMICS: MU'+KU=F

%https://www.wias-berlin.de/people/peschka/teaching/fem_present.pdf



%% Define the rectangular domain
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



%% Discretize 2D rectangular domain into triangular elements
% And choose a basis of tent functions for interior nodes

myelement =@(x,y,a,b) max(1 - abs(x-a) - abs(y-b),0);

[x,y] = meshgrid(1:15,1:15);
T = delaunay(x,y);
% [x,y] = meshgrid(0:dx:Lx,0:dy:Ly);
% T = delaunay(x,y);

z1 = myelement(x,y,10,10);
z2 = myelement(x,y,5,5);
% z1 = myelement(x,y,3*dx,4*dy);
% z2 = myelement(x,y,1*dx,6*dy);

figure (1)
trisurf(T,x,y,z1,'facecolor',[0 0.4470 0.7410])
alpha 0.75; hold on;
trisurf(T,x,y,z2,'facecolor',[0.8500 0.3250 0.0980])
alpha 0.75; hold on;
trisurf(T,x,y,0*z2,'facecolor','w')
alpha 0.75;
axis equal;
colormap jet;

title('Triangulation of domain','Interpreter','latex','FontSize',fs)
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)
legend('Example shape function','Example shape function',...
    'Interpreter','latex','FontSize',fs-4)




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

% Initial condition (t = 0)

% Time stepping (RK-4?)
t0 = 0;  % initial time
tf = 10; % final time

dt = 0.0015; % time discretization
Nt = ceil((tf-t0)/dt);
% Stability conidition is apparently dt <= 1/(2*((1/dx)^2+(1/dy)^2))

%Turk = questdlg('Cook the turkey?','God calling:','Oh yeah','Nay','Nay');



% DUMPATRON
% %      ----------------------       NODE NUMBERING
% %      |/\/\/\/\/\/\/\/\/\/\|       1    -> 2    -> ... -> Nx
% %      |/\/\/\/\/\/\/\/\/\/\|       Nx+1 -> Nx+2 -> ... -> 2*Nx
% %      |/\/\/\/\/\/\/\/\/\/\|                  ...
% %      ----------------------       ... ->  ...  -> ... -> Nx*Ny
% 
% % Triangular elements are defined by three boundary nodes, e.g.:
% %     1 -- 2                 2
% %     |  /                  /|
% %     | /         or       / |      etc.
% %     |/                  /  |
% %   (Nx+1)           (Nx+1)--(Nx+2)
% 
% nodes = zeros(Ny,Nx); % nodes for discretized domain
% 
% % Define the Cartesian coords of the nodes, e.g.:
% % (x,y)=(0dx,0dy)       (x,y)=(1dx,0dy)             (x,y)=(Lx,0dy)
% %     (1) ----------------- (2) ----- ...  ------------ (Nx)
% % (x,y)=(0dx,1dy)       (x,y)=(1dx,1dy)             (x,y)=(Lx,1dy)
% %    (Nx+1) ------------- (Nx+2) ----- ...  ---------- (2*Nx)
% %                           ...
% % (x,y)=(0dx,Ly)        (x,y)=(1dx,Ly)              (x,y)=(Lx,Ly)
% %    ([Nx-1]*Ny) ------ ([Nx-1]*Ny + 1) --- ...  ----- (Nx*Ny)
% 
% x = 0:dx:Lx;
% y = 0:dy:Ly;
% [X,Y] = meshgrid(x,y); % Cartesian coordinates at 'nodes'
% 
% % Define the Cartesian coords of the vertices for each element
% % (Listed in element order, 2*(Nx-1)*(Ny-1) total elements)
% 
% % Initialize (3 vertices) x (2*(Nx-1)*(Ny-1) element) array
% num_elements = 2*(Nx-1)*(Ny-1);
% vertices = zeros(num_elements,6);
% 
% % F
% counter = 0;
% for i = 1:(Nx-1)
%     for j = 1:(Ny-1)
%         counter = counter + 1;
%         fprintf('Vertex %.0f (up)\n',counter)
%         fprintf('(%.2f,%.2f), (%.2f,%.2f), (%.2f,%.2f)\n\n',...
%             X(j,i),Y(j,i),X(j,i+1),Y(j,i+1),X(j+1,i),Y(j+1,i))
%         vertices(counter,:) = [X(j,i),Y(j,i),X(j,i+1),Y(j,i+1),X(j+1,i),Y(j+1,i)];
%         
%         counter = counter + 1;
%         fprintf('Vertex %.0f (down)\n',counter)
%         fprintf('(%.2f,%.2f), (%.2f,%.2f), (%.2f,%.2f)\n\n',...
%             X(j+1,i),Y(j+1,i),X(j,i+1),Y(j,i+1),X(j+1,i+1),Y(j+1,i+1))
%         vertices(counter,:) = [X(j+1,i),Y(j+1,i),X(j,i+1),Y(j,i+1),X(j+1,i+1),Y(j+1,i+1)];
%     end
% end
% 
% figure (1)
% for j = 1:num_elements
%     plot([vertices(j,1);vertices(j,3);vertices(j,5);vertices(j,1)],...
%     [vertices(j,2);vertices(j,4);vertices(j,6);vertices(j,2)],'-k')
%     hold on
% end
% grid on
% xlim([-0.1,Lx+0.1])
% ylim([-0.1,Ly+0.1])
% title('Triangulation of domain','Interpreter','latex','FontSize',fs)
% xlabel('$x$','Interpreter','latex','FontSize',fs)
% ylabel('$y$','Interpreter','latex','FontSize',fs)
