function results = waves_2d(Nx,Ny,x0,y0,xi0,eta0,T1,T2,ht,ind)

%%
% Solution the the 2d wave equation 
% 
% rho(x)u_tt  - div(sigma(x)grad u) = 0,  (x,t) in [(-1,1)^2]x(0,T)
% u(-L,t) = u(L,t) = 0,                   t in (0,T)
% u(x,0) = u0(x), u_t(x,0) = u1(x)        x in (-1,1)^2
% 
% with a highly concentrated and oscillating initial data u0, using 
% finite differences on a uniform or non-uniform mesh and explicit Euler 
% for the time integration.
%
% Nx, Ny = number of points in the space mesh.
% L = size of the space interval.
% T1, T2 = size of the time interval.
% ind = choiche of the mesh: ind == 0 --> uniform mesh;
%                            ind == 1 --> refined mesh g(x) = L*tan((pi/(4*L)).*x);
%                            otherwise --> refined mesh g(x) = 2*L*sin((pi/(6*L)).*x);
%                                         
% (x0,y0,xi0,eta0) = postion and frequency of the initial datum.


%% Definiton of the mesh

hx = 1/(Nx+1); % uniform mesh size in the x direction
hy = 1/(Ny+1); % uniform mesh size in the y direction
x = -1:hx:1;   % mesh in the x direction
y = -1:hy:1;   % mesh in the x direction
lxx = length(x);
lyy = length(y); 
lx = lxx-2;
ly = lyy-2; 

% Choice of uniform/non-uniform mesh
if ind == 0
    g = @(x) x;
elseif ind == 1
    g = @(x) tan(0.25*pi*x);
elseif ind == 2
    g = @(x) 2*sin(pi*x/6);
end

gx = g(x);
gy = g(y);

%% Construction of the stiffness matrix and definition of the initial datum

dpx = 2*ones(1,lx)./((gx(3:lx+2)-gx(2:lx+1)).*(gx(2:lx+1)-gx(1:lx)));
dix = -2*ones(1,lx-1)./((gx(3:lx+1)-gx(2:lx)).*(gx(4:lx+2)-gx(2:lx)));
dsx = -2*ones(1,lx-1)./((gx(3:lx+1)-gx(2:lx)).*(gx(3:lx+1)-gx(1:lx-1)));
Ax = diag(dpx,0)+diag(dix,-1)+diag(dsx,1); 

dpy = 2*ones(1,ly)./((gy(3:ly+2)-gy(2:ly+1)).*(gy(2:ly+1)-gy(1:ly)));
diy = -2*ones(1,ly-1)./((gy(3:ly+1)-gy(2:ly)).*(gy(4:ly+2)-gy(2:ly)));
dsy = -2*ones(1,ly-1)./((gy(3:ly+1)-gy(2:ly)).*(gy(3:ly+1)-gy(1:ly-1)));
Ay = diag(dpy,0)+diag(diy,-1)+diag(dsy,1); 

[Vx,Dx] = eig(Ax);
[Vy,Dy] = eig(Ay');
dx = diag(Dx); 
dy = diag(Dy);
[DY,DX] = meshgrid(dy,dx); 
D = DY+DX;
[Y,X] = meshgrid(y,x);
[GY,GX] = meshgrid(gy,gx);
gamma = hx^(-0.9);

Data = Vx\(exp(-gamma*(X(2:lxx-1,2:lxx-1)-4*atan(x0)/pi).^2 ...
       -gamma*(Y(2:lyy-1,2:lyy-1)-4*atan(y0)/pi).^2).*...
       exp(1i*X(2:lxx-1,2:lxx-1)*xi0/hx+1i*Y(2:lyy-1,2:lyy-1)*eta0/hy))*Vy;

lambda = @(x,xi) (8/pi)*sin(0.5*xi)*(cos(0.25*pi*x))^2;
Lambda = @(x,y,xi,eta) lambda(x,xi)^2+lambda(y,eta)^2;
r0 = sqrt(Lambda(x0,y0,xi0,eta0));

raysx = @(t,x) [-4*sin(x(2))/(pi*r0*(x(1)^2+1));...
                -32*x(1)*(sin(0.5*x(2))^2)/(pi*r0*(x(1)^2+1)^2)];

raysy = @(t,y) [-4*sin(y(2))/(pi*r0*(y(1)^2+1));...
                -32*y(1)*(sin(0.5*y(2))^2)/(pi*r0*(y(1)^2+1)^2)];
            
in_datax = [x0,xi0];          
in_datay = [y0,eta0];          
tham = 0:ht:T1; 

[~,solhamx] = ode45(raysx,tham,in_datax);
[~,solhamy] = ode45(raysy,tham,in_datay);

t = 0:ht:T2;
lt = length(t);
Sol = zeros(lxx,lyy);

matr = zeros(lxx,lyy);
data = cell(0);

% figure('units','normalized','outerposition',[0 0 1 1])
% plot3(solhamx(:,1),solhamy(:,1),tham,'LineWidth',2)
% axis([-1 1 -1 1 0 T1])
% axis off
% view(0,90)
% print(gcf,'-dtiff','rays2D_top.tiff')
% % print('rays2D_top','-dpng')

%% SOlution of the problem

data2 = {};
for m = 1:lt-1
    Sol(2:lxx-1,2:lyy-1) = Vx*(Data.*exp(1i*t(m)*sqrt(D)))/Vy;
    data2{m} = abs(Sol);

    matr = max(Sol,matr);
    data{m} = matr;
    
end
%%

results.data2 = data2;
results.solhamx = solhamx;
results.solhamy = solhamy;
results.tham = tham;
results.lt = lt;
return

isurf = surf(abs(data2{1}));

shading interp
view(0,90)
zlim([0,3])
%caxis([0 1])
hold on 
plot3(solhamx(:,1),solhamy(:,1),tham,'w','LineWidth',3)


for m = 1:lt-1
    isurf.ZData = abs(data2{m});
    pause(0.1)
end

%% Plots
 



Pabs = abs(data{end});
Mabs = ceil(10*max(max(Pabs)))/10;
iind = num2str(ind);
string2 = strcat('waves2D_npv_mesh_',iind,'.png');

% Top view

figure('units','normalized','outerposition',[0 0 1 1])
surface(GX,GY,Pabs)
colormap jet
set(gca,'XTick',[-1,0,1],'XTickLabel',[-1,0,1],'fontsize',20,'fontweight','bold')
set(gca,'YTick',[-1,0,1],'YTickLabel',[-1,0,1],'fontsize',20,'fontweight','bold')
axis([-1 1 -1 1 0 T1])
caxis([0 1.2]);
box on
set(gca,'linewidth',2)
shading interp;
hold on
plot3(solhamx(:,1),solhamy(:,1),tham,'w','LineWidth',3)
view(0,90)
% print(string2,'-dpng')

% Side view

figure('units','normalized','outerposition',[0 0 1 1])
surface(GX,GY,Pabs)
colormap jet
set(gca,'XTick',[-1,0,1],'XTickLabel',[-1,0,1],'fontsize',20,'fontweight','bold')
set(gca,'YTick',[-1,0,1],'YTickLabel',[-1,0,1],'fontsize',20,'fontweight','bold')
set(gca,'ZTick',[0,Mabs],'ZTickLabel',[0,Mabs],'fontsize',20,'fontweight','bold')
axis([-1 1 -1 1 0 Mabs+0.1])
caxis([0 1.2]);
box on
set(gca,'linewidth',2)
shading interp;
view(45,65)
print(gcf,'-dtiff','waves2D_abs_side.tiff')
% % print('waves2D_abs_side','-dpng')
end