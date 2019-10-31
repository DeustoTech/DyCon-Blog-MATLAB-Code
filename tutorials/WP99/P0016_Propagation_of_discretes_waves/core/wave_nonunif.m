function result = wave_nonunif(Nx,L,T,y0,xi0,rho,sigma,ind)

%%
% Solution the the 1d wave equation 
% 
% rho(x)u_tt  - (sigma(x)u_x)_x = 0,  (x,t) in (-L,L)x(0,T)
% u(-L,t) = u(L,t) = 0,                t in (0,T)
% u(x,0) = u0(x), u_t(x,0) = u1(x)     x in (-L,L)
% 
% with a highly concentrated and oscillating initial data u0, using 
% finite differences on a uniform or non-uniform mesh and explicit Euler 
% for the time integration.
%
% Nx = number of points in the space mesh.
% L = size of the space interval.
% T = size of the time interval.
% ind = choiche of the mesh: ind == 0 --> uniform mesh;
%                            ind == 1 --> refined mesh g(x) = L*tan((pi/(4*L)).*x);
%                            otherwise --> refined mesh g(x) = 2*L*sin((pi/(6*L)).*x);
%                                         
% (y0,xi0) = postion and frequency of the initial data.
% 


%% Definiton of the mesh

hx = (2*L)/(Nx+1);
mu = 0.1;
ht = mu*hx; % CFL condition
t = 0:ht:T;
lt = length(t);
mesh = mesh_nonunif(L,Nx,ind);
ly = length(mesh.yi);

%% Coefficients and stiffness matrix

rho = @(x) rho(x);
sigma = @(x) sigma(x);

A = matrix(mesh,sigma);

%% Initial datum 

gamma = hx^(-0.9);

if ind == 0
    data = exp(-0.5*gamma*(mesh.yi'-y0).^2).*exp(1i*xi0*mesh.yi'/hx);
else
    g = mesh.g;
    g = sym(g);
    ginv = finverse(g);
    ginv = matlabFunction(ginv);
    data = exp(-0.5*gamma*(ginv(L,mesh.yi')-ginv(L,y0)).^2) ...
                                         .*exp(1i*xi0*ginv(L,mesh.yi')/hx);
end

[V,D] = eigs(A,Nx);
d = diag(D);
four = zeros(1,length(d));

for j = 1:length(d)
    four(j) = data(1:length(d))*V(:,j);
end

u0 = zeros(ly,1);
u1 = zeros(ly,1);
for j = 1:length(d)
    u0 = u0 + four(j)*V(:,j);
    u1 = u1 + 1i*sqrt(d(j))*four(j)*V(:,j);
end

%% Solution of the equation 

Rho = diag(rho(mesh.yi));

sol = zeros(ly,lt);
sol(:,1) = u0;
sol(:,2) = u0 + ht*u1;

for j = 3:lt
    sol(:,j) = (2*eye(ly)-ht^2*(Rho\A))*sol(:,j-1)-sol(:,j-2);
end

sol = abs(sol);

%% Plots.
result.L = L;
result.T = T;
result.mesh = mesh;
result.t = t;
result.ind = ind;
result.sol = sol;

%
end

%% Auxiliary functions employed in the program.

function [mesh] = mesh_nonunif(L,N,ind)

g1 = @(x,L) x;
g2 = @(x,L) L*tan((pi/(4*L)).*x);
g3 = @(x,L) 2*L*sin((pi/(6*L)).*x);

if ind == 0
    g = g1;
    y = g(linspace(-L,L,N+2)');
elseif ind == 1
    g = g2;
    y = g(linspace(-L,L,N+2)',L);
else
    g = g3;
    y = g(linspace(-L,L,N+2)',L);
end

yimd(1:N) = (y(2:N+1) + y(1:N))/2 ;   % y_{i-1/2} points on the left of the point i
yipd(1:N) = (y(2:N+1) + y(3:N+2))/2;  % y_{i+1/2} points on the right of the point i

himd = y(2:end-1) - y(1:end-2);
hipd = y(3:end) - y(2:end-1);
hi = 0.5*(hipd+himd);
y = y(2:end-1);

mesh = struct('name','refined mesh',...
    'dim',1,...
    'yi',y,...
    'yipd',yipd',...
    'yimd',yimd',...
    'himd',himd,...
    'hipd',hipd,...
    'hi',hi,...
    'step',max(hi),...
    'N',N,...
    'g',g);
end

function A = matrix(mesh,sigma)

if mesh.dim==1 
    A = zeros(mesh.N,mesh.N);
    bloc = zeros(mesh.N,mesh.N);
    coeff = sigma;
    coeffipd = feval(coeff,mesh.yipd);
    coeffimd = feval(coeff,mesh.yimd);
    i = mesh.N;
    bloc(i,i-1) = -coeffimd(i)/(mesh.himd(i)*mesh.hi(i));
    bloc(i,i) = coeffipd(i)/(mesh.hipd(i)*mesh.hi(i)) + ...
                                     coeffimd(i)/(mesh.himd(i)*mesh.hi(i));
    
    i=1;
    bloc(i,i+1) = -coeffipd(i)/(mesh.hipd(i)*mesh.hi(i));
    bloc(i,i) = coeffipd(i)/(mesh.hipd(i)*mesh.hi(i)) + ...
                                     coeffimd(i)/(mesh.himd(i)*mesh.hi(i));
    for i = 2:mesh.N-1
        bloc(i,i+1) = -coeffipd(i)/(mesh.hipd(i)*mesh.hi(i));
        bloc(i,i-1) = -coeffimd(i)/(mesh.himd(i)*mesh.hi(i));
        bloc(i,i) = coeffipd(i)/(mesh.hipd(i)*mesh.hi(i)) +...
                                     coeffimd(i)/(mesh.himd(i)*mesh.hi(i));
    end    
    A(1:mesh.N,1:mesh.N) = bloc;
end
end