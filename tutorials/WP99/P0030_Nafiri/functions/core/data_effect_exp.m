%--------------------------------------------------------------------------
%Created by: Salem Nafiri (FSSM - Faculty of Sciences Semlalia Marrakesh)
%Problem: 1d Thermoelastic Problem,
%Method: Finite element method
%Version date: 06/07/2014
% \partial_tt u = \partial_xx u - \Gamma*\teta_x , 
% \partial_t\teta   = \partial_xx\teta +\Gamma*\partial_tx u
% u|0,L=0 && \teta|0,L=0 (Dirichlet/Dirichlet conditions)
% I.C: u0,v0=\partial u_t0, teta0
%--------------------------------------------------------------------------
%              decay of energy of beam thermoelastic equation
%--------------------------------------------------------------------------
% 
% parameters
%
%  n        : size of matrix
% Gamma: is a positive constant
%
function [u,v,theta,Et] = data_effect_exp(T,dt,n,k,Gamma)
%
%Data of the system
%
% col=['b','r','g'];
L=pi;
%
h=L./(n+1);  %Space step
x=h:h:n*h;   %Space variable
t=0:dt:T;    %time variable 
%Initial data
u0=zeros(1,n);
v0=100*ones(1,n);
%v0=sqrt(2/pi)*sin(k*x); % sin(x) or sin(j*pi*x/L) which is less regular
teta0=zeros(1,n);
u0=u0';
v0=v0';
teta0=teta0';
U0=[u0;v0;teta0];
%
% construction de la matrice D, A et I
%
% Build a vector of ones
e1=ones(n,1);
% Construct matrix Mn
Mn1=spdiags([-1/2*e1 e1 -1/2*e1],[-1 0 1],n,n);
Mn23=spdiags([1/4*e1 e1 1/4*e1],[-1 0 1],n,n);
Mn=[Mn1      zeros(n) zeros(n);
    zeros(n)  Mn23     zeros(n);
    zeros(n) zeros(n) Mn23];
% Build matrix Bn
Dn=1/h*spdiags([-sqrt(3)/2*e1 sqrt(3)*e1 -sqrt(3)/2*e1],[-1 0 1],n,n);
Fn=1/h*spdiags([3/4*e1 0.*e1 -3/4*e1],[-1 0 1],n,n);
Gn=1/h/h*spdiags([-3/2*e1 3*e1 -3/2*e1],[-1 0 1],n,n);
Bn=[zeros(n) Dn  zeros(n);
    -Dn zeros(n) -Gamma*Fn;
    zeros(n) Gamma*Fn' -Gn];
% Cholesky decomposition Mn=L'*L
L1=decltl(1,-1/2,n);
L2=decltl(1,1/4,n);
L3=L2;
% Generate An (D/D BC)
An=[zeros(n) inv(L1')*Dn*inv(L2)  zeros(n);
    -inv(L2')*Dn*inv(L1) zeros(n) -Gamma*inv(L2')*Fn*inv(L3);
    zeros(n) Gamma*inv(L3')*Fn'*inv(L2) -inv(L3')*Gn*inv(L3)];
A=full(An);
%------------------------------------
%
% Allocation de l'espace memoire de e
e=[];
%u=[];
%A=full(A);
E0=1/2*(norm(U0(1:n))^2 + norm(U0(n+1:2*n))^2 + norm(U0(2*n+1:3*n))^2);
%%
Ut = zeros(length(U0),length(t));
Et = zeros(1,length(t));

for i=1:length(t)
    U=expm(t(i)*A)*U0;
    E=1/2*(norm(U(1:n))^2 + norm(U(n+1:2*n))^2 + norm(U(2*n+1:3*n))^2);
    E=E/E0;
    Ut(:,i) = U;
    Et(i)   = E;
end
%%
dim = length(U0)/3;
u     = Ut(1         : dim    , :).';
v     = Ut(dim+1     : 2*dim  , :).';
theta = Ut(2*dim + 1 : end    , :).';
end