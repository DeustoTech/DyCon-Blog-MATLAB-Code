% Stabilité exponentielle du système thermoelastique
%--------------------------------------------------------------------------
%Created by: Salem Nafiri (FSSM - Faculty of Sciences Semlalia Marrakesh)
%Problem: 1d Thermoelastic Problem, Spectral method
%Method: Modal method
%Version date:05/04/2013
% \partial_tt u = \partial_tt u - Gamma*\teta , 
% \partial_t    = Gamma*\partial_t u + k*\partial_xx \teta
% u|0,pi=0 && \teta|0,pi=0 (Dirichlet conditions)
% I.C: u0,v0=\partial u_t0, teta0
%--------------------------------------------------------------------------
%              Resolution of thermoelastic equation
%--------------------------------------------------------------------------


function [d]=stabexpsem(Nm)
%[e] = stabexp(Nm)
%This function computes the eigenvalues of the matrix An
%clc;
Gamma=0.1;
% Build bloc matrix An
for k=1:4
D=diag(1:k*Nm,0);
I=eye(k*Nm);
for i=1:k*Nm
    for j=1:k*Nm
        if mod(abs(i-j),2)==0 F(i,j)=0;
        else   F(i,j)=(-4/pi)*((i*j)/(i*i-j*j));
        end
    end
end
% Generate An (Dirichlet-Dirichlet BC)
%-------------------------------------
An=[zeros(k*Nm) D zeros(k*Nm);
     -D      zeros(k*Nm) -Gamma*F ;
    zeros(k*Nm)  Gamma*F'  -D^2];
%------------------------------------
% Get the eigenvalues of An
e = eig(An);
disp(['pour n=',num2str(k*Nm) ', d=' num2str(min(-real(e)),'%e')]);
%close all  % Closes all currently open figures.
%figure (k)
%plot real and imaginary parts
subplot(2,2,k)
plot(real(e),imag(e),'b.')
axis([-0.005 0 -10*k 10*k])
%axis square
hold on
%xlabel('Real')
%ylabel('Imaginary')
t1=['n=', num2str(k*Nm)];
title(t1)
end
end
 