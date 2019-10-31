%--------------------------------------------------------------------------
%Created by: Salem Nafiri (FSSM - Faculty of Sciences Semlalia Marrakesh)
%Problem: 1d Thermoelastic Problem, SM
%Method: Spectral method
%Version date:06/07/2014
% \partial_tt u = \partial_tt u - Gamma*\teta , 
% \partial_t    = Gamma*\partial_t u + k*\partial_xx \teta
% u|0,pi=0 && \teta|0,pi=0 (Dirichlet/Dirichlet B.C)
% I.C: u0,v0=\partial u_t0, teta0
%--------------------------------------------------------------------------
%              Resolution of thermoelastic equation
%--------------------------------------------------------------------------

function stabpolysem(Nm)
%[e] = stabpoly(Nm)
%This function computes the eigenvalues of the matrix An

Gamma=0.1;
% Build bloc matrix An
for k=1:4
Dn=diag(1:k*Nm,0);
In=eye(k*Nm);
% Generate An
An=[zeros(k*Nm) Dn zeros(k*Nm);
      -Dn      zeros(k*Nm) -Gamma*In ;
     zeros(k*Nm)  Gamma*In  -Dn^2];
% An*An'
% An'*An
% Get the eigenvalues of A
e = eig(An);

%close all  % Closes all currently open figures.
%figure (k)
%plot real and imaginary parts
subplot(2,2,k)
plot(real(e),imag(e),'b.'), hold on

axis([-0.005 0 -20*k 20*k])
%axis square
%xlabel('Real')
%ylabel('Imaginary')
t1=['n=', num2str(k*Nm)];
title(t1)
end
end
 