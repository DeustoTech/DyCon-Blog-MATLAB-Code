%% Movil Control
% In this tutorial, we show the numerical implementation of movil control
% strategy. For this, we consider the following model:
%%
% $$ 
% \begin{cases} 
% \dot{u}   = \Delta u + f\chi_{\omega(\textbf{d})} \\
% \dot{\textbf{d}} = \text{v} \\
% \dot{\textbf{v}} = \text{g} \\
% \end{cases}
% $$
%%
%
% where $u$

%%
% Create the mesh variables
clear;
Ns = 6;
Nt = 30;
xline = linspace(-1,1,Ns);
yline = linspace(-1,1,Ns);
[xms,yms] = meshgrid(xline,yline);

%%
% We create the B() function 
xwidth = 0.5;
ywidth = 0.5;
B = @(xms,yms,xs,ys) WinWP05(xms,xs,xwidth).*WinWP05(yms,ys,ywidth);
Bmatrix =  @(xs,ys) [diag(reshape(B(xms,yms,xs,ys),1,Ns^2)) ;zeros(Ns^2)];

%%
A  = FDLaplacial2D(xline,yline);

Atotal = zeros(2*Ns^2+4,2*Ns^2+4);
%
Atotal( 1:Ns^2  , 1:Ns^2 ) = A;
%
Atotal( Ns^2+1 : 2*Ns^2   ,   1    :  Ns^2   )  =  eye(Ns^2);
Atotal(    1   :  Ns^2    , Ns^2+1 : 2*Ns^2  )  =  50*eye(Ns^2); % z = 50*y

RumbaMatrixDynamics = [0 0 1 0; ...
                       0 0 0 1; ...
                       0 0 0 0; ...
                       0 0 0 0 ];

             
Atotal(2*Ns^2+1:end,2*Ns^2+1:end) = RumbaMatrixDynamics;
Atotal = sparse(Atotal);
%%
% gaussian function in 2D
gs = @(x,y,x0,y0,alpha) exp(-(x-x0).^2/alpha^2  - (y-y0).^2/alpha^2);

alpha = 0.1;
alphamid = 0.2;
% We build the initial condition

Y0 = + 2*gs(xms,yms,+0.25,+0.25,alpha) ...
     + 2*gs(xms,yms,+0.25,-0.25,alpha) ...
     + 2*gs(xms,yms,-0.25,-0.25,alpha);
%%
surf(xms,yms,Y0);
%%

Y0 = [Y0(:) ;Y0(:)*0];
%%
T = 0.5;
tspan = linspace(0,T,Nt+1);

%%
Ysym = sym('Y',[2*Ns^2+4 1]);
Usym = sym('U',[Ns^2+2 1]);
FDT = @(t,Y,U,Params) Atotal*Y+ [Bmatrix(Y(end-3),Y(end-2))*U(1:end-2) ;0; 0; U(end-1:end)];

idynammics = pde(FDT,Ysym,Usym);
idynammics.mesh = {xline,yline};
idynammics.InitialCondition = [Y0;0;0;0;0];
idynammics.Nt = Nt+1;
idynammics.FinalTime = T;
idynammics.Solver = @eulere;
%
[~ , Xnum_free] = solve(idynammics);

Xnum_free = Xnum_free';

%%
figure
ax1 = subplot(1,2,1); ax2 = subplot(1,2,2); 
isurf = surf(reshape(Xnum_free(1:Ns^2,1),Ns,Ns),'Parent',ax1);
jsurf = surf(reshape(Xnum_free(Ns^2+1:2*Ns^2,1),Ns,Ns),'Parent',ax2);

zlim([ax1,ax2],[0 0.3])

for it = 1:Nt
    isurf.ZData = reshape(Xnum_free(1:Ns^2,it),Ns,Ns);
    jsurf.ZData = reshape(Xnum_free(Ns^2+1:2*Ns^2,it),Ns,Ns);

    pause(0.05)
end
%%
opti = casadi.Opti();  % CasADi optimization structure

% ---- Input variables ---------
Ycas = opti.variable(2*Ns^2+4,Nt+1); % state trajectory
Ucas = opti.variable(Ns^2+2,Nt+1);   % control

% ---- Dynamic constraints --------
Fcas = @(y,u) Atotal*y+ [Bmatrix(y(end-3),y(end-2))*u(1:end-2) ;0;0; u(end-1:end)]; % dx/dt = f(x,u)

for k=1:Nt % loop over control intervals
   % Euler forward method
   x_next = Ycas(:,k) + (T/Nt)*Fcas(Ycas(:,k),Ucas(:,k)); 
   opti.subject_to(Ycas(:,k+1)==x_next); % close the gaps
end

% ---- State constraints --------
opti.subject_to(Ycas(:,1)==[Y0;0.2;0.2; 1 ;1]);

% HeatCas = Ucas(1:end-4,:);
% Max = 1e2;
% opti.subject_to(HeatCas(:) <= Max)
% opti.subject_to(HeatCas(:) >= -Max)
% opti.subject_to(HeatCas(:) <= 10)

% ---- Optimization objective  ----------
beta = 1e5;
Cost = (Ycas(1:end-4,Nt+1))'*(Ycas(1:end-4,Nt+1)) + beta*sum(sum((Ucas(1:end-2,:))'*(Ucas(1:end-2,:))));
%Cost = (Xcas(1:end-4,Nt+1))'*(Xcas(1:end-4,Nt+1));

opti.minimize(Cost); % minimizing L2 at the final time

% ---- initial guesses for solver ---
opti.set_initial(Ycas, Xnum_free);
opti.set_initial(Ucas, 0);


% ---- solve NLP              ------
p_opts = struct('expand',false);
s_opts = struct('acceptable_tol',1e-3,'constr_viol_tol',1e-3,'compl_inf_tol',1e-3);
opti.solver('ipopt',p_opts,s_opts); % set numerical backend
tic
sol = opti.solve();   % actual solve
toc

%%
Xnum_with_control = sol.value(Ycas);
%% interpolation in time

AnimationMovilControl(Xnum_with_control,Xnum_free,tspan,xline,yline,Bmatrix)
%%
