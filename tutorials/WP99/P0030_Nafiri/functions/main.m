
mode = 1;
dt = 0.1;
FinalTime = 500;
Nx = 50;
Gamma = 0.1;

%[u,v,theta,Et] = data_effect_exp(FinalTime,dt,Nx,mode,Gamma);
[u,v,theta,Et] = p1energyfem1(FinalTime,dt,Nx,mode,Gamma);


AniThermalDisplacement(u,v,theta,Et)