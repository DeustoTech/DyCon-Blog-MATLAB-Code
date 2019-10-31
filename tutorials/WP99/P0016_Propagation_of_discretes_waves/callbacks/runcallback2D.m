function runcallback(obj,event,h)
%RUNCALLBACK Summary of this function goes here
%   Detailed explanation goes here
N = 100;
T = 8;

x0 = h.DA.twoD.slidex.Value;
y0 = h.DA.twoD.slidey.Value;

eta0 = str2num(h.DA.twoD.editfreceta0.String)*pi;
xi0 = str2num(h.DA.twoD.editfrecxi0.String)*pi;
%% fig 1

x0 = 0;
y0 = 1/2;

eta0 = 0.25*pi;
xi0  = 0.25*pi;
%% fig 2

x0 = 0;
y0 = 0;

eta0 = pi;
xi0  = pi;
%% fig 2

y0 = 0;
x0 = tan(acos(0.5^(0.25)));

eta0 = 0.5*pi;
xi0  = pi;
%%
hwait = waitbar(0);

axes = {h.DA.twoD.ax1,h.DA.twoD.ax2,h.DA.twoD.ax3};


for i = 1:3
    ind = i-1;
    result2D(i) =  waves_2d(N,N,x0,y0,xi0,eta0,15,15,0.1,ind);
    waitbar(i/3,hwait)
end
delete(hwait)

 show2D(result2D,axes,h.DA.irect)
 

end

