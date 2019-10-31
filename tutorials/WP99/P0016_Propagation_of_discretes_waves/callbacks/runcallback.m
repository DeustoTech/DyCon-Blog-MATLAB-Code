function runcallback(obj,event,h)
%RUNCALLBACK Summary of this function goes here
%   Detailed explanation goes here
N = 100;
T = 5;
y0 = h.DA.oneD.slide.Value;

w = str2num(h.DA.oneD.editfrec.String)*pi;

% %% fig 1
% y0 = 0;
% w  = pi/4;
% %% fig 2
% y0 = 0;
% w  = pi;
% %% fig 3
% y0 = 0.5;
% w  = pi;
%%
hwait = waitbar(0);

axes = {h.DA.oneD.ax1,h.DA.oneD.ax2,h.DA.oneD.ax3};


for i = 1:3
    ind = i-1;
    result1D =  wave_nonunif(N,1,T,y0,w,@(x) 0*x + 1,@(x) 0*x + 1,ind);
    show1D(result1D,axes{i})
    waitbar(i/3,hwait)
end
delete(hwait)
end

