function animationsparce(Uall,N,n)
%ANIMATIONSPARCE Summary of this function goes here
%   Detailed explanation goes here
figure('unit','norm','pos',[0.1 0.1 0.8 0.8])

tspan = 1:n;
new_tspan = 1:0.01:n;
Uall = interp1(tspan,Uall,new_tspan);
isurf = surf(reshape(Uall(1,:),N,length(Uall(1,:))/N),'FaceLighting','gouraud');
title('Evolution of sparse sources identified ','FontSize',18)

lightangle(62,62)
lightangle(62,-62)

shading interp
colormap jet
axis('off')
zlim([-1 10])
view(-39,30)
%caxis([0 20])
gif('evolution.gif','frame',isurf.Parent.Parent,'nodither')
for it = 1:length(new_tspan)
isurf.ZData = (reshape(Uall(it,:),N,length(Uall(it,:))/N));
view(-158+0.1*it,20+0.05*it)
gif
pause(0.01)
end

for it = 1:50
isurf.ZData = (reshape(Uall(end,:),N,length(Uall(end,:))/N));
gif('frame',isurf.Parent.Parent)
pause(0.01)
end

end

