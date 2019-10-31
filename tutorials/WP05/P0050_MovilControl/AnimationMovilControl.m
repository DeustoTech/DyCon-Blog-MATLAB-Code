function AnimationMovilControl(Xnum_with_control,Xnum_free,tspan,xline,yline,Bmatrix)
%ANIMATIONMOVILCONTROL Summary of this function goes here
%   Detailed explanation goes here

Ns = length(xline);
[xms,yms] = meshgrid(xline,yline);

tspan_fine = linspace(tspan(1),tspan(end),8*length(tspan));

Xnum_with_control    = interp1(tspan,Xnum_with_control(:,1:end)',tspan_fine,'linear')';
Xnum_free            = interp1(tspan,Xnum_free',tspan_fine,'linear')';
%%
xline_fine = linspace(xline(1),xline(end),4*length(xline));
yline_fine = linspace(yline(1),yline(end),4*length(yline));
[xms_fine,yms_fine] = meshgrid(xline_fine,yline_fine);
%%


figure('unit','norm','pos',[0 0 1 1],'Color','k');

ax1 = subplot(2,2,1,'Color','none');
hold on
Z = reshape(Xnum_with_control(1:Ns^2,1),Ns,Ns);
isurf = surf(xms_fine,yms_fine,interp2(xms,yms,Z,xms_fine,yms_fine,'linear'),'Parent',ax1);
ax2 = subplot(2,2,2,'Color','none');
Z = reshape(Xnum_free(1:Ns^2,1),Ns,Ns);
jsurf = surf(xms_fine,yms_fine,interp2(xms,yms,Z,xms_fine,yms_fine,'linear'),'Parent',ax2);
%%
x = Xnum_with_control(end-3,1);
y = Xnum_with_control(end-2,1);
radius = 0.1;
z =  1.5;
% Create the mesh of shere
[Xs,Ys,Zs] = sphere;
X = Xs*radius+x ; Y = Ys*radius+y; Z = Zs*radius*0.1+z;
% create the sphere object
%
%

%%
ksurf =surf(X,Y,Z,'Parent',ax1);

%%
ax3 = subplot(2,2,3,'Color','none');
Z = reshape(Xnum_free(1+Ns^2:2*Ns^2,1),Ns,Ns);
lsurf = surf(xms_fine,yms_fine,interp2(xms,yms,Z,xms_fine,yms_fine,'linear'),'Parent',ax3);
ax4 = subplot(2,2,4,'Color','none');
Z = reshape(Xnum_with_control(1+Ns^2:2*Ns^2,1),Ns,Ns);
msurf = surf(xms_fine,yms_fine,interp2(xms,yms,Z,xms_fine,yms_fine,'linear'),'Parent',ax4);


zlim(  [ax1 ax2],[-0.5 2])
zlim(  [ax3 ax4],[-0.05 0.3])

xlim(  [ax1 ax2 ax3 ax4],[-1.5 1.5])
ylim(  [ax1 ax2 ax3 ax4],[-1.5 1.5])
axis(  [ax1 ax2 ax3 ax4],'off')

az = 30;el = 30;
view(  [ax1 ] ,az,el)
view(  [ax2] , az,el)
view(  [ax3 ] ,az,el)
view(  [ax4] ,az,el)

caxis( [ax1 ],[-0.1 0.5])
caxis( [ax2 ],[-0.1 0.5])
caxis( [ax3 ],[-0.01 0.2])
caxis( [ ax4],[-0.01 0.2])
%
shading(   [ax1 ],'interp')
shading(   [ax2],'interp')
shading(   [ax3 ],'interp')
shading(   [ax4],'interp')

lighting(  [ax1 ]','gouraud')
lighting(  [ ax2  ]','gouraud')
lighting(  [  ax3 ]','gouraud')
lighting(  [   ax4]','gouraud')

lightangle([ax1 ],40,40)
lightangle([ ax2],40,40)
lightangle([ax3 ],40,40)
lightangle([ ax4],40,40)

%daspect([ax1 ],[1 1 0.5])
%
title(ax1,'Control Dynamics','Color','w','FontSize',15)
title(ax2,'Free Dynamics','Color','w','FontSize',15)
title(ax3,'Memory Control Dynamics','Color','w','FontSize',15)
title(ax4,'MemoryFree Dynamics','Color','w','FontSize',15)
%%
ax2.Color = 'none';
for it = 1:length(tspan_fine)
    % controled state
    Z = reshape(Xnum_with_control(1:Ns^2,it),Ns,Ns);
    isurf.ZData = interp2(xms,yms,Z,xms_fine,yms_fine,'spline');
    %
    Z = reshape(Xnum_free(1:Ns^2,it),Ns,Ns);
    jsurf.ZData = interp2(xms,yms,Z,xms_fine,yms_fine,'spline');
    %
    Z = reshape(Xnum_with_control(1+Ns^2:2*Ns^2,it),Ns,Ns);
    lsurf.ZData = interp2(xms,yms,Z,xms_fine,yms_fine,'spline');

    Z = reshape(Xnum_free(1+Ns^2:2*Ns^2,it),Ns,Ns);
    msurf.ZData = interp2(xms,yms,Z,xms_fine,yms_fine,'spline');

    %
    %%
    x = Xnum_with_control(end-3,it);
    y = Xnum_with_control(end-2,it);

    % Create the mesh of shere
    X = x+(Xs)*radius ; Y = y+(Ys)*radius; Z = z+Zs*radius*0.5;

    ksurf.XData = X;
    ksurf.YData = Y;
    ksurf.ZData = Z;

    pause(0.05)
end

end

