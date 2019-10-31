function show1D(result,ax)
%SHOW1D Summary of this function goes here
%   Detailed explanation goes here

L = result.L;
T = result.T;
mesh = result.mesh;
t = result.t;
ind = result.ind ;
sol = result.sol ;



surf(mesh.yi,t',sol','EdgeColor','none','LineStyle','none','FaceLighting','gouraud','Parent',ax)
lightangle(ax,40,40)
colormap(ax,'jet')
shading(ax,'interp')
axis(ax,[mesh.yi(1) mesh.yi(end) 0 T 0 2])
xticks(ax,[mesh.yi(1) -0.5*L 0 0.5*L mesh.yi(end)])
xticklabels(ax,{'-1','-0.5','0','0.5','1'})
yticks(ax,[0 1 2 3 4 5])
view(ax,0,90)
ax.FontSize = 12;
xlabel(ax,'Space')
ylabel(ax,'Time')


end

