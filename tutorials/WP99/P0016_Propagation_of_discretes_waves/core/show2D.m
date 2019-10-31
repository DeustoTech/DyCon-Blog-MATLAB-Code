function show2D(results,ax,irect)
%SHOW2D Summary of this function goes here
%   Detailed explanation goes here

for i=1:3
delete(ax{i}.Children)
data2 = results(i).data2;
solhamx = results(i).solhamx;
solhamy = results(i).solhamy;
tham    = results(i).tham;
lt      = results(i).lt;

isurf(i) = surf(abs(data2{1}),'Parent',ax{i},'FaceLighting','gouraud');
N = length(data2{1});
xlim(ax{i},[0 N]);
ylim(ax{i},[0 N]);

shading(ax{i},'interp')
axis(ax{i},'off')
view(ax{i},90,-90)
zlim(ax{i},[0,3])
lightangle(ax{i},-50,-50)
%caxis([0 1])
hold(ax{i},'on')
plot3(solhamx(:,1),solhamy(:,1),tham,'w','LineWidth',3,'Parent',ax{i})
end
fig = ax{1}.Parent.Parent.Parent.Parent;

gif('fig.3.gif','frame',fig)
for m = 1:lt-1
    for i=1:3
        data2 = results(i).data2;
        isurf(i).ZData = abs(data2{m});
    end
    
    irect.Position(3) = m/(lt-1);
    pause(0.05)
    if strcmp(ax{1}.Tag,'stop')
       ax{1}.Tag = 'run';
       irect.Position(3) = 0;
       return 
    end
    pause(0.05)
    gif('frame',fig)

end
irect.Position(3) = 0;


end

