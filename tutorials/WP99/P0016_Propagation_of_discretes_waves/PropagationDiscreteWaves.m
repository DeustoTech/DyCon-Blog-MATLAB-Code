function PropagationDiscreteWaves

pathfile = which('PropagationDiscreteWaves');
path = replace(pathfile,'PropagationDiscreteWaves.m','');
addpath(genpath(path))
h = fighandle;
h.fig = figure('unit','norm','pos',[0.1 0.1 0.8 0.8],'MenuBar','none','Toolbar','none');
%% 
uig = uitabgroup(h.fig);
uitab1 = uitab('Title','1D','Parent',uig);
uitab2 = uitab('Title','2D','Parent',uig);
%%
L = 1;
N = 15;
g1 = @(x) x;
y1 = g1(linspace(-L,L,N+2)');

g2 = @(x) L*tan((pi/(4*L)).*x);
y2 = g2(linspace(-L,L,N+2)');

   
g3 = @(x) 2*L*sin((pi/(6*L)).*x);
y3 = g3(linspace(-L,L,N+2)');


% ind = choiche of the mesh: ind == 0 --> uniform mesh;
%                            ind == 1 --> refined mesh g(x) = L*tan((pi/(4*L)).*x);
%                            otherwise --> refined mesh g(x) = 2*L*sin((pi/(6*L)).*x);
xline = y1;


yy = 0.9;
hh = 0.1;
axtop1 = panelaxis('',   uitab1,[0.00 yy 1/3 hh],'noaxis');
plot(y1,g1(y1),'.','Parent',axtop1,'MarkerSize',15)
xticks(axtop1,y1)
xticklabels(axtop1,[])
axtop1.Title.String = 'Uniform';
axtop1.Title.FontSize = 12;

axtop2 = panelaxis('',uitab1,[1/3 yy 1/3 hh],'noaxis');
plot(y2,g2(y2),'.','Parent',axtop2,'MarkerSize',15)
axtop2.Title.String = 'L*tan((\pi/(4*L)).*x)';
axtop2.Title.FontSize = 12;

xticks(axtop2,y2)

xticklabels(axtop2,[])

axtop3 = panelaxis('',uitab1,[2/3 yy 1/3 hh],'noaxis');
plot(y3,g3(y3),'.','Parent',axtop3,'MarkerSize',15)
axtop3.Title.String = '2*L*sin((\pi/(6*L)).*x';
axtop3.Title.FontSize = 12;

xticks(axtop3,y3)

xticklabels(axtop3,[])

ax1 = panelaxis('',   uitab1,[0.00 0.1 1/3 0.8]);
ax2 = panelaxis('',uitab1,[1/3 0.1 1/3 0.8]);
ax3 = panelaxis('',uitab1,[2/3 0.1 1/3 0.8]);

%%
uip = uipanel('pos',[0.0 0.0 1.0 0.1],'Parent',uitab1);
btn = uicontrol('string','run','Parent',uip,'Callback',{@runcallback,h},'unit','norm','pos',[0.85 0.1 0.1 0.3]);

uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.3 0.7 0.3 0.2],'String','Ray''s Initial Position');
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.3 0.1 0.3 0.3],'String','frequency');
slide = uicontrol('Style','slider','Parent',uip,'unit','norm','pos',[0.3 0.15 0.3 0.25]);
slide.Min = -1;
slide.Max = 1;

uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.7 0.5 0.05 0.3],'String','frequency');
editfrec = uicontrol('Style','edit','Parent',uip,'unit','norm','pos',[0.7 0.1 0.05 0.3],'String',1);
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.75 0.075 0.025 0.3],'String','pi');


h.DA.oneD.ax1 = ax1;
h.DA.oneD.ax2 = ax2;
h.DA.oneD.ax3 = ax3;
%
h.DA.oneD.slide = slide;
h.DA.oneD.editfrec = editfrec;




%%
%%
%%
%%
%%
FontSize = 12;
axtop1 = panelaxis('',   uitab2,[0.00 yy 1/3 hh],'noaxis');
plot(y1,g1(y1),'.','Parent',axtop1,'MarkerSize',25)
xticks(axtop1,y1)
xticklabels(axtop1,[])
axtop1.Title.String = 'Uniform';
axtop1.Title.FontSize = FontSize;

axtop2 = panelaxis('',uitab2,[1/3 yy 1/3 hh],'noaxis');
plot(y2,g2(y2),'.','Parent',axtop2,'MarkerSize',25)
axtop2.Title.String = 'L*tan((\pi/(4*L))*x)';
axtop2.Title.FontSize = FontSize;
xticks(axtop2,y2)

xticklabels(axtop2,[])

axtop3 = panelaxis('',uitab2,[2/3 yy 1/3 hh],'noaxis');
plot(y3,g3(y3),'.','Parent',axtop3,'MarkerSize',25)
axtop3.Title.String = '2*L*sin((\pi/(6*L))*x';
axtop3.Title.FontSize = FontSize;
xticks(axtop3,y3)

xticklabels(axtop3,[])

ax1 = panelaxis('',   uitab2,[0.00 0.15  1/3 0.75]);
ax2 = panelaxis('',uitab2,[1/3  0.15  1/3 0.75]);
ax3 = panelaxis('',uitab2,[2/3  0.15  1/3 0.75]);

% time 
axtime = panelaxis('',uitab2,[0  0.1  1 0.05]);
xlim(axtime,[0 1])
axis(axtime,'off')
axtime.Units = 'norm';

axtime.Position = [0.025 0.25 0.95 0.5];
rectangle('parent',axtime);
irect = rectangle('parent',axtime,'FaceColor','red');
irect.Position(3) = 0.0;
h.DA.irect = irect;
%%
uip = uipanel('pos',[0.0 0.0 1.0 0.1],'Parent',uitab2);
btn = uicontrol('string','run','Parent',uip,'Callback',{@runcallback2D,h},'unit','norm','pos',[0.825 0.15 0.05 0.4]);
%
stopbtn  = uicontrol('style','pushbutton','Parent',uip,'unit','norm','pos',[0.9 0.15 0.05 0.4],'String','stop','Callback',{@stopani,h});
%%
yy = 0.2;
%
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.05 0.7 0.2 0.2],'String','Ray''s x - Initial Position');
slidex = uicontrol('Style','slider','Parent',uip,'unit','norm','pos',[0.05 yy 0.2 0.25]);
slidex.Min = -1;
slidex.Max = 1;
%
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.35 0.7 0.2 0.2],'String','Ray''s y - Initial Position');
slidey = uicontrol('Style','slider','Parent',uip,'unit','norm','pos',[0.35 yy 0.2 0.25]);
slidey.Min = -1;
slidey.Max = 1;
% 
yy = 0.2;
% Frecuency 1
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.6 0.5 0.05 0.3],'String','Frecuency');
editfrecxi0 = uicontrol('Style','edit','Parent',uip,'unit','norm','pos',[0.6 yy 0.05 0.3],'String',1);
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.65 yy-0.05 0.025 0.3],'String','pi');
% Frecuency 2
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.7 0.5 0.05 0.3],'String','Frecuency');
editfreceta0 = uicontrol('Style','edit','Parent',uip,'unit','norm','pos',[0.7 yy 0.05 0.3],'String',1);
uicontrol('Style','text','Parent',uip,'unit','norm','pos',[0.75 yy-0.05 0.025 0.3],'String','pi');
% 
h.DA.twoD.ax1 = ax1;
h.DA.twoD.ax2 = ax2;
h.DA.twoD.ax3 = ax3;
%
h.DA.twoD.slidex = slidex;
h.DA.twoD.slidey = slidey;

h.DA.twoD.editfreceta0 = editfreceta0;
h.DA.twoD.editfrecxi0  = editfrecxi0;

%%

function ax= panelaxis(String,Parent,Pos,noaxis)
    uip = uipanel('Parent',Parent,'Title',String,'Unit','norm','Pos',Pos);
    ax = axes('Parent',uip,'Tag',String);
    if exist('noaxis')
       axis(ax,'off') 
       uip.Title = '';
       ax.YLim = [-1 1];
    end
end
end