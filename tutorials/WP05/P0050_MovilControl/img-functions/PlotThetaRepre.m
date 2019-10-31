clear;
fig = figure

xl = linspace(0,1,500);

plot(xl,thetaWP05(xl - (0.4)) - thetaWP05(xl - (0.8)),'-','LineWidth',3)
hold on

plot(xl,thetaWP05(xl - (0.4)),'--','LineWidth',2)

plot(xl,thetaWP05(xl - (0.8)),'--','LineWidth',2)


yticks([0 1])
ylim([-0.1 1.1])
title('\chi_{(a,b)}=\chi_{(0.4,0.8)}')
ylabel('')
xlabel('x')
xline(0.4)
xline(0.8)
legend({'W(x,a,b)','\Theta(x-a)','\Theta(x-b)'},'Location','northeastoutside')

saveas(fig,'img/WTHETA.png')