%PLOT_STOCH_TRAJ produces illustrative plot of system and controller
%trajectories for the feedback-cooling model, highlighting the role of
%time-scale separation
%
% OUTPUTS:
%  outputs eps figure of example trajectories with insets
%
% author:  JEhrich
% version: 1.4 (2022-05-06)
% changes: removed potential energy labels
clear
close all
clc
% set font size, line width, and marker size
fS = 18;
lW = 2.0;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% parameters
% seed RNG
rng(1)
% controller mobility during relaxation
nu_low = 3E-3;
% controller mobility during control operation
nu_high = 8E0;
% measurement error
s2 = 0.005;
% total time interval
ts = 0.095;
% measurement time interval
tau = 0.005;
% integration time-step
dt = 1E-6;
% number of time steps
K = 4E0;

% number of steps per time steps
n = round(ts/dt);
% number of steps during measurement-feedback
n_meas = round(tau/dt);
n_relax = n - n_meas;

% measurement: initial stiffness
k0 = 1/((s2 - 1)*exp(-2*ts) + 1);
% measurement: final stiffness
k1 = 1/s2;

%% simulate process
% initialize with equilibrium distribution and measurement drawn from
% pre-measurement distribution
x = randn;
z = x + sqrt(1/k0)*randn;
% data structures for trajectories
x_traj = nan(K*n+1,1);
z_traj = nan(K*n+1,1);
x_traj(1) = x;
z_traj(1) = z;
for ii = 1:K
    % control
    for jj = 1:n_meas
        % stiffness
        k = k0 + jj/n_meas*(k1 - k0);
        dx = -k*(x-z)*dt + sqrt(2*dt)*randn;
        dz = nu_high*k*(x-z)*dt + sqrt(2*nu_high*dt)*randn;
        x = x + dx;
        z = z + dz;
        x_traj((ii-1)*n+jj+1) = x;
        z_traj((ii-1)*n+jj+1) = z;
    end
    % relaxation
    for jj = 1:n_relax
        dx = -(x-z)*dt + sqrt(2*dt)*randn;
        dz = nu_low*(x-z)*dt + sqrt(2*nu_low*dt)*randn;
        x = x + dx;
        z = z + dz;
        x_traj((ii-1)*n+n_meas+jj+1) = x;
        z_traj((ii-1)*n+n_meas+jj+1) = z;
    end
end





%% main figure
% set inset axes limits
inset1 = [3*ts-1*tau,3*ts+2*tau,-0.32,0.5];

figure('Position',[1300,1000,560,420]);
axes('Position',[0.09 .14 .85 .82])
plot(nan,nan,'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'b','LineWidth',lW,'MarkerSize',mS);
% mark different regions
for ii = 0:K
    % mark control region
    patch([ii*ts, ii*ts+tau, ii*ts+tau, ii*ts],[-9 -9 9 9],...
        [0,0,1],'EdgeColor','none','FaceAlpha',0.15);
    % mark relaxation region
    patch([(ii-1)*ts+tau, ii*ts, ii*ts, (ii-1)*ts+tau],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.15);
end
plot(0:dt*1E2:K*ts,x_traj(1:1E2:end),'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(0:dt*1E2:K*ts,z_traj(1:1E2:end),'b','LineWidth',lW,'MarkerSize',mS);
axis([0,K*ts,-0.85,1.8]);
set(gca,'XTick',(0:4)*ts,...
    'XTicklabel',{'$0$','$t_\mathrm{s}$','$2t_\mathrm{s}$',...
    '$3t_\mathrm{s}$','$4t_\mathrm{s}$'},'YTick',[]);
ylabel('state','interpreter','latex');
xlabel('time','interpreter','latex');

box off
%axis off
% mark inset 1
plot([inset1(1) inset1(2)],[inset1(3) inset1(3)],'k');
plot([inset1(2) inset1(2)],[inset1(3) inset1(4)],'k');
plot([inset1(1) inset1(2)],[inset1(4) inset1(4)],'k');
plot([inset1(1) inset1(1)],[inset1(3) inset1(4)],'k');

% mark relaxation
annotation('doublearrow','Position',[0.102,0.7,0.197,0],'Head1Style','plain',...
    'Head2Style','plain','Linewidth',lW,'color','r');
annotation('textbox','Position',[0.1495,0.68, 0.1, 0.1],'String','relaxation',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','r');
annotation('textbox','Position',[0.1495,0.62, 0.1, 0.1],'String','$\nu_\mathrm{low}$',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','r');

% mark feedback
plot([ts+tau,ts+tau],[0.41 0.54],'b');
plot([ts,ts],[0.41 0.54],'b');
annotation('textbox','Position',[0.3075,0.54,0.0,0.1],'String','control',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','b');
annotation('arrow',[0.258 .298],0.55*[1 1],'Linewidth',lW,'HeadStyle','plain','color','b');
annotation('arrow',[0.357 0.317],0.55*[1 1],'Linewidth',lW,'HeadStyle','plain','color','b');
annotation('textbox','Position',[0.3075,0.45,0.0,0.1],'String','$\nu_\mathrm{high}$',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','b');

set(gca,'FontSize',fS);

% mark tau
plot([2*ts,2*ts],[-0.635 -0.765],'k');
plot([2*ts+tau,2*ts+tau],[-0.635 -0.765],'k');
annotation('arrow',[0.47 0.51],0.188*[1 1],'Linewidth',lW,'HeadStyle','plain');
annotation('arrow',[0.568 0.528],0.188*[1 1],'Linewidth',lW,'HeadStyle','plain');
annotation('textbox','Position',[0.519,0.175, 0.0, 0.1],'String','$\tau$',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none');


legend({'system $X$','controller $Z$'},'Location','NorthWest');

%% inset
axes('Position',[0.55 .66 .37 .27])
box on
% mark regions
patch([3*ts-1*tau, 3*ts, 3*ts, 3*ts-1*tau],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
patch([3*ts, 3*ts+tau, 3*ts+tau, 3*ts],[-9 -9 9 9],...
        [0,0,1],'EdgeColor','none','FaceAlpha',0.2);
patch([3*ts+tau, 3*ts+2*tau, 3*ts+2*tau, 3*ts+tau],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
hold on;
plot(3*ts-1*tau:dt*1E1:3*ts+2*tau,x_traj(3*n-1*n_meas:1E1:3*n+2*n_meas),'r','LineWidth',lW,'MarkerSize',mS);
axis(inset1);
set(gca,'XTick',[3*ts, 3*ts+tau],'XTicklabel',{'$3t_\mathrm{s}$','$3t_\mathrm{s}\!+\!\tau$'},'YTick',[]);
set(gca,'FontSize',fS);

% mark states 1
plot(3*ts,x_traj(3*n+1),'.r','LineWidth',lW,'MarkerSize',45);
text(3*ts-0.0007,x_traj(3*n+1)+0.1,'$x_3$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(3*ts,z_traj(3*n+1),'.b','LineWidth',lW,'MarkerSize',45);
text(3*ts-0.0007,z_traj(3*n+1)-0.1,'$\approx\!z_2$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');

% mark states 2
%plot(3*ts+tau,x_traj(3*n+n_meas+1),'.r','LineWidth',lW,'MarkerSize',45);
%text(3*ts+tau+0.0005,x_traj(3*n+n_meas+1)+0.15,'$x_3$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(3*ts+tau,z_traj(3*n+n_meas+1),'.b','LineWidth',lW,'MarkerSize',45);
text(3*ts+tau-0.001,z_traj(3*n+n_meas+1)-0.05,'$z_3$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');

plot(3*ts-1*tau:dt*1E1:3*ts+2*tau,z_traj(3*n-1*n_meas:1E1:3*n+2*n_meas),'b','LineWidth',lW,'MarkerSize',mS);

%% connect insets with main plot
annotation('line',[0.5505 0.716],[0.66 0.557]);
annotation('line',[0.92 0.75],[0.66 0.557]);

%% export
saveas(gcf, '../../doc/feedback_process.eps','epsc')

