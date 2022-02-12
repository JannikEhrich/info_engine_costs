%SCRIPTNAME does...
%
% OUTPUTS:
%  outputs XXX
%
% author:  JEhrich
% version: 0.0 (2022-02-10)
% changes: -
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
nu_relax = 3E-3;
% controller mobility during measurement-feedback
nu_meas = 8E0;
% measurement error
s2 = 0.005;
% total time interval
ts = 0.095;
% measurement time interval
tau = 0.005;
% integration time-step
dt = 1E-6;
% number of time steps
K = 5E0;

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
% initialize with equilibrium distribution
x = randn;
z = x + sqrt(s2)*randn;
% data structures for rajectories
x_traj = nan(K*n+1,1);
z_traj = nan(K*n+1,1);
x_traj(1) = x;
z_traj(1) = z;
for ii = 1:K
    % relaxation
    for jj = 1:n_relax
        dx = -(x-z)*dt + sqrt(2*dt)*randn;
        dz = nu_relax*(x-z)*dt + sqrt(2*nu_relax*dt)*randn;
        x = x + dx;
        z = z + dz;
        x_traj((ii-1)*n+jj+1) = x;
        z_traj((ii-1)*n+jj+1) = z;
    end
    % measurement-feedback
    for jj = 1:n_meas
        % stiffness
        k = k0 + jj/n_meas*(k1 - k0);
        dx = -k*(x-z)*dt + sqrt(2*dt)*randn;
        dz = nu_meas*k*(x-z)*dt + sqrt(2*nu_meas*dt)*randn;
        x = x + dx;
        z = z + dz;
        x_traj((ii-1)*n+n_relax+jj+1) = x;
        z_traj((ii-1)*n+n_relax+jj+1) = z;
    end
end





%% main figure
% set inset axes limits
inset1 = [2*ts-2*tau,2*ts+1*tau,-0.75,0.47];
inset2 = [3*ts-2*tau,3*ts+1*tau,-0.75,0.47];

figure('Position',[1000,1000,560,420]);
axes('Position',[.05 .07 .9 .48])
plot(nan,nan,'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'b','LineWidth',lW,'MarkerSize',mS);
% mark different regions
for ii = 1:K
    % mark relaxation region
    patch([(ii-1)*ts, ii*ts-tau, ii*ts-tau, (ii-1)*ts],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
    % mark measurement-feedback region
    patch([ii*ts-tau, ii*ts, ii*ts, ii*ts-tau],[-9 -9 9 9],...
        [0,0,1],'EdgeColor','none','FaceAlpha',0.2);
    %plot(ii*ts-tau*[1,1],[-2,2],'LineWidth',lW,'MarkerSize',mS,'color',0.5*[1,1,1]);
    %hold on;
    %plot(ii*ts*[1,1],[-2,2],'LineWidth',lW,'MarkerSize',mS,'color',0.5*[1,1,1]);
end
plot(0:dt*1E2:K*ts,x_traj(1:1E2:end),'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(0:dt*1E2:K*ts,z_traj(1:1E2:end),'b','LineWidth',lW,'MarkerSize',mS);
axis([0,K*ts,-0.9,0.6]);
set(gca,'XTick',[0,ts-tau,ts,2*ts-tau,2*ts,3*ts-tau,3*ts,4*ts-tau,4*ts,5*ts-tau,5*ts],...
    'XTicklabel',[],'YTick',[]);

% axes tick labels
text(0,-1,'$0$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(ts-tau-0.015,-1,'$t_\mathrm{s}\!\!-\!\tau$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(ts+0.009,-1,'$t_\mathrm{s}$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(2*ts-tau-0.02,-1,'$2t_\mathrm{s}\!\!-\!\tau$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(2*ts+0.015,-1,'$2t_\mathrm{s}$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(3*ts-tau-0.02,-1,'$3t_\mathrm{s}\!\!-\!\tau$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(3*ts+0.015,-1,'$3t_\mathrm{s}$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(4*ts-tau-0.02,-1,'$4t_\mathrm{s}\!\!-\!\tau$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(4*ts+0.015,-1,'$4t_\mathrm{s}$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(5*ts-tau-0.02,-1,'$5t_\mathrm{s}\!\!-\!\tau$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(5*ts+0.015,-1,'$5t_\mathrm{s}$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');


box off
%axis off
% mark inset 1
plot([inset1(1) inset1(2)],[inset1(3) inset1(3)],'k');
plot([inset1(2) inset1(2)],[inset1(3) inset1(4)],'k');
plot([inset1(1) inset1(2)],[inset1(4) inset1(4)],'k');
plot([inset1(1) inset1(1)],[inset1(3) inset1(4)],'k');
% mark inset 2
plot([inset2(1) inset2(2)],[inset2(3) inset2(3)],'k');
plot([inset2(2) inset2(2)],[inset2(3) inset2(4)],'k');
plot([inset2(1) inset2(2)],[inset2(4) inset2(4)],'k');
plot([inset2(1) inset2(1)],[inset2(3) inset2(4)],'k');

% mark relaxation
annotation('doublearrow','Position',[0.05,0.34,0.17,0],'Head1Style','plain',...
    'Head2Style','plain','Linewidth',lW,'color','r');
annotation('textbox','Position',[0.09,0.32, 0.1, 0.1],'String','relaxation',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','r');
annotation('textbox','Position',[0.09,0.25, 0.1, 0.1],'String','V(x,z)',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','r');

annotation('textarrow',[0.26 0.225],[0.545 0.41],'String','feedback',...
    'Interpreter','latex','FontSize',fS,'Linewidth',lW,'HeadStyle','plain','color','b');
annotation('doublearrow','Position',[0.1925,0.62,0.14,0],'Head1Style','plain',...
    'Head2Style','plain','Linewidth',lW,'color','b');
annotation('textbox','Position',[0.1925,0.6,0.14,0.1],'String','U(x,z;t)',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','b');
legend({'system $X$','controller $Z$'},'Location','SouthEast');
set(gca,'FontSize',fS);

%% inset 1
axes('Position',[0.05 .6 .425 .35])
box on
% mark regions
patch([2*ts-2*tau, 2*ts-tau, 2*ts-tau, 2*ts-2*tau],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
patch([2*ts-tau, 2*ts, 2*ts, 2*ts-tau],[-9 -9 9 9],...
        [0,0,1],'EdgeColor','none','FaceAlpha',0.2);
patch([2*ts, 2*ts+tau, 2*ts+tau, 2*ts],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
hold on;
plot(2*ts-2*tau:dt*1E1:2*ts+tau,x_traj(2*n-2*n_meas:1E1:2*n+n_meas),'r','LineWidth',lW,'MarkerSize',mS);
plot(2*ts-2*tau:dt*1E1:2*ts+tau,z_traj(2*n-2*n_meas:1E1:2*n+n_meas),'b','LineWidth',lW,'MarkerSize',mS);
axis(inset1);
set(gca,'XTick',[],'YTick',[]);

% mark states 1
plot(2*ts-tau,x_traj(2*n-n_meas+1),'.r','LineWidth',lW,'MarkerSize',45);
text(2*ts-tau,x_traj(2*n-n_meas+1)+0.16,'$x_1^+$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(2*ts-tau,z_traj(2*n-n_meas+1),'.b','LineWidth',lW,'MarkerSize',45);
text(2*ts-tau-0.0005,z_traj(2*n-n_meas+1)+0.16,'$z_1^+$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');

% mark states 2
plot(2*ts,x_traj(2*n+1),'.r','LineWidth',lW,'MarkerSize',45);
text(2*ts-0,x_traj(2*n+1)+0.16,'$x_2$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(2*ts,z_traj(2*n+1),'.b','LineWidth',lW,'MarkerSize',45);
text(2*ts-0,z_traj(2*n+1)-0.16,'$z_2$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');


%% inset 2
axes('Position',[0.525 .6 .425 .35])
box on
% mark regions
patch([3*ts-2*tau, 3*ts-tau, 3*ts-tau, 3*ts-2*tau],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
patch([3*ts-tau, 3*ts, 3*ts, 3*ts-tau],[-9 -9 9 9],...
        [0,0,1],'EdgeColor','none','FaceAlpha',0.2);
patch([3*ts, 3*ts+tau, 3*ts+tau, 3*ts],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
hold on;
plot(3*ts-2*tau:dt*1E1:3*ts+tau,x_traj(3*n-2*n_meas:1E1:3*n+n_meas),'r','LineWidth',lW,'MarkerSize',mS);
plot(3*ts-2*tau:dt*1E1:3*ts+tau,z_traj(3*n-2*n_meas:1E1:3*n+n_meas),'b','LineWidth',lW,'MarkerSize',mS);
axis(inset2);
set(gca,'XTick',[3*ts-tau,3*ts],'XTickLabel',[],'YTick',[]);
text(3*ts-tau,-0.85,'$3t_\mathrm{s}\!\!-\!\tau$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
text(3*ts,-0.85,'$3t_\mathrm{s}$','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');


% mark states 2
plot(3*ts-tau,x_traj(3*n-n_meas+1),'.r','LineWidth',lW,'MarkerSize',45);
text(3*ts-tau,x_traj(3*n-n_meas+1)-0.16,'$x_2^+$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(3*ts-tau,z_traj(3*n-n_meas+1),'.b','LineWidth',lW,'MarkerSize',45);
text(3*ts-tau,z_traj(3*n-n_meas+1)+0.16,'$z_2^+$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');

% mark states 3
plot(3*ts,x_traj(3*n+1),'.r','LineWidth',lW,'MarkerSize',45);
text(3*ts-0,x_traj(3*n+1)-0.16,'$x_3$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(3*ts,z_traj(3*n+1),'.b','LineWidth',lW,'MarkerSize',45);
text(3*ts-0,z_traj(3*n+1)+0.16,'$z_3$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');

%% connect insets with main plot
annotation('line',[0.405 0.405],[0.507 0.60]);
annotation('line',[0.585 0.585],[0.507 0.60]);

%% export
saveas(gcf, 'feedback_process.eps','epsc')

