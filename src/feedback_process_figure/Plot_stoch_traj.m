%PLOT_STOCH_TRAJ produces illustrative plot of system and controller
%trajectories for the feedback-cooling model, highlighting the role of
%time-scale separation
%
% OUTPUTS:
%  outputs eps figure of example trajectories with insets
%
% author:  JEhrich
% version: 1.2 (2022-03-01)
% changes: renamed nu_high and nu_low
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
        dz = nu_low*(x-z)*dt + sqrt(2*nu_low*dt)*randn;
        x = x + dx;
        z = z + dz;
        x_traj((ii-1)*n+jj+1) = x;
        z_traj((ii-1)*n+jj+1) = z;
    end
    % control
    for jj = 1:n_meas
        % stiffness
        k = k0 + jj/n_meas*(k1 - k0);
        dx = -k*(x-z)*dt + sqrt(2*dt)*randn;
        dz = nu_high*k*(x-z)*dt + sqrt(2*nu_high*dt)*randn;
        x = x + dx;
        z = z + dz;
        x_traj((ii-1)*n+n_relax+jj+1) = x;
        z_traj((ii-1)*n+n_relax+jj+1) = z;
    end
end





%% main figure
% set inset axes limits
inset1 = [3*ts-2*tau,3*ts+1*tau,-0.42,0.4];

figure('Position',[1000,1000,560,420]);
axes('Position',[0.09 .14 .85 .82])
plot(nan,nan,'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'b','LineWidth',lW,'MarkerSize',mS);
% mark different regions
for ii = 1:K
    % mark relaxation region
    patch([(ii-1)*ts, ii*ts-tau, ii*ts-tau, (ii-1)*ts],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.15);
    % mark measurement-feedback region
    patch([ii*ts-tau, ii*ts, ii*ts, ii*ts-tau],[-9 -9 9 9],...
        [0,0,1],'EdgeColor','none','FaceAlpha',0.15);
    %plot(ii*ts-tau*[1,1],[-2,2],'LineWidth',lW,'MarkerSize',mS,'color',0.5*[1,1,1]);
    %hold on;
    %plot(ii*ts*[1,1],[-2,2],'LineWidth',lW,'MarkerSize',mS,'color',0.5*[1,1,1]);
end
plot(0:dt*1E2:K*ts,x_traj(1:1E2:end),'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(0:dt*1E2:K*ts,z_traj(1:1E2:end),'b','LineWidth',lW,'MarkerSize',mS);
axis([0,K*ts,-0.9,1.7]);
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
annotation('doublearrow','Position',[0.091,0.7,0.195,0],'Head1Style','plain',...
    'Head2Style','plain','Linewidth',lW,'color','r');
annotation('textbox','Position',[0.143,0.68, 0.1, 0.1],'String','relaxation',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','r');
annotation('textbox','Position',[0.143,0.61, 0.1, 0.1],'String','$V_\mathrm{r}(x,z), \nu_\mathrm{low}$',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','r');

% mark feedback
plot([ts-tau,ts-tau],[0.35 0.45],'b');
plot([ts,ts],[0.35 0.45],'b');
annotation('textbox','Position',[0.2965,0.54,0.0,0.1],'String','feedback',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','b');
annotation('arrow',[0.248 0.288],0.55*[1 1],'Linewidth',lW,'HeadStyle','plain','color','b');
annotation('arrow',[0.345 0.305],0.55*[1 1],'Linewidth',lW,'HeadStyle','plain','color','b');
annotation('textbox','Position',[0.2965,0.45,0.0,0.1],'String','$V_\mathrm{c}(x,z;t), \nu_\mathrm{high}$',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none','color','b');

set(gca,'FontSize',fS);

% mark tau
plot([2*ts-tau,2*ts-tau],[-0.7 -0.8],'k');
plot([2*ts,2*ts],[-0.7 -0.8],'k');
annotation('arrow',[0.46 0.5],0.188*[1 1],'Linewidth',lW,'HeadStyle','plain');
annotation('arrow',[0.558 0.518],0.188*[1 1],'Linewidth',lW,'HeadStyle','plain');
annotation('textbox','Position',[0.509,0.175, 0.0, 0.1],'String','$\tau$',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none');


legend({'system $X$','controller $Z$'},'Location','NorthWest');
%% inset
axes('Position',[0.54 .65 .37 .29])
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
axis(inset1);
set(gca,'XTick',[3*ts-tau, 3*ts],'XTicklabel',{'$3t_\mathrm{s}\!\!-\!\tau$','$3t_\mathrm{s}$'},'YTick',[]);
set(gca,'FontSize',fS);

% mark states 1
plot(3*ts-tau,x_traj(3*n-n_meas+1),'.r','LineWidth',lW,'MarkerSize',45);
text(3*ts-tau-0.0007,x_traj(3*n-n_meas+1)+0.12,'$x_3^-$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(3*ts-tau,z_traj(3*n-n_meas+1),'.b','LineWidth',lW,'MarkerSize',45);
text(3*ts-tau+0.0007,z_traj(3*n-n_meas+1)+0.13,'$z_3^-$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');

% mark states 2
plot(3*ts,x_traj(3*n+1),'.r','LineWidth',lW,'MarkerSize',45);
text(3*ts+0.0007,x_traj(3*n+1)-0.1,'$x_3$','color','r','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');
plot(3*ts,z_traj(3*n+1),'.b','LineWidth',lW,'MarkerSize',45);
text(3*ts+0.0007,z_traj(3*n+1)+0.1,'$z_3$','color','b','FontSize',fS,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'middle');


%% connect insets with main plot
annotation('line',[0.54 0.705],[0.65 0.55]);
annotation('line',[0.91 0.739],[0.65 0.55]);

%% export
saveas(gcf, '../../doc/feedback_process.eps','epsc')

