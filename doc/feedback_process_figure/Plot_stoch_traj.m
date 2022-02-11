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
inset1 = [2*ts-2*tau,2*ts+1*tau,-0.6,0.4];
inset2 = [3*ts-2*tau,3*ts+1*tau,-0.4,0.33];

figure('Position',[1000,1000,560,420]);
axes('Position',[.05 .05 .9 .5])
% mark different regions
for ii = 1:K
    % mark relaxation region
    patch([(ii-1)*ts, ii*ts-tau, ii*ts-tau, (ii-1)*ts],[-9 -9 9 9],...
        [1,0,0],'EdgeColor','none','FaceAlpha',0.2);
    hold on;
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
set(gca,'XTick',[],'YTick',[]);
box on
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
annotation('doublearrow','Position',[0.05,0.31,0.17,0],'Head1Style','plain',...
    'Head2Style','plain','Linewidth',lW);
annotation('textbox','Position',[0.09,0.30, 0.1, 0.1],'String','relaxation',...
    'Interpreter','latex','FontSize',fS,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Linestyle','none');
annotation('textarrow',[0.18 0.225],[0.48 0.37],'String','feedback',...
    'Interpreter','latex','FontSize',fS,'Linewidth',lW,'HeadStyle','plain');

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
set(gca,'XTick',[],'YTick',[]);



