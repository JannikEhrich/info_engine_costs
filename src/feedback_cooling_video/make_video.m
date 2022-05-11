%SCRIPTNAME does...
%
% OUTPUTS:
%  outputs XXX
%
% author:  JEhrich
% version: 0.0 (2022-05-09)
% changes: -
clear
close all
clc
% set font size, line width, and marker size
fS = 60;
lW = 8;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% parameters
% measurement error
s2 = 0.1;
% total time interval
ts = 0.4;
% measurement time interval
tau = 0.001;
% video length for control
t_vid_c = 2;
% video length for relaxation
t_vid_r = 4;
% framerate
fps = 30;
% set resolution
res = [1920 1080];
% second for fade
t_fade = 1.5;

% measurement: initial stiffness
k0 = 1/((s2 - 1)*exp(-2*ts) + 1)
% measurement: final stiffness
k1 = 1/s2

% x_axis
x = linspace(-5,5,1E3);
ax = [min(x),max(x),0,50];

% factor for probability scaling
f = 35;

v = VideoWriter('test','MPEG-4');
open(v)
tic

%% fade to control step
n_fade = t_fade*fps;
V0 = 1/2*x.^2;
V1 = k0/2*x.^2;
p = 1/sqrt(2*pi*1/k0)*exp(-x.^2/2*k0);

for ii = 1:n_fade
    ii/n_fade
    figure('Position',[200,100,res]);
    %set(gcf,'Visible', 'off')
    p0 = plot(x,V0,'color',[1,0,0,1-(ii-1)/n_fade],'linewidth',lW);
    hold on;
    p1 = plot(x,V1,'color',[1,0,0,(ii-1)/n_fade],'linewidth',lW);
    inBetween = [f*p, zeros(size(p))];
    pF = patch([x,fliplr(x)], inBetween, 'b','EdgeColor','none');
    alpha(pF,.2)
    plot(x,p*f,'b','linewidth',lW)
    xlabel('$x-z$','Interpreter','latex')
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    axis(ax)
    text(0,47,'relaxation $\to$ control','HorizontalAlignment', 'center','Interpreter','latex','FontSize',fS)
    set(gca,'FontSize',fS)
    F = getframe(gcf);
    [A, ~] = frame2im(F);
    writeVideo(v,A);
end


%% control step
t_c = linspace(0,tau,t_vid_c*fps);
for ii = 1:length(t_c)
    t = t_c(ii)
    k = k0 + (ii-1)*(k1-k0)/(t_vid_c*fps);
    V = k/2*x.^2;

    c = 1/k;
    p = 1/sqrt(2*pi*c)*exp(-x.^2/2/c);

    figure('Position',[200,100,res]);
    set(gcf,'Visible', 'off')
    plot(x,V,'r','linewidth',lW);
    hold on;
    inBetween = [f*p, zeros(size(p))];
    pF = patch([x,fliplr(x)], inBetween, 'b','EdgeColor','none');
    alpha(pF,.2)
    plot(x,p*f,'b','linewidth',lW)
    xlabel('$x-z$','Interpreter','latex')
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    axis(ax)
    set(gca,'FontSize',fS)
    text(0,47,'control','HorizontalAlignment', 'center', 'Interpreter','latex','FontSize',fS,'color','b')
    F = getframe(gcf);
    [A, ~] = frame2im(F);
    writeVideo(v,A);
    
end


%% fade to relaxation step
n_fade = t_fade*fps;
V0 = k1/2*x.^2;
V1 = 1/2*x.^2;
p = 1/sqrt(2*pi*1/k1)*exp(-x.^2/2*k1);

for ii = 1:n_fade
    ii/n_fade
    figure('Position',[200,100,res]);
    set(gcf,'Visible', 'off')
    p0 = plot(x,V0,'color',[1,0,0,1-(ii-1)/n_fade],'linewidth',lW);
    hold on;
    p1 = plot(x,V1,'color',[1,0,0,(ii-1)/n_fade],'linewidth',lW);
    inBetween = [f*p, zeros(size(p))];
    pF = patch([x,fliplr(x)], inBetween, 'b','EdgeColor','none');
    alpha(pF,.2)
    plot(x,p*f,'b','linewidth',3)
    xlabel('$x-z$','Interpreter','latex')
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    axis(ax)
    text(0,47,'control $\to$ relaxation','HorizontalAlignment', 'center','Interpreter','latex','FontSize',fS)
    set(gca,'FontSize',fS)
    F = getframe(gcf);
    [A, ~] = frame2im(F);
    writeVideo(v,A);
end

%% relaxation step
t_r = linspace(0,ts,t_vid_r*fps);
for ii = 1:length(t_r)
    t = t_r(ii)
    k = k0 + (ii-1)*(k1-k0)/(t_vid_c*fps);
    V = 1/2*x.^2;

    c = s2 + (1-s2)*(1 - exp(-2*t));
    p = 1/sqrt(2*pi*c)*exp(-x.^2/2/c);

    figure('Position',[200,100,res]);
    set(gcf,'Visible', 'off')
    plot(x,V,'r','linewidth',lW);
    hold on;
    inBetween = [f*p, zeros(size(p))];
    pF = patch([x,fliplr(x)], inBetween, 'b','EdgeColor','none');
    alpha(pF,.2)
    plot(x,p*f,'b','linewidth',lW)
    xlabel('$x-z$','Interpreter','latex')
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    axis(ax)
    set(gca,'FontSize',fS)
    text(0,47,'relaxation','HorizontalAlignment', 'center','Interpreter','latex','FontSize',fS,'color','r')
    F = getframe(gcf);
    [A, ~] = frame2im(F);
    writeVideo(v,A);
end

close(v)
toc





