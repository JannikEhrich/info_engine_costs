%PLOT_FEEDBACK_COOLING_PERFORMANCE plots output work and information flow
%(= input work) for teh feedback cooling model. It also plots total work
%(= entropy production) per time step and efficiency
%
% OUTPUTS:
%  vertical figure with 4 panels plotting output work, input work, total
%  work, and efficiency
%
% author:  JEhrich
% version: 0.0 (2022-02-25)
% changes: -
clear
close all
clc
% set font size, line width, and marker size
fS = 16;
lW = 2.7;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% model parameters
% time step intervals
ts_vec = [1 0.3 0.01];
% measurement precision
s2_vec = logspace(-4,3,1E3);

%% output work and information flow
W_fb = nan(length(s2_vec),length(ts_vec));
I_np_Z = nan(length(s2_vec),length(ts_vec));

for ii = 1:length(ts_vec)
    ts = ts_vec(ii);
    for jj = 1:length(s2_vec)
        s2 = s2_vec(jj);
        
        W_fb(ii,jj) = -((exp(-2*ts) - 1)*(s2 - 1))/2;
        I_np_Z(ii,jj) = -log(s2)/2 + log(exp(2*ts) + s2 - 1)/2 - ts;
    end
end

%% plot output work
colors = hsv(length(ts_vec));
figure('Position',[400,1000,560,900]);
ax1 = axes('Position',[0.13 0.77 0.77 0.2]);
for ii = 1:length(ts_vec)
    semilogx(nan,nan,'Color',colors(ii,:),'linewidth',lW);
    hold on;
end
semilogx([1E-12,1E12],[0,0],'--','color',[1,1,1]*0.5,'linewidth',1);
semilogx([1,1],[-10,10],'--','color',[1,1,1]*0.5,'linewidth',1);
for ii = 1:length(ts_vec)
    semilogx(s2_vec,W_fb(ii,:)/ts_vec(ii),'Color',colors(ii,:),'linewidth',lW);
end
set(gca,'FontSize',fS);
% xlabel('$\sigma^2$','Interpreter','latex');
% ylabel('$\left\langle \delta w_\mathrm{fb}\right\rangle/t_s$','Interpreter','latex')
title('feedback work $\left\langle \delta w_\mathrm{fb}\right\rangle/t_s$','FontWeight','normal','Interpreter','latex');
lgd = legend({...
    ['$' num2str(ts_vec(1)) '$'],...
    ['$' num2str(ts_vec(2)) '$'],...
    ['$' num2str(ts_vec(3)) '$']},'Location','NorthWest');
legend boxoff 
lgd.Title.String = '$t_s$';
lgd.FontSize = fS;
set(gca,'XTickLabels',[]);
axis([1E-2,1E1,-1.3,2.5]);

%% plot input work
ax2 = axes('Position',[0.13 0.53 0.77 0.2]);
for ii = 1:length(ts_vec)
    semilogx(nan,nan,'Color',colors(ii,:),'linewidth',lW);
    hold on;
end
semilogx([1E-12,1E12],[0,0],'--','color',[1,1,1]*0.5,'linewidth',1);
semilogx([1,1],[-10,10],'--','color',[1,1,1]*0.5,'linewidth',1);
for ii = 1:length(ts_vec)
    semilogx(s2_vec,I_np_Z(ii,:)/ts_vec(ii),'Color',colors(ii,:),'linewidth',lW);
end
set(gca,'FontSize',fS);
% xlabel('$\sigma^2$','Interpreter','latex');
% ylabel('$\lim\limits_{k \to \infty} \delta I^\mathrm{nonpred}_Z(k)/t_s$','Interpreter','latex')
title('control work $\lim_{k \to \infty} \delta I^\mathrm{nonpred}_Z(k)/t_s$','FontWeight','normal','Interpreter','latex');
% lgd = legend({...
%     ['$' num2str(ts_vec(1)) '$'],...
%     ['$' num2str(ts_vec(2)) '$'],...
%     ['$' num2str(ts_vec(3)) '$']},'Location','NorthWest');
% legend boxoff 
% lgd.Title.String = '$t_s$';
% lgd.FontSize = fS;
set(gca,'XTickLabels',[]);
axis([1E-2,1E1,-1.3,6.4]);

%% plot total work
ax3 = axes('Position',[0.13 0.29 0.77 0.2]);
for ii = 1:length(ts_vec)
    semilogx(nan,nan,'Color',colors(ii,:),'linewidth',lW);
    hold on;
end
% semilogx([1E-12,1E12],[0,0],'--','color',[1,1,1]*0.5,'linewidth',1);
semilogx([1,1],[-10,100],'--','color',[1,1,1]*0.5,'linewidth',1);
for ii = 1:length(ts_vec)
    semilogx(s2_vec,(W_fb(ii,:)+I_np_Z(ii,:))/ts_vec(ii),'Color',colors(ii,:),'linewidth',lW);
end
set(gca,'FontSize',fS);
% xlabel('$\sigma^2$','Interpreter','latex');
% ylabel('$\lim\limits_{k \to \infty} \left[\left\langle \delta w_\mathrm{fb}\right\rangle +\delta I^\mathrm{nonpred}_Z(k)\right]/t_s$','Interpreter','latex')
title('total work','FontWeight','normal','Interpreter','latex');
% lgd = legend({...
%     ['$' num2str(ts_vec(1)) '$'],...
%     ['$' num2str(ts_vec(2)) '$'],...
%     ['$' num2str(ts_vec(3)) '$']},'Location','NorthWest');
% legend boxoff 
% lgd.Title.String = '$t_s$';
% lgd.FontSize = fS;
set(gca,'XTickLabels',[]);
axis([1E-2,1E1,0,8]);




%% plot efficiency
ax4 = axes('Position',[0.13 0.05 0.77 0.2]);
for ii = 1:length(ts_vec)
    semilogx(nan,nan,'Color',colors(ii,:),'linewidth',lW);
    hold on;
end
semilogx([1E-12,1E12],[0,0],'--','color',[1,1,1]*0.5,'linewidth',1);
semilogx([1,1],[-10,10],'--','color',[1,1,1]*0.5,'linewidth',1);
for ii = 1:length(ts_vec)
    semilogx(s2_vec,-W_fb(ii,:)./I_np_Z(ii,:),'Color',colors(ii,:),'linewidth',lW);
end
set(gca,'FontSize',fS);
xlabel('$\sigma^2$','Interpreter','latex');
% ylabel('$\lim\limits_{k \to \infty} \left\langle \delta w_\mathrm{fb}\right\rangle/\delta I^\mathrm{nonpred}_Z(k)$','Interpreter','latex')
title('efficiency $\eta = -\lim_{k \to \infty} \left\langle \delta w_\mathrm{fb}\right\rangle/\delta I^\mathrm{nonpred}_Z(k)$','FontWeight','normal','Interpreter','latex');
% lgd = legend({...
%     ['$' num2str(ts_vec(1)) '$'],...
%     ['$' num2str(ts_vec(2)) '$'],...
%     ['$' num2str(ts_vec(3)) '$']},'Location','NorthWest');
% legend boxoff 
% lgd.Title.String = '$t_s$';
% lgd.FontSize = fS;
axis([1E-2,1E1,0,1]);

%% export
saveas(gcf, '../../doc/feedback_cooling_performance.eps','epsc')
