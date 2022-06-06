%PLOT_FEEDBACK_COOLING_PERFORMANCE plots output work and information flow
%(= input work) for teh feedback cooling model. It also plots total work
%(= entropy production) per time step and efficiency
%
% OUTPUTS:
%  vertical figure with 4 panels plotting output work, input work, total
%  work, and efficiency
%
% author:  JEhrich
% version: 1.6 (2022-06-06)
% changes: changed output figure to a 2x2 grid
% power"
clear
close all
clc
% set font size, line width, and marker size
fS = 20;
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
W_c_min = nan(length(s2_vec),length(ts_vec));

for ii = 1:length(ts_vec)
    ts = ts_vec(ii);
    for jj = 1:length(s2_vec)
        s2 = s2_vec(jj);
        
        W_fb(ii,jj) = -((exp(-2*ts) - 1)*(s2 - 1))/2;
        W_c_min(ii,jj) = -log(s2)/2 + log(exp(2*ts) + s2 - 1)/2 - ts;
    end
end

%% plot output work
colors = hsv(length(ts_vec));
lStyle{1} = '-';
lStyle{2} = '-.';
lStyle{3} = '--';
figure('Position',[400,1000,900,550]);
ax1 = axes('Position',[0.07 0.57 0.41 0.405]);
for ii = 1:length(ts_vec)
    semilogx(nan,nan,'Color',colors(ii,:),'LineStyle',lStyle{ii},'linewidth',lW);
    hold on;
end
semilogx([1E-12,1E12],[0,0],'--','color',[1,1,1]*0.5,'linewidth',1);
semilogx([1,1],[-10,10],'--','color',[1,1,1]*0.5,'linewidth',1);
for ii = 1:length(ts_vec)
    semilogx(s2_vec,W_fb(ii,:)/ts_vec(ii),'LineStyle',lStyle{ii},'Color',colors(ii,:),'linewidth',lW);
end

% add marker
ts_m = 0.3;
s2_m = 0.1;
W_fb_m = -((exp(-2*ts_m) - 1)*(s2_m - 1))/2;
W_c_min_m = -log(s2_m)/2 + log(exp(2*ts_m) + s2_m - 1)/2 - ts_m;
plot(s2_m,W_fb_m/ts_m,'xk','linewidth',lW,'markerSize',mS);

set(gca,'FontSize',fS);
ylabel('$\left\langle w^\mathrm{fb}_k\right\rangle/t_s$','FontWeight','normal','Interpreter','latex');
lgd = legend({...
    ['$' num2str(ts_vec(1)) '$'],...
    ['$' num2str(ts_vec(2)) '$'],...
    ['$' num2str(ts_vec(3)) '$']},'Location','NorthWest');
legend boxoff 
lgd.Title.String = '$t_s$';
lgd.FontSize = fS;
set(gca,'XTickLabels',[]);
axis([1E-2,1E1,-1.3,2.5]);
text(3.2E-3,2.5,'(a)','interpreter','latex','FontSize',fS+2);

%% plot control work
ax2 = axes('Position',[0.57 0.57 0.41 0.405]);
semilogx([1E-12,1E12],[0,0],'--','color',[1,1,1]*0.5,'linewidth',1);
hold on;
semilogx([1,1],[-10,10],'--','color',[1,1,1]*0.5,'linewidth',1);
for ii = 1:length(ts_vec)
    semilogx(s2_vec,W_c_min(ii,:)/ts_vec(ii),'Color',colors(ii,:),'LineStyle',lStyle{ii},'linewidth',lW);
end
set(gca,'FontSize',fS);
ylabel('$\left\langle w^\mathrm{c}_k\right\rangle_\mathrm{min}/t_s$','FontWeight','normal','Interpreter','latex');
set(gca,'XTickLabels',[]);
axis([1E-2,1E1,-1.3,6.4]);
text(3.2E-3,6.4,'(b)','interpreter','latex','FontSize',fS+2);
% add marker
plot(s2_m,W_c_min_m/ts_m,'xk','linewidth',lW,'markerSize',mS);

%% plot total work
ax3 = axes('Position',[0.07 0.1 0.41 0.405]);
semilogx([1,1],[-10,100],'--','color',[1,1,1]*0.5,'linewidth',1);
hold on;
for ii = 1:length(ts_vec)
    semilogx(s2_vec,(W_fb(ii,:)+W_c_min(ii,:))/ts_vec(ii),'Color',colors(ii,:),'LineStyle',lStyle{ii},'linewidth',lW);
end
set(gca,'FontSize',fS);
xlabel('$\sigma^2$','Interpreter','latex');
ylabel('total power','FontWeight','normal','Interpreter','latex');
axis([1E-2,1E1,0,8]);
text(3.2E-3,8,'(c)','interpreter','latex','FontSize',fS+2);
% add marker
plot(s2_m,(W_fb_m+W_c_min_m)/ts_m,'xk','linewidth',lW,'markerSize',mS);


%% plot efficiency
ax4 = axes('Position',[0.57 0.1 0.41 0.405]);
semilogx([1,1],[-10,10],'--','color',[1,1,1]*0.5,'linewidth',1);
hold on;
for ii = 1:length(ts_vec)
    semilogx(s2_vec,-W_fb(ii,:)./W_c_min(ii,:),'Color',colors(ii,:),'LineStyle',lStyle{ii},'linewidth',lW,'MarkerSize',mS);
end
set(gca,'FontSize',fS);
xlabel('$\sigma^2$','Interpreter','latex');
ylabel('efficiency $\eta_\mathrm{inf}$','FontWeight','normal','Interpreter','latex');
axis([1E-2,1E1,0,1]);
text(3.2E-3,1,'(d)','interpreter','latex','FontSize',fS+2);
% add marker
plot(s2_m,-W_fb_m/W_c_min_m,'xk','linewidth',lW,'markerSize',mS);

%% export
saveas(gcf, '../../doc/feedback_cooling_performance.eps','epsc')
