%SIM_FEEDBACK_COOLING_EXPLICIT simulates explicit controller model for
%feedback cooling information engine
%
% OUTPUTS:
%  outputs eps figure showing feedback and control work and efficiency
%
% author:  JEhrich
% version: 1.2 (2022-05-20)
% changes: removed legend in panel (a), label curves directly instead
clear
close all
clc
% set font size, line width, and marker size
fS = 20;
lW = 2.0;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% system parameters
% time-scale separation factors
nu_h_vec = logspace(0,5,20);
nu_l = 1E-2;
% measurement error
s2 = 0.1;
% total time interval
ts = 0.3;
% measurement time interval
tau = 0.001;

%% simulation parameters
% integration time-step
dt = 1E-7;
% number of time steps
K = 1E4;
% number of initial steps
K_ini = 10;

%% main loop
% data structures for work mean and variance
W_fb_mean = nan(length(nu_h_vec),1);
W_c_mean = nan(length(nu_h_vec),1);
W_fb_var = nan(length(nu_h_vec),1);
W_c_var = nan(length(nu_h_vec),1);

s2_meas = nan(length(nu_h_vec),1);
s2_meas_err = nan(length(nu_h_vec),1);


tic
parfor ii = 1:length(nu_h_vec)
    ii
    nu_h = nu_h_vec(ii);
    % relaxation step propagator covariance matrix
    cxx = (-exp(2*(nu_l + 1)*(-ts + tau)) + 1 + (-2*tau + 2*ts)*nu_l^2 + (-2*tau + 2*ts)*nu_l)/(nu_l + 1)^2;
    cxz = 2*(exp(2*(nu_l + 1)*(-ts + tau))/2 + (ts - tau)*nu_l + ts - tau - 1/2)*nu_l/(nu_l + 1)^2;
    czz = 2*(-nu_l*exp(2*(nu_l + 1)*(-ts + tau))/2 + (ts - tau + 1/2)*nu_l + ts - tau)*nu_l/(nu_l + 1)^2;
    C_prop = [cxx, cxz; cxz, czz];
    x = 0;
    z = randn*sqrt(s2);
    % initial equilibration
    for kk = 1:K_ini
        [x,z,~,~,~] = sim_one_step(x,z,ts,s2,tau,nu_l,nu_h,C_prop,dt);
    end
    
    % simulate K steps
    W_fb = nan(K,1);
    W_c = nan(K,1);
    % data structures for post-measurement deviation
    d_meas = nan(K,1);
    for kk = 1:K
        [x,z,W_fb(kk),W_c(kk),d_meas(kk)] = sim_one_step(x,z,ts,s2,tau,nu_l,nu_h,C_prop,dt);
    end
    
    % calculate mean and variance of works
    W_fb_mean(ii) = mean(W_fb);
    W_c_mean(ii) = mean(W_c);
    W_fb_var(ii) = var(W_fb);
    W_c_var(ii) = var(W_c);

    % calculate variance of z around x after measurement
    s2_meas(ii) = mean(d_meas.^2);
    s2_meas_err(ii) = sqrt(var(d_meas.^2)/K)

end
toc

%% error calculation
W_fb_err = sqrt(W_fb_var/K);
W_c_err = sqrt(W_c_var/K);
eta_err = sqrt((1./W_c_mean).^2.*W_fb_err.^2 + (W_fb_mean./W_c_mean.^2).^2.*W_c_err.^2);

%% analytical expectations
W_fb_ana = -((exp(-2*ts) - 1)*(s2 - 1))/2;
W_c_ana = -log(s2)/2 + log((s2 - 1)*exp(-2*ts) + 1)/2;

%% plot convergence
figure('Position',[400,1000,560,850]);
ax1 = axes('Position',[0.13 0.69 0.77 0.29]);
%errorbar(nu_h_vec,W_c_mean/ts,W_c_err/ts,'bs','MarkerSize',mS,'lineWidth',lW);
plot(nu_h_vec,W_c_mean/ts,'bs','MarkerSize',mS,'lineWidth',lW);
hold on;
plot(nu_h_vec,W_fb_mean/ts,'rs','MarkerSize',mS,'lineWidth',lW);
plot(nu_h_vec,ones(size(nu_h_vec))*W_c_ana/ts,'b','MarkerSize',mS,'lineWidth',lW);
%errorbar(nu_h_vec,W_fb_mean/ts,W_fb_err/ts,'rs','MarkerSize',mS,'lineWidth',lW);
plot(nu_h_vec,ones(size(nu_h_vec))*W_fb_ana/ts,'r','MarkerSize',mS,'lineWidth',lW);
set(gca,'XScale','log','FontSize',fS);
set(gca,'XTick',10.^[0,1,2,3,4,5]);
set(gca,'XTicklabels',[]);
%xlabel('$\nu_\mathrm{high}$','Interpreter','latex');
ylabel('rate of work','Interpreter','latex');
%legend({'$\left\langle w^\mathrm{c} \right\rangle/t_\mathrm{s}$',...
%    '$\left\langle w^\mathrm{fb} \right\rangle/t_\mathrm{s}$'},...
%    'Location','NorthWest');
%legend boxoff 
axis([min(nu_h_vec),max(nu_h_vec),-1,4.5]);
text(1E2, 0 , '$\left\langle w^\mathrm{fb} \right\rangle/t_\mathrm{s}$',...
    'Interpreter','latex','FontSize',fS,'Color','r');
text(8E2, 3.8 , '$\left\langle w^\mathrm{c} \right\rangle/t_\mathrm{s}$',...
    'Interpreter','latex','FontSize',fS,'Color','b');
text(2E-1,4.5,'(a)','interpreter','latex','FontSize',fS+2);


% plot deviation from measurement
ax2 = axes('Position',[0.13 0.38 0.77 0.29]);
%errorbar(nu_h_vec,s2_meas,s2_meas_err,'ks','MarkerSize',mS,'lineWidth',lW);
plot(nu_h_vec,s2_meas,'ks','MarkerSize',mS,'lineWidth',lW);
hold on;
plot(nu_h_vec,ones(size(nu_h_vec))*s2,'k','MarkerSize',mS,'lineWidth',lW);
%xlabel('$\nu_\mathrm{high}$','Interpreter','latex');
ylabel('$\left\langle (x_{k} - z_{k})^2 \right\rangle$','Interpreter','latex');
set(gca,'XScale','log','FontSize',fS);
set(gca,'XTick',10.^[0,1,2,3,4,5]);
set(gca,'XTicklabels',[]);
text(3E-2, 0.14 , '$\sigma^2$',...
    'Interpreter','latex','FontSize',fS);
text(2E-1,1,'(b)','interpreter','latex','FontSize',fS+2);


% plot efficiency
ax3 = axes('Position',[0.13 0.07 0.77 0.29]);
errorbar(nu_h_vec,-W_fb_mean./W_c_mean,eta_err,'ks','MarkerSize',mS,'lineWidth',lW);
hold on;
plot(nu_h_vec,-ones(size(nu_h_vec))*W_fb_ana/W_c_ana,'k','MarkerSize',mS,'lineWidth',lW);
xlabel('$\nu_\mathrm{high}$','Interpreter','latex');
ylabel('efficiency','Interpreter','latex');
set(gca,'XScale','log','FontSize',fS);
set(gca,'XTick',10.^[0,1,2,3,4,5]);
%set(gca,'XTicklabels',[]);
axis([min(nu_h_vec),max(nu_h_vec),0,0.28]);
text(2E-1,0.28,'(c)','interpreter','latex','FontSize',fS+2);
text(1E2, 0.13 , '$-\left\langle w^\mathrm{fb} \right\rangle/\left\langle w^\mathrm{c} \right\rangle$',...
    'Interpreter','latex','FontSize',fS);
text(1E1, 0.24 , '$\eta_\mathrm{inf}$',...
    'Interpreter','latex','FontSize',fS);

% export
saveas(gcf, '../../doc/feedback_cooling_explicit_sim.eps','epsc')




%% support function simulating one time step
function [x,z,W_fb,W_c,d_meas] = sim_one_step(x0,z0,ts,s2,tau,nu_l,nu_h,C_prop,dt)
    x = x0;
    z = z0;
    % control operation, initial energy
    E0 = 1/2*(x-z)^2;
    % initial and final stiffness
    kappa0 = 1/(exp(-2*ts)*s2 - exp(-2*ts) + 1);
    kappa1 = 1/s2;
    % energy after switch to control potential
    E1 = kappa0/2*(x-z)^2;
    E2 = E1;
    % linearly interpolate between initial and final stiffness
    W_c = 0;
    for ii = 1:ceil(tau/dt)
        t = ii*dt;
        % current stiffness
        kappa = kappa0 + t/tau*(kappa1-kappa0) ;
        % new energy and control work
        E3 = kappa/2*(x-z)^2;
        W_c = W_c + E3 - E2;
        % integrate x and z
        dx = -kappa*(x-z)*dt + sqrt(2*dt)*randn;
        dz = nu_h*kappa*(x-z)*dt + sqrt(2*dt*nu_h)*randn;
        x = x + dx;
        z = z + dz;
        % energy after system movement
        E2 = kappa/2*(x-z)^2;
    end
    % boundary terms of control work
    W_c = W_c + E1 - E2;
    % energy after switch to feedback potential
    E4 = 1/2*(x-z)^2;
    % feedback work
    W_fb = E4 - E0;

    % calculate post-measurement deviation of system and controller
    d_meas = x-z;

    % relaxation, draw next states
    mu_x = ((x - z)*exp((nu_l + 1)*(-ts + tau)) + nu_l*x + z)/(nu_l + 1);
    mu_z = (-(x - z)*exp((nu_l + 1)*(-ts + tau))*nu_l + nu_l*x + z)/(nu_l + 1);
    x = mu_x + randn*sqrt(C_prop(1,1));
    z = mu_z + C_prop(1,2)/C_prop(1,1)*(x-mu_x) + randn*sqrt(C_prop(2,2)-C_prop(1,2)^2/C_prop(1,1));
    
end




