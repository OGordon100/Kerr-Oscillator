%% Oli Gordon & Matthew Edmondson QCD Project
clearvars
close all
set(0,'defaulttextInterpreter','latex');

%% Settings

% Pick thing to vary
varying = 'ep';
in_range = 2:-0.1:0;

% Pick other constants
dp = 0;
ep = 1;
kappa = 1;
K = -0.5;

% Preallocate matrices
trace_long = zeros(1,length(in_range));
trace2_long = zeros(1,length(in_range));
adag_a_long = zeros(1,length(in_range));
semi_class = zeros(1,length(in_range));
g20_long = zeros(1,length(in_range));
in_leg = cellstr(num2str(in_range'))';

%% Plot
hold all
for in_loop = 1:length(in_range)
    % Read in
    assignin('base',varying,in_range(in_loop))
    filename = ['Data/Outputs_K=',num2str(K),'_dp=',num2str(dp),...
        '_ep=',num2str(ep),'_kappa=',num2str(kappa),'.mat'];
    load(filename);
    
    % Store in vectors
    trace_long(in_loop) = abs(traces(end));
    trace2_long(in_loop) = abs(traces2(end));
    adag_a_long(in_loop) = abs(adag_a(end));
    g20_long(in_loop) = abs(g20(end));
    %semi_class(in_loop) = semiclass;
    
    plotter(:,in_loop) = abs(g20c);
    
    % Add to plot
    figure(1)
    hold all
    plot(timeaxis_plot,abs(traces))
    figure(2)
    hold all
    plot(timeaxis_plot,abs(traces2))
    figure(3)
    hold all
    plot(timeaxis_plot,abs(adag_a))
    figure(4)
    hold all
    plot(timeaxis_plot,abs(g20))
end

%% Niceties

set(0,'defaulttextInterpreter','latex');

if strcmp(varying,'K') == 1
    varying_lgd = 'K/\kappa';
elseif strcmp(varying,'dp') == 1
    varying_lgd = '\Delta_\rho/\kappa';
elseif strcmp(varying,'ep') == 1
    varying_lgd = '\epsilon_\rho/\kappa';
elseif strcmp(varying,'kappa') == 1
    varying_lgd = '\kappa';
end

% % Plot Tr[p]
% figure(1)
% lg1 = legend(in_leg);
% title(lg1,varying_lgd)
% xlabel('Time (s)')
% ylabel('Tr$[\rho]$')
% ylim([0,1.2])
% 
% % Plot Tr[p^2]
% figure(2)
% lg2 = legend(in_leg);
% title(lg2,varying_lgd)
% xlabel('Time (s)')
% ylabel('Tr$[\rho^2]$')
% ylim([0,1.2])
% 
% % Plot Tr[adag a]
% figure(3)
% lg3 = legend(in_leg);
% title(lg3,varying_lgd)
% xlabel('Time (s)')
% ylabel('$\langle a^\dagger a\rangle$','FontSize',12)
% 
% % Plot Tr[g^2]
% figure(4)
% lg4 = legend(in_leg);
% title(lg4,varying_lgd)
% xlabel('Time (s)')
% ylabel(['$\frac{\langle a^\dagger a^\dagger aa \rangle}'...
%     '{{\langle a^\dagger a \rangle}^2}$ '],'FontSize',16)

% Plot long term behaviors
figure(5)
plot(in_range,trace_long)
xlabel(varying_lgd,'Interpreter','tex')
ylabel('Tr$[\rho]$')
ylim([0,1.2])

figure(6)
plot(in_range,trace2_long)
xlabel(varying_lgd,'Interpreter','tex')
ylabel('Tr$[\rho^2]$')
ylim([0,1.2])
set(gca,'fontsize',12)

figure(7)
set(gca,'fontsize',12)
hold on
plot(in_range,adag_a_long,'b')
plot(in_range,semi_class,'r--')
xlabel(varying_lgd,'Interpreter','tex')
ylabel('$\langle a^\dagger a\rangle$','FontSize',16)
%legend({'Numerical';'Semi-Classical'},...
%    'Location','best')
% 
% 
% subplot(1,2,2)
% set(gca,'fontsize',12)
% plot(in_range,g20_long)
% set(gca,'fontsize',12)
% xlabel(varying_lgd,'Interpreter','tex')
% ylabel(['$\frac{\langle a^\dagger a^\dagger aa \rangle}'...
%     '{{\langle a^\dagger a \rangle}^2}$ '],'FontSize',18)
% 
% 
% Save figures
figs = [figure(1),figure(2),figure(3),figure(4)];
%savefig(figs,'Data/K_vary')

set(0,'defaulttextInterpreter','tex');

figure(6)
figure(7)