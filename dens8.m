%% Oli Gordon & Matthew Edmondson QCD Project
clearvars
close all

%% Settings
maxsize = 25;                     % Grid size
startn = 0;                       % Starting points
startn_ = 0;

maxtime = 15;                     % Maximum time
ode_inc = 1;                      % Timestep for ODE to calculate over
odesteps = 100;                   % Number of calculations per ODE
saving = 0;                       % Choose to save 0|1 (no|yes)

K = -0.5;                         % Other variables
dp = 0;
ep = 1;
kappa = 1;

%% Setup
timesteps = maxtime/ode_inc;

% Preallocate arrays, including +2 of padding
p_pad = zeros(maxsize+2,maxsize+2);
pdot_pad = zeros(maxsize+2,maxsize+2);
traces = zeros(1,(odesteps*maxtime)+1);
traces2 = zeros(1,(odesteps*maxtime)+1);
adag_a = zeros(1,(odesteps*maxtime)+1);
g20 = zeros(1,(odesteps*maxtime)+1);
amp = zeros(1,(odesteps*maxtime)+1);
timeaxis_calc = 0:ode_inc:(ode_inc*timesteps);
timeaxis_plot = zeros(1,(odesteps*maxtime)+1);
timeaxis_plot(end) = maxtime;

% Calculate raising and lowering operators in matrix form
a = diag(sqrt(1:maxsize-1),1);
a_dag = diag(sqrt(1:maxsize-1),-1);

% Create list of <n,n'> to call through
row = 0:maxsize-1;
col = 0:maxsize-1;
[X,Y] = meshgrid(row,col);
nall = reshape(Y,[1,maxsize^2]);
n_all = reshape(X,[1,maxsize^2]);

% Create and explore initial p
p_pad(startn_+2,startn+2) = 1;
p = p_pad(2:end-1,2:end-1);

%% Main Calculations
% Split calls of the ODE into parts, limiting matrix size to save RAM and
% therefore allow us to investigate longer timescales and higher
% resolutions :)
for timers = 1:timesteps
    %% Integration
    % Determine time period ODE is performed over
    tspan = timeaxis_calc(timers):ode_inc/odesteps:...
        timeaxis_calc(timers+1);
    
    % Perform ODE
    p = reshape(p,[1,maxsize^2]);
    [t,p_new] = ode45(@(t,p) ...
        p_dot_cal(p,maxsize,nall,n_all,dp,ep,K,kappa),tspan, p);
    
    % Store outputs of ODE
    timeaxis_plot((((timers-1)*odesteps)+1):(((timers)*odesteps))) ...
        = t(1:end-1);
    
    %% Calculations
    for diag_loop = 1:size(p_new,1)-1
        % p (should be =1), p^2 (0-1, expect to decrease over time)
        p_segment = reshape(p_new(diag_loop,:),[maxsize,maxsize])';
        traces(((timers-1)*odesteps)+diag_loop) = ...
            sum(diag(p_segment));
        traces2(((timers-1)*odesteps)+diag_loop) = ...
            sum(diag((p_segment^2)));
        adag_a(((timers-1)*odesteps)+diag_loop) = ...
            sum(diag(a_dag*a*p_segment));
        g20(((timers-1)*odesteps)+diag_loop) = ...
            sum(diag((a_dag*a_dag*a*a*p_segment)/...
            (adag_a(((timers-1)*odesteps)+diag_loop))^2));
        amp(((timers-1)*odesteps)+diag_loop) = ...
            sum(diag(a*p_segment));
    end
    
    % Reshape last value of p to a 2D matrix for the next time around
    p = reshape(p_new(end,:),[maxsize,maxsize])';
    
end

% Get calculations for very last result
traces(end) = sum(diag(p));
traces2(end) = sum(diag((p^2)));
adag_a(end) = sum(diag(a_dag*a*p));
g20(end) = sum(diag(a_dag*a_dag*a*a*p))/((adag_a(end))^2);
amp(end) = sum(diag(a*p));

% Make semiclassical approximation
A = K^2;
B = 2*dp*K;
C = (kappa^2)/4 + dp^2;
D = -(ep^2);
p = [A B C D];
semiclass_all = roots(p);
[~, idx] = sort(abs(imag(semiclass_all)));
semiclass_sorted = semiclass_all(idx);
semiclass = real(semiclass_sorted(1));

%% Plot Graphs
set(0,'defaulttextInterpreter','latex');

% Plot Figures
figure()
subplot(2,2,1)
plot(timeaxis_plot,abs(traces))
ylim([0,1.2])
xlabel('Time (s)')
ylabel('Tr$[\rho]$')
%set(gca,'fontsize',12)

subplot(2,2,2)
plot(timeaxis_plot,abs(traces2))
ylim([0,1.2])
xlabel('Time (s)')
ylabel('Tr$[\rho^2]$')
%set(gca,'fontsize',12)

subplot(2,2,3)
hold on
plot(timeaxis_plot,abs(adag_a))
plot(timeaxis_plot,semiclass.*ones(1,length(timeaxis_plot)),'r--')
xlabel('Time (s)')
ylabel('$\langle a^\dagger a\rangle$','FontSize',12)
legend({'Numerical Integration';'Semi-Classical Approximation'},...
    'Location','best')
set(gca,'fontsize',12)

subplot(2,2,4)
plot(timeaxis_plot,abs(g20))
xlabel('Time (s)')
ylabel(['$\frac{\langle a^\dagger a^\dagger aa \rangle}'...
    '{{\langle a^\dagger a \rangle}^2}$ '],'FontSize',16)
ylim([0 inf])
%set(gca,'fontsize',12)

% subplot(2,3,5)
% plot(timeaxis_plot,abs(amp))
% xlabel('Time (Arbitrary Units)')
% ylabel('$\langle a\rangle$','FontSize',12)


set(0,'defaulttextInterpreter','tex');

% Save Output
if saving == 1
    save(['Data/Outputs_K=',num2str(K),'_dp=',num2str(dp),...
        '_ep=',num2str(ep),'_kappa=',num2str(kappa),'.mat'],...
        'traces','traces2','adag_a','g20','timeaxis_plot','semiclass')
end