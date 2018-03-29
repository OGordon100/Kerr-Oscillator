clearvars;
close all;

%% Settings
                        % Other variables
dp = 0;
ep = 1;
kappa = 1;
K = -0.5;
ep_vary = 0:0.1:2;

parfor ep_all = 1:length(ep_vary)
    ep = ep_vary(ep_all);
    dens8_function(K,dp,ep,kappa);
    disp(['done ',num2str(dp)])
end



beep