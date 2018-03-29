function [pdot] = p_dot_cal(p_1d,maxsize,nall,n_all,dp,ep,K,kappa)

% Set up arrays
pdot_pad = zeros(maxsize+2,maxsize+2);

p = reshape(p_1d,[maxsize,maxsize]);
p_pad = padarray(p,[1,1]);

% For all n,n'
for posi = 1:maxsize^2
    % Extract position
    n = nall(posi);
    n_ = n_all(posi);
    
    % Calculate pdot from equation
    pdot_pad(n+2,n_+2) = (p_pad(n+2,n_+2))*( (1i * dp * (n_ - n))...
        + 1i * K/2 * (n_ *(n_ -1) - n * (n-1)) - 0.5*kappa * (n + n_)) ...
        + p_pad(n+3,n_+3)*( kappa*sqrt(n+1)*sqrt(n_+1))...
        + (p_pad(n+1,n_+2)*( ep*sqrt(n)))...
        - p_pad(n+3,n_+2)*( ep*sqrt(n+1))...
        - p_pad(n+2,n_+3)*( ep*sqrt(n_+1))...
        +p_pad(n+2,n_+1)*( ep*sqrt(n_));
end

% Strip padding
pdot = pdot_pad(2:end-1,2:end-1);
pdot = reshape(pdot,[maxsize^2,1]);
end