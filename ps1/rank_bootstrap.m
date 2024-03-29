function [percent, efron] = rank_bootstrap(thetahat, Rhat, B)
% Set parameters
n = size(thetahat,1); % dimension of theta vector

% Initialize matrix to hold B roots
Root_dist = zeros(n, B);

% Run the boostrap
for b = 1:B
    % Draw theta in bootstrap world
    thetastar = randn(n,1) + thetahat;
        
    % Construct rankings
    [~,rankposition] = sort(thetastar);
    rankstar = (1:n)';
    rankstar(rankposition) = rankstar;

    % Compute root
    rootstar = sqrt(n).* (rankstar - Rhat);
    Root_dist(:,b) = rootstar;
end

% Compute quantiles of root bootstrap distribution
q = prctile(Root_dist,[2.5 97.5],2);

% Compute bootstrap percentile and Efron confidence intervals
percent_lower = Rhat - (1/sqrt(n)) * q(:,2);
percent_upper = Rhat - (1/sqrt(n)) * q(:,1);
percent = [percent_lower percent_upper];

efron_lower = Rhat + (1/sqrt(n)) * q(:,1);
efron_upper = Rhat + (1/sqrt(n)) * q(:,2);
efron = [efron_lower efron_upper];
end
