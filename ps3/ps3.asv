%% ECO 539b: Problem Set 3

%% Regression 1 (a)

% Import data
data = readtable('reg1.csv');
N = size(data,1); % Number of observations
M = length(unique(data{:,'prov'})); % Number of clusters

% Run OLS
y1 = data{:, 'ldeaths'};
X1 = [ones(N,1) data{:,4:end}];
beta1 = X1\y1;
eps1 = y1 - X1*beta1;
se1 = se_ehw(X1,eps1,2:3);


% Run Bootstrap for part (a) (non-clustered standard errors)
B = 50000; % Number of bootstrap iterations
beta1a_dist = zeros(B,2);
t1a_dist = zeros(B,2);
for b = 1:B
    % Draw bootstrap data and run OLS
    datastar = data(ceil(rand(N,1)*(N)),:);
    y1star = datastar{:,'ldeaths'};
    X1star = [ones(N,1) datastar{:,4:end}];
    beta1star = X1star\y1star;
    beta1a_dist(b,:) = beta1star(2:3)';
    eps1star = y1star - X1star*beta1star;

    % Calculate EHW standard error and t-stat
    stderr_ehw_star = se_ehw(X1star,eps1star, 2:3);
    tstat_ehw_star = (beta1star(2:3) - beta1(2:3))./stderr_ehw_star;
    t1a_dist(b,:) = tstat_ehw_star';
end

% Calculate bootstrap standard errors and percentile-t CI 
bootstrapse_1a = std(beta1a_dist)';
ptile = prctile(t1a_dist,[2.5 97.5]);
percentilet_1a = zeros(2);
percentilet_1a(:,1) = beta1(2:3) - (se1 .* ptile(2,:)');
percentilet_1a(:,2) = beta1(2:3) - (se1 .* ptile(1,:)');

% Display results
disp('bootstrap standard errors, reg 1, no cluster')
disp(bootstrapse_1a)
disp('bootstrap 95% CI, reg 1, no cluster')
disp(percentilet_1a)


%% Regression 1 (b)

% Run Bootstrap for part (b) (clustered standard errors)
    % Set up cluster index - stores indices of each cluster in real data
[~,~,clusterids1] = unique(data{:,'prov'});
se1_cluster = se_cluster(X1,eps1,2:3,clusterids1); % cluster SEs on real data
clusterindex = cell(M,1);
for m = 1:M
clusterindex{m} = find(clusterids1 == m);
end

    % Run bootstrap
beta1b_dist = zeros(B,2);
t1b_dist = zeros(B,2);
for b = 1:B
        % Draw M clusters to form bootstrap data
    clusterstar = ceil(rand(M,1)*M);
    datastar = [];
    for m = 1:M
        si = size(data(clusterindex{clusterstar(m)},:),1);
        clusteridstar = m*ones(si,1);
        newcluster = [data(clusterindex{clusterstar(m)},:) table(clusteridstar)];
        datastar = [datastar; newcluster];
    end
    clusteridstar = datastar{:,end};
    Nstar = size(datastar,1);
    y1star = datastar{:,'ldeaths'};
    X1star = [ones(Nstar,1) datastar{:,4:end-1}];

        % Run OLS
    beta1star = X1star\y1star;
    beta1b_dist(b,:) = beta1star(2:3)';
    eps1star = y1star - X1star*beta1star;

        % Calculate cluster-robust standard errors and t-stat
    stderr_cluster_star = se_cluster(X1star, eps1star, 2:3, clusteridstar);
    tstat_cluster_star = (beta1star(2:3) - beta1(2:3))./stderr_cluster_star;
    t1b_dist(b,:) = tstat_cluster_star';
    disp(b)
end

% Calculate bootstrap standard errors and percentile-t CI 
bootstrapse_1b = std(beta1b_dist)';
ptile = prctile(t1b_dist,[2.5 97.5]);
percentilet_1b = zeros(2);
percentilet_1b(:,1) = beta1(2:3) - (se1_cluster .* ptile(2,:)');
percentilet_1b(:,2) = beta1(2:3) - (se1_cluster .* ptile(1,:)');

% Display results
disp('bootstrap standard errors, reg 1, cluster')
disp(bootstrapse_1b)
disp('bootstrap 95% CI, reg 1, cluster')
disp(percentilet_1b)


%% Regression 2(a)
% Import data
data = readtable('reg2.csv');
N = size(data,1); % Number of observations
M = length(unique(data{:,'prov'})); % Number of clusters

% Run OLS
y2 = data{:, 'ldeaths'};
X2 = [ones(N,1) data{:,4:end}];
beta2 = X2\y2;
eps2 = y2 - X2*beta2;
se2 = se_ehw(X2,eps2,2:3);


% Run Bootstrap for part (a) (non-clustered standard errors)
B = 50000; % Number of bootstrap iterations
beta2a_dist = zeros(B,2);
t2a_dist = zeros(B,2);
for b = 1:B
    % Draw bootstrap data and run OLS
    datastar = data(ceil(rand(N,1)*(N)),:);
    y2star = datastar{:,'ldeaths'};
    X2star = [ones(N,1) datastar{:,4:end}];
    beta2star = X2star\y2star;
    beta2a_dist(b,:) = beta2star(2:3)';
    eps2star = y2star - X2star*beta2star;

    % Calculate EHW standard error and t-stat
    stderr_ehw_star = se_ehw(X2star,eps2star, 2:3);
    tstat_ehw_star = (beta2star(2:3) - beta2(2:3))./stderr_ehw_star;
    t2a_dist(b,:) = tstat_ehw_star';
end

% Calculate bootstrap standard errors and percentile-t CI 
bootstrapse_2a = std(beta2a_dist)';
ptile = prctile(t2a_dist,[2.5 97.5]);
percentilet_2a = zeros(2);
percentilet_2a(:,1) = beta2(2:3) - (se2 .* ptile(2,:)');
percentilet_2a(:,2) = beta2(2:3) - (se2 .* ptile(1,:)');

% Display results
disp('bootstrap standard errors, reg 2, no cluster')
disp(bootstrapse_2a)
disp('bootstrap 95% CI, reg 2, no cluster')
disp(percentilet_2a)

%% Regression 1 (b)

% Run Bootstrap for part (b) (clustered standard errors)
    % Set up cluster index - stores indices of each cluster in real data
[~,~,clusterids2] = unique(data{:,'prov'});
se2_cluster = se_cluster(X2,eps2,2:3,clusterids2); % cluster SEs on real data
clusterindex = cell(M,1);
for m = 1:M
clusterindex{m} = find(clusterids2 == m);
end

    % Run bootstrap
beta2b_dist = zeros(B,2);
t2b_dist = zeros(B,2);
for b = 1:B
        % Draw M clusters to form bootstrap data
    clusterstar = ceil(rand(M,1)*M);
    datastar = [];
    for m = 1:M
        si = size(data(clusterindex{clusterstar(m)},:),1);
        clusteridstar = m*ones(si,1);
        newcluster = [data(clusterindex{clusterstar(m)},:) table(clusteridstar)];
        datastar = [datastar; newcluster];
    end
    clusteridstar = datastar{:,end};
    Nstar = size(datastar,1);
    y2star = datastar{:,'ldeaths'};
    X2star = [ones(Nstar,1) datastar{:,4:end-1}];

        % Run OLS
    beta2star = X2star\y2star;
    beta2b_dist(b,:) = beta2star(2:3)';
    eps2star = y2star - X2star*beta2star;

        % Calculate cluster-robust standard errors and t-stat
    stderr_cluster_star = se_cluster(X2star, eps2star, 2:3, clusteridstar);
    tstat_cluster_star = (beta2star(2:3) - beta2(2:3))./stderr_cluster_star;
    t2b_dist(b,:) = tstat_cluster_star';
end

% Calculate bootstrap standard errors and percentile-t CI 
bootstrapse_2b = std(beta2b_dist)';
ptile = prctile(t2b_dist,[2.5 97.5]);
percentilet_2b = zeros(2);
percentilet_2b(:,1) = beta2(2:3) - (se2_cluster .* ptile(2,:)');
percentilet_2b(:,2) = beta2(2:3) - (se2_cluster .* ptile(1,:)');

% Display results
disp('bootstrap standard errors, reg 2, cluster')
disp(bootstrapse_2b)
disp('bootstrap 95% CI, reg 2, cluster')
disp(percentilet_2b)








