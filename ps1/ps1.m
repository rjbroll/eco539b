%% Q1

% Set parameters
M = 5000; % Number of Monte Carlo simulations
B = 1000; % Number of Bootstrap iterations
n = 10; % Dimension of theta


%% Design 1: theta = (1:10)'
    % Construct true theta
theta = (1:n)';

    % Construct true ranks
[~,Rposition] = sort(theta);
R = (1:n)';
R(Rposition) = R;

    % Initialize vector of results
design1_results = zeros(n,2,M);

    % Run the Montel Carlo simulations
for m = 1:M
        % Draw thetahat, construct Rhat
    thetahat = randn(n,1) + theta;
    [~, Rhatposition] = sort(thetahat);
    Rhat = (1:n)';
    Rhat(Rhatposition) = Rhat;

        % Run the bootstrap
    [percent, efron] = rank_bootstrap(thetahat, Rhat, B);
    design1_results(:,1,m) = (R >= percent(:,1)) & (R <= percent(:,2));
    design1_results(:,2,m) = (R >= efron(:,1)) & (R <= efron(:,2));
end

    % Get results
design1_results = mean(design1_results,3);

%% Design 2: theta = 10*(1:10)'
    % Construct true theta
theta = 10*(1:n)';

    % Construct true ranks
[~,Rposition] = sort(theta);
R = (1:n)';
R(Rposition) = R;

    % Initialize vector of results
design2_results = zeros(n,2,M);

    % Run the Monte Carlo simulations
for m = 1:M
        % Draw thetahat, construct Rhat
    thetahat = randn(n,1) + theta;
    [~, Rhatposition] = sort(thetahat);
    Rhat = (1:n)';
    Rhat(Rhatposition) = Rhat;

        % Run the bootstrap
    [percent, efron] = rank_bootstrap(thetahat, Rhat, B);
    design2_results(:,1,m) = (R >= percent(:,1)) & (R <= percent(:,2));
    design2_results(:,2,m) = (R >= efron(:,1)) & (R <= efron(:,2));
end

    % Get results
design2_results = mean(design2_results,3);

%% Design 3: theta = .1*(1:10)'
    % Construct true theta
theta = .1*(1:n)';

    % Construct true ranks
[~,Rposition] = sort(theta);
R = (1:n)';
R(Rposition) = R;

    % Initialize vector of results
design3_results = zeros(n,2,M);

    % Run the Montel Carlo simulations
for m = 1:M
        % Draw thetahat, construct Rhat
    thetahat = randn(n,1) + theta;
    [~, Rhatposition] = sort(thetahat);
    Rhat = (1:n)';
    Rhat(Rhatposition) = Rhat;

        % Run the bootstrap
    [percent, efron] = rank_bootstrap(thetahat, Rhat, B);
    design3_results(:,1,m) = (R >= percent(:,1)) & (R <= percent(:,2));
    design3_results(:,2,m) = (R >= efron(:,1)) & (R <= efron(:,2));
end

    % Get results
design3_results = mean(design3_results,3);

