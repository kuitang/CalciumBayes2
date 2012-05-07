clear
run_parallel = 1;

%% use the same random number sequence for debugging purposes
%RandStream.setDefaultStream ...
%     (RandStream('mt19937ar','seed',2));

%% Load data
% load('../data/truncdata.mat')
% truncate the data even more!
% truncdata = truncdata(:,1:100);
% n = truncdata;
data = '../data/25n_noise';
load([data '.mat'])
n = sim.n(1:12,1:15000);
%load('../data/good_sim_data_01.mat')
%n = sim.n(1:10,:);

%% Set optimization options
optim_options = optimset('LargeScale','on','Algorithm', ...
    'trust-region-reflective','GradObj','on','Hessian','user-supplied', 'MaxIter',30);

%% Set physical parameters
% Physical parameters (TODO: Set up and figure out scale!)
% Currently using unit (discrete) time and bullshit values
%sigma = 0.01;
sigma = .05;
tau = .020; %set to match simulator
delta = .010;

beta_bound = 5;
w_bound = 5;
w_reg = 10;
beta_reg = 5;



[N T] = size(n);
S = 20; % Indirect time window
M = 50; % size of particle sampler

% Parameter matrices
% TODO: Set up priors!
iters = 1;
%% Set codistributed arrays
if run_parallel
    spmd(N)
        codist = codistributor1d(1);
        beta = zeros(N, N, S-1, codist);
        % lambda = ones(N, S);
        b = ones(N, 1, codist);
        w = zeros(N, N, codist) .* .1*exprnd(.9,N);
        p_weights = zeros(N,T,M,codist);
        % for i=1:N/5
        %     rand_inh
        %     w(i,:) = -exprnd(2.3,N,1);
        % end
        h = zeros(N,N,T,M, codist);

        w = w .* binornd(1,.1,N,N);%second arg is "sparesness"

        llv = zeros(N,1,codist);
%        for i = drange(1:N)
%            w(i,i) = -abs(normrnd(.6,.2));
%        end
    end
else
    %% Set noncodistributed arrays
    beta = zeros(N, N, S-1);
    % lambda = ones(N, S);
    b = zeros(N, 1);
    w = ones(N, N) .* .1*exprnd(.9,N);
    p_weights = zeros(N,T,M);
    % for i=1:N/5
    %     rand_inh
    %     w(i,:) = -exprnd(2.3,N,1);
    % end

    w = w .* binornd(1,.1,N,N);%second arg is "sparesness"

    h = zeros(N,N,T,M);

    for i = 1:N
        w(i,i) = -abs(normrnd(.6,.2));
    end
end
% load('first_e_step_complete.mat');
% first = 1;

%% Main loop
%until the change in connectivity matrix w is below threshold change

w_prev = ones(size(w)) * 500;
thresh_w = .001;

ll = -Inf;

while(norm(w - w_prev) > thresh_w)    
    
    
    disp(w);    
    disp(['NORM: ' num2str(norm(w - w_prev))]); 
    disp('********BEGINNING OF LOOP********');    
    disp('*********************************');
    w_prev = w; 
    
   disp(['iters: ' num2str(iters)]); 
    spmd(N)  
        for i = drange(1:N)
            disp(['Neuron ' num2str(i) '/' num2str(N)]);            

            %% Initialize the intrinsic parameters
            theta = [b(i) w(i,:) reshape(beta(i, :, :),1,N*S-N)];
                
                
            %% E step (SMC) for one neuron                
            beta_subset = reshape(beta(i,:,:), N, S - 1);
            [p_weights(i,:,:) h(i,:,:,:)] = e_step_smc(i, M, tau, delta, sigma, beta_subset, b(i), w(i,:), n);

            %% M step for the intrinsic parameters for one neuron
            theta = m_step_full(theta, optim_options, N, beta_bound, w_bound, beta_subset, w(i,:), squeeze(h(i,:,:,:)), n, i,...
                delta, tau, sigma, squeeze(p_weights(i,:,:)), w_reg, beta_reg);
            b(i,1) = theta(1);
            w(i,:) = reshape(theta(2:N+1),1,N);
            beta(i,:,:) = reshape(theta(N+2:end), 1, N, (S - 1));
            disp('new params:');
            disp(b(i,1));
            disp(w(i,:));
        end
           disp(['NEURON ' num2str(i) ' DONE!']);
    end

    disp('OVER!');
    disp(b);
    disp(w);
    iters = iters + 1;
    w_gathered = gather(w);
    beta_gathered = gather(beta);
    b_gathered = gather(b);
    %% Log likelihood for whole model'
    spmd
        for i = drange(1:N)
           llv(i,1) = log_likelihood(i, reshape(beta(i,:,:),N,S-1), b(i), w(i,:), squeeze(h(i,:,:,:)),n,delta,squeeze(p_weights(i,:,:)));
        end
    end
    nll = sum(gather(llv));
        ll = [ll nll];
    disp('ll =');
    disp(ll);
    disp('diff = ');
    disp(ll(iters) - ll(iters-1));
    save([data '_biggerreg.mat'], 'iters','sigma', 'tau', 'delta', 'w_gathered', 'beta_gathered', 'b_gathered','data','ll');


end



