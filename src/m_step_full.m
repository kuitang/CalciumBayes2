function [ theta ] = m_step_full(theta, options, N, beta_bound, w_bound, varargin)
%M_STEP_FULL Summary of this function goes here
%   Detailed explanation goes here

    %options = optimset(options, 'DerivativeCheck', 'on');
    %options = optimset(options,'FinDiffType', 'central');
    lb = ones(size(theta)) * -beta_bound;
    ub = ones(size(theta)) * beta_bound;
    lb(1) = 0; lb(2:N+1) = -w_bound;
    ub(1) = 5; ub(2:N+1) = w_bound;
    [theta,fval,exitflag,output] = fmincon('q_single_neuron_full',theta, [], [], [], [], lb, ub,[],options, varargin{:});

end

