function [ theta ] = m_step_full_j(theta, options, N, varargin)
%M_STEP_FULL Summary of this function goes here
%   Detailed explanation goes here

    %options = optimset(options, 'DerivativeCheck', 'on');
    %options = optimset(options,'FinDiffType', 'central');
    lb = ones(size(theta)) * -.1;
    ub = ones(size(theta)) * .1;
    lb(1) = 0; lb(2:N+1) = -5;
    ub(1) = 5; ub(2:N+1) = 5;
    [theta,fval,exitflag,output] = fmincon('q_single_neuron_full_j',theta, [], [], [], [], lb, ub,[],options, varargin{:});

end

