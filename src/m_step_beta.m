function [ theta ] = m_step_pos_beta(theta, options, N, beta_bound, varargin)
%M_STEP_FULL Summary of this function goes here
%   Detailed explanation goes here

    %options = optimset(options, 'DerivativeCheck', 'on');
    %options = optimset(options,'FinDiffType', 'central');
    lb = ones(size(theta)) * -beta_bound;
    ub = ones(size(theta)) * beta_bound;
    [theta,fval,exitflag,output] = fmincon('q_single_neuron_beta',theta, [], [], [], [], lb, ub,[],options, varargin{:});

end

