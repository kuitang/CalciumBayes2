function [ theta ] = m_step_pos_beta(theta, options, w_bound, varargin)
%M_STEP_FULL Summary of this function goes here
%   Detailed explanation goes here

    %options = optimset(options, 'DerivativeCheck', 'on');
    %options = optimset(options,'FinDiffType', 'central');
    lb = ones(size(theta)) * -w_bound;
    ub = ones(size(theta)) * w_bound;
    lb(1) = 0;
    ub(1) = 5;
    [theta,fval,exitflag,output] = fmincon('q_single_neuron_bw',theta, [], [], [], [], lb, ub,[],options, varargin{:});

end

