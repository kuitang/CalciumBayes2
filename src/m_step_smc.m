function [theta_intrinsic] = m_step_smc(theta_intrinsic, options, varargin)

    %options = optimset(options,'DerivativeCheck', 'on');
    %options = optimset(options,'FinDiffType', 'central');
    lb = ones(size(theta_intrinsic)) * -.1;
    ub = ones(size(theta_intrinsic)) * .1;
    lb(1) = 0; lb(2) = -5;
    ub(1) = 5; ub(2) = 5;
    [theta_intrinsic,fval,exitflag,output] = fmincon('q_single_neuron',theta_intrinsic, [], [], [], [], lb, ub,[],options, varargin{:});

end

