function ll = log_likelihood(beta, b, w, h, n, delta, p_weights)

% error('TODO: Fix the log_likelihood interface to only use the history
% terms of one neuron at a time')
% isn't this already the case??

%log(P(data | params))

S  = size(beta,2) + 1;
[N, T] = size(n);
M = size(h,3);



ll = 0;
        for t = S+1:T
                        

            I_terms = beta .* n(:,(t-2):-1:(t-S));
            I = sum(I_terms(:));
            
            for m = 1:M                
                J = b + I + w * h(:,t,m)';
                eJd = exp(J)*delta;                
                if n(t)
                    eeJd = exp(-eJd);
                    Qm = log(1 - eeJd);                    
                else
                    Qm = -eJd;                    
                end
                ll = ll + p_weights(t,m) * Qm;
            end
        end                


end
