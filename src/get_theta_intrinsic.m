function [ b_i w_ii beta_iis ] = get_theta_intrinsic( i, b, w, beta )
%GET_THETA_INTRINSIC Summary of this function goes here
%   Detailed explanation goes here

S = size(beta, 3) + 1;
b_i = b;
w_ii = w(i,i);
beta_iis = reshape(beta(i, :),1,S-1);

end

