function [trans_shape,error] = shape_phase_trans(est_shape,true_shape) 
% match two shapes with one tranlation constant


n = length(est_shape);
error = norm(est_shape-true_shape);
trans_shape = est_shape;

for i = 1:n-1
    trans_est_shape = circshift(est_shape.',i).';
    trans_error = norm(trans_est_shape-true_shape)/sqrt(n);
    if trans_error<error
        error = trans_error;
        trans_shape= trans_est_shape;
    end
end


