function e = obj_fn(p_var, p_fix,data,tspan,tu, sigma_nu,sigma_mu,plt)
%4 unknown model parameters:
% p_var = [Sg k3 Si G0];
%3 known parameters:
% p_fix = [Gb Ib x0(2)];
%p = [Sg, Gb, k3, Si, Ib];
p = [p_var(1), p_fix(1), p_var(2), p_var(3), p_fix(2)];
%x0 = [G0, X0]
x0 = [p_var(4), p_fix(3)];
gluc = gluc_sim(tspan,x0,tu, p, sigma_nu,sigma_mu,0);
%LSQNONLIN: objective function should return the model error
e = gluc-data;
if plt==1 %fast update during estimation process
 N=length(tspan);
 figure(2)
 plot(1:N,gluc,'b',1:N,data,'r'); drawnow
end