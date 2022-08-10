function gluc = gluc_sim(tspan,x0,tu, p, sigma_nu,sigma_mu,plt)
%SIM_INPUT_SIM Simulation of glucose minimal model.
%x0: initial conditions for glucose conc. and state variable denoting
% insulin action
%gluc: model output, glucose conc. in plasma [mg/dL]
if tspan(1) < tu(1,1)
 error(' no input defined at start of simulation (compare tspan(1))vs. tu(1,1)')
end
if 1 %Forward Euler integration
 td = .2*min(diff(tspan)); %simulation time step is taken as 1/5 of
the
 %minimal time interval between the samples of the experimental time
 %vector
 tspan_res = [];
 for i = 1:length(tspan)-1
 tspan_res = [tspan_res, tspan(i) :td: tspan(i+1)-td];
%resample while assuring that all time
 %samples of tspan also occur
in tspan_res
 end
 tspan_res = [tspan_res, tspan(end)];
 TD=diff(tspan_res); %vector with integration time steps (most equal
to td) 
20

 u = interp1(tu(:,1),tu(:,2), tspan_res);
 N = length(tspan_res);
 x = x0;
 for k = 2:N
 dx(:,k) = gluc_ode([],x(k-1,:),u(k), p,sigma_nu); %column
 x(k,:) = x(k-1,:) + dx(:,k)'*TD(k-1); %different
time samples in rows
 end
 for i = 1:length(tspan)
 id(i) = find(tspan_res==tspan(i)); %select samples
corresponding to experiment
 end
 x = x(id,:);
elseif 0 %
 ode_options = [];
 [t,x] = ode45(@gluc_ode,tspan,x0,ode_options, tu, p,sigma_nu);
 %[t,x] = ode15s(@gluc_ode,tspan,x0,ode_options, tu, p,sigma_nu);
end
%Output
 gluc = x(:,1);
if plt==1
 gluc_plt(tspan,x,tu,p,sigma_nu,sigma_mu)
end 