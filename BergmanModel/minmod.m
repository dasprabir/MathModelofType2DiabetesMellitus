function minmod
%MINMOD Implement the glucose-insulin minimal model
%   MINMOD solves the glucose-insulin minimal model that assesses
%   insulin sensitivity based on observed insulin concentrations
%   and using a forcing function. 
clc;
clear all;
p.ib     = 10;     % uU/mL
p.gb     = 100;    % mg/dL
p.SG     = 0.03;   % 1/min  
p.SI     = 0.0008; % 1/min * mL/uU
p.p2     = 0.04;   % 1/min
p.doseIV = 21000;  % mg
p.V      = 70;     % Vd glucose dL 
[t,y] = ode15s(@derivatives,[0 180],[p.doseIV/p.V 0],[],p);
subplot(3,1,1)
plot(t,y(:,1),'LineWidth',2);
legend('glucose')
ylabel('Glucose [mg/dL]');           % xlabel('Time [min]');
subplot(3,1,2)
plot(t,y(:,2),'LineWidth',2);
legend('X')
ylabel('Insulin Action [1/min]');   % xlabel('Time [min]');
subplot(3,1,3,'LineWidth',2)
plot(t,i(t),'r')
legend('insulin')
ylabel('Insulin [uU/mL]'); xlabel('Time [min]'); 
print('-dtiff','-r900','minmod')

end


function dydt = derivatives(t, y, p)
%DERIVATIVES Compute the right-hand side of the ODE.
% DYDT = DERIVATIVES(T, Y, P) calculates |DYDT|, the right-hand
% side of the ODE model, at points defined by the vector of
% dependent variables |Y|, time |T|, and with parameters |P|.

glu  = y(1);     % glucose
x    = y(2);     % insulin action as forcing function
dglu = -(p.SG+x)*glu + p.SG*p.gb;
dx   = - p.p2*(x - p.SI*(i(t)-p.ib));
dydt = [dglu; dx];

end


function ins = i(t)
%I interpolate an insulin value
%   INS = I(T) computes an insulin value |INS| at time |T|,
%   using the linear interpolation given observation times and 
%   values.

tObs   = [ 0  1   2  5 10 20 29 30  31  35 40 50 100 180]';
insObs = [10 25 100 50 50 50 50 450 200 50 30 10  10  10]';
interpolation = fit(tObs, insObs, 'linearinterp');
ins = feval(interpolation, t);

end