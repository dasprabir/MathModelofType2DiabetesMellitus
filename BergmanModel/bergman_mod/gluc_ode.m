function dxout = gluc_ode(t,xin,tu, p,sigma_nu)
%GLUC_ODE ODE's of glucose minimal model
% xin: state values at previous time sample x(k-1) = [G(k-1); X(k-1)]
% dxout: new state derivatives dx(k) = [dG(k); dX(k)], column vector
% tu: u(k) (input signal at sample k) OR matrix tu
% p = [Sg, Gb, k3, Si, Ib];
%History
%17-Jun-03, Natal van Riel, TU/e
idG = 1; idX = 2; %glucose conc. [mg/mL] and state variable denoting
insulin action [uU/mL]
Sg = p(1); %[1/min] glucose effectiveness
Gb = p(2); %[mg/mL] baseline glucose conc.
k3 = p(3); %[1/min]
Si = p(4); %[mL/uU*min] insulin sensitivity
Ib = p(5); %[uU/mL] baseline insulin conc. in plasma
if length(tu)==1 %u(k) is provided as input to this function
 u = tu;
else %calculate u(k) by interpolation of tu
 u = interp1(tu(:,1),tu(:,2), t);
end
%[t u]
%ode's
dG = Sg*(Gb - xin(idG)) - xin(idX)*xin(idG);
dX = k3*( Si*(u-Ib) - xin(idX) );
dxout = [dG; dX];
end