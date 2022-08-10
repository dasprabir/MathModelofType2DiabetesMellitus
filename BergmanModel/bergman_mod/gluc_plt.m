function gluc_plt(tspan,x,tu,p,sigma_nu,sigma_mu)
Gb = p(2); 
21
Ib = p(5);
figure
subplot(221); h = plot(tspan,x(:,1), '-','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)
%baseline level
ylabel('glucose level [mg/dL]')
legend(h, 'model simulation')
u = interp1(tu(:,1),tu(:,2), tspan); %reconstruct used input signal
subplot(222); plot(tspan,u, '--or','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)
%baseline level
ylabel('insulin level [\muU/mL]'); xlabel('time [min]')
legend('interpolated measured test data')
%figure
subplot(223); plot(tspan, x(:,2), 'Linewidth',2)
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')
end