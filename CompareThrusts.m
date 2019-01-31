percentage_error = ((max_thrust-opt_thrust)./opt_thrust)*100;
avg_percentage_error = sum(percentage_error)/length(percentage_error)


phi = 0:1:15;
figure
plot(phi, max_thrust, phi, opt_thrust)
title('Calculated Maximum Thrust and Thrust at P_e = P_b at Different Phi Values')
xlabel('Phi')
ylabel('Thrust')
xlim([0 15])
ylim([0 6.5e+04])
grid on