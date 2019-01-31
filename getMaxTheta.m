function [theta_max] = getMaxTheta(M1, gamma)

global Gamma M
Gamma = gamma;
M = M1;

beta_theta_max = fminbnd(@theta,0,pi/2);
theta_max = acot((tan(beta_theta_max))*((((Gamma+1)*(M^2))/(2*(((M^2)*((sin(beta_theta_max))^2))-1)))-1));
theta_max = theta_max*180/pi;
end

