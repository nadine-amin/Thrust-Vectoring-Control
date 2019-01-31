function [theta] = theta(beta)

global Gamma M
theta = -acot((tan(beta))*((((Gamma+1)*(M^2))/(2*(((M^2)*((sin(beta))^2))-1)))-1));

end

