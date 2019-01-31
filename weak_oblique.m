function [beta, M2, pressure_ratio, tot_pressure_ratio, temp_ratio] = weak_oblique(M1, theta, gamma)

    flag = 0;
    if (theta == 0)
        theta = 0.001;
        flag = 1;
    end
    lambda = ((((M1^2)-1)^2)-((3)*(1+(0.5*(M1^2)*(gamma-1)))*(1+(0.5*(M1^2)*(gamma+1)))*((tand(theta))^2)))^0.5;
    X = ((((M1^2)-1)^3)-((9)*(1+(0.5*(M1^2)*(gamma-1)))*(1+(0.5*(M1^2)*(gamma-1))+(0.25*(M1^4)*(gamma+1)))*((tand(theta))^2)))/(lambda^3);
    betaW = atan(((M1^2)-(1)+(2*(lambda)*(cos(((4*pi)+(acos(X)))/3))))/((3)*(1+(0.5*(M1^2)*(gamma-1)))*(tand(theta))));
    beta = betaW*180/pi;
    Mn1W = M1*sin(betaW);
    P21W = 1 + (((2*gamma)/(gamma+1))*(((Mn1W)^2)-1));
    pressure_ratio = P21W;
    Mn2W = sqrt(((Mn1W^2)+(2/(gamma-1)))/((((2*gamma)/(gamma-1))*(Mn1W^2))-1));
    M2 = (Mn2W)/(sin(betaW - (theta*pi/180)));
    tot_pressure_ratio = P21W*(((1+((M2^2)*(gamma-1)/2))/(1+((M1^2)*(gamma-1)/2)))^(gamma/(gamma-1)));
    temp_ratio = P21W*((2+((gamma-1)*(Mn1W^2)))/((gamma+1)*(Mn1W^2)));
    
    if (flag)
        M2 = M1;
        pressure_ratio = 1;
        tot_pressure_ratio = 1;
        temp_ratio = 1;
    end
end

