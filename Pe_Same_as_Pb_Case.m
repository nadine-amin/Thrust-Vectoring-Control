clear
clc

gamma = 1.4;
R = 287.05;
height = 0.5;
width = 0.5;
flap = 0.25;
Pb = 100000;

% Region 1
M1 = 3;
P1 = 12500;
T0 = 2000;
theta_max1 = getMaxTheta(M1, gamma);
Pb1 = Pb/P1;

P01 = P1*((1+((M1^2)*(gamma-1)/2))^(gamma/(gamma-1)));
T1 = T0/(1+((M1^2)*(gamma-1)/2));

%%
theta_optimum_thrust = zeros(2,16);
opt_thrust = zeros(1,16);
for phi = 0:1:15;
theta3 = phi
theta4 = phi + theta3;

while (1)
        if (theta3 > theta_max1)
            theta3 = theta3 - 1;
        elseif (theta3 < 0)
            theta3 = theta3 + 1;
        else
            break;
        end
end
    
[beta3, M3, P31, P031, T31] = weak_oblique(M1, theta3, gamma);
while (1)
        if ((theta3 > theta_max1) || (theta4 > getMaxTheta(M3, gamma)))
            theta3 = theta3 - 1;
            theta4 = phi + theta3;
            [beta3, M3, P31, P031, T31] = weak_oblique(M1, theta3, gamma);
        elseif ((theta3 < 0) || (theta4 < 0))
            theta3 = theta3 + 1;
            theta4 = phi + theta3;
            [beta3, M3, P31, P031, T31] = weak_oblique(M1, theta3, gamma);
        else
            break;
        end
end

flag = 0;
while (1)
    [beta3, M3, P31, P031, T31] = weak_oblique(M1, theta3, gamma);
    
    theta4 = phi + theta3;
    [beta4, M4, P43, P043, T43] = weak_oblique(M3, theta4, gamma);
    
    T4 = T43*T31*T1;
    P41 = P43*P31;
    P4 = P41*P1;
    P04 = P043*P031*P01;
    
    if ((P41 - Pb1)> 1e-4)
        if ((P41 - Pb1)> 1e-2)
            theta3 = theta3 - 1e-4;
            theta4 = phi + theta3;
        else 
            theta3 = theta3 - 1e-5;
            theta4 = phi + theta3;
        end
        if ((theta3 < 0) || (theta4 < 0))
            break;
        end
    elseif ((Pb1 - P41)> 1e-4)
        if ((Pb1 - P41)> 0.1)
            theta3 = theta3 + 1e-4;
            theta4 = phi + theta3;
        else 
            theta3 = theta3 + 1e-5;
            theta4 = phi + theta3;
        end
        [beta3, M3, P31, P031, T31] = weak_oblique(M1, theta3, gamma);
        if ((theta3 > theta_max1) || (theta4 > getMaxTheta(M3, gamma)))
            break;
        end
    else
        flag = 1; 
        break;
    end

end

flag2 = 0;
if (flag)
    theta2 = theta4;
    while (1)
        if (theta2 > theta_max1)
            theta2 = theta2 - 1;
        elseif (theta2 < 0)
            theta2 = theta2 + 1;
        else
            break;
        end
    end
    
    [beta2, M2, P21, P021, T21] = weak_oblique(M1, theta2, gamma);
    theta5 = theta2 - phi;
    
    while (1)
        if ((theta2 > theta_max1) || (theta5 > getMaxTheta(M2, gamma)))
            theta2 = theta2 - 1;
            theta5 = theta2 - phi;
            [beta2, M2, P21, P021, T21] = weak_oblique(M1, theta2, gamma);
        elseif ((theta2 < 0) || (theta5 < 0))
            theta2 = theta2 + 1;
            theta5 = theta2 - phi;
            [beta2, M2, P21, P021, T21] = weak_oblique(M1, theta2, gamma);
        else
            break;
        end
    end

    while (1)
                [beta2, M2, P21, P021, T21] = weak_oblique(M1, theta2, gamma); 
                
                % Region 5
                theta5 = theta2 - phi;
                
                [beta5, M5, P52, P052, T52] = weak_oblique(M2, theta5, gamma);
                T5 = T52*T21*T1;
                P51 = P52*P21;
                P5 = P51*P1;
                P05 = P052*P021*P01;

                if ((P51 - P41)> 1e-4)
                    if ((P51 - P41)> 0.1)
                        theta2 = theta2 - 1e-4;
                        theta5 = theta2 - phi;
                    else 
                        theta2 = theta2 - 1e-5;
                        theta5 = theta2 - phi;
                    end
                    if ((theta2 < 0) || (theta5 < 0))
                        break;
                    end
                elseif ((P41-P51)> 1e-4)
                    if ((P41 - P51)> 0.1)
                        theta2 = theta2 + 1e-4;
                        theta5 = theta2 - phi;
                    else 
                        theta2 = theta2 + 1e-5;
                        theta5 = theta2 - phi;
                    end
                    [beta2, M2, P21, P021, T21] = weak_oblique(M1, theta2, gamma);
                    if ((theta2 > theta_max1) || (theta5 > getMaxTheta(M2, gamma)))
                        break;
                    end
                else
                    flag2 = 1; 
                    break;
                end
    end
end

if (flag2)
    theta_optimum_thrust(:, phi+1) = [theta3;theta2];
    
    v_exit_up = M4*sqrt(gamma*R*T4);
    v_exit_down = M5*sqrt(gamma*R*T5);

    % Calculating Exit Areas
    h3 = height/((sind(beta3))+((cosd(beta3))*(tand(beta2))));
    h2 = h3*(cosd(beta3))/(cosd(beta2));
    l3 = sqrt((flap^2)+(h3^2)-(2*flap*h3*cosd(beta3-theta3)));
    l2 = sqrt((flap^2)+(h2^2)-(2*flap*h2*cosd(beta2-theta2)));
    zeta3 = 180 - phi - beta3 - asind((flap)*(sind(beta3-theta3))/(l3));
    zeta2 = 180 + phi - beta2 - asind((flap)*(sind(beta2-theta2))/(l2));
    area_up = width*l3*sind(zeta3);
    area_down = width*l2*sind(zeta2);

    m_dot_up = (gamma*P04*area_up*M4)/((sqrt(gamma*R*T0))*((1+((gamma-1)*(M4^2)/2))^((gamma+1)/(2*(gamma-1)))));
    m_dot_down = (gamma*P05*area_down*M5)/((sqrt(gamma*R*T0))*((1+((gamma-1)*(M5^2)/2))^((gamma+1)/(2*(gamma-1)))));

    optimum_thrust = m_dot_up*v_exit_up + m_dot_down*v_exit_down + (P4 - Pb)*area_up + (P5 - Pb)*area_down;
    opt_thrust(1, phi+1) = optimum_thrust;
end
end
theta_optimum_thrust
opt_thrust