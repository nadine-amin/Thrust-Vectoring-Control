clear
clc

gamma = 1.4;
R = 287.05;
nozzle_height = 0.5;
nozzle_width = 0.5;
flap_length = 0.25;
Pb = 100000;

% Region 1
M1 = 3;
P1 = 12500;
T0 = 2000;
theta_max1 = getMaxTheta(M1, gamma);

P01 = P1*((1+((M1^2)*(gamma-1)/2))^(gamma/(gamma-1)));
T1 = T0/(1+((M1^2)*(gamma-1)/2));

%%
theta_max_thrust = zeros(2,16);
max_thrust = zeros(1,16);
for phi=0:1:15
theta_combinations = [];
thrust = [];
for theta3 = 0:0.1:theta_max1
    flag = 0;
    % Region 3
    [beta3, M3, P31, P031, T31] = weak_oblique(M1, theta3, gamma);

    % Region 4
    theta4 = phi + theta3;
    if ((theta4 > getMaxTheta(M3, gamma)) || (theta4 < 0))
        continue;
    end
    [beta4, M4, P43, P043, T43] = weak_oblique(M3, theta4, gamma);

    T4 = T43*T31*T1;
    P41 = P43*P31;
    P4 = P41*P1;
    P04 = P043*P031*P01;

    % Region 2
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
    % Region 5
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
                    if ((P51 - P41)> 1e-2)
                        theta2 = theta2 - 1e-4;
                        theta5 = theta2 - phi;
                    elseif ((P51 - P41)> 1e-3)
                        theta2 = theta2 - 1e-5;
                        theta5 = theta2 - phi;
                    else 
                        theta2 = theta2 - 1e-6;
                        theta5 = theta2 - phi;
                    end
                    if ((theta2 < 0) || (theta5 < 0))
                        break;
                    end
                elseif ((P41-P51)> 1e-4)
                    if ((P41 - P51)> 1e-2)
                        theta2 = theta2 + 1e-4;
                        theta5 = theta2 - phi;
                    elseif ((P41 - P51)> 1e-3)
                        theta2 = theta2 + 1e-5;
                        theta5 = theta2 - phi;
                    else 
                        theta2 = theta2 + 1e-6;
                        theta5 = theta2 - phi;
                    end
                    [beta2, M2, P21, P021, T21] = weak_oblique(M1, theta2, gamma);
                    if ((theta2 > theta_max1) || (theta5 > getMaxTheta(M2, gamma)))
                        break;
                    end
                else
                    flag = 1; 
                    break;
                end
    end
    
    if (flag)
        theta_combinations = [theta_combinations [theta3;theta2]];
        
        v_exit_up = M4*sqrt(gamma*R*T4);
        v_exit_down = M5*sqrt(gamma*R*T5);

        % Calculating Exit Areas
        h3 = nozzle_height/((sind(beta3))+((cosd(beta3))*(tand(beta2))));
        h2 = h3*(cosd(beta3))/(cosd(beta2));
        l3 = sqrt((flap_length^2)+(h3^2)-(2*flap_length*h3*cosd(beta3-theta3)));
        l2 = sqrt((flap_length^2)+(h2^2)-(2*flap_length*h2*cosd(beta2-theta2)));
        zeta3 = 180 - phi - beta3 - asind((flap_length)*(sind(beta3-theta3))/(l3));
        zeta2 = 180 + phi - beta2 - asind((flap_length)*(sind(beta2-theta2))/(l2));
        area_up = nozzle_width*l3*sind(zeta3);
        area_down = nozzle_width*l2*sind(zeta2);

        m_dot_up = (v_exit_up*area_up*P4)/(R*T4);
        m_dot_down = (v_exit_down*area_down*P5)/(R*T5);

        thrust_value = m_dot_up*v_exit_up + m_dot_down*v_exit_down + (P4 - Pb)*area_up + (P5 - Pb)*area_down;
        thrust = [thrust thrust_value];
    end
end
theta_max_thrust(:, phi+1) = theta_combinations(:, find(thrust==max(thrust)));
max_thrust(1, phi+1) = max(thrust);

figure
plot((0:length(thrust)-1),thrust)
title('Thrust at Different Combinations of Theta 2 and 3')
xlabel('Number of Combination')
ylabel('Thrust')
xlim([0 length(thrust)])
grid on
end
theta_max_thrust
max_thrust

Pe_Same_as_Pb_Case
CompareThrusts