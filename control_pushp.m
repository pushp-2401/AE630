function [U] = control_pushp(X, pd, pd_prime, pd_double_prime, ad, ad_prime, ad_double_prime)
global m g J zeta_pos omega_pos zeta_att omega_att
U = [0; 0; 0; 0];

%------Outer loop--------
pos = X(1:3);
vel = X(4:6);
att = X(7:9);
omega = X(10:12);

pos_prime = utils_pushp.get_X_prime(att, vel); 
pos_double_prime = (pd_double_prime - (2*zeta_pos*omega_pos)*(pos_prime - pd_prime) - (omega_pos^2)*(pos-pd));

U(1) = m*sqrt(pos_double_prime(1)^2 + pos_double_prime(2)^2 + (g - pos_double_prime(3))^2); % Thrust req
ad(1) = asin((-m*pos_double_prime(1)*sin(ad(3)) + m*pos_double_prime(2)*cos(ad(3)))/U(1)); % Phi req
ad(2) = asin((-m*pos_double_prime(1)*cos(ad(3)) - m*pos_double_prime(2)*sin(ad(3)))/(U(1)*cos(ad(1)))); % Theta req

%-------Inner Loop--------
phi = att(1);
theta = att(2);

A = [1  sin(phi)*tan(theta)     cos(phi)*tan(theta);
     0  cos(phi)                -sin(phi);
     0  sin(phi)/cos(theta)     cos(phi)/cos(theta)];

att_prime = utils_pushp.get_E_prime(att, omega);
phi_prime = att_prime(1);
theta_prime = att_prime(2);

A_prime = [0  cos(phi)*phi_prime*tan(theta)+sin(phi)*theta_prime*sec(theta)^2                        -sin(phi)*phi_prime*tan(theta)+cos(phi)*theta_prime*sec(theta)^2;
        0  -sin(phi)*phi_prime                                                                 -cos(phi)*phi_prime;
        0  (cos(theta)*cos(phi)*phi_prime + sin(phi)*sin(theta)*theta_prime)/(cos(theta)^2)      (-cos(theta)*sin(phi)*phi_prime + sin(theta)*cos(phi)*theta_prime)/(cos(theta)^2)];

U(2:4) = -J/A*((A_prime*omega) + (2*zeta_att*omega_att)*(A*omega - ad_prime) + ...
         (omega_att^2)*(att-ad) - ad_double_prime)...
         + cross(omega, J*omega);


end

