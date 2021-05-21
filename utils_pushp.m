classdef utils_pushp
    methods(Static)
        function [vel_prime] = get_vel_prime(X, U)
            global m g
            vel = X(4:6);
            att = X(7:9);
            omega = X(10:12);
            
            Rib = utils_pushp.InertialToBody(att);
            vel_prime = [0; 0; -U(1)]/m + Rib*[0; 0; g] - (cross(omega, vel));
        end
        
        function [omega_prime] = get_omega_prime(X, U)
            global J
            omega = X(10:12);
            Mb = U(2:4);
            omega_prime = inv(J)*(Mb - cross(omega,J*omega));
        end
        
        function [E_prime] = get_E_prime(E, omega)

            phi = E(1);
            theta = E(2);
            A = [1     sin(phi)*sin(theta)/cos(theta)       cos(phi)*sin(theta)/cos(theta);
                 0     cos(phi)                             -sin(phi);
                 0     sin(phi)/cos(theta)                  cos(phi)/cos(theta)];
            E_prime = A*omega;
        end
        
        function [R] = InertialToBody(E)

            phi = E(1);
            theta = E(2);
            psi = E(3);
            R = [cos(theta)*cos(psi)                                      cos(theta)*sin(psi)                               -sin(theta)   ;
                -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)     cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)      sin(phi)*cos(theta);
                 sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)     -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)     cos(phi)*cos(theta)];

        end
        
        function [X_prime] = get_X_prime(E, vel)
            phi = E(1);
            theta = E(2);
            psi = E(3);
            R = [cos(theta)*cos(psi)                                      cos(theta)*sin(psi)                               -sin(theta)   ;
                -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)     cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)      sin(phi)*cos(theta);
                 sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)     -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)     cos(phi)*cos(theta)];

            Rbi = R';
             
            X_prime = Rbi*vel;
        end
        
    end
end