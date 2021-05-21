function [X_new] = state_update_pushp(X, U)
    global dt
    X_new = X;
    k1_1 = [utils_pushp.get_X_prime(X(7:9),X(4:6)); utils_pushp.get_vel_prime(X, U)];
    k1_2 = [utils_pushp.get_E_prime(X(7:9),X(10:12)); utils_pushp.get_omega_prime(X, U)];
    X(1:6)= X(1:6)+ 0.5*k1_1*dt;
    X(7:12) = X(7:12) + 0.5*k1_2*dt;
    k2_1 = [utils_pushp.get_X_prime(X(7:9),X(4:6)); utils_pushp.get_vel_prime(X, U)];
    k2_2 = [utils_pushp.get_E_prime(X(7:9),X(10:12)); utils_pushp.get_omega_prime(X, U)];
    X(1:6)= X(1:6)+ 0.5*k2_1*dt;
    X(7:12) = X(7:12) + 0.5*k2_2*dt;
    k3_1 = [utils_pushp.get_X_prime(X(7:9),X(4:6)); utils_pushp.get_vel_prime(X, U)];
    k3_2 = [utils_pushp.get_E_prime(X(7:9),X(10:12)); utils_pushp.get_omega_prime(X, U)];
    X(1:6)= X(1:6)+ 0.5*k3_1*dt;
    X(7:12) = X(7:12) + 0.5*k3_2*dt;
    k4_1 = [utils_pushp.get_X_prime(X(7:9),X(4:6)); utils_pushp.get_vel_prime(X, U)];
    k4_2 = [utils_pushp.get_E_prime(X(7:9),X(10:12)); utils_pushp.get_omega_prime(X, U)];
    X(1:6)= X(1:6)+ 0.5*k4_1*dt;
    X(7:12) = X(7:12) + 0.5*k4_2*dt;
    X_new(1:6) = X(1:6) + (1/6)*dt*(k1_1+2*k2_1+2*k3_1+k4_1);
    X_new(7:12) = X(7:12) + (1/6)*dt*(k1_2+2*k2_2+2*k3_2+k4_2);
    
end

