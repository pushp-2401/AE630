clc;
clear all;

global dt m g J X zeta_pos omega_pos zeta_att omega_att dt

dt = 0.01;
t_final = 10;
t = 0.01:dt:t_final;
N = ceil(t_final/dt)+1 ;

X = zeros(12,N); % x,y,z,u,v,w,phi,theta,psi,p,q,r
Xd = zeros(12,N);
pd = zeros(3,N-1);
pd_prime = zeros(3,N-1);
pd_double_prime = zeros(3,N-1);
ad = zeros(3,N-1);
ad_prime = zeros(3,N-1);
ad_double_prime = zeros(3,N-1);

U = zeros(4, N-1); % Thrust,l,m,n

header_pushp();

for k = 1:N-1
    
    %trajectory
%     pos_des(:,k) = [0.1*sin(pi*k); 0; 0];
%     pos_des_dot(:,k) = [0.1*pi*cos(pi*k); 0; 0];
%     pos_des_dot_dot(:,k) = [-0.1*pi*pi*sin(pi*k); 0; 0];
    pd(:,k) = [1; 2; 10];
    pd_prime(:,k) = [0; 0; 0];
    pd_double_prime(:,k) = [0; 0; 0];
    ad(:,k) = [0; 0; 0];
    ad_prime(:,k) = [0; 0; 0];
    ad_double_prime(:,k) = [0; 0; 0];
    
    %Control
    U(:,k)= control_pushp(X(:,k), pd(:,k), pd_prime(:,k), pd_double_prime(:,k), ad(:,k), ad_prime(:,k), ad_double_prime(:,k));
    
    %State Update
    X(:,k+1) = state_update_pushp(X(:,k), U(:,k));
end

figure(1)
subplot(3,1,1)
plot(t,(X(1,1:N-1)),'r',t,(pd(1,:)))
xlabel('t')
ylabel('x (m)')
subplot(3,1,2)
plot(t,(X(2,1:N-1)),'r',t,(pd(2,:)))
xlabel('t')
ylabel('y (m)')
subplot(3,1,3)
plot(t,(X(3,1:N-1)),'r',t,(pd(3,:)))
xlabel('t')
ylabel('z (m)')

figure(2)
subplot(3,1,1)
plot(t,rad2deg(X(7,1:N-1)),'r',t,rad2deg(ad(1,:)))
xlabel('t')
ylabel('Roll (^o)')
subplot(3,1,2)
plot(t,rad2deg(X(8,1:N-1)),'r',t,rad2deg(ad(2,:)))
xlabel('t')
ylabel('Pitch (^o)')
subplot(3,1,3)
plot(t,rad2deg(X(9,1:N-1)),'r',t,rad2deg(ad(3,:)))
xlabel('t')
ylabel('Yaw (^o)')

figure(3)
subplot(2,2,1)
plot(t,U(2,:))
xlabel('t')
ylabel('Rolling Moment (N-m)')
subplot(2,2,2)
plot(t,U(3,:))
xlabel('t')
ylabel('Pitching Moment (N-m)')
subplot(2,2,3)
plot(t,U(4,:))
xlabel('t')
ylabel('Yawing Moment (N-m)')
subplot(2,2,4)
plot(t,U(1,:))
xlabel('t')
ylabel('Thrust (N)')