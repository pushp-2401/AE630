clear all
clc
data = [-16.2500   -0.6742
  -16.0000   -0.9991
  -15.7500   -1.0281
  -15.5000   -1.0555
  -15.2500   -1.0806
  -15.0000   -1.0994
  -14.7500   -1.1090
  -14.5000   -1.1130
  -14.2500   -1.1162
  -14.0000   -1.1178
  -13.7500   -1.1192
  -13.5000   -1.1182
  -13.2500   -1.1181
  -13.0000   -1.1180
  -12.7500   -1.1174
  -12.5000   -1.1162
  -12.2500   -1.1150
  -12.0000   -1.1150
  -11.7500   -1.1047
  -11.5000   -1.0933
  -11.2500   -1.0824
  -11.0000   -1.0715
  -10.7500   -1.0600
  -10.5000   -1.0485
  -10.2500   -1.0371
  -10.0000   -1.0256
   -9.7500   -1.0121
   -9.5000   -0.9989
   -9.2500   -0.9857
   -9.0000   -0.9722
   -8.7500   -0.9582
   -8.5000   -0.9438
   -8.2500   -0.9289
   -8.0000   -0.9135
   -7.7500   -0.8988
   -7.5000   -0.8739
   -7.2500   -0.8430
   -7.0000   -0.8110
   -6.7500   -0.7783
   -6.5000   -0.7480
   -6.2500   -0.7188
   -6.0000   -0.6868
   -5.7500   -0.6530
   -5.5000   -0.6184
   -5.2500   -0.5914
   -5.0000   -0.5603
   -4.7500   -0.5272
   -4.5000   -0.4932
   -4.2500   -0.4586
   -4.0000   -0.4335
   -3.7500   -0.4036
   -3.5000   -0.3741
   -3.2500   -0.3462
   -3.0000   -0.3249
   -2.7500   -0.3003
   -2.5000   -0.2716
   -2.2500   -0.2413
   -2.0000   -0.2188
   -1.7500   -0.1898
   -1.5000   -0.1595
   -1.2500   -0.1352
   -1.0000   -0.1072
   -0.7500   -0.0794
   -0.5000   -0.0544
   -0.2500   -0.0258
         0         0
    0.2500    0.0259
    0.7500    0.0794
    1.0000    0.1072
    1.2500    0.1352
    1.5000    0.1595
    1.7500    0.1898
    2.0000    0.2188
    2.2500    0.2413
    2.5000    0.2716
    2.7500    0.3002
    3.0000    0.3249
    3.2500    0.3462
    3.5000    0.3741
    3.7500    0.4037
    4.0000    0.4334
    4.2500    0.4586
    4.5000    0.4933
    4.7500    0.5272
    5.0000    0.5603
    5.2500    0.5914
    5.5000    0.6184
    5.7500    0.6530
    6.0000    0.6868
    6.2500    0.7188
    6.5000    0.7481
    6.7500    0.7783
    7.0000    0.8111
    7.2500    0.8431
    7.5000    0.8741
    7.7500    0.8987
    8.0000    0.9133
    8.2500    0.9288
    8.5000    0.9436
    8.7500    0.9581
    9.0000    0.9721
    9.2500    0.9856
    9.5000    0.9987
    9.7500    1.0120
   10.0000    1.0255
   10.2500    1.0370
   10.5000    1.0485
   10.7500    1.0600
   11.0000    1.0715
   11.2500    1.0825
   11.5000    1.0935
   11.7500    1.1050
   12.0000    1.1154
   12.2500    1.1154
   12.5000    1.1168
   12.7500    1.1180
   13.0000    1.1187
   13.2500    1.1189
   13.5000    1.1192
   13.7500    1.1204
   14.0000    1.1190
   14.2500    1.1177
   14.5000    1.1151
   14.7500    1.1104
   15.0000    1.1008
   15.2500    1.0814
   15.5000    1.0554
   15.7500    1.0285
   16.0000    0.9995];
rho = 1.225; %density
nb = 2;      %number of blade
R = 0.355;       %radius

C = 0.032;     %blade chord length
%c_tip = 0.5*C;

CDo = 0.01;  %profile drag coefficient
%Cla = 0.0876;  %lift curve slope obtained from data given for NACA0012 airfoil lookup table
omega = 157; %rotor angular rate
%Rct = 0.2;     %root cut out
%Cth = 0.006;   % thrust coefficient for hover

sigma = nb*C/(pi*R);
%sigma1 = nb*c_tip/ (pi*R); %ideal twist with ideal taper

% theta_tip = 4.86*pi/180;  %theta_tip decided by gaussaian quadrature
% theta_tip1 = 8*pi/180;  %ideal twist with ideal 

r = linspace ( 0.0000000000000001, 1, 129);
W = [0.171324492,0.360761573,0.467913935,0.467913935,0.360761573,0.171324492];
m = [0.932469514, 0.661209386, 0.02386191860,-0.02386191860,-0.661209386,-0.932469514];
%ct = 0
%theta_tip_closed = (4*Cth/(sigma*Cla)) + (Cth/2)^0.5;    %closed form of theta tip

% Use of gaussian quaderature to decide theta_tip
% for i = 1:length(r)-1
%     a = r(i);
%     b = r(i+1);
%     sum = 0;
%     for j = 1:6
%         ti = (((b-a)/2)* m(j))+ ((a+b)/2);
%         lambda = sigma *(Cla/16)* ((1+(32* theta_tip/(sigma* Cla)))^0.5 -1) ;
%         dct = 0.5* sigma* Cla* ((theta_tip*ti) - lambda*ti^2)*W(j)* ((b-a)/2);
%         sum = sum+ dct
%     end
%     ct = sum + ct;
% end
% error = ((ct - Cth)/Cth)*100 
% ct1 = 0;
% Use of gaussian quaderature to decide theta_tip for ideal twist with  ideal taper case
% for g = 1:length(r)-1
%     a = r(g);
%     b = r(g+1);
%     sum1 = 0;
%     for h = 1:6
%         ti1 = (((b-a)/2)* m(h))+ ((a+b)/2);
%         lambda1 = sigma1 *(Cla/16)* ((1+(32* theta_tip1/(sigma1* Cla)))^0.5 -1) ;
%         dct1 = 0.5* sigma1* Cla* ((theta_tip1*ti1) - lambda1*ti1^2)*W(h)* ((b-a)/2);
%         sum1 = sum1+ dct1
%     end
%     ct1 = sum1 + ct1;
% end
% error1 = ((ct1 - Cth)/Cth)*100     
% 
% lambda_hover = (Cth/2)^0.5;
% 
% theta_75 = ((2*Cth/(sigma*Cla))+ (lambda_hover/2))*3;

%lambda = sigma *(Cla/16)* ((1+(32* theta_tip/(sigma* Cla)))^0.5 -1) %constant lambda for theta tip
%lambda1 =sigma1 *(Cla/16)* ((1+(32* theta_tip/(sigma1* Cla)))^0.5 -1); %ideal twist with ideal taper


% dCt_by_dr = rand (1, 10);
% dCt_by_dr_var = rand (1, 10);
% dCt_by_dr_var1 = rand (1, 10);
% dCt_by_dr_idtap = rand (1,10);
% 
% phi = rand (1,10);
% phi_var = rand (1,10);
% phi_var1 = rand (1,10);
% phi_idtap = rand (1,10);
% 
% %theta = rand(1,10);
% theta_var = rand(1,10);
% theta_var1 = rand(1,10);
% theta_idtap = rand (1,10);
% 
% lambda_var = rand(1,10);
% lambda_var1 = rand(1,10);
% 
% alpha_idtap = rand(1,10);
% alpha_var = rand(1,10);
% alpha = rand(1,10);
% alpha_var1 = rand(1,10);

% dCqi_by_dr = rand(1,10);
% dCqi_by_dr_var = rand(1,10);
% dCqi_by_dr_var1 = rand(1,10);
% dCqi_by_dr_idtap = rand(1,10);
% 
% dCq0_by_dr = rand(1,10);
% dCq0_by_dr_var = rand(1,10);
% dCq0_by_dr_var1 = rand(1,10);
% dCq0_by_dr_idtap = rand(1,10);
% 
% dCq_by_dr = rand (1,10);
% dCq_by_dr_var = rand (1,10);
% dCq_by_dr_var1 = rand (1,10);
% dCq_by_dr_idtap = rand (1,10);
% 
% dCt_by_dr_tiploss = rand(1,10);
% dCt_by_dr_var_tiploss = rand(1,10);
% dCt_by_dr_var1_tiploss  = rand(1,10);
% dCt_by_dr_idtap_tiploss= rand(1,10);

% lambda_tiploss = rand(1,10);
% lambda_var_tiploss = rand(1,10);
% lambda_var1_tiploss = rand(1,10);
% lambda_idtap_tiploss = rand(1,10);
% 
% dCq_by_dr_tiploss = rand (1,10);
% dCq_by_dr_var_tiploss = rand (1,10);
% dCq_by_dr_var1_tiploss= rand (1,10);
% dCq_by_dr_idtap_tiploss = rand (1,10);
% 
% dCqi_by_dr_tiploss = rand(1,10);
% dCqi_by_dr_var_tiploss = rand(1,10);
% dCqi_by_dr_var1_tiploss = rand(1,10);
% dCqi_by_dr_idtap_tiploss = rand(1,10);
% 
% dCq0_by_dr_tiploss = rand(1,10);
% dCq0_by_dr_var_tiploss = rand(1,10);
% dCq0_by_dr_var1_tiploss = rand(1,10);
% dCq0_by_dr_idtap_tiploss = rand(1,10);

% f = rand(1,10);
% F = rand(1,10);
lambda = rand(1,128);
phi = rand(1,128);
alpha = rand(1,128);
CT = rand(1,17); 
Cla = rand(1,28);
cl = data(:, 2);
Theta = 0:1:16;
theta = rand(1,17);
P = zeros(1,17)
for j = 1:length(Theta)
        %for ideal twist with theta_tip
        %theta(i) = theta_tip/ r(i);  %pitch angle for theta tip 
        %dCt_by_dr(i) =  0.5* sigma*Cla*((theta^2* r(i)) - (lambda(i)*r(i)) );   %thrust distribution equation for theta tip
        
        theta(j) = Theta(j)*(pi/180);
        for i = 1:length(r)-1
            lambda(i) = sqrt(sigma* abs(cl(i))* r(i)/8);
            phi(i)= lambda(i) / r(i);   %aoa for theta tip
            alpha(i) = theta(j) - phi(i);
            Cla(i) = cl(i)/alpha(i);
            a = r(i);
            b = r(i+1);
            sum1 = 0;
            for h = 1:6
                ti1 = (((b-a)/2)* m(h))+ ((a+b)/2);
                %lambda(i) = sigma1 *(Cla(i)/16)* ((1+(32* theta/(sigma1* Cla)))^0.5 -1) ;
                dct1 = 0.5* sigma* Cla(i)* ((theta(j)^2*ti1) - lambda(i)*ti1)*W(h)* ((b-a)/2);
                sum1 = sum1+ dct1;
            end
            P(j) = sum1 + P(j)
        end
        CT(j) = P(j);
%         dCqi_by_dr(i) = 0.5*sigma* Cla* alpha(i)* lambda* (r(i)^2);
%         dCq0_by_dr(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr(i) = dCqi_by_dr(i) +  dCq0_by_dr(i);
%         
%         f(i) = (nb/2)* ((1-r(i))/(r(i)*sin(phi(i))));
%         F(i) = (2/pi)* acos(exp(-f(i)));
%         
%         lambda_tiploss(i) = (sigma*Cla/(16*F(i)))*((1+ (32*F(i)*theta(i)*r(i)/(sigma* Cla)))^0.5 - 1);
%         dCt_by_dr_tiploss(i) =  0.5* sigma*Cla*((theta_tip* r(i)) - (lambda_tiploss(i)*r(i)) );
%         dCqi_by_dr_tiploss(i) = 0.5*sigma* Cla* alpha(i)* lambda_tiploss(i)* (r(i)^2);
%         dCq0_by_dr_tiploss(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr_tiploss(i) = dCqi_by_dr_tiploss(i) +  dCq0_by_dr_tiploss(i);
%         
%         % for linear twist
%         theta_tw = 0;
%         theta_var(i) = theta_75 + (theta_tw*(r(i)-0.75));      %variable theta
%         lambda_var(i) = sigma *(Cla/16)* ((1+(32* theta_var(i)* r(i)/(sigma* Cla)))^0.5 -1);
%         phi_var (i) = lambda_var(i)/r(i);
%         alpha_var(i) = theta_var(i) - phi_var(i);
%         dCt_by_dr_var(i) =  0.5* sigma*Cla*((theta_var(i)* r(i)*r(i)) - (lambda_var(i)*r(i)) );
%         dCqi_by_dr_var(i) = 0.5*sigma* Cla* alpha_var(i)* lambda_var(i)* (r(i)^2);
%         dCq0_by_dr_var(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr_var(i) = dCqi_by_dr_var(i) +  dCq0_by_dr_var(i);
%         
%         lambda_var_tiploss(i) = sigma*Cla/(16*F(i))*((1+ (32*F(i)*theta_var(i)*r(i)/(sigma* Cla)))^0.5 - 1);
%         dCt_by_dr_var_tiploss(i) =  0.5* sigma*Cla*((theta_tip* r(i)) - (lambda_var_tiploss(i)*r(i)) );
%         dCqi_by_dr_var_tiploss(i) = 0.5*sigma* Cla* alpha_var(i)* lambda_var_tiploss(i)* (r(i)^2);
%         dCq0_by_dr_var_tiploss(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr_var_tiploss(i) = dCqi_by_dr_var_tiploss(i) +  dCq0_by_dr_var_tiploss(i);
%         
%         
%         %for linear twist with theta_tw = -15deg
%         theta_tw1 = -15 *pi/180;
%         theta_var1(i) = theta_75 + (theta_tw1*(r(i)-0.75));      %variable theta
%         lambda_var1(i) = sigma *(Cla/16)* ((1+(32* theta_var1(i)* r(i)/(sigma* Cla)))^0.5 -1);
%         phi_var1 (i) = lambda_var1(i)/r(i);
%         alpha_var1(i) = theta_var1(i) - phi_var1(i);
%         dCt_by_dr_var1(i) =  0.5* sigma*Cla*((theta_var1(i)* r(i)*r(i)) - (lambda_var1(i)*r(i)) );
%         dCqi_by_dr_var1(i) = 0.5*sigma* Cla* alpha_var1(i)* lambda_var1(i)* (r(i)^2);
%         dCq0_by_dr_var1(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr_var1(i) = dCqi_by_dr_var1(i) +  dCq0_by_dr_var1(i);
%         
%         lambda_var1_tiploss(i) = sigma*Cla/(16*F(i))*((1+ (32*F(i)*theta_var1(i)*r(i)/(sigma* Cla)))^0.5 - 1);
%         dCt_by_dr_var1_tiploss(i) =  0.5* sigma*Cla*((theta_tip* r(i)) - (lambda_var1_tiploss(i)*r(i)) );
%         dCqi_by_dr_var1_tiploss(i) = 0.5*sigma* Cla* alpha_var1(i)* lambda_var1_tiploss(i)* (r(i)^2);
%         dCq0_by_dr_var1_tiploss(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr_var1_tiploss(i) = dCqi_by_dr_var1_tiploss(i) +  dCq0_by_dr_var1_tiploss(i);
%         
%         
%         % for ldeal twist with ideal taper
%         theta_idtap(i) = theta_tip1/r(i);
%         phi_idtap(i) = lambda1/ r(i);
%         alpha_idtap(i) = theta_idtap(i) - phi_idtap(i);
%         dCt_by_dr_idtap(i) =  0.5* sigma1*Cla*((theta_tip1* r(i)) - (lambda1*r(i)) );
%         dCqi_by_dr_idtap(i) = 0.5*sigma* Cla* alpha_idtap(i)* lambda1* (r(i)^2);
%         dCq0_by_dr_idtap(i) = 0.5*sigma* CDo* (r(i)^3);
%         dCq_by_dr_idtap(i) = dCqi_by_dr_idtap(i) +  dCq0_by_dr_idtap(i);
%         
%         lambda_idtap_tiploss(i) = sigma1*Cla/(16*F(i))*((1+ (32*F(i)*theta_idtap(i)*r(i)/(sigma1* Cla)))^0.5 - 1);
%         dCt_by_dr_idtap_tiploss(i) =  0.5* sigma1*Cla*((theta_tip1* r(i)) - (lambda_idtap_tiploss(i)*r(i)) );
%         dCqi_by_dr_var1_tiploss(i) = 0.5*sigma1* Cla* alpha_idtap(i)* lambda_idtap_tiploss(i)* (r(i)^2);
%         dCq0_by_dr_var1_tiploss(i) = 0.5*sigma1* CDo* (r(i)^3);
%         dCq_by_dr_var1_tiploss(i) = dCqi_by_dr_idtap_tiploss(i) +  dCq0_by_dr_idtap_tiploss(i);
        
        
        
    
end

figure(1)
plot(theta, P, 'b--o')
% 
% %figure
% % plot(r, dCt_by_dr, r, dCt_by_dr_var, 'r',  r, dCt_by_dr_var1, 'g', r,  dCt_by_dr_idtap, 'k')
% % title ('variation of thrust distribution with non dimensional radial location')
% % xlabel ('non dimensional radial location(r)')
% % ylabel ('Thrust distribution')
% % legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg' , 'ideal twist with ideal taper')
% 
% % figure
% % plot(r, alpha, r , alpha_var, 'r', r, alpha_var1 ,'g', r, alpha_idtap , 'k')
% % title ('plot of AOA for constant lambda and closed form sloution')
% % xlabel ('non dimensional radial location(r)')
% % ylabel('AOA')
% % legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, theta, r, theta_var, 'r' , r, theta_var1, 'g', r , theta_idtap, 'k')
% title('variation of pitch angle with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('pitch angle')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, lambda, r, lambda_var, 'r', r, lambda_var1, 'g', r, lambda1, 'k')
% title('variation of induced inflow with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('induced inflow ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, dCqi_by_dr, 'b', r, dCqi_by_dr_var, 'r', r , dCqi_by_dr_var1, 'g', r, dCqi_by_dr_idtap, 'k')
% title('variation of induced torque distribution with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('induced torque distribution ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, dCq0_by_dr, 'b', r, dCq0_by_dr_var, 'r', r , dCq0_by_dr_var1, 'g', r, dCq0_by_dr_idtap, 'k')
% title('variation of profile torque distribution with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('profile torque distribution ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, dCq_by_dr, 'b', r, dCq_by_dr_var, 'r', r , dCq_by_dr_var1, 'g', r, dCq_by_dr_idtap, 'k')
% title('variation of total torque distribution with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('total torque distribution ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, lambda_tiploss, 'b', r, lambda_var_tiploss, 'r', r ,lambda_var1_tiploss , 'g', r, lambda_idtap_tiploss, 'k')
% title('variation of inflow with tip loss with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('inflow with tiploss ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% figure
% plot (r, dCt_by_dr_tiploss, 'b',r, dCt_by_dr_var_tiploss, 'r', r ,dCt_by_dr_var1_tiploss, 'g', r, dCt_by_dr_idtap_tiploss, 'k')
% title('variation of thrust distribution with tip loss with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('thrust distribution with tip loss ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')
% 
% 
% 
% figure
% plot (r, dCq_by_dr_tiploss, 'b',r, dCq_by_dr_var_tiploss, 'r', r ,dCq_by_dr_var1_tiploss, 'g', r, dCq_by_dr_idtap_tiploss, 'k')
% title('variation of torque distribution with tip loss with non dimensional radial location')
% xlabel ('non dimensional radial location(r)')
% ylabel('torque distribution with tip loss ')
% legend ('ideal twist', 'linear twist with 0 deg'  , 'twist with -15deg', 'ideal twist with ideal taper')

