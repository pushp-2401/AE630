%%%% AE630 HW-2 submission
clear
clc
Payload = 0.200;% payload of quadrotor in kg 
gtow = 0.600; %gross take off weight of quadrotor in kg
rho = 1.225;    %density of air
E = 1/3;       % endurance in HOURS
S = 3;        % number of cells in battery
%theta = that_0 + 17

%inputs assumed
nr = 4;   %number of rotor
c = 0.09;  %chord length
nb = 2;   %since it is a quadrotor
eta_motor = 0.8;   % given in question,small quadrotor has motor efficiency in this range
dl = 55;  %small quadrotor has disc loading in the range between 40-90
R =  0.01; % rotor radius is assumed to be 50mm
sigma = nb * c/ pi*R;  %soldity value
theta = 50* pi/180;
cl_alpha = 5.73;
cd0 = 0.01;
capacity = 3; %capacity of battery, since endurance is 20 mins thus capacity of battery will be 1/(1/3hour) = 3

r = linspace ( 0, 1, 100);
W = [0.171324492,0.360761573,0.467913935,0.467913935,0.360761573,0.171324492];
m = [0.932469514, 0.661209386, 0.02386191860,-0.02386191860,-0.661209386,-0.932469514];
ct = 0;
prev_weight = 0.600;
current_weight = 0;



while (  abs (prev_weight -  current_weight)  >= 0.005)
    if (current_weight ~= 0)
        prev_weight = est_weight;
    end
    for i = 1:length(r)-1
        
        
        lambda= sigma* cl_alpha*abs(sqrt((1+32*theta* r(i)/sigma* cl_alpha)-1)/16);
        %phi(i)= lambda(i) / r(i);   %aoa for theta tip
        %alpha1(i) = theta(j) - phi(i);
        %Cla(i) = cl(i)/alpha(i);
        
        a = r(i);
        b = r(i+1);
        sum1 = 0;
        sum2=0;
        for h = 1:6
            ti = (((b-a)/2)* m(h))+ ((a+b)/2);
            dct = 0.5* sigma* cl_alpha* ((theta^2*ti) - lambda*ti)*W(h)* ((b-a)/2);
            %dT = dct1* rho* pi*R^2*omega^2*R^2;
            %dcp = 0.5* sigma* (phi(i)* cl(i)+ (cd0+k*cl(i)^2))*ti^2;
            %dP =  dcp* rho * pi* R^2* omega^3* R^3;
            %sum1 = sum1+ dT;
            sum2 = sum2 + dct;
        end
        ct = sum2 + ct;
        %T = sum1 + T;
        
    end
    
    %------------ calculation of rotational speed---------------------------------------------------------------------------
    T = prev_weight*9.8/(0.8 * nr);
    omega_for_hover = abs (sqrt (T/ct*rho* pi* R^4));
    omega_for_design = abs (sqrt (2*T/ct*rho* pi* R^4));
    P = T* omega_for_hover;
    I = P/(3.7*S);
    Pmax = 2* T* omega_for_design;
    Imax = Pmax / (3.7* S);
    %capacity = P/3.7*S*3;
    
    rotormass = 0.0195 * R^2.0859 * sigma^-0.02038 * nb^0.5344;         % mass of rotor for design
    batterymass = 0.0418* capacity^0.9327*S^1.0725;                     %mass of battery
    
    %------------calculation for motor-----------------------------------------------------------------------
    kv = omega_for_hover / (3.7*S) ;  %speed constant for motor;
    l_bl = 4.8910* I^0.01751* P^0.2476;  %motor casing length;
    d_bl =  41.45* kv^-0.1919 * P^0.1935;
    motormass = 0.0109* kv^0.5122 * P^-0.01902 * (log10(l_bl))^2.5582  * (log10(d_bl))^12.8502;   %mass of battery
    escmass = 0.977 * Imax^0.8483;
    airframemass = 1.3119* R^1.2767* (batterymass)^0.4587;
    est_weight = rotormass + batterymass + motormass + escmass + airframemass + Payload;
    current_weight = est_weight;
end


table ( est_weight ,prev_weight , current_weight,rotormass, batterymass, motormass, escmass, airframemass )
table ( Payload , R ,dl , nr , nb , c , theta , cl_alpha , P , I , Pmax , Imax , kv , capacity  )

