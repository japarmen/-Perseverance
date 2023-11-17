clear all;
initPlots;

% Simulation parameters
g = 3.7;
m = 1250 ;
beta = 0.1132 ;
Jx = 1250 * (1.4^2) ;
Jy = 1250 * (1.5^2) ;
Jz = 1250 * (1.25^2) ;
J = diag([ Jx , Jy , Jz ]);
zI = [0;0;1];
nx = 12;
ny = 12;
x0 = [0; 0; -2100; 0 ; 0 ; 0 ; 0; 0; 89; 0; 0; 0];
Dt = 0.01;
t = 0:Dt:35;
T1 =0.4*m*g ;% small constant thrust
T2 = T1 ;
T3 = T1 ;
T4 = T1 ;
ang = 20*pi/180 ;
r1x = 1.4 ; r1y = 1.5 ; r1z =1.25 ;
r2x = 1.4 ; r2y = -1.5 ; r2z =1.25 ;
r3x = -1.4 ; r3y = 1.5 ; r3z =1.25 ;
r4x = -1.4 ; r4y = -1.5 ; r4z =1.25 ;

fp1 = T1*[  0 ; -sin(ang) ; -cos(ang)] ;
fp2 = T2*[ 0 ; sin(ang) ; -cos(ang)] ;
fp3 = T3*[0 ; -sin(ang) ; -cos(ang)] ;
fp4 = T4*[0 ; sin(ang) ; -cos(ang)] ;

fp = fp1 + fp2 + fp3 + fp4 ;

p1 = [ r1x , r1y , r1z ]  ;
p2 = [ r2x , r2y , r2z ]  ;
p3 = [ r3x , r3y , r3z ]  ;
p4 = [ r4x , r4y , r4z ]  ;

np1 = skew(p1)*fp1 ;
np2 = skew(p2)*fp2 ;
np3 = skew(p3)*fp3 ;
np4 = skew(p4)*fp4 ;

np = np1 + np2 + np3 + np4 ;

u_NL = [ fp ; np ] * ones(size(t));

% simulate nonlinear system
Nsim = length(t);
x = zeros(nx,Nsim);
y = zeros(ny,Nsim);
x(:,1) = x0;

for k = 1:Nsim
    % prepare variables:
    p   = x(1:3,k);
    lbd  = x(4:6,k);
    v = x(7:9,k);
    om  = x(10:12,k);
    R = Euler2R(lbd);
    T = u_NL(3,k);
    np = u_NL(4:6,k);
    
    % compute state derivative:
    p_dot = Euler2R(lbd)*v;
    v_dot = -skew(om)*v + g*R'*zI + (1/m)*beta * v.^2 + (1/m)*fp ;
    lbd_dot = Euler2Q(lbd)*om;
    om_dot = -inv(J)*skew(om)*J*om + inv(J)*np;
    x_dot = [p_dot;lbd_dot;v_dot;om_dot];

    % integrate state
    x(:,k+1) = x(:,k) + Dt*x_dot;

    % compute current output:
    y(:,k) = x(:,k);
end

figure (1)
plot(t,u_NL);
grid on;


figure(2);
plot(t,y(1:2,:),'-.');
grid on;
legend('$$p_x$$ (nlin)','$$p_y$$ (nlin)')

figure(3);
plot(t,y(3,:),'-.',t,y(6,:), '-.');
grid on;
legend( '$$p_z$$ (nlin)' , '$$\psi$$ (nlin)');

vx = 2 ;
vy = 4 ;
vz = 0 ;
theta = 10*pi/180 ;
fi = 5*pi/180 ;
psi = 20*pi/180 ;
px = 0 ;
py = 0 ;
pz = - 21.3 ; 
wx =0 ;
wy =0 ;
wz = 0 ;
xop1 = [px ; py ; pz ; fi ; theta ; psi ; vx ; vy ; vz ; wx ; wy ; wz ] ;

A1 = zeros (3,3) ; %f1/p

A2 = [ vy*(cos (fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi))+vz*(-sin(fi)*sin(theta)*cos(psi)+cos(fi)*sin(psi)) , vx*(-sin(theta)*cos(psi))+vy*(sin(fi)*cos(theta)*cos(psi))+vz*(cos(fi)*cos(theta)*cos(psi)) , vx*(-cos(theta)*sin(psi))+vy*(-sin(fi)*sin(theta)*sin(psi)-cos(fi)*cos(psi))+vz*(-cos(fi)*sin(theta)*sin(psi)+sin(fi)*cos(psi)) ;
       vy*(cos(fi)*sin(theta)*sin(psi)-sin(fi)-cos(psi))+vz*(-sin(fi)*sin(theta)*sin(psi)-cos(fi)*cos(psi)) , vx*(-sin(theta)*sin(psi))+vy*(sin(fi)*cos(theta)*sin(psi))+vz*(cos(fi)*cos(theta)*sin(psi)) , vx*(cos(theta)*cos(psi))+vy*(sin(fi)*sin(theta)*cos(psi)-cos(fi)*sin(psi))+vz*(cos(fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi)) ;
       cos(theta)*(vy*cos(fi)-vz*sin(fi)) , -vx*cos(theta)+vy*(sin(fi)*sin(theta))+vz*(cos(fi)*sin(theta)) , 0 ] ; %f1/lambda

A3 = [ cos(theta)*cos(psi) , sin(fi)*sin(theta)*cos(psi)-cos(fi)*sin(psi) , cos(fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi) ;
       cos(theta)*sin(psi) , sin(fi)*sin(theta)*sin(psi)+cos(fi)*cos(psi) , cos(fi)*sin(theta)*sin(psi)-sin(fi)*cos(psi) ;
       -sin(theta) , sin(fi)*cos(theta) , cos(fi)*cos(theta) ] ; %f1/v

A4 = zeros (3,3) ; %f1/w

A5 =  zeros (3,3) ; %f2/p

A6 = [wy*tan(theta)*cos(fi) - wz*tan(theta)*sin(fi) , wy*sin(fi)* (sec (theta))^2 + wz*cos(fi)*(sec (theta))^2 , 0 ;
      -wy*(sin(fi)-wz*cos(theta)) , 0 , 0 ;
      1/cos(theta)*(wy*cos(fi)-wz*sin(theta)) , (sec(theta)*tan(theta))*(wy*sin(fi)+wz*cos(theta)) , 0 ] ; %f2/lambda

A7 = zeros (3,3) ; %f2/v

A8 = [ 1 , sin(fi)*tan(theta) , cos(fi)*tan(theta) ;
       0 , cos(fi) , -sin(theta) ;
       0 , sin(fi)/cos(theta) , cos(fi)/cos(theta) ] ; % f2/w

A9 = zeros (3,3) ; %f3/p

A10 = [0 , -cos(theta)*g, 0;
       cos(fi)*cos(theta)*g, -sin(theta)*sin(fi)*g, 0;
       -sin(fi)*cos(theta)*g, -sin(theta)*cos(fi)*g, 0]; %f3/lambda

A11 = [2*beta*vx*(1/m), wz, -wy;
       -wz, 2*beta*vy*(1/m), wx;
       wy, -wx, 2*beta*vz*(1/m)]; %f3/v

A12 = [0, -vz, vy;
       vz, 0, -vx;
       -vy, vx, 0]; %f3/w


A13 = zeros(3,3);

A14 = zeros(3,3);

A15 = zeros(3,3);

A16 = [0, wz*(Jy-Jz)*(1/Jx), wy*(Jy-Jz)*(1/Jx);
wz*(Jz-Jx)*(1/Jy), 0, wx*(Jz-Jx)*(1/Jy);
wy*(Jx-Jy)*(1/Jz), wx*(Jx-Jy)*(1/Jz), 0] ;

A = [A1 , A2 , A3 , A4 ;
     A5 , A6 , A7 , A8 ;
     A9 , A10 , A11 , A12 ;
     A13 , A14 , A15 , A16 ] 

 B1 = zeros (3,1) ; 
 
 B2 = zeros (3,1) ;
 
 B3 = zeros (3,1) ;
 
 B4 = zeros (3,1) ;
 
 B5 =  zeros (3,1) ;
 
 B6 =  zeros (3,1) ;
 
 B7 =  zeros (3,1) ;
 
 B8 =  zeros (3,1) ;
 
 B9 = 1/m*[0 ; -T1*sin(ang) ; -T1*cos(ang)] ;
 
 B10 = 1/m*[0 ; T2*sin(ang) ; -T2*cos(ang)] ;
 
 B11 = 1/m*[0 ; T3*sin(ang) ; -T3*cos(ang)] ;
 
 B12 =  1/m*[0 ; -T4*sin(ang) ; -T4*cos(ang)] ;
 
 B13 = [-T1/Jx*(r1z*sin(ang)*(-1)+r1y*cos(ang)) ; T1/Jy * (r1x*cos(ang)) ; T2/Jz * (r2x*sin(ang)*(-1)) ] ;
 
 B14 = [-T2/Jx*(r2z*sin(ang)+r2y*cos(ang)) ; T2/Jy * (r2x*cos(ang)) ; T2/Jz * (r2x*sin(ang)) ] ;
 
 B15 = [-T3/Jx*(r3z*sin(ang)*+r3y*cos(ang)) ; T3/Jy * (r1x*cos(ang)) ; T3/Jz * (r3x*sin(ang)) ] ;
 
 B16 =[-T4/Jx*(r4z*sin(ang)*(-1)+r4y*cos(ang)) ; T4/Jy * (r4x*cos(ang)) ; T4/Jz * (r4x*sin(ang)*(-1)) ] ;
 
B = [B1 , B2 , B3 , B4 ;
     B5 , B6 , B7 , B8 ;
     B9 , B10 , B11 , B12 ;
     B13 , B14 , B15 , B16 ] 

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);

u_L = [ 2.4*T1 ; 2.4*T2 ; 2.4*T3 ; 2.8*T4 ]*(t>=0);

% simulate linear system:
y_L = lsim(sys,u_L,t,xop1)';



[Vj,Jor] = jordan(A)
[V,DL,W] = eig(A);
mode_obs = C*V
mode_ctrl = W'*B


