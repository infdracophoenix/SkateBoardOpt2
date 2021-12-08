clc;close all;clear;

syms thetD thetB thetBdot thet1 thet1dot thet2 thet2dot thet3 thet3dot V Vdot ...
    thetBdotdot thet1dotdot thet2dotdot thet3dotdot ...
    tau1 tau2 tau3 K C B m g Ib t P L1 L2 L3 alph real
%%

rinv=C*thetB;
Iy=m*(sin(thetB+thet1-thet2+thet3)*L3/2)^2;
Iz=m*(L3/2)^2;
rx=sin(thetB+thet1)*L1+sin(thetB+thet1-thet2)*L2+sin(thetB+thet1-thet2+thet3)*L3/2;
ry=cos(thetB+thet1)*L1+cos(thetB+thet1-thet2)*L2+cos(thetB+thet1-thet2+thet3)*L3/2;
grad=@(y,x) gradient(y,x);
thet=[thetB thet1 thet2 thet3];
thetdot=[thetBdot thet1dot thet2dot thet3dot];
thetdotdot=[thetBdotdot thet1dotdot thet2dotdot thet3dotdot];
state=[thet V thetdot];
statedot=[thetdot Vdot thetdotdot];
tau=[tau1 tau2 tau3];

Vp=[thetdot*grad(rx,thet);...
    thetdot*grad(ry,thet)];

KE=(m*(V-V*rinv*rx)^2 + Iy*(V*rinv)^2 + Iz*(sum(thetdot)-2*thet2dot)^2 ...
    + m*sum(Vp.^2) + Ib*thetBdot^2)/2;
PE=m*g*(-sin(alph)*(cos(thetD)*P - sin(thetD)*rx)+cos(alph)*ry) ... This line has the gravitational PE
    - sum(tau.*thet(2:4))+ (1/2)*K*thetB^2 + V*B*P;% the partials of this line give the forces/torques from friction, the board spring, and the torques at the joints
Lag=KE-PE;

% eq for P and V (P is the arc length)
dLagdP=grad(Lag,P)+grad(Lag,thetD)*rinv;% dthetD/dP=rinv
dLagdP=subs(dLagdP,P,0);%P is a virtual displacement along the arc, so we set it to zero
dLagdV=grad(Lag,V);
dLagdV=subs(dLagdV,P,0);
ddLagdVdt=statedot*grad(dLagdV,state);

EQ(1)=simplify(simplify(ddLagdVdt-dLagdP));

% eqs for all thetas except thetD

dLagdthet=grad(Lag,thet);
dLagdthetdot=grad(Lag,thetdot);
ddLagdthetdotdt=jacobian(dLagdthetdot,state)*statedot';

EQ(2:5,1)=simplify(simplify(ddLagdthetdotdt-dLagdthet));

%%
EoMs=solve(EQ,[Vdot thetdotdot]);
EoMs.thetDdot=V*rinv;
EoMs.Xdot=V*sin(thetD);
EoMs.Ydot=-V*cos(thetD);



