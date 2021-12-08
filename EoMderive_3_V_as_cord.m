clc;close all;clear;

syms thetD thetDdot thetDdotdot thetB thetBdot thetBdotdot thet1 thet1dot thet1dotdot...
    thet2 thet2dot thet2dotdot thet3 thet3dot thet3dotdot ...
    P V Vdot thetDv thetBv thet1v thet2v thet3v Pv  ...
    tau1 tau2 tau3 Fc Tc K C B m_t m_k g Ib t L1 L2 L3 L3cross alph real
assume(g>0 & K>0 & m_t>0 & m_k>0 & Ib>0 & C>0 & B>0 & L1>0 & L2>0 & L3>0 & L3cross>0);
%%



Cords=[P thetD thetB thet1 thet2 thet3];
Cordsdot=[V thetDdot thetBdot thet1dot thet2dot thet3dot];
Cordsdotdot=[Vdot thetDdotdot thetBdotdot thet1dotdot thet2dotdot thet3dotdot];
state=[Cords Cordsdot];
statedot=[Cordsdot Cordsdotdot];
virtualDisp=[Pv thetDv thetBv thet1v thet2v thet3v];
Vars=[Cordsdotdot Tc];

rinv=C*thetB;
Iy=m_t*((sin(thetB+thet1-thet2+thet3)*L3/2)^2 + (cos(thetB+thet1-thet2+thet3)*L3cross/2)^2)/2;
Iz=m_t*((L3/2)^2 + (L3cross/2)^2)/2;
rx_knee=sin(thetB+thet1)*L1;
ry_knee=cos(thetB+thet1)*L1;
rx_torso=rx_knee+sin(thetB+thet1-thet2)*L2+sin(thetB+thet1-thet2+thet3)*L3/2;
ry_torso=ry_knee+cos(thetB+thet1-thet2)*L2+cos(thetB+thet1-thet2+thet3)*L3/2;


grad=@(y,x) gradient(y,x);

thetdot=thetBdot+thet1dot-thet2dot+thet3dot;

V_knee=[statedot*grad(rx_knee,state);...
               statedot*grad(ry_knee,state);
                        V - thetDdot*rx_knee];

V_torso=[statedot*grad(rx_torso,state);...
               statedot*grad(ry_torso,state);
                            V - thetDdot*rx_torso];

KE_torso=(m_t*sum(V_torso.^2) + Iy*thetDdot^2 + Iz*thetdot^2)/2;
KE_knee=m_k*sum(V_knee.^2)/2;
KE_board=(Ib*thetBdot^2)/2;

KE=KE_torso + KE_knee + KE_board;

PE_torso=m_t*g*(sin(alph)*(-cos(thetD)*P + sin(thetD)*rx_torso)+cos(alph)*ry_torso);
PE_knee=m_k*g*(sin(alph)*(-cos(thetD)*P + sin(thetD)*rx_knee)+cos(alph)*ry_knee);
PE_board=(1/2)*K*thetB^2;

PE=PE_torso + PE_knee + PE_board;

W = Tc*thetDv + tau1*thet1v + tau2*thet2v + tau3*thet3v ... 
     - B*(V*Pv) ;

Lag=KE-PE;

% constraint EQs
EQ(1,1)= thetDdot-C*thetB*V;
EQ(1,1)= statedot*grad(EQ(1,1),state);

% euler lagrange EQs

dLagdCords=grad(Lag,Cords);
dLagdCordsdot=grad(Lag,Cordsdot);
ddLagdCordsdotdt=jacobian(dLagdCordsdot,state)*statedot';
dWdvdisp=grad(W,virtualDisp);

EQ(2:7,1)=ddLagdCordsdotdt-dLagdCords-dWdvdisp;
EQ=subs(EQ,P,0);

M=jacobian(EQ,Vars);
F=EQ-M*Vars';



%%

EoMs=simplify(M\(-F));



%%

save('EoMs_V_as_cord','EoMs');



