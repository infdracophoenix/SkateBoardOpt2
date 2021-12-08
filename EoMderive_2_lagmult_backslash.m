clc;close all;clear;

syms thetD thetDdot thetDdotdot thetB thetBdot thetBdotdot thet1 thet1dot thet1dotdot...
    thet2 thet2dot thet2dotdot thet3 thet3dot thet3dotdot ...
    x xdot xdotdot y ydot ydotdot thetDv thetBv thet1v thet2v thet3v yv xv  ...
    tau1 tau2 tau3 Fc Tc K C B m_t m_k g Ib t L1 L2 L3 L3cross alph V Vdot real
assume(g>0 & K>0 & m_t>0 & m_k>0 & Ib>0 & C>0 & B>0 & L1>0 & L2>0 & L3>0 & L3cross>0);
%%
Cords=[x y thetD thetB thet1 thet2 thet3];
Cordsdot=[xdot ydot thetDdot thetBdot thet1dot thet2dot thet3dot];
Cordsdotdot=[xdotdot ydotdot thetDdotdot thetBdotdot thet1dotdot thet2dotdot thet3dotdot];
state=[Cords Cordsdot];
statedot=[Cordsdot Cordsdotdot];
virtualDisp=[xv yv thetDv thetBv thet1v thet2v thet3v];
Vars=[Cordsdotdot Fc Tc];

rinv=C*thetB;
Iy=m_t*((sin(thetB+thet1-thet2+thet3)*L3/2)^2 + (cos(thetB+thet1-thet2+thet3)*L3cross/2)^2)/2;
Iz=m_t*((L3/2)^2 + (L3cross/2)^2)/2;
rx_knee=sin(thetB+thet1)*L1;
ry_knee=cos(thetB+thet1)*L1;
rx_torso=rx_knee+sin(thetB+thet1-thet2)*L2+sin(thetB+thet1-thet2+thet3)*L3/2;
ry_torso=ry_knee+cos(thetB+thet1-thet2)*L2+cos(thetB+thet1-thet2+thet3)*L3/2;


grad=@(y,x) gradient(y,x);

thetdot=thetBdot+thet1dot-thet2dot+thet3dot;

Vplanar_knee=[statedot*grad(rx_knee,state);...
               statedot*grad(ry_knee,state);
                            0];
Vtot_knee=Vplanar_knee+Rot(2,-thetD)*[xdot; 0; -ydot];

Vplanar_torso=[statedot*grad(rx_torso,state);...
               statedot*grad(ry_torso,state);
                            0];
Vtot_torso=Vplanar_torso+Rot(2,-thetD)*[xdot; 0; -ydot];

KE_torso=(m_t*sum(Vtot_torso.^2) + Iy*thetDdot^2 + Iz*thetdot^2)/2;
KE_knee=m_k*sum(Vtot_knee.^2)/2;
KE_board=(Ib*thetBdot^2)/2;

KE=KE_torso + KE_knee + KE_board;

PE_torso=m_t*g*(sin(alph)*(y + sin(thetD)*rx_torso)+cos(alph)*ry_torso);
PE_knee=m_k*g*(sin(alph)*(y + sin(thetD)*rx_knee)+cos(alph)*ry_knee);
PE_board=(1/2)*K*thetB^2;

PE=PE_torso + PE_knee + PE_board;

W = Tc*thetDv + tau1*thet1v + tau2*thet2v + tau3*thet3v ... 
     + Fc*(cos(thetD)*xv+sin(thetD)*yv) - B*(xdot*xv + ydot*yv) ;

Lag=KE-PE;

% constraint EQs
EQ(1,1)= thetDdot-C*thetB*([xdot ydot]*[sin(thetD) -cos(thetD)]');
EQ(1,1)= statedot*grad(EQ(1,1),state);
EQ(2,1)= xdot+tan(thetD)*ydot;
EQ(2,1)= statedot*grad(EQ(2,1),state);

% euler lagrange EQs

dLagdCords=grad(Lag,Cords);
dLagdCordsdot=grad(Lag,Cordsdot);
ddLagdCordsdotdt=jacobian(dLagdCordsdot,state)*statedot';
dWdvdisp=grad(W,virtualDisp);

EQ(3:9,1)=ddLagdCordsdotdt-dLagdCords-dWdvdisp;


M=jacobian(EQ,Vars);
F=EQ-M*Vars';

%M(1,:)=simplify(subs(M(1,:),[xdot ydot],V*[sin(thetD) -cos(thetD)]));

%%

EoMs=M\(-F);



%%
% EoMs_ZeroVel=limit(subs(EoMs,[xdot ydot],V*[sin(thetD) -cos(thetD)]),V,0);
% EoMs=simplify(EoMs);
%%

save('EoMs_no_singularities','EoMs');
