clear;close all; clc;

syms alph thetD delthetD thetDdot thetDdotdot thetB delthetB thetBdot thetBdotdot ...
    x delx xdot xdotdot y dely ydot ydotdot m Ib Id g Tc Fc C V real

state=[x y thetD thetB xdot ydot thetDdot thetBdot];
statedot=[xdot ydot thetDdot thetBdot xdotdot ydotdot thetDdotdot thetBdotdot];
Vars=[xdotdot ydotdot thetDdotdot thetBdotdot Fc Tc];


KE=(Id*thetDdot^2 + Ib*thetBdot^2 +m*(xdot^2 + ydot^2))/2;
PE=m*g*sin(alph)*y;
L=KE-PE;
W=Tc*delthetD + Fc*(cos(thetD)*delx+sin(thetD)*dely);

grad=@(y,x) gradient(y,x);

% constraint eqs
EQ(1,1)= thetDdot-C*thetB*((xdot^2+ydot^2)^(1/2));
EQ(1,1)= statedot*grad(EQ(1,1),state);
EQ(2,1)= xdot/ydot+tan(thetD);
EQ(2,1)= statedot*grad(EQ(2,1),state);

%euler lagrange eqs
%x
dLdx=grad(L,x);  dLdxdot=grad(L,xdot);
ddLdxdotdt=statedot*grad(dLdxdot,state);
dWddelx=grad(W,delx);

EQ(3,1)=ddLdxdotdt-dLdx-dWddelx;

%y
dLdy=grad(L,y);  dLdydot=grad(L,ydot);
ddLdydotdt=statedot*grad(dLdydot,state);
dWddely=grad(W,dely);

EQ(4,1)=ddLdydotdt-dLdy-dWddely;

%thetD
dLdthetD=grad(L,thetD);  dLdthetDdot=grad(L,thetDdot);
ddLdthetDdotdt=statedot*grad(dLdthetDdot,state);
dWddelthetD=grad(W,delthetD);

EQ(5,1)=ddLdthetDdotdt-dLdthetD-dWddelthetD;

%thetB
dLdthetB=grad(L,thetB);  dLdthetBdot=grad(L,thetBdot);
ddLdthetBdotdt=statedot*grad(dLdthetBdot,state);
dWddelthetB=grad(W,delthetB);

EQ(6,1)=ddLdthetBdotdt-dLdthetB-dWddelthetB;


M=jacobian(EQ,Vars);

F=EQ-M*Vars';

EoMs=M\(-F);


EoMsZeroVel=simplify(limit(subs(EoMs,[xdot ydot],V*[sin(thetD) -cos(thetD)]),V,0));




