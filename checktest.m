clc;close all;clear;

syms thet thetdot thetdotdot m g real 
state=[thet,thetdot];
statedot=[thetdot,thetdotdot];
KE=(1/2)*thetdot^2;
PE=-cos(thet)*m*g;
Lag=KE-PE;
dLagdthet=gradient(Lag,thet);
dLagdthetdot=gradient(Lag,thetdot);
ddLagdthetdotdt=statedot*gradient(dLagdthetdot,state);

EQ=ddLagdthetdotdt-dLagdthet;

soln.thetdotdot=solve(EQ,thetdotdot);
soln.gump=thet^2;
check=subs(EQ,soln);