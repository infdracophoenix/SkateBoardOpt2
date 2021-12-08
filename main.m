%% ME 7385 Course Project
%
% Inspiration for multiple shooting taken from Matthew Kelly's trajectory 
% optimization tutorials: 
% https://github.com/MatthewPeterKelly/dscTutorials/blob/master/TrajectoryOptimization/Example_2_CartPole/MAIN. m
%% Initialize Code
% 
clc; clear; close;

%% Nonlinear Programming Problem
%

body_weight = 75; % body weight in kg taken from Archit's measurements
body_height = 1.65; % body height in m taken from Archit's measurements
Param.m_t = 0.551*body_weight; % mass in kg of total trunk taken from https://exrx.net/Kinesiology/Segments
Param.m_k = 0.167*body_weight; % mass in kg of total leg taken from https://exrx.net/Kinesiology/Segments
Param.alph = pi/6; % ramp angle in radians
Param.C = 6/(5*pi); % 
Param.g = 9.81; % gravitational acceleration in m/s^2
Param.K = 10.0; % board spring constant
Param.B = 10.0; % damping constant
Param.L1 = 0.63*body_height; % length of trunk in m taken from https://exrx.net/Kinesiology/Segments
Param.L2 = 0.433*body_height; % length of thigh in m taken from https://exrx.net/Kinesiology/Segments
Param.L3 = 0.434*body_height; % length of leg in m taken from https://exrx.net/Kinesiology/Segments
Param.L3cross = .15*Param.L3; % length in m cross link
Param.Ib = 1; % moment of inertia of the board  in kg*m^2
Param.avgVel = -0.6*(Param.m_t + Param.m_k)*sin(Param.alph)/Param.B;
Param.taugridNum=20;
Param.numShoot= 10;
% lower and upper bounds
Param.tauLim=200;
Param.xUB=12; Param.xLB=-12;
Param.xendUB=12; Param.xendLB=-12;
Param.yUB = 0; Param.yLB = -50;
Param.yendUB = -10; Param.yendLB = -50;
Param.thet1UB = pi/4; Param.thet1LB = -pi/4;
%Param.thet1endUB = pi/4; Param.thet1LB = -pi/4;
Param.thet2UB = 11*pi/12; Param.thet2LB = 0;
%Param.thet2endUB = 11*pi/12; Param.thet2endLB = 0;
Param.thet3UB = pi-atan(Param.L3cross/Param.L3); Param.thet3LB = -pi/4;
%Param.thet3endUB = pi-atan(L3cross/L3); Param.thet3endLB = -pi/4;

UB = []; LB = []; Aeq = []; Beq = []; Aineq = []; Bineq = [];

% initial guesses
taugrid0 = zeros(Param.taugridNum*3,1);
tdur0=Param.yendUB/Param.avgVel;
slackVars0= zeros(Param.numShoot*12,1);
decisionVars0=[taugrid0; tdur0; slackVars0];

options = optimset('display','iter','MaxFunEvals',200000,'MaxIter',20000,'diffmaxchange',1.1*1e5, ...
    'diffminchange',1e-5,'UseParallel',true);

% optimize 
tic;
[decisionVarsResult,optFVal] = fmincon(@skateboardObjective,decisionVars0,Aineq,Bineq,Aeq,Beq,LB,UB,@multipleShootingConstraints,options,Param);
optTime0 = toc;






