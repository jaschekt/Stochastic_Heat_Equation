% heat.m
%clear; close all; clc;

% load and mesh the domain
%msh = load('pi.mat','G');
%[msh.P,msh.E,msh.T] = initmesh(msh.G,'Hmax',1e-2);
msh.P=p;
msh.E=e;
msh.T=t;
% data
a = 1e-1;
tSpan = [0 0.5] ;
u0 = (sin(pi*msh.P(1,:)).*sin(pi*msh.P(2,:)))';
GammaD = @(x1,x2) true(size(x1));

%__________________________________________________________________________
%backwardEuler 10 steps
dt = 1e-2/2;
tSteps = (tSpan(2)-tSpan(1))/dt;

fh = 50*randn(tSteps,length(msh.P(:,msh.T(1,:))));
%f = @(t) fh(round(t/dt),:);
f = @(t) zeros(length(fh(round(t/dt),:)),1)';
% set up the Heat equation
[F,M,Pf,PD,process] = discretiseHeat(a,f,GammaD,msh,true,'tim2');
u0 = Pf*u0;
% time stepping
backwardsEuler(F,M,tSpan,tSteps,u0,1e-9,process);


