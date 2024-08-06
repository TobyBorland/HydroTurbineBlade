//    Copyright 2012 Manolo Venturin, EnginSoft S.P.A.
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

// Exercice 1: Least square fitting

// Close all opened figures and clear workspace
xdel(winsid());
clear;
clc;

pathdir = get_absolute_file_path('ex1.sce');
exec(pathdir + "polyfit.sci");

// Exponential decay
deff('[Nt]=exp_dacay(t,N0,lambda)','Nt = N0.*exp(-lambda*t)');

np = 100;
noise = 0.5*(rand(np,1)-0.5);

t = linspace(0,1,np)';
N0 = 10;
lambda = 1;
Nt = exp_dacay(t,N0,lambda)
Ntnoise = Nt+noise;

// Estimate coefficient
exec("polyfit.sci",-1);
pcoeff = coeff(polyfit(t, log(Ntnoise), 1));

N0est = exp(pcoeff(1));
lambdaest = - pcoeff(2);
Ntest = N0est.*exp(-lambdaest*t);

// Plot
scf(1);
clf(1);
plot(t,Ntnoise,'k-');
plot(t,Nt,'b-');
// e = gce();
// e.children.thickness=3;
plot(t,Ntest,'r-');
xlabel("t");
ylabel("N(t)");
title("Exponential decay");
legend(["Data";"True";"Estimated"]);

