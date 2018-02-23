%Project :Numerical Approximation using Adam-Bashforth fourth order method 
%Engineer : Nikhil.P.Lokhande
%email:nikhil.l.1@aol.com

clc
clf
clear all

%    AB4 method

% Solve dy/dt = f(t,y). In general it can be a function of both 
% variables t and y. If your function is only a function of t then
% you will need to add a 0*y to your function.

%  define f(t,y)
   fcnstr='1+(y/t)' ;
   f=inline(fcnstr,'t','y') ;

% t0, initial time
   t0=1;
% y0, corresponding value of y at t0
   y0=2;
% tf, upper boundary of time t (end time). 
   tf=2;
% n, number of steps to take
   n= 10;
   h=(tf-t0)/n;
%
ta(1)= t0;
w(1) = y0;
% using RK4
for i=1:n
  
% Adding Step Size
  ta(i+1)=ta(i)+h ;

% Calculating k1, k2, k3, and k4
  rk1 = f(ta(i),w(i)) ;
  rk2 = f(ta(i)+0.5*h,w(i)+0.5*rk1*h) ;
  rk3 = f(ta(i)+0.5*h,w(i)+0.5*rk2*h) ;
  rk4 = f(ta(i)+h,w(i)+rk3*h) ;
     
% Using 4th Order Runge-Kutta formula
  w(i+1)=w(i)+1/6*(rk1+2*rk2+2*rk3+rk4)*h ;
end

for i=1:n
  
% Adding Step Size
  ta(i+1)=ta(i)+h;
  if(i>3)
% Calculating k1, k2, k3, and k4
  k1 = f(ta(i),w(i));
  k2 = f(ta(i-1),w(i-1));
  k3 = f(ta(i-2),w(i-2));
  k4 = f(ta(i-3),w(i-3));
 
     
% Using 4th Order Adam-Bashforth formula
   w(i+1)=(w(i)+h/24*(55*k1-59*k2+ 37*k3-9*k4));
  end
end



% % The following finds what is called the 'Exact' solution by ode45
% tspan = [t0 tf];
% [t,y]=ode45(f,tspan,y0);
% [yfi dummy]=size(y);
% yf=y(yfi);
 t=ta;
 y=(t .* (log(t)+2)); 

% Plotting the Exact and Approximate solution of the ODE.
hold on
xlabel('x');ylabel('y');
title('Exact and Approximate Solution of the ODE by the 4th Order Adam-Bashforth Method');
plot(t,y,'-','LineWidth',2,'Color',[1 0 0]);            
plot(ta,w,'--','LineWidth',2,'Color',[0 0 1]);
legend('Exact','Approximation');

for j = 1:n+1, 
    disp(sprintf('t=%f,  w=%0.11f,  y(t)=%0.11f',t(j),w(j),y(j)));
end

