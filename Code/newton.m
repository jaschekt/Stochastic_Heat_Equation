function x = newton(fun,x,tol)
%NEWTON this function does the following:
%Given an initial guess for a root of a function F, it uses Newton's method
%to find a better approximation of the analytic root. (Newton's method is
%known from basic calculus and will not be described here)

%input:
%   fun -   is a function handle. For given x it returns the value of F and DF at x. 
%   x -     x represents initial guess.
%   tol -   float number that is used as variable for tolerance. 
%           Our programm will consider a float as root of F if F(x)<tol (stopps the
%           while slope)
%output:
%   x -     new approximation for root of F

[F,DF] = fun(x); %compute F(X),DF(x) for initial guess of root

while norm(F,inf)>tol
    x = x-DF\F; %find root of tangent according to Newton method
    [F,DF] = fun(x); %compute F(x),DF(x) for the new approximation
end

