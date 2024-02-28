function p = t_to_p(t,n)
% two-sided
t = abs(t);
p = (1-tcdf(t,n-1))*2;