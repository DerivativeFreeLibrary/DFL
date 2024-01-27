clear all
close all
clc

format longeng;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROBLEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n          = 4;
lb         = zeros(n,1);
ub         = 10.0*ones(n,1);
index_int  = [false; true; false; true];
is_integer = index_int;
step       = [0.0; 0.5; 0.0; 0.5];
scale_int  = step;
x0         = (ub+lb)./2.0;
fob        = funex1(x0);
finiz      = fob;
num_funct  = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET ALGORITHM PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alfa_stop  = 1.e-3;
nf_max     = 20000;
iprint     = -1;
   
fprintf(' -------------------------------------------------\n');
fprintf(' objective function at xo:\n');
fprintf(' %26.16e\n',fob);
fprintf(' -------------------------------------------------\n');
for i = 1:n
    fprintf(' x(%3d) = %26.16e %26.16e %1d\n',i,x0(i),step(i),is_integer(i));
end
fprintf(' -------------------------------------------------\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% NOW CALL THE OPTIMIZER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,f,num_iter,nf,istop] = sd_box(@funex1,x0,index_int,scale_int,lb,ub,alfa_stop,nf_max,iprint);

fprintf(' -------------------------------------------------\n');
if (istop == 1)
	fprintf(' END - stopping criterion satisfied \n');
end
if (istop == 2)
    fprintf(' END - maximum number of function calculation reached  = %d\n',nf_max);
end
if (istop == 3)
    fprintf(' END -  maximum number of iterations reached = %d\n',nf_max);
end

fprintf(' -------------------------------------------------\n');
fprintf(' objective function at x*:\n');
fprintf(' %26.16e\n',f);
fprintf(' -------------------------------------------------\n');
for i = 1:n
    fprintf('x*(%3d) = %26.16e \n',i,x(i));
end
fprintf(' -------------------------------------------------\n');

