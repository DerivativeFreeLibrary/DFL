clear all
close all
clc

format longeng;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROBLEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n          = 4;
lbint      = zeros(n,1);
ubint      = 100.0*ones(n,1);
index_int  = [true; false; true; true];
is_integer = index_int;
step       = [1.0; 0.0; 1.0; 1.0];
scale_int  = step;
x0         = (ubint+lbint)./2;
startp     = [0.25; 0.39; 0.415; 0.39];
lb         = startp-10.0;
ub         = startp+10.0;
lbint(2)   = lb(2);
ubint(2)   = ub(2);
fhandle    = @(x)funex2(x,lb,ub,lbint,ubint,index_int);
fob        = feval(fhandle,x0);
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

[x,f,num_iter,nf,istop] = sd_box(@funex1,x0,index_int,scale_int,lbint,ubint,alfa_stop,nf_max,iprint);

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

