-----------------------------------------------------------
 How to use the derivative-free optimizer DFL for MINLP
-----------------------------------------------------------

1- Gunzip and untar the archive in a folder on your computer.

-----------------------------------------------------------
   FORTRAN SOURCE
-----------------------------------------------------------

2- Edit file parameter.f to set the number of variables.

3- Edit file wrap_mixed.f90 to define your own MINLP problem.
   In particular, subroutines:
   INIZPAR(n,x,lb,ub)  
	sets initial point and lower and upper bounds

   which_integer(n,lb,ub,is_integer,step,x)
	defines which variables are integers and the 
	allowed step size for those variables
	also adjust initial point to respect integrality
	constraint

   funct(n,x,f)
	defines the objective function

4- At command prompt execute 

     $> make
 
   which will create the executable 'dfl'

5- execute

     $> ./dfl

-----------------------------------------------------------
   MATLAB SOURCE
-----------------------------------------------------------

You can find the matlab version of the package in folder MATLAB.
The folder contains the following files:

- sd_box.m		(optimizer routine)
- example1.m    (main script to run example 1)
- funex1.m 		(function used in example 1)
- example2.m	(main script to run example 2)
- funex2.m		(function used in example 2)

