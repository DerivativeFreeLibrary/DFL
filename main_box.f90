!============================================================================================
!    DFL - Derivative-Free Linesearch program for Mixed Integer Nonlinear Programming 
!    Copyright (C) 2011  G.Liuzzi, S.Lucidi, F.Rinaldi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G.Liuzzi, S.Lucidi, F.Rinaldi. Derivative-free methods for bound constrained 
!    mixed-integer optimization, Computational Optimization and Applications, 
!    53:505-526 (2012). DOI: 10.1007/s10589-011-9405-3
!============================================================================================

!-----------------------------------------------------------------------------!
!          Program DFL for bound-constrained mixed integer problems           !
!-----------------------------------------------------------------------------! 
!
!
!We solve the following problem:
!                                
!                            min f(x)
!                             l<=x<=u            
!                            
!                            l>-\infty
!                            u<\infty
!
!The file "parameter.f" contains the number of variables n
!
!
!-----------------------------------------------------------------------------!

program main_box_discr
      implicit none

!parameter.f: file containing the number of variables of the problem
include 'parameter.f'
	  
!-----------------------------------------------------------------------------	  
	  integer	:: i, istop, icheck, index_int(n)
	  real		:: tbegin, tend
	  logical :: is_integer(n)

      real*8 :: x(n),bl(n),bu(n),scale_int(n),step(n)


!-----------------------------------------------------------------------------

	  integer ::            num_funct,num_iter
	  real*8             :: f,alfamax,delta 
	  real*8			   :: fob,fapp
	  real*8			   :: violiniz, violint, finiz, fex
      real*8             :: alfa_stop
      integer            :: nf_max,iprint

	  common /num/f
  	  common /calfamax/alfamax
!------------------------------------------------------------------------------

 	  call cpu_time(tbegin)




!-----------------------------------------------------------------------
!      Starting point and bound calculation
!-----------------------------------------------------------------------

	  call inizpar(n,x,bl,bu)

	  do i=1,n

         if((x(i).lt.bl(i)).or.(x(i).gt.bu(i))) then
		   write(*,*) ' punto iniziale non in box'
		   stop
		 endif
     
	  enddo

2002  format(2d20.10)
198   format(1x,a3,i3,a4,es26.16,1x,es26.16,1x,l1)
199   format(1xes26.16)






!-----------------------------------------------------------------------
!     show starting point info
!-----------------------------------------------------------------------


      index_int=0
	  scale_int= 0.d0

	if(.true.) then  ! ATTENZIONE RIMUOVERE DOPO DEBUG
	  call which_integer(n,bl,bu,is_integer,step,x)
	 
      do i = 1,n
		if(is_integer(i)) then
		  index_int(i)=1
		  scale_int(i) = step(i) 
		endif
      enddo
	endif

!-----------------------------------------------------------------------
!     calculate starting point violation 
!-----------------------------------------------------------------------

    violiniz=0.0d0
    
	do i=1,n
       violiniz=max(violiniz,x(i)-bu(i),bl(i)-x(i)) 
	end do
    

	  	call funct(n,x,fob)


        finiz=fob
   
        write(*,*) '------------------------------------------------- '

        write(*,*) 'objective function at xo:'
	    write(*,199) fob

!       ---- objective function value ----

        write(*,*) '------------------------------------------------- '


!      ---- x(i) =  i-th variable value----

	  do i=1,n
		write(*,198) 'xo(',i,')=',x(i),step(i),is_integer(i)
	  enddo

        write(*,*) '------------------------------------------------- '
!------------------------------------

	   num_funct   = 1 
	   num_iter    = 0
	   alfa_stop   = 1.d-3
  	   nf_max      = 20000
       iprint      = 0

       call sd_box(n,index_int,scale_int,x,f,bl,bu,alfa_stop,nf_max,num_iter,num_funct,iprint,istop)

	   call funct(n,x,fob) 

      !-----------------------------------------------------------------------
      !     integrality constraints violation for x* (last solution found) 
      !-----------------------------------------------------------------------

		violint=0.0d0
    
		do i=1,n

		   if(is_integer(i)) then
		   		if((bl(i) > -1.d+10).and.(bu(i)) < 1.d+10) then
					     violint=max( violint,abs(x(i)-bl(i)-( floor( ( x(i)-bl(i))/step(i)+0.5d0 )*step(i) ) )  )
		  
				elseif(bl(i) > -1.d+10) then
					 violint=max( violint,abs(x(i)-bl(i)-( floor( ( x(i)-bl(i))/step(i)+0.5d0 )*step(i) ) )  )
		  
				elseif(bu(i) <  1.d+10) then
					violint=max( violint,abs(bu(i)-x(i)-( floor( (bu(i)-x(i))/step(i)+0.5d0 )*step(i) ) )  )
		  
				else
				   violint=max( violint,abs(x(i)-( floor( x(i)/step(i)+0.5d0 )*step(i) ) )  )
		  
				endif
			   
			endif   	   
		  
		end do


      !-----------------------------------------------------------------------

	   call cpu_time(tend)

!       write(*,987) n,violiniz,finiz,fob,violint,num_funct,num_iter
! 
! 987 format(' & ', i3,' & ',es9.2,' & ',es14.7,' & ',es14.7,' & ',es9.2,' & ',i5,' & ',i5,'\\')

      
	   write(*,*) '------------------------------------------------------------------------------'     
	   if(istop.eq.1) then
         write(*,*)  ' END - stopping criterion satisfied '
	   endif
       if(istop.eq.2) then
         write(*,*)  ' END - maximum number of function calculation reached  =',nf_max
	   endif

       if(istop.eq.3) then
         write(*,*)  ' END -  maximum number of iterations reached =',nf_max
	   endif

	   write(*,*) '------------------------------------------------------------------------------'  


        write(*,*) ' objective function = ',fob

!      ---- x(i) = i-th variable ----

        write(*,*) ' ------------------------------------------------- '
	    do i=1,n
		   write(*,*) 'x(',i,') =',x(i)
        enddo

!      ---- nftot = number of function evaluations ----

        write(*,*) ' ------------------------------------------------- '
		write(*,*)
	   
end program main_box_discr

      subroutine sd_box(n,index_int,scale_int,x,f,bl,bu,alfa_stop,nf_max,ni,nf,iprint,istop)
      implicit none
	  logical :: cambio_eps
      integer :: n,i,j,i_corr,nf,ni,nf_max,index_int(n)
      integer :: num_fal,istop
      integer :: iprint,i_corr_fall
	  integer :: flag_fail(n)

      real*8 :: x(n),z(n),d(n)
      real*8 :: alfa_d(n),alfa,alfa_max, alfa_d_old
      real*8 :: f,fz , eta
	  real*8 :: bl(n),bu(n),alfa_stop,maxeps,scale_int(n) 
	  logical:: discr_change

!     values of f calculated on a n+1 simplex

      real*8 :: fstop(n+1)


!     num_fal number of failures

!     i_corr is the index of the current direction


!     initialization

	  discr_change = .false. 

	  eta = 1.d-6

      flag_fail=0

	  num_fal=0

      istop = 0

      fstop=0.d0

!     ---- choice of the starting stepsizes along the directions --------
	  alfa_max = 0.d0
      do i=1,n

        if(index_int(i).eq.0) then
        
           alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
		   alfa_max  = max(alfa_max,alfa_d(i))
      
		else

		   alfa_d(i)=dmax1(scale_int(i),dmin1(2.d0*scale_int(i),dabs(x(i))))
      
		end if
        if(iprint.ge.1) then
          write(*,*) ' alfainiz(',i,')=',alfa_d(i)
        endif
      end do

      do i=1,n      
        d(i)=1.d0 
      end do
!     ---------------------------------------------------------------------  
     
      
      call funct(n,x,f)
	  nf=nf+1

	  i_corr=1

      fstop(i_corr)=f

      do i=1,n
	    z(i)=x(i)
      end do

      if(iprint.ge.2) then
        write(*,*) ' ----------------------------------'
        write(*,*) ' finiz =',f
        do i=1,n
          write(*,*) ' xiniz(',i,')=',x(i)
        enddo
        write(*,*) ' ----------------------------------'
      endif

!---------------------------   
!     main loop
!---------------------------

      do 

         if(iprint.ge.0) then
         
           write(*,100) ni,nf,f,alfa_max
100        format(' ni=',i4,'  nf=',i5,'   f=',d12.5,'   alfamax=',d12.5)
         endif
         if(iprint.ge.2) then
	       do i=1,n
                write(*,*) ' x(',i,')=',x(i)
            enddo
         endif
!-------------------------------------
!    sampling along coordinate i_corr
!-------------------------------------
         if(index_int(i_corr).eq.0) then 
 
                call linesearchbox_cont(n,x,f,d,alfa,alfa_d,z,fz,i_corr,&
                           alfa_max,i_corr_fall,iprint,bl,bu,nf)

         else
                alfa_d_old=alfa_d(i_corr)
                call linesearchbox_discr(n,eta,index_int,scale_int,x,f,d,alfa,alfa_d,z,fz,i_corr,&
                           alfa_max,i_corr_fall,iprint,bl,bu,nf,discr_change,flag_fail)

         endif

         if(dabs(alfa).ge.1.d-12) then
		    
			flag_fail(i_corr)=0
		               
            x(i_corr) = x(i_corr)+alfa*d(i_corr)
            f=fz
 	        fstop(i_corr)=f
			     
            num_fal=0
      
         else

			flag_fail(i_corr)=1

			if ((index_int(i_corr).eq.1).and.(alfa_d_old.gt.scale_int(i_corr)))  flag_fail(i_corr)=0

	        if(i_corr_fall.lt.2) then 

		      fstop(i_corr)=fz         

              num_fal=num_fal+1

	        endif

	     end if

         ni=ni+1
		 z(i_corr) = x(i_corr)

         if(i_corr.lt.n) then
            i_corr=i_corr+1
         else
		    if(.not.discr_change) then
				do i = 1,n
					if((index_int(i).eq.1).and.(alfa_d(i) > 1)) then
						discr_change = .true.
						exit
					endif
				enddo
				if(.not.discr_change) then
					eta = eta/2.d0
				endif
			endif
            i_corr=1
	   	    discr_change = .false. 
         end if 

         call stop(n,index_int,scale_int,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max,flag_fail)

         if (istop.ge.1) exit


      enddo
      return
    


      end
        


!     #######################################################

      subroutine stop(n,index_int,scale_int, alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max, flag_fail)
      implicit none
      
      integer :: n,istop,i,nf,ni,nf_max
	  integer :: index_int(n), flag_fail(n)

      real*8 :: alfa_d(n),alfa_max,fstop(n+1),ffstop,ffm,f,alfa_stop
	  real*8 :: scale_int(n) 

	  logical :: test

      istop=0

      alfa_max=0.0d0


      do i=1,n				 
	    if (index_int(i).eq.0) then 
          if(alfa_d(i).gt.alfa_max) then
            alfa_max=alfa_d(i)
          end if
		end if
      end do
     

      if(ni.ge.(n+1)) then
        ffm=f
        do i=1,n
          ffm=ffm+fstop(i)
        enddo
        ffm=ffm/dfloat((n+1))

        ffstop=(f-ffm)*(f-ffm)
        do i=1,n
           ffstop=ffstop+(fstop(i)-ffm)*(fstop(i)-ffm)
        enddo
 
        ffstop=dsqrt(ffstop/dfloat(n+1))



	  endif


        
      if(alfa_max.le.alfa_stop) then
	    test=.true.
		do i=1,n
		
         if (index_int(i).eq.1) then 
        
		  if((alfa_d(i).ne.scale_int(i)).or.(flag_fail(i).ne.1)) then
		    test=.false.
	      end if
        
		 end if
	  

		end do
        if (test.eqv..true.) then
		   istop = 1
		end if
        
	  end if
      


      if(nf.gt.nf_max) then
        istop = 2
      end if

      if(ni.gt.nf_max) then
        istop = 3
      end if

      return

      end




!     *********************************************************
!     *         
!     *                 Continuous Linesearch
!     *
!     *********************************************************
           
 
      subroutine linesearchbox_cont(n,x,f,d,alfa,alfa_d,z,fz,i_corr,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,nf)
      
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int
      real*8 :: delta,delta1,fpar,fzdelta

	  
	  gamma=1.d-6      !-6

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0

!     index of current direction

      j=i_corr

	  if(iprint.ge.1) then
			write(*,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif


	  if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
			alfa=0.d0
			if(iprint.ge.1) then
				 write(*,*) '  alfa piccolo'
				 write(*,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
			endif
			return
	  endif
      
!     choice of the direction

	  do ielle=1,2

		 if(d(j).gt.0.d0) then

		     if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
   			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
			    alfa=bu(j)-x(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' point on the boundary. *'
				endif
			 endif

		  else

			 if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
				alfa=x(j)-bl(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *'
				endif
			 endif

		  endif

		  if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
  
			 d(j)=-d(j)
			 i_corr_fall=i_corr_fall+1
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta per alfa piccolo'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
			  endif
			  alfa=0.d0
			  cycle

		  endif

		  alfaex=alfa

		  z(j) = x(j)+alfa*d(j)
    
		 
	      call funct(n,z,fz)
		  

		  nf=nf+1

		  if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
		  endif
		  if(iprint.ge.2) then
			  do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
			  enddo
		  endif

		  fpar= f-gamma*alfa*alfa

!         test on the direction

		  if(fz.lt.fpar) then

!         expansion step

			 do

		   	   if((ifront.eq.1)) then

			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif
				     alfa_d(j)=delta*alfa

				     return

				 end if

				 if(d(j).gt.0.d0) then
							
					 if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=bu(j)-x(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
						 endif
					 end if

				 else

					 if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=x(j)-bl(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
						 endif
					 end if

				 endif
						 
				 z(j) = x(j)+alfaex*d(j) 
				   
     
				
			     call funct(n,z,fzdelta)
							      
				
				 nf=nf+1

				 if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
				 endif
				 if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
					  enddo
				 endif

				 fpar= f-gamma*alfaex*alfaex

				 if(fzdelta.lt.fpar) then

					 fz=fzdelta
					 alfa=alfaex

				 else               
					 alfa_d(j)=delta*alfa
			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto fz =',fz,'   alfa =',alfa
			         endif
					 return
				 end if

		     enddo 

		  else   !opposite direction    

			 d(j)=-d(j)
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
			 endif

		  endif       ! test on the direction
			  
	  enddo       

	  if(i_corr_fall.eq.2) then
			 alfa_d(j)=alfa_d(j)
	  else
			 alfa_d(j)=delta*alfa_d(j)
	  end if

	  alfa=0.d0

	  if(iprint.ge.1) then
			write(*,*) ' failure along the direction'
	  endif

	  return      
	  
      end

!     *********************************************************
!     *         
!     *         Linesearch along the discrete variables
!     *
!     ********************************************************
           
 
      subroutine linesearchbox_discr(n,eta,index_int,scale_int,x,f,d,alfa,alfa_d,z,fz,i_corr,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,nf,discr_change,flag_fail)
      
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
	  integer :: index_int(n),flag_fail(n)
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n),scale_int(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int,eta
      real*8 :: delta,delta1,fpar,fzdelta
	  logical:: discr_change, test

	  
      gamma_int=1.d-0 

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0

!     index of the current direction

      j=i_corr


	  if(iprint.ge.1) then
			   write(*,*) 'variabile discreta  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif

	  test=.true.

      if(alfa_d(i_corr).eq.scale_int(i_corr)) then  
        
		  do i=1,n
		  
		    if((flag_fail(i).eq.0)) then

		      test=.false.
			  exit

	        end if

		  enddo

          if(test) then

		     alfa=0.d0

		     if(iprint.ge.1) then
			    write(*,*) ' direzione gia analizzata'
		     endif

             return
          endif
                  
	   end if
      
	   do ielle=1,2

		   if(d(j).gt.0.d0) then

				if( ( ( bu(j)-x(j)-alfa_d(j) ) ).lt.0.d0 ) then  
				   alfa=      bu(j)-x(j)
				   ifront=1
				   if (alfa.eq.0.d0) then 
				      d(j)=-d(j)
					  ifront=0
				      cycle
				   endif
                else
				   alfa=alfa_d(j)
				end if		

		   else

				if( ((x(j)-alfa_d(j)-bl(j))).lt.0.0d0 ) then
				   alfa=      x(j)-bl(j)
				   if(iprint .gt. 1) then
					   write(*,*) 'alfa =',alfa
				   endif
				   ifront=1
				   if (alfa.eq.0.d0) then
				      d(j)=-d(j)
					  ifront=0
				      cycle
				   endif
				 else
				   alfa=alfa_d(j)
				endif

		   endif

           alfaex=alfa

		   z(j) = x(j)+alfa*d(j)
    
		   
		   call funct(n,z,fz)
		   

		   nf=nf+1

		   if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
		   endif
		   if(iprint.ge.2) then
			   do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
			   enddo
		   endif

		   fpar= f-gamma_int*eta

!          test on the direction

		   if(fz.lt.fpar) then
			  
			  discr_change = .true.

!             expansion step

			  do 
                  if(ifront.eq.1) then 

                     if(iprint.ge.1) then
				              write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif

				     return

				  endif

				  if(d(j).gt.0.d0) then
							
					 if((bu(j)-x(j)-2.0d0*alfa ).lt.(0.0d0)) then

					    alfaex=      bu(j)-x(j)
				        ifront=1

				        if (alfaex.le.alfa) then
						 
						   alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))

			               if(iprint.ge.1) then
				              write(*,*) ' accetta punto quasi frontiera fz =',fz,'   alfa =',alfa
			               endif

						   return

					    endif
						 
					  else

					     alfaex=alfa*2.0d0						
					
					  end if

				   else

					  if(( x(j)-2.0d0*alfa-bl(j) ).lt.(0.0d0)) then

					      alfaex=      x(j)-bl(j)
						  ifront=1

						  if (alfaex.le.alfa) then
						 
						     alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))

			                 if(iprint.ge.1) then
				               write(*,*) ' accetta punto quasi frontiera fz =',fz,'   alfa =',alfa
			                 endif

						     return

						  endif
						 
					  else

					      alfaex=alfa*2.0d0						
					
					  end if

				   endif
						 
				   z(j) = x(j)+alfaex*d(j) 
				   
     
				  
				   call funct(n,z,fzdelta)
				   			      
				
				   nf=nf+1

				   if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
				   endif
				   if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
					  enddo
				   endif

				   fpar= f-gamma_int*eta

				   if(fzdelta.lt.fpar) then

					  fz=fzdelta
					  alfa=alfaex

				   else               
					   alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))
			           if(iprint.ge.1) then
				         write(*,*) ' accetta punto  fz =',fz,'   alfa =',alfa
			           endif

					  return
				   end if

				enddo

			 else 

				d(j)=-d(j)
				ifront=0

				if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				endif

			 endif
			  
		  enddo

		  alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))

		  alfa=0.d0
		  if(iprint.ge.1) then
			  write(*,*) ' fallimento direzione'
		  endif

		  return
	  	     
      end


