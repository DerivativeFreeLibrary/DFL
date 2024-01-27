SUBROUTINE INIZPAR(N,X,A,B)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
	DIMENSION X(N),A(N),B(N)
	DO I=1,N
	A(I)=0.d0
	B(I)=10.d0
	X(I)=(A(I)+B(I))/2.D0
	ENDDO
	RETURN
END SUBROUTINE INIZPAR

subroutine which_integer(n,A,B,is_integer,step,x)

!----------------------------------------------------
! set even variables to integer
! 
! integer variable can assume 51 integer values 
! between A and B
! step   = (B-A) / 50
! tildex = A + x * step
! set starting value of the variables 
! integer variables at A+B/2
!----------------------------------------------------

	implicit none
	integer, intent(IN)	:: n
	real*8,  intent(INOUT)	:: A(n), B(n)
	logical, intent(OUT)	:: is_integer(n)
	real*8,  intent(OUT)	:: step(n)
	real*8,  intent(INOUT)	:: x(n)
	integer			:: i
	
	is_integer = .false.
	step       = 0.d0
    
	do i = 2,n,2
		is_integer(i) = .true.
		if((A(i) > -1.d+10).and.(B(i)) < 1.d+10) then
			step(i)       = ( B(i) - A(i) ) / 20.0d0
			x(i)          = ( A(i) + B(i) ) /  2.0d0
		
		elseif(A(i) > -1.d+10) then
			step(i)       = 1.d0
			x(i)          = A(i)
		
		elseif(B(i) <  1.d+10) then
			step(i)       = 1.d0
			x(i)          = B(i)
		
		else
			step(i)       = 1.d0
			x(i)          = 0.d0
		endif
	enddo
	
	return

end subroutine which_integer


SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION :: A(10,4), C(10), FA
      INTEGER          :: I, J

      DO I=1,4
         A(1,I)=4.D0
         A(2,I)=1.D0
         A(3,I)=8.D0
         A(4,I)=6.D0
      END DO
      DO I=1,2
         A(5,2*(I-1)+1)=3.D0
         A(5,2*I)=7.D0
         A(6,2*(I-1)+1)=2.D0
         A(6,2*I)=9.D0
         A(7,I)=5.D0
         A(7,I+2)=3.D0
         A(8,2*(I-1)+1)=8.D0
         A(8,2*I)=1.D0
         A(9,2*(I-1)+1)=6.D0
         A(9,2*I)=2.D0
         A(10,2*(I-1)+1)=7.D0
         A(10,2*I)=3.6D0
      END DO

      C(1)=0.1D0
      C(2)=0.2D0
      C(3)=0.2D0
      C(4)=0.4D0
      C(5)=0.4D0
      C(6)=0.6D0
      C(7)=0.3D0
      C(8)=0.7D0
      C(9)=0.5D0
      C(10)=0.5D0

      F  = 0.0D0
      FA = 0.0D0

      DO I=1,10
         DO J=1,4
            FA = FA +(X(J)-A(I,J))**2
         END DO
         IF ((FA+C(I)).EQ.0.D0) THEN
            F=1.D25
            RETURN
         ENDIF
         F  = F -1.0D0/(FA+C(I))
         FA = 0.0D0
      END DO

      RETURN

END SUBROUTINE FUNCT

