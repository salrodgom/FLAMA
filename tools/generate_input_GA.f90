PROGRAM generate_GA_input 
!
! This program generates an input file to use the Genetic Algorithm.
!
IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: n_param = 26 
INTEGER, PARAMETER :: n_refit = 1 !800 
INTEGER, PARAMETER :: seed = 2123
REAL, PARAMETER :: max_change = 0.1 ! 0.1 is is a change of 10%  
! Variables
INTEGER :: i,j    
real :: phenotype(n_param),dum 

OPEN(110,file="phenotype.dat")
OPEN(111,file="generated_input")

write(111,'(a)') 'physically_constrained'
write(111,'(a,i3)') 'n_par', n_param
DO i=1,n_param 
 read(110,*) phenotype(i)
 if (i<10) then
  write(111,'(a,i1,a)') 'xxxxxxxxx',i,'  Buck'
 else
  write(111,'(a,i2,a)') 'xxxxxxxx',i,'  Buck'
 end if
END DO
CLOSE(110)
write(111,'(a)') 'ffit?  .true.'
write(111,'(a)') 'refit? .true.'
write(111,*) n_refit

call srand(seed)
DO i=1,n_refit
 DO j=1,n_param
  if (i==1) then
   dum=phenotype(j)
  else
   dum=phenotype(j)*(rand()*2.0*max_change+1.0-max_change) 
    if (dum<0.0) then
     dum=0.0
    end if
  end if
   write(111,'(b32.32,a,e20.15)') dum, ' # ', dum
 END DO
END DO

end program
