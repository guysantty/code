!
module mod_data
!
contains
!-----------------------------------------------------------------------
!
subroutine write_data (u_field, NSTEP, flout)
  !
  implicit none
  !
  ! external variables
  integer, parameter :: dp=kind(1.0d0)
  integer :: NX, NSTEP, i, j
  real(dp), dimension(:,:) :: u_field
  character *300 :: flout
  !
  NX = size(u_field, 1)  

  !
  open(unit = 1, file = flout, status = 'unknown', &
       form = 'formatted', access = 'sequential', action = 'write')
  !
     do j =1000,1001
       do i = 1,NX
            write (unit=1, fmt='(E30.20)') (u_field(i,j))
       end do
     end do
  !
  close(unit=1)
  !  
  
 
 
 ! write 2D matrix in text file for plotting in any programming language (MATLAB, Pytton, etc).
 !open (unit = 16, file = flout, status = 'unknown', form = 'formatted',&
   !      access='sequential', action='write')

  !write (unit = 16, fmt = '(1600E30.20/)') ((u_field(j,i), i = 1,1600), j = 1,201)
 !close(unit = 16) 
 !

 end subroutine write_data
 !
 

end module mod_data
