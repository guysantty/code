program main
!
!                       1d FWI code based on :
!                 - optimally accurate operators O(2,2)
!                 - optimally accurate operators O(2,4)
!                 - conventional operators O(2,2)
!                 - conventional operators O(2,4)
!                 - staggered grid O(2,2)
!                           - Guy Mabyalaht                    
!                              12-March-2016
!
!---------------------------------------------------------------------------
!  	                  MODELLING SETTINGS
!---------------------------------------------------------------------------
  use mod_FW
  use mod_data
  use mod_pml
  !
  implicit none
  !
  !
  ! variables to store the options used in the modelling 
  logical :: CONV2, OPT2, CONV4, OPT4, HOM, HET, BACK
  integer :: choice

  ! file
  character *300 :: fichier
 
  ! variable that stores double precision declaration
  integer, parameter :: dp=kind(1.0d0)
  
  ! displacement field data
  real(dp), dimension(:,:), allocatable :: u_field, dobs, u_field1, u_field2, residuals, &
  & adjoint_source
  
  ! pml data
  real(dp), dimension(:), allocatable :: sig, cpe 
  
  ! reference and starting models
  real(dp), dimension(:), allocatable :: rho, mu, u, cp, cs, rho1, mu1, cp1, cs1
  
  ! other data
  integer ::  i, NSTEP, ISOURCE, npml, k=1, NX,NX1, mg, NREC
  real(dp):: time=0.08d0, DELTAT, XMAX=100.0d0, DELTAX, XMIN=0,&
  & xsource, mindist, courant_number1, vm, a=2.5, b=.5, vmin, &
  & courant_number2, dx2, dist
  integer, dimension(:), allocatable :: ZREC

  !--------------------------------------------------------------------------------

  ! Displays the Modelling operators availabe in the FWI code
  write(*,*) 'Welcome to the 1D FWI code ' 
  write(*,*) 'Here are the options for the modelling operators to be used'
  write(*,*) '1 for the Opt(2,2)'
  write(*,*) '2 for the Opt(2,2)'
  write(*,*) '3 for the Conv(2,2)'
  write(*,*) '4 for the Conv(2,4)'
  write(*,*) '5 for the staggered grid'
  write (*,*) 'the choice must be defined in the input_data.don file'
  !
   
  ! option decided by the user
  read(*,*) choice
  !
  ! input data
  read(*,*) NX     ! model size
  read(*,*) NX1    ! model size
  read(*,*) vmin   ! vmin
  read(*,*) DELTAT ! time step
  read(*,*)  mg
  read (*,*) NREC  !number of receivers
  !
  !
  allocate(rho(NX+1), mu(NX+1), u(NX+1), cp(NX+1), cs(NX+1), &
   rho1(NX+1), mu1(NX+1), cp1(NX+1), cs1(NX+1))
  !
  
  ! selects the operator to be used
  select case (choice)
  ! Opt(2,2)
  case (1)
    write(*,*) 'opt(2,2) has been selected'
    CONV2 = .false.
    CONV4 = .false.
    OPT2  = .true.
    OPT4  = .false.

  ! Opt(2,4)
  case(2) 
   write(*,*) 'opt(2,4) has been selected'
    CONV2 = .false.
    CONV4 = .false.
    OPT2  = .false.
    OPT4  = .true.
  
  ! Conv(2,2)
  case(3) 
   write(*,*) 'conv(2,2) has been selected'
    CONV2 = .true.
    CONV4 = .false.
    OPT2  = .false.
    OPT4  = .false.

  ! Conv(2,4)
  case(4)  
   write(*,*) 'conv(2,4) has been selected'
    CONV2 = .false.
    CONV4 = .true.
    OPT2  = .false.
    OPT4  = .false.
  
  ! No choice has been done
  case default
    write(*,*) 'Choice out of range, the code will stop'
    stop
  end select
 !
 

 !----------------------------------------------------------------------------
 
 ! time step in seconds
 DELTAT=DELTAT/mg

 ! total number of time steps
   NSTEP= nint(time/DELTAT)

 ! X SPACING
 DELTAX = (XMAX -XMIN)/NX
  dx2=DELTAX**2
  
 allocate(dobs(NX+1, NSTEP), ZREC(NREC) ) 

 ! define i for the source
   xsource = 0.1_dp*XMAX    ! arbitrary source position in [m]
   mindist = XMAX
   do i = 2, NX
   dist = xsource - DELTAX*i
      if (abs(dist) < mindist) then
         mindist = dist
         ISOURCE = i
      end if
   end do

  ! receiver location
   ZREC=(/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
 
  ! reference velocity models (cp, cs)
   cp=(/((vmin+k*a), k=0,NX)/)
   cs=cp/1.732_dp
  
 ! reference rho (density)
   rho=(/((2400_dp+k*b), k=0,NX)/)
 
 ! reference mu (stiffness/rigidity matrix) 
   do i=1,NX+1
       mu(i)=rho(i)*cs(i)*cs(i)
   end do

 ! STARTING MODEL
  cp1 = 2500
  cs1 = cp1/1.732_dp
  rho1 = 1500
   do i=1,NX+1
       mu1(i)=rho1(i)*cs1(i)*cs1(i)
   end do

 
 ! stability 
    courant_number1=2500_dp*DELTAT*DELTAX
   write(*,*)'courant number1 = ', courant_number1
    courant_number2 = 3000_dp*DELTAT*DELTAX
   write(*,*)'courant number2 = ', courant_number2
   if((courant_number1>1_dp) .or.( courant_number2>1_dp)) then
      write(*,*) ' time step is too large, simulation will be unstable'
   end if
 
 !
  

 ! predictor corrector data vectors
  allocate(u_field(NX+1, NSTEP), u_field1(NX+1, NSTEP), u_field2(NX+1, NSTEP),&
       residuals(NX+1, NSTEP), adjoint_source(sizeof(ZREC), NX+1)) ! save the field
  !
  allocate(sig(NX+1+2*npml))
  !
  npml=10
  !
  allocate(cpe(NX+1+2*npml))  
  !
  vm=sum(cp)/size(cp,1)
  !
  !call pml (DELTAX, vm, npml, sig, NX) 
  !
  !call extend (cp, npml, cpe) 
  write(*,*) DELTAT, DELTAX, isource
  !stop
 
 !---------------------- CALLS TO MODELLING TOOLS--------------------------------

  ! FORWARD MODELLING O(2,2) and Conv(2,2)
 if(((OPT2).eqv.(.true.)).or.((CONV2).eqv.(.true.))) then
   BACK=.false.

    ! REFERENCE MODEL
   HET=.true. ! hetereogeneous exact model
   HOM=.false. 
   call FW_modelling_heterogeneous (u_field,rho,mu,NSTEP,NX,dx2,DELTAT,ISOURCE, &
         CONV2, OPT2, CONV4, OPT4, HET, HOM, BACK, ZREC, adjoint_source)
      
     fichier = 'output_dobs_22.txt'
     call write_data(u_field, NSTEP, fichier) ! ASCII measured data
  
    u_field1 = u_field
    u_field=0

    ! STARTING MODEL
    HET=.false. !
    HOM=.true. 
   call FW_modelling_heterogeneous (u_field,rho1,mu1,NSTEP,NX,dx2,DELTAT,ISOURCE, &
         CONV2, OPT2, CONV4, OPT4, HET, HOM, BACK, ZREC, adjoint_source)
 
     fichier = 'output_calc_22.txt'
     call write_data(u_field, NSTEP, fichier) ! ASCII measured data
   
   u_field2 = u_field

   residuals = u_field1-u_field2
 
   fichier = 'residuals_22.txt'
     call write_data(residuals, NSTEP, fichier) ! ASCII measured data
 
   adjoint_source = residuals(ZREC, :)
   fichier = 'adjoint1_22.txt'
     call write_data(adjoint_source, NSTEP, fichier) ! ASCII measured data
  
 
   adjoint_source = adjoint_source(:, NSTEP:1:-1)
   fichier = 'adjoint2_22.txt'
     call write_data(adjoint_source, NSTEP, fichier) ! ASCII measured data

 
 end if
 !

  u_field=0
 ! FORWARD MODELLING O(2,4) and Conv(2,4)
  if(((OPT4).eqv.(.true.)).or.((CONV4).eqv.(.true.))) then

   BACK=.false.
    ! REFERENCE MODEL
    HET=.true. ! hetereogeneous exact model
    HOM=.false. 
    call FW_modelling_heterogeneous (u_field,rho,mu,NSTEP,NX,dx2,DELTAT,ISOURCE, &
         CONV2, OPT2, CONV4, OPT4, HOM, HET, BACK, ZREC, adjoint_source)

     fichier = 'output_dobs_24.txt'
     call write_data(u_field, NSTEP, fichier) ! ASCII measured data
     
    ! STARTING MODEL
      HET=.false. ! hetereogeneous exact model
      HOM=.true.
    call FW_modelling_heterogeneous (u_field,rho1,mu1,NSTEP,NX,dx2,DELTAT,ISOURCE, &
         CONV2, OPT2, CONV4, OPT4, HOM, HET, BACK, ZREC, adjoint_source)

     fichier = 'output_calc_24.txt'
     call write_data(u_field, NSTEP, fichier) ! ASCII measured data
  end if
  !
  

  ! FORWARD MODELLING WITH AN ADJOINT SOURCE O(2,2) and Conv(2,2)
    BACK=.true.
    
    if(((OPT2).eqv.(.true.)).or.((CONV2).eqv.(.true.)).and.(BACK.eqv.(.true.))) then
   
     HET=.false. 
     HOM=.true. 
     call FW_modelling_heterogeneous (u_field,rho1,mu1,NSTEP,NX,dx2,DELTAT,ISOURCE, &
         CONV2, OPT2, CONV4, OPT4, HOM, HET, BACK, ZREC, adjoint_source)

    fichier = 'output_back_22.txt'
     call write_data(u_field, NSTEP, fichier) ! ASCII measured data
    end if
  !


  ! FORWARD MODELLING WITH AN ADJOINT SOURCE O(2,4) and Conv(2,4)
    BACK=.false.
   if(((OPT4).eqv.(.true.)).or.((CONV4).eqv.(.true.)).and.BACK.eqv.(.true.)) then
    HET=.false. 
     HOM=.true. 
    call FW_modelling_heterogeneous (u_field,rho1,mu1,NSTEP,NX,dx2,DELTAT,ISOURCE, &
         CONV2, OPT2, CONV4, OPT4, HOM, HET, BACK, ZREC, adjoint_source)

    fichier = 'output_back_24.txt'
     call write_data(u_field, NSTEP, fichier) ! ASCII measured data 

   end if

!----------------END OF CALLS-------------------------------------------------------------
    
  ! Write data 
  !call write_data(u_field, dobs, NSTEP) ! ASCII 
  !call write_data_bin(u_field, dobs, NSTEP) ! BINARY
  
 
end program main


