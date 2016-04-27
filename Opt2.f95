
module mod_FW
!
contains
!------------------------------------------------------------------------
!
subroutine FW_modelling_heterogeneous (u_field,rho,mu,NSTEP,NX,dx2,DELTAT,ISOURCE, &
 CONV2, OPT2, CONV4, OPT4, HET, HOM, BACK, ZREC, adjoint_source)
  !
  implicit none
  !
  ! external variables
  integer, parameter :: dp=kind(1.0d0)
  real(dp) :: DELTAT
  integer :: NSTEP, ISOURCE, NX
  real(dp), dimension(:) :: rho, mu
  real(dp), dimension(:,:) :: u_field, adjoint_source
 
  !  
  
  !internal variables
  integer :: i, it
  real(dp), parameter ::pi=3.141592653589793238462643d0
  real(dp) ::umaxv, to, fo, a, ampl, dx2, dt2, t, source_term, &
  & value_du_dxx
  real(dp), dimension(1:NX+1) :: unm1, unm2, u, g
  real(dp), dimension(3,5,NX-1) :: A_con4, A_opt4, K_con4, K_opt4, dA4, dK4, &
  & UI4
  real(dp), dimension(3,3,NX-1) :: A_con2, A_opt2, K_con2, K_opt2, &
  & dA2, dK2, UI2
  logical :: CONV2, OPT2, CONV4, OPT4, HET, HOM, BACK
  integer, dimension(:) :: ZREC
  real(dp) :: source
  
  !


  ! 
  dt2= DELTAT**2
  !
  
  ! operators as well as their corresponding residuals
  !
  ! SECOND ORDER OPERATORS
  !
  do i = 2, NX-1
   
   ! HETEROGENEOUS CASE
   if (HET.eqv.(.true.)) then
    if ((CONV2 .eqv. (.true.)).or.((OPT2).eqv.(.true.))) then

     A_con2(:,:,i)=(reshape((/real(dp)::0,0,0 ,rho(i), -2_dp*rho(i), rho(i), 0,0,0/), (/3,3/)))*(1._dp/dt2)

     k_con2(:,:,i)=(reshape((/real(dp)::0,(mu(i-1)+mu(i)),0,0,-(mu(i-1)+2_dp*mu(i)+mu(i+1)), 0, 0,(mu(i)+mu(i+1)),0/), &
                   (/3,3/)))*(1._dp/(2_dp*dx2))
    end if
     
   
    if ((OPT2).eqv.(.true.)) then

     A_opt2(:,:,i) = (reshape((/real(dp)::rho(i), -2*rho(i), rho(i),10*rho(i),-20*rho(i),10*rho(i),rho(i),&
                  -2*rho(i), rho(i)/), (/3,3/)))*(1._dp/(12_dp*dt2))

     K_opt2(:,:,i) = (reshape((/real(dp)::(mu(i-1)+mu(i)),10_dp*(mu(i-1)+mu(i)), (mu(i-1)+mu(i)),-(mu(i-1)+2_dp*mu(i)+mu(i+1)),&
                    -10_dp*(mu(i-1)+2_dp*mu(i)+mu(i+1)),-(mu(i-1)+2_dp*mu(i)+mu(i+1)), (mu(i)+mu(i+1)),&
                     10_dp*(mu(i)+mu(i+1)),(mu(i)+mu(i+1))/), (/3,3/)))*(1._dp/(24_dp*dx2))

     dA2(:,:,i) = A_opt2(:,:,i)-A_con2(:,:,i)

     dK2(:,:,i) = K_opt2(:,:,i)-K_con2(:,:,i)
    end if
   end if
  !
 
  
  ! HOMOGENEOUS CASE

   if (HOM.eqv.(.true.)) then
    if ((CONV2 .eqv. (.true.)).or.((OPT2).eqv.(.true.))) then
     A_con2(:,:,i)=(reshape((/real(dp)::0,0,0 ,rho(i), -2_dp*rho(i), rho(i), 0,0,0/), (/3,3/)))*(1._dp/dt2)

     k_con2(:,:,i)=(reshape((/real(dp)::0,mu(i),0,0,-2_dp*mu(i), 0, 0, mu(i),0/), &
                   (/3,3/)))*(1._dp/dx2)
    end if
   
    if ((OPT2).eqv.(.true.)) then

     A_opt2(:,:,i) = (reshape((/real(dp)::rho(i), -2_dp*rho(i), rho(i),10_dp*rho(i),-20_dp*rho(i),10_dp*rho(i),rho(i),&
                    -2_dp*rho(i), rho(i)/), (/3,3/)))*(1_dp/(12_dp*dt2))

     K_opt2(:,:,i) = (reshape((/real(dp)::(mu(i)),10_dp*mu(i), mu(i),-2_dp*mu(i),&
                    -20_dp*mu(i),-2_dp*mu(i), mu(i), 10_dp*mu(i),mu(i)/), (/3,3/)))*(1_dp/(24_dp*dx2))

     dA2(:,:,i) = A_opt2(:,:,i) - A_con2(:,:,i)

     dK2(:,:,i) = K_opt2(:,:,i) - K_con2(:,:,i)
    end if
   end if
   
  end do 
   !


   ! FOURTH ORDER ACCURACY
   !
  do i = 3, NX-2
       
   !HETEROGENEOUS CASE
   if (HET.eqv.(.true.))then
    if ((CONV4 .eqv. (.true.)).or.((OPT4).eqv.(.true.))) then
  
     A_con4(:,:,i)=(reshape((/real(dp):: 0,0,0 ,0,0,0, rho(i), -2._dp*rho(i), rho(i), 0,0,0, 0,0,0/), &
                    (/3,5/)))*(1._dp/dt2)

     k_con4(:,:,i)=(reshape((/real(dp):: 0, -(mu(i-2)+mu(i)), 0, 0, 16._dp*(mu(i)+mu(i-1)), 0, 0, &
                                   -16*(mu(i-1)+2*mu(i)+mu(i+1))+(mu(i-2)+2*mu(i)+mu(i+2)),&
                     0, 0, 16*(mu(i)+mu(i+1)), 0, 0,-(mu(i)+mu(i+2)),0/), (/3,5/)))*(1._dp/(24_dp*dx2))
    end if


    if ((OPT4).eqv.(.true.)) then

     A_opt4(:,:,i)= (reshape((/real(dp):: -rho(i), 2._dp*rho(i), -rho(i), 4._dp*rho(i), -8._dp*rho(i), 4._dp*rho(i), 84._dp*rho(i),&
                  -168._dp*rho(i), 84._dp*rho(i), 4._dp*rho(i), -8._dp*rho(i), 4._dp*rho(i), -rho(i), 2._dp*rho(i), -rho(i)/), &
                   (/3,5/)))*(1._dp/(90._dp*dt2))

     K_opt4(:,:,i)=(reshape((/real(dp):: -(mu(i-2)+mu(i)), -10._dp*(mu(i-2)+mu(i)), -(mu(i-2)+mu(i)), 16._dp*(mu(i)+mu(i-1)), &
              16._dp*10._dp*(mu(i)+mu(i-1)), 16._dp*(mu(i)+mu(i-1)), -16._dp*(mu(i-1)+2*mu(i)+mu(i+1))+(mu(i-2)+2*mu(i)+mu(i+2)), &
         -16*10*(mu(i-1)+2*mu(i)+mu(i+1))+(mu(i-2)+2*mu(i)+mu(i+2)), -16*(mu(i-1)+2*mu(i)+mu(i+1))+10*(mu(i-2)+2*mu(i)+mu(i+2)), &
         16*(mu(i)+mu(i+1)), 16*10*(mu(i)+mu(i+1)), 16*(mu(i)+mu(i+1)), -(mu(i)+mu(i+2)), -10*(mu(i)+mu(i+2)), -(mu(i)+mu(i+2))/),&
              (/3,5/)))*(1/(288*dx2))

     dA4(:,:,i)=A_opt4(:,:,i)-A_con4(:,:,i)

     dK4(:,:,i)=K_opt4(:,:,i)-K_con4(:,:,i)

    end if
   end if
  

   !HOMOGENEOUS CASE

   if (HOM.eqv.(.true.)) then
    if ((CONV4 .eqv. (.true.)).or.((OPT4).eqv.(.true.))) then
  
     A_con4(:,:,i)=(reshape((/real(dp):: 0,0,0 ,0,0,0, rho(i), -2._dp*rho(i), rho(i), 0,0,0, 0,0,0/), &
                    (/3,5/)))*(1._dp/dt2)

     k_con4(:,:,i)=(reshape((/real(dp):: 0, -mu(i), 0, 0, 16._dp*mu(i), 0, 0, &
                    -30*mu(i),0, 0, 16*mu(i), 0, 0,-mu(i),0/), (/3,5/)))*(1._dp/(12._dp*dx2))
        
    end if
    
  
    if ((OPT4).eqv.(.true.)) then
     A_opt4(:,:,i)= (reshape((/real(dp):: -rho(i), 2._dp*rho(i), -rho(i), 4._dp*rho(i), -8._dp*rho(i), 4._dp*rho(i), 84._dp*rho(i),&
                  -168._dp*rho(i), 84._dp*rho(i), 4._dp*rho(i), -8._dp*rho(i), 4._dp*rho(i), -rho(i), 2._dp*rho(i), -rho(i)/), &
                   (/3,5/)))*(1/(90._dp*dt2))

     K_opt4(:,:,i)=(reshape((/real(dp):: -mu(i), -10._dp*mu(i), -mu(i), 16._dp*mu(i),160._dp*mu(i), 16._dp*mu(i), -30._dp*mu(i),&
         -300*mu(i), -30*mu(i),16*mu(i), 160*mu(i), 16*mu(i), -mu(i), -10*mu(i), -mu(i)/), (/3,5/)))*(1/(144*dx2))
 
     dA4(:,:,i) = A_opt4(:,:,i) - A_con4(:,:,i)

     dK4(:,:,i) = K_opt4(:,:,i) - K_con4(:,:,i)

    end if
   end if
  end do 
  !

  !

  !parameters for the source
  fo=100._dp    ! frequency
  to=1.20_dp/fo ! time of source activation
  ampl=1.0d9
  a=pi*pi*fo*fo ! amplitude of the signal
  umaxv=(ampl*exp(-a*to**2))/2._dp*DELTAT

  !ADJOINT SOURCE
  if (BACK .eqv. (.true.)) then
   u_field = 0
  end if
      


!--------------------------------------------------------------------------------------
!                           PREDICTOR-CORRECTOR SCHEMES
!--------------------------------------------------------------------------------------

  do it = 1, NSTEP
     u= 0.
     g = 0.
     ! ricker source time function (second derivative of a Gaussian)
       t = (it-1)*DELTAT
       source_term=ampl*exp(-a*(t-to)**2)*DELTAT**2/rho(ISOURCE)

     
   ! SECOND ORDER ACCURACY ------------------------------------------------------------------
    if ((CONV2 .eqv. (.true.)).or.((OPT2).eqv.(.true.))) then
      
   ! PREDICTOR SCHEME               
      do i = 2, NX-1
      
        if (BACK .eqv. (.true.)) then
          if (any (i .eq. ZREC)) then
            source = adjoint_source(i, it) 
          else 
            source = 0
          end if
        end if   
     
        if (BACK.eqv.(.false.)) then
          source = 0
        end if
          
         value_du_dxx = unm1(i-1)*K_con2(2,1,i)+unm1(i)*K_con2(2,2,i)+unm1(i+1)*K_con2(2,3,i)
      u(i) = (value_du_dxx-unm1(i)*A_con2(2,2,i)-unm2(i)*A_con2(3,2,i))/A_con2(1,2,i)+source*dt2/rho(i)
   
      end do

       
     
    ! SOURCE TERM 
     if ((BACK.eqv.(.false.)).and.(CONV2.eqv.(.true.)))  then
       u(ISOURCE) = source_term + u(ISOURCE)
      end if


       ! UPDATE
      if ((CONV2).eqv.(.true.)) then
      unm2 = unm1 
      unm1 = u  
      u_field(:,it)=u
      end if
   
    end if
    ! 


   ! CORRECTOR SCHEME 
    if ((OPT2).eqv.(.true.)) then
      !UI2(:,:,:)= 0._dp

       ! SECONDARY SOURCE
      do i = 2, NX-1
         UI2(1:3,1:3,i) = reshape((/real(kind=dp):: u(i-1),unm1(i-1),unm2(i-1),u(i),unm1(i),&
            unm2(i),u(i+1),unm1(i+1),unm2(i+1)/),(/3,3/))
       
        g(i) = -sum(sum((dA2(:,:,i)-dK2(:,:,i))*UI2(:,:,i),dim=1),dim=1)/A_con2(1,2,i)
      end do

      ! CORRECTED DISPLACEMENT FIELD
      do i = 2, NX-1
            u(i) = u(i) + (g(i)/A_con2(1,2,i))
      end do


     ! SOURCE TERM 
      if ((BACK.eqv.(.false.)).and.(OPT2.eqv.(.true.)))  then
        u(ISOURCE) = source_term + u(ISOURCE)
      end if

      
       !UPDATE
     if ((OPT2).eqv.(.true.)) then
     unm2 = unm1
     unm1 = u
     u_field(:, it) = u 
     end if 
     
    end if
  ! END OF CORRECTOR SCHEME

  ! END OF SECOND ORDER ACCURACY--------------------------------------------------------------

   

  ! FOURTH ORDER ACCURACY----------------------------------------------------------------------
   
    if ((CONV4 .eqv. (.true.)).or.((OPT4).eqv.(.true.))) then

      if (BACK .eqv. (.true.)) then
        if (any (i .eq. (ZREC))) then
          source = adjoint_source(i, it) 
        else 
         source = 0
        end if
      end if   
     
      if (BACK.eqv.(.false.)) then
        source = 0
      end if
    
   ! PREDICTOR SCHEME
      do i = 3, NX-2
         value_du_dxx = unm1(i-2)*K_con4(2,1,i) + unm1(i-1)*K_con4(2,2,i) + unm1(i)*K_con4(2,3,i) &
                       +unm1(i+1)*K_con4(2,4,i)+unm1(i+2)*K_con4(2,5,i)
       u(i) = (value_du_dxx -unm1(i)*A_con4(2,3,i)-unm2(i)*A_con4(3,3,i))/A_con4(1,3,i)+(source*dt2/rho(i))
      end do

    ! SOURCE TERM
      if ((BACK.eqv.(.false.)).and.(CONV4.eqv.(.true.))) then
         u(ISOURCE) = u(ISOURCE) + source_term
      end if

       ! update
    if ((CONV4).eqv.(.true.)) then
     unm2 = unm1 
     unm1 = u  
     u_field(:,it)=u
    end if

    
    end if
  !END OF PREDICTOR SCHEME
   

  !CORRECTOR SCHEME
    if ((OPT4).eqv.(.true.)) then
     UI4(:,:,:) = 0._dp

      do i = 3, NX-2
        UI4(1:3,1:5,i) = reshape((/real(kind=dp):: u(i-2), unm1(i-2), unm2(i-2), u(i-1), unm1(i-1), unm2(i-1), &
          u(i), unm1(i), unm2(i), u(i+1), unm1(i+1), unm2(i+1), u(i+2), unm1(i+2), unm2(i+2)/),(/3,5/))
       
        g(i) = -sum(sum((dA4(:,:,i)-dK4(:,:,i))*UI4(:,:,i),dim=1),dim=1)*dt2/rho(i)
      end do
   !

   !CORRECTED DISPLACEMENT FIELD
      do i = 3, NX-2
          u(i) = u(i) + g(i)/A_con4(1,3,i)
      end do

   ! SOURCE TERM
    if ((BACK.eqv. (.false.) ).and. (OPT4.eqv.(.true.))) then
       u(ISOURCE) = u(ISOURCE) + source_term
      end if

       !UPDATE
    if ((OPT4).eqv.(.true.)) then
     unm2 = unm1
     unm1 = u
     u_field(:, it) = u 
    end if 
      
    end if
   !
  
  end do

 end subroutine FW_modelling_heterogeneous 
 !
 
end module mod_FW
!--------------------------------------------------------------------------------------------------------------------------------------------


