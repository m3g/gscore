!
! Subroutine that performs finite difference and analytical gradient
! comparision. Used only for testing purpouses
!
! L. Martinez, Dec, 4, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine test_grad(n,x,g,computef,computeg)

  implicit none
  integer ::  n, i, iworst
  double precision :: x(n), g(n), fx, step, gcomp, gbest, eworst, error, steperror, stepbest
  real :: time0, tarray(2), etime
  external :: computef, computeg

  write(*,*)
  write(*,*) ' Comparing analytical and finite-difference gradients... may take a while. '
  write(*,*)
  write(*,*) ' Five first coordinates variables: '
  do i = 1, min(5,n)
    write(*,"( 3(tr2,f12.5) )" ) x((i-1)*3+1), x((i-1)*3+2), x((i-1)*3+3)
  end do
  write(*,*) 
  write(*,*) ' Computing gradient ... ' 
  
  call computef(n,x,fx) 
  write(*,*) ' Function value at test point: ', fx
  open(98, file = 'chkgrad.log',status='unknown') 
  write(98, *)'Function Value = ', fx 
  call computeg(n,x,g) 
  
  write(98,"( t2,'Component',t16,'Analytical',t33,'Discrete',t51,'Error',t62,'Best step')")
  write(*,"( t2,'Component',t16,'Analytical',t33,'Discrete',t51,'Error',t62,'Best step')")
  time0 = etime(tarray)
  iworst = 0
  eworst = 0.d0
  do i = 1, n
    if(etime(tarray)-time0.gt.10.) then
      time0 = etime(tarray)
      write(*,*) ' Computing the ',i,'th of ',n,' components. '
    end if
    error = 1.d20
    step = 1.d-2
    do while(error.gt.1.d-6.and.step.ge.1.d-20)
      call discret(n,i,x,gcomp,step,computef)
      if(dmin1(abs(g(i)),abs(gcomp)).gt.1.d-10) then
        steperror = abs( ( gcomp - g(i) ) / g(i) )
      else
        steperror = abs( gcomp - g(i) )
      end if
      if( steperror .lt. error ) then
        error = steperror
        gbest = gcomp
        stepbest = step
      end if
      step = step / 10.d0
    end do
    write(98,"(i10,4(tr2,d13.6))") i, g(i), gbest, error, stepbest 
    flush(98)
    if ( i <= 9 ) then
      write(*,"(i10,4(tr2,d13.6))") i, g(i), gbest, error, stepbest 
    end if
    if(error.gt.eworst) then
      iworst = i
      eworst = error
      write(*,*) ' Greatest relative error up to now: ', eworst
    end if
  end do
  write(98,*) 'Maximum difference = ', iworst,' Error= ', eworst
  
  write(*,*) ' Done. '

end     

subroutine discret(n,icomp,x,gcomp,step,computef)
  implicit none
  integer :: icomp, n
  double precision :: save, step, x(*), fplus, fminus, gcomp
  external :: computef
  save = x(icomp) 
  x(icomp) = save + step 
  call computef(n,x,fplus) 
  x(icomp) = save - step 
  call computef(n,x,fminus) 
  gcomp = (fplus - fminus) / (2.d0 * step) 
  x(icomp) = save 
  return      
end

