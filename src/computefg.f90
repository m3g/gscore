!
! Subroutine that computes the complete function value
!

subroutine computef(n,x,f)

  use functionpars
  implicit none
  integer :: i, n
  double precision :: x(4), f, oneovers

  f = 0.d0
  do i = 1, nbins
    oneovers = 1.d0 / s(i) - 1.d0
    f = f + (probs(i) - x(1)* ( 1.d0/(1.d0+exp(-x(2)*(s(i)-x(3)))) ) * oneovers * exp(-x(4)*oneovers ))**2
  end do

end subroutine computef

!
! Subroutine that computes the complete gradient
!

subroutine computeg(n,x,g)

  use functionpars
  implicit none
  integer :: n, i
  double precision :: x(4), g(4), gall
  double precision :: oneovers

  do i = 1, n
    g(i) = 0.d0
  end do
  
  do i = 1, nbins

    oneovers = 1.d0 / s(i) - 1.d0

    gall = 2.d0*(probs(i) - x(1)* ( 1.d0/(1.d0+exp(-x(2)*(s(i)-x(3)))) ) * oneovers * exp(-x(4)*oneovers ))

    g(1) = g(1) - gall * ( 1.d0/(1.d0+exp(-x(2)*(s(i)-x(3)))) ) * oneovers * exp(-x(4)*oneovers )

    g(2) = g(2) - gall * x(1)*( 1.d0/(1.d0+exp(-x(2)*(s(i)-x(3))))**2 )*exp(-x(2)*(s(i)-x(3) ))*(s(i)-x(3))*&
                         oneovers * exp(-x(4)*oneovers )

    g(3) = g(3) - gall * x(1)*( -1.d0/(1.d0+exp(-x(2)*(s(i)-x(3))))**2 )*exp(-x(2)*(s(i)-x(3) ))*(x(2))*&
                        oneovers * exp(-x(4)*oneovers )

    g(4) = g(4) - gall * x(1)*( 1.d0/(1.d0+exp(-x(2)*(s(i)-x(3)))) ) * &
                         oneovers * exp(-x(4)*oneovers )*(-1.d0*oneovers)

  end do

end subroutine computeg

