!
! Subroutine that calls the cgnewton optimization code
!

subroutine callcgnewton(n,x,f,outputunit,evalf,evalg,optpars)

  implicit none
  integer :: optpars(10)
  integer :: n, maxkon, outputunit, maxnef, nopro, ipri, maxcg, ier
  double precision :: x(n), f, eps, frac, epsnop, g(n), xmax, &
                      xn(n), d(n), aux(3*n)
  external :: evalf, evalg

  if ( outputunit == 0 ) then
    ipri = 0
  else if( outputunit < 0 ) then
    ipri = 2
  else
    ipri = 1
  end if

  maxnef = optpars(1)
  maxcg = optpars(2)

  eps = 1.d-3
  xmax = 1000.d10
  maxkon = 1000000
  nopro = 1000
  frac = 0.d0
  epsnop = 1.d-8
  if ( ipri == 1 ) then
    write(outputunit,*) ' maxcg = ', maxcg
  else if ( ipri == 2 ) then
    write(*,*) ' maxcg = ', maxcg
  end if
  call cgnewton(n, x, f, g, maxcg, ier, eps, frac, ipri, xmax, &
                maxkon, maxnef, nopro, epsnop, d, aux, xn, &
                outputunit, evalf, evalg)    

end subroutine callcgnewton

!
! Interface for funcgn
!

subroutine funcgn(n, x, f, warn, evalf)
  implicit none
  integer :: n
  double precision :: x(n), f
  logical :: warn
  external :: evalf

  warn = .false.

  call evalf(n,x,f)

  return
end subroutine funcgn

!
! Interface for gracgn
!

subroutine gracgn(n, x, g, warn, evalg)

  implicit none
  integer :: n
  double precision :: x(n), g(n)
  logical :: warn
  external :: evalg

  warn = .false.

  call evalg(n,x,g)

return
end

!
! Subroutine that computes the product of the hessian by a vector
!

subroutine hesvec(n, x, xnor, xn, b, p, ap, evalg)

  implicit none
  integer :: n, i
  double precision :: x(n), b(n), p(n), ap(n), xn(n), h, pnor, xnor
  logical :: warn
  external :: evalg

  pnor = 0.d0
  do i = 1, n
    pnor = dmax1(pnor, dabs(p(i)))
  end do

  h = 1.d-6*dmax1(0.001, xnor)/dmax1(1.d0, pnor)

  do i = 1, n
    xn(i) = x(i) + h*p(i)
  end do

  call gracgn(n, xn, ap, warn, evalg)
  do i = 1, n
    ap(i) = (ap(i) + b(i))/h 
  end do     

return
end

subroutine cgnewton(n, x, f, g, maxcg, ier, eps, frac, ipri, xmax,&
                    maxkon, maxnef, nopro, epsnop, d, aux, xn,& 
                    outputunit, evalf, evalg)
 
  implicit none
  integer :: n, kon, i, j, maxcg, ipri, maxkon, maxnef, nopro, konop, iercg, kcg, ier, nef
  integer :: outputunit
  double precision :: x(n), g(n), d(n), aux(3*n), f, gnor, xmax, epsnop
  double precision :: xn(n), pesca, fn, t, delta, xnor, eps, dnor, frac
  double precision :: epscg
  logical :: warn, negra
  external :: evalf, evalg

  ! This subroutine tries to minimize the unconstrained function f(x)
  ! using the Inexact-Newton method that employs Conjugate-Gradients
  ! for solving the Newtonian Linear systems. 
  ! n is the number of variables
  ! x is the initial point and the best point obtained.
  ! maxcg is the maximal number of iterations allowed at each call of cg.
  ! eps is the convergence criterion.
  ! frac (between 0 and 1) (generally 0.1) is the fraction of the gradient
  ! that is required at each call to the CG-method.
  ! ipri should be negative if you dont want printing
  ! xmax is a maximal value for the norm of x 
  ! nopro is the number of iterations you tolerate with progress smaller than epsnop
  ! f, g, in output are the function values and gradient.
  ! aux and xn are auxiliar vectors.
  ! maxkon is the maximal number of iterations.
  ! maxnef is the maximal number of evaluations.
  ! The user needs to provide the subroutine funcgn, that computes the
  ! objective function and the subroutine gradcgn that computes the gradient.
  ! The subroutine hesvec, that should compute hessian-vector products is optional,
  ! since it can be replaced by the one provided together with the subroutine cgnewton.
  ! The subroutine cgnewton is, essentially, a simplification of Gencan, that only
  ! deals with unconstrained problems. 
  !
  ! J. M. Martinez, IMECC-UNICAMP
  ! http://www.ime.unicamp.br/~martinez
  !

  ier = 0
  nef = 1
  kon = 0
  konop = 0
  call funcgn(n, x, f, warn, evalf)
1 continue
  call gracgn(n, x, g, warn, evalg)
  konop = 0
  if(warn) then
    ier = 10
    if ( ipri == 1 ) then
      write(outputunit,*)' Function or gradient not possible to evaluate'
      write(outputunit,*)' at cgnewton. Return with ier = 10'
      write(outputunit,*)' Function or gradient not possible to evaluate'
      write(outputunit,*)' at cgnewton. Return with ier = 10'
    else if(ipri.eq.2) then
      write(*,*)' Function or gradient not possible to evaluate'
      write(*,*)' at cgnewton. Return with ier = 10'
      write(*,*)' Function or gradient not possible to evaluate'
      write(*,*)' at cgnewton. Return with ier = 10'
    end if
    return
  endif

  gnor = 0.d0
  do i = 1, n
    gnor = dmax1(gnor, dabs(g(i)))
  end do

  if(ipri.eq.1) then
    write(outputunit,*)
    write(outputunit,*)' cgnewton iteration:', kon
    write(outputunit,*)' x(1), x(2), x(n) =', x(1), x(2), x(n)
    write(outputunit,*)' f(x) = ', f
    write(outputunit,*)' g(1), g(2), g(n) =', g(1), g(2), g(n)
    write(outputunit,*)' sup norm of gradient:', gnor
    write(outputunit,*)' function evaluations:', nef
  else if (ipri == 2) then
    write(*,*)
    write(*,*)' cgnewton iteration:', kon
    write(*,*)' x(1), x(2), x(n) =', x(1), x(2), x(n)
    write(*,*)' f(x) = ', f
    write(*,*)' g(1), g(2), g(n) =', g(1), g(2), g(n)
    write(*,*)' sup norm of gradient:', gnor
    write(*,*)' function evaluations:', nef
  end if

  if(gnor.le.eps) then
    ier = 0
    if(ipri.eq.1) then
      write(outputunit,*)' sup-norm of gradient small enough', gnor
      write(outputunit,*)' cgnewton returns with ier = 0'
    else if ( ipri == 2 ) then
      write(*,*)' sup-norm of gradient small enough', gnor
      write(*,*)' cgnewton returns with ier = 0'
    endif
    return
  endif

  xnor = 0.d0
  do i = 1, n
    xnor = dmax1(xnor, dabs(x(i)))
  end do

  if(kon.eq.0) delta = xnor 

  if(xnor.gt.xmax) then
    ier = 1
    if(ipri.eq.1) then
      write(outputunit,*)' sup-norm of x very big (>', xmax
      write(outputunit,*)' cgnewton returns with ier = 1'
    else if ( ipri .eq. 2 ) then
      write(*, *)' sup-norm of x very big (>', xmax
      write(*, *)' cgnewton returns with ier = 1'
    end if
    return
  endif
 
  if(kon.ge.maxkon) then
    ier = 2 
    if(ipri.eq.1) then
      write(outputunit,*)' Excessive cgnewton iterations:', kon
      write(outputunit,*)' cgnewton returns with ier = 2'
    else if ( ipri .eq. 2 ) then
      write(*, *)' Excessive cgnewton iterations:', kon
      write(*, *)' cgnewton returns with ier = 2'
    endif
    return
  endif
 
  if(nef.ge.maxnef) then
    ier = 3 
    if(ipri.eq.1) then
      write(outputunit,*)' Excessive cgnewton evaluations:', nef
      write(outputunit,*)' cgnewton returns with ier = 3'
    else if ( ipri .eq. 2 ) then
      write(*, *)' Excessive cgnewton evaluations:', nef
      write(*, *)' cgnewton returns with ier = 3'
    endif
    return
  endif
  
  if(konop.ge.nopro) then
    ier = 4 
    if(ipri.eq.1) then
      write(outputunit,*)konop,' cgnewton iterations without progress'
      write(outputunit,*)' cgnewton returns with ier = 4'
    else if ( ipri .eq. 2 ) then
      write(*, *)konop,' cgnewton iterations without progress'
      write(*, *)' cgnewton returns with ier = 4'
    endif
    return
  endif

3 if(ier.eq.5) then
    if(ipri.eq.1) then
      write(outputunit,*)' Stagnation of cgnewton ier = 5, nef =', nef
      write(outputunit,*)' cgnewton returns with ier = 5'
    else if ( ipri .eq. 2 ) then
      write(*,*)' Stagnation of cgnewton ier = 5, nef =', nef
      write(*,*)' cgnewton returns with ier = 5'
    endif
    return
  endif

  ! A new cgnewton iteration begins here

  kon = kon + 1

  ! Compute cgnewton direction 

  do i = 1, n
    g(i) = - g(i)
  end do

  negra = .false.
  epscg = dmax1(eps, frac*gnor)

  call cg(n, g, x, xnor, xn, d, aux, aux(n+1), aux(2*n+1),&
          iercg, kcg, maxcg, epscg, outputunit, evalg) 
 
  if(ipri.ge.0) then
    if ( ipri.eq.1 ) then
      write(outputunit,*)' cg returned using ',kcg-1,' leq ', maxcg,' iterations'  
    else if (ipri.eq.2) then
      write(*,*)' cg returned using ', kcg-1,' leq ', maxcg,' iterations'
    end if
    if(iercg.eq.1) then
      if ( ipri.eq.1 ) then
        write(outputunit,*)' having found a negative-curvature direction'
      else if (ipri.eq.2) then
        write(*,*)' having found a negative-curvature direction'
      end if
    endif
  endif

  !  Line search
  !  We test whether we obtained a descent direction

  pesca = 0.d0
  do i = 1, n
    pesca = pesca -  g(i)*d(i)
  end do
  if(pesca.ge.0.d0) then
    if(ipri.eq.1) then
      write(outputunit,*)' The direction obtained by cg is not descent'
    else if ( ipri.eq.2 ) then
      write(*,*)' The direction obtained by cg is not descent'
    endif
    !  We will use the gradient as descent direction
    j = 1
    call cg(n, g, x, xnor, xn, d, aux, aux(n+1),aux(2*n+1), iercg, kcg, j, eps,&
            outputunit, evalg)
  endif

  dnor = 0.d0
  do i = 1, n
    dnor = dmax1(dnor, dabs(d(i)))
  end do

  ! Test whether cg stopped with negative-curvature

  if(iercg.eq.1) then
    if(ipri.eq.1) then
      write(outputunit, *)' cg stopped by negative curvature'
    else if ( ipri.eq.2 ) then
      write(*, *)' cg stopped by negative curvature'
    endif
    if(kcg.eq.1) then
      if(ipri.eq.1) then
        write(outputunit, *)' The gradient of f is negative-curvature'
      else if ( ipri.eq.2 ) then
        write(*,*)' The gradient of f is negative-curvature'
      endif
      do i = 1, n
        d(i) = g(i)
      end do
      pesca = 0.d0
      do i = 1, n
        d(i) = d(i)*delta/gnor
        pesca = pesca -  d(i)*g(i)
      end do
      dnor = delta
      negra = .true.
    endif
  endif

  ! Test whether the function decreases enough at d

  t = 1.d0
  ier = 0
2 do i = 1, n
    xn(i) = x(i) + t*d(i)
  end do
  call funcgn(n, xn, fn, warn, evalf)
  nef = nef + 1

  if(fn.le.f + 1.d-4*t*pesca) then
    if(t.ge.1.d0.and.negra) delta = dmin1(100.d0*xnor, 3.d0*delta)
    do i = 1, n
      x(i) = xn(i)
    end do
    if(fn.ge. f - epsnop*(dabs(f)+1.d0)) then
      konop = konop + 1
    else
      konop = 0
  endif
  f = fn
  if(t.lt.1.d0.and.negra) delta = delta/2.d0
  go to 1
  endif
  if(t*dnor.le.dmax1(1.d-30, 1.d-10*xnor)) then
    if(kcg.eq.1) ier = 5
    ! goto3 means Stagnation
    if(ier.eq.5) go to 3
    ier = 5
    t = 1.d0
    do i = 1, n
      d(i) = g(i)*xnor/gnor
    end do
    go to 2     
  endif
  call interp(t, f, fn, pesca)
  go to 2

  return
  end

subroutine interp (amb, f, fn, deriv)
  implicit double precision (a-h,o-z)
  a = (fn - f - deriv * amb) / (amb * amb)
  ambn = - deriv / (2.d0 * a)
  if(ambn .le. 0.05d0 * amb .or. ambn .ge. 0.5d0 * amb) then
    ambn = amb / 2.d0
  endif
  amb = ambn
return
end

subroutine cg(n, b, x, xnor,  xn, s, r, p, ap, ier, k, max, eps, &
              outputunit, evalg)

  implicit double precision (a-h, o-z)
  dimension :: b(n), x(n), r(n), p(n), ap(n), s(n), xn(n)
  integer :: outputunit
  external :: evalg

  !  This subroutine tries to solve A s = b
  !  using conjugate gradients.
  !  A maximum of max iterations is allowed.
  !  The matrix A is assumed to be the Hessian of a function
  !  computed at x.
  !  The matrix A is not given explicitly. Instead, products
  !  A * p must be coded by the user in a subroutine that must
  !  be called hesvec (hessian times vector).
  !  x and b must be auxiliar parameters of hesvec since they
  !  possibly are useful for computing the product. For example,
  !  when cg is called from a nonlinear minimization code 
  !  b is minus-gradient. 
  !  This allows the user to possibly compute Hessian-vector
  !  products by means of difference of gradients:
  !  A * p  \approx [g(x+ h p) + b]/h
  !  When the code finishes having detected a negative or null
  !  curvature direction, it returns with ier=1. In any other
  !  situation it returns with ier=0.
  !  The negative curvature direction is returned in p
  !  ap is an auxiliary vector.  

  ier = 0
  do i = 1, n
    s(i) = 0.d0
    r(i) = b(i)
  end do
  r2 = 0.d0

  ! Here is the loop

  do k = 1, max
    r2a = r2
    r2 = 0.d0
    rnor = 0.d0
    do i = 1, n
      rnor = dmax1(rnor, dabs(r(i)))
      r2 = r2 + r(i)**2
    end do
    if(rnor.le.eps) then
      if ( outputunit > 0 ) then
        write(outputunit,*)' cg finishes with rnor = ', rnor
      else if ( outputunit < 0 ) then
        write(*,*)' cg finishes with rnor = ', rnor
      end if
      return
    endif
    beta = 0.d0
    if(k.gt.1) beta = r2/r2a
    if(k.eq.1) then
      do i = 1, n
        p(i) = r(i)
      end do
    else
      do i = 1, n
        p(i) = r(i) + beta * p(i)
      end do
    end if
    call hesvec(n, x, xnor, xn, b, p, ap, evalg)
    pap = 0.d0
    do i = 1, n
      pap = pap + p(i) * ap(i)
    end do
    if(pap.le.0.d0) then
      ier = 1
      return
    else
      alfa = r2/pap
      do i = 1, n
        s(i) = s(i) + alfa * p(i)
        r(i) = r(i) - alfa * ap(i)
      end do
    endif
  end do


return
end


