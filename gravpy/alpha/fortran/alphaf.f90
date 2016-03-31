module alphaf
  implicit NONE

contains
  subroutine general(res,x,y,modelargs,len)
    integer, intent(in) :: len
    real(8), intent(in) :: x(len)
    real(8), intent(in) :: y(len)
    real(8), intent(in) :: modelargs(7)
    real(8), intent(out):: res(6,len)
    !!f2py integer intent(hide), depend(x) :: len = shape(x,0)

    integer :: i

    do i=1,len
       call single_eval(res(:,i),x(i),y(i),modelargs)
    enddo

    
  end subroutine general

  subroutine single_eval(res,x,y,modelargs)
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: modelargs(7)
    real(8), intent(out):: res(6)

    real(8) :: b,x0,y0,e,te,s,a
    real(8) :: r,cost,sint,front,pot,phir,phir_r,phirr,phix,phiy,phixx,phiyy,phixy,temparr(6)

    external digamma, hyp
    real(8) :: digamma, ifault, hypr, hypim
    
    b  = modelargs(1)
    x0 = modelargs(2)
    y0 = modelargs(3)
    e  = modelargs(4)
    te = modelargs(5)
    s  = modelargs(6)
    a  = modelargs(7)

    r = sqrt(x*x+y*y)
    if (r>0) then
       cost = x/r
       sint = y/r
    else
       res = (/ 0.,0.,0.,0.,0.,0. /) / 0
       return
    endif

    front = (1-e)/2.0*b**(2.0-a)

    if (abs(a) > 0.01) then
       front = front * a
    endif

    if (abs(e) < 1.0e-6) then
       phir   = front*alpha0phir(r,s,a)

       if (r==0) then
          phir_r = front*s**(a-2.0)
       else
          phir_r = phir/r
       endif

       phirr  = -phir_r + 2.0*front*(s*s+r*r)**(0.5*a-1.0)
       phix   = phir_r*x
       phiy   = phir_r*y
       phixx  = phir_r*sint*sint + phirr*cost*cost
       phiyy  = phir_r*cost*cost + phirr*sint*sint
       phixy  = (phirr-phir_r)*sint*cost

       if (abs(s) <= 2.0*r) then 
          call hyp(-(s*s)/(r*r),-0.5*a,-0.5*a,1.0-0.5*a,hypr,hypim)
          pot = hypr* b*b*(r/b)**a/(a*a)
          pot = pot - b*b*(s/b)**a/a*(log(r/s)+0.5*(0.577216-digamma(-0.5*a,ifault)))

          if (abs(a) > 0.01) then
             pot = pot* a
          endif
       else
          if (r==0) then
             pot = 0.
          else
             pot = front*alpha0phir_integral(0.0d0,r,s,a) ! alpha0phir(r,s,a)
          endif
       endif
    else
       temparr = alpha_integral(0.0d0,1.0d0,r,e,sint,cost,s,a)

       pot  = front*(temparr(1)/2.)
       phix = front*(temparr(2)*x)
       phiy = front*(temparr(3)*y)
       phixx= front*(temparr(2)+2.*x*x*temparr(4))
       phiyy= front*(temparr(3)+2.*y*y*temparr(5))
       phixy= front*(2.*x*y*temparr(6))


    endif

    res = (/ pot,phix,phiy,phixx,phiyy,phixy /)

  end subroutine single_eval

  function alpha0phir(r,s,a) result(res)
    real(8), intent(in) :: r,s,a
    real(8) :: res
    
    if (r/s < 1.0e-4) then
       res = r*s**(a-2.0)
    else if (a==0.0) then
       res = log(1.0+r*r/(s*s))/r
    else
       res = 2.0/(a*r)*((s*s+r*r)**(a/2.0) - s**a)
    endif

  end function alpha0phir

  
  function alpha0phir_integral(lower,upper,s,a) result(integral)
    use cui
    real(8), intent(in) :: lower,upper,s,a
    real(8)  :: integral
    integer  :: ndim,maxeval,key,neval,fail
    integer  :: rgtype
    real(8)  :: epsrel,epsabs
    real(8)  :: limits(1,2),abserr

    ndim = 1
    epsrel = 1d-3
    epsabs = 1d-9
    maxeval = 100000
    key = 0
    limits(1,1) = lower
    limits(1,2) = upper
    rgtype = 1
    fail = -1

    call cubatr(ndim,f,limits,rgtype,integral,abserr,&
         key=key, maxpts=maxeval,neval=neval, ifail=fail)
    

  contains
    function f(ncomp,x) result(res)
      real(8), intent(in) :: x(:)
      integer, intent(in) :: ncomp
      real(8) :: res(ncomp),r
      r = x(1)
      
      if (r/s < 1.0e-4) then
         res = r*s**(a-2.0)
      else if (a==0.0) then
         res = log(1.0+r*r/(s*s))/r
      else
         res = 2.0/(a*r)*((s*s+r*r)**(a/2.0) - s**a)
      endif
    end function f

  end function alpha0phir_integral


  function alpha_integral(lower,upper,r,e,sint,cost,s,a) result(integral)
    use cui
    real(8), intent(in) :: lower,upper,r,e,sint,cost,s,a
    real(8)  :: integral(6)
    integer  :: ndim,ncomp,maxeval,key,nregions,neval,fail
    integer  :: rgtype(1)
    real(8)  :: epsrel,epsabs
    real(8)  :: limits(1,2,1),abserr(6)

    ndim = 1
    ncomp = 6
    nregions = 1
    epsrel = 1d-3
    epsabs = 1d-9
    maxeval = 100000
    key = 0
    limits(1,1,1) = lower
    limits(1,2,1) = upper
    rgtype(1) = 1
    fail = -1

    call cubatr(ndim,ncomp,f,nregions,limits,rgtype,integral,abserr,&
         key=key, maxpts=maxeval,neval=neval, ifail=fail)

  contains
    function f(ncomp,x) result(res)
      real(8), intent(in) :: x(:)
      integer, intent(in) :: ncomp
      real(8) ::  res(ncomp)
      
      real(8) :: u,q,t0,t1,t3,t5,t6,t7,mphiu,k,kp
      u = x(1)
      
      q  = 1.0-e
      t0 = 1.0-(1.0-q*q)*u
      t1 = 1.0/sqrt(t0)
      t3 = t1/t0
      t5 = t3/t0
      t6 = cost*cost+sint*sint/(1.0 - (1.0-q*q)*u)*r*r
      t7 = t6*u

      if (u==0.0) then
         mphiu = t6*s**(a-2.0)
      else
         mphiu = sqrt(t7)*alpha0phir(sqrt(t7),s,a)/u
      endif
      
      k  =(s*s+t7)**(a/2.0-1.0)
      kp = k*(a/2.0-1.0)/(s*s+t7)

      res = (/ t1*mphiu,t1*k,t3*k,t1*kp*u,t5*kp*u,t3*kp*u /)
    end function f
  end function alpha_integral
  
end module alphaf
