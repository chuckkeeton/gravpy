module sief
  implicit NONE

contains

  subroutine elliptical(z,x,y,modelargs,xlen)

    !arguments
    integer, intent(in) :: xlen
    real(8), intent(in) :: x(xlen)
    real(8), intent(in) :: y(xlen)
    real(8), intent(in) :: modelargs(6)
    real(8), intent(out):: z(6,xlen)
    !integer, intent(in) :: xlen
    !!f2py integer intent(hide), depend(x) :: xlen = shape(x,0)

    !locals
    real(8) :: b,x0,y0,e,te,s
    real(8), dimension(xlen) :: x2,y2,s2,q,q2,om,rt,psi,psis
    real(8), dimension(xlen) :: phix,phiy,invdenom,phixx,phiyy,phixy,pot

    b  = modelargs(1)
    x0 = modelargs(2)
    y0 = modelargs(3)
    e  = modelargs(4)
    te = modelargs(5)
    s  = modelargs(6)

    if (s == 0.0) then
       s = 1e-4
    endif

    x2  = x*x
    y2  = y*y
    s2  = s*s
    q   = 1.0-e
    q2  = q*q
    om  = 1.0-q2
    rt  = sqrt(om)
    psi = sqrt(q2*(s2+x2)+y2)
    psis= psi + s

    phix = b*q/rt *atan(rt*x/psis)
    phiy = b*q/rt *atanh(rt*y/(psi+s*q2))

    invdenom = 1/(psi*(om*x2+psis*psis))
    phixx = b*q*(psi*psis-q2*x2)*invdenom
    phiyy = b*q*(x2+s*psis)*invdenom
    phixy = -b*q*x*y*invdenom

    pot = b*q*s*(-0.5*log(psis*psis+om*x2) + log(s*(1.0+q)) ) + x*phix+y*phiy


    z(1,:) = pot
    z(2,:) = phix
    z(3,:) = phiy
    z(4,:) = phixx
    z(5,:) = phiyy
    z(6,:) = phixy

  end subroutine elliptical

  subroutine spherical(z,x,y,modelargs,xlen)

    !arguments
    integer, intent(in) :: xlen
    real(8), intent(in) :: x(xlen)
    real(8), intent(in) :: y(xlen)
    real(8), intent(in) :: modelargs(6)
    real(8), intent(out):: z(6,xlen)

    !locals
    real(8) :: b,x0,y0,e,te,s
    real(8), dimension(xlen) :: rad,sprad,invdenom
    real(8), dimension(xlen) :: pot,phix,phiy,phixx,phiyy,phixy    

    b  = modelargs(1)
    x0 = modelargs(2)
    y0 = modelargs(3)
    e  = modelargs(4)
    te = modelargs(5)
    s  = modelargs(6)

    if (s == 0.0) then
       s = 1e-4
    endif


    rad = sqrt(x*x+y*y+s*s)
    sprad = s + rad
    invdenom = 1/(rad*sprad*sprad)

    pot = b * (rad-s*(1+log(sprad/(2*s))))
    phix = b * x / sprad
    phiy = b * y / sprad
    phixx = b * (s*sprad + y*y) * invdenom
    phiyy = b * (s*sprad + x*x) * invdenom
    phixy = -b*x*y *invdenom

    z(1,:) = pot
    z(2,:) = phix
    z(3,:) = phiy
    z(4,:) = phixx
    z(5,:) = phiyy
    z(6,:) = phixy


  end subroutine spherical

  subroutine which_function(z,x,y,modelargs)
    !arguments
    real(8), intent(out) :: z(:,:)
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(:)
    real(8), intent(in) :: modelargs(6)

    !locals
    real(8) :: e 
    e = modelargs(4)


    if (e==0.0) then
       call spherical(z,x,y,modelargs,size(x))
    else
       call elliptical(z,x,y,modelargs,size(x))
    endif

  end subroutine which_function


  subroutine phiarray(z,x,y,modelargs,vec,xlen) 
    !arguments
    real(8), intent(in) :: x(:) !x(xlen)
    real(8), intent(in) :: y(:) !y(xlen)
    real(8), intent(in) :: modelargs(6)
    logical, intent(in) :: vec
    !f2py logical optional,intent(in) :: vec
    integer, intent(in) :: xlen
    !f2py integer, intent(in), depend(x) :: xlen = size(x)
    real(8), intent(out):: z(6,xlen)

    call which_function(z,x,y,modelargs)

  end subroutine phiarray

  subroutine phiarray_vec(z,x,y,modelargs,vec,xlen,masses)
    !arguments
    real(8), intent(in) :: x(:,:) !x(xlen)
    real(8), intent(in) :: y(:,:) !y(xlen)
    real(8), intent(in) :: modelargs(:,:)
    logical, intent(in) :: vec
    !f2py logical optional,intent(in) :: vec
    integer, intent(in) :: xlen
    !f2py integer, intent(in), depend(x) :: xlen = shape(x,0)
    integer, intent(in) :: masses
    !f2py integer, intent(in), depend(modelargs) :: masses = shape(modelargs,1)
    real(8), intent(out):: z(masses,6,xlen)
    integer :: i

    do i = 1,masses
       call which_function(z(i,:,:),x(i,:),y(i,:),modelargs)
    end do
  end subroutine phiarray_vec

end module sief
