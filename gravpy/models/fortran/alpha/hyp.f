c***********************************************************************
c
c   program :  hypergeometric function
c
c   notation :  F(a,b;c;z)
c
c   reference:  see article `Computing the hypergeometric function'
c               by R.C. Forrey, J. Comp. Phys. 137, 79-100 (1997).
c
c   send comments to:
c
c        Robert C. Forrey
c        Institute for Theoretical Atomic and Molecular Physics
c        Harvard-Smithsonian Center for Astrophysics
c        60 Garden Street Cambridge, MA 02138
c        rforrey@cfa.harvard.edu
c
c***********************************************************************
c
c  subroutine name    - hyp
c
c  computation
c  performed          - calculates the hypergeometric function
c
c  usage              - call hyp(z,a,b,c,re,im)
c
c  arguments     
c                  z  - the independent variable of the hypergeometric
c                       function (must be real).
c
c               a,b,c - real parameters of the hypergeometric function.
c
c               re,im - the real and imaginary parts of the
c                       hypergeometric function.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine hyp(z,a,b,c,re,im)

      real*8  zero,one,two,half
      parameter (zero=0.d0,one=1.d0,two=2.d0,half=0.5d0)
      integer flag,flag2,neps
      real*8   a,b,c,z,w,f,f1,f2,gamm,tol,test,pi,machep,re2,
     #         alpha0,alpha1,rn,binom,eps,re,im,x1,x2,x3,x4,
     #         coeff1,coeff2,coeff3,coeff4,temp1,temp2,term,
     #         a1,b1,c1,a2,b2,c2,alpha2
      logical fix
      common /bcoeff/binom(5151)

c  tabulate the binomial coefficients and set the defaults

      fix=.false.
      call binomc
      tol=.1d0
      im=zero
      nmax=100
      n=5

      call geteps(machep,neps)
      pi=dacos(-1.d0)
c     write(20,'(/a,i3/)') ' machine epsilon is machep = (1/2)**',neps
c     write(20,'(a,d23.15//)') ' to machine accuracy, pi = ',pi

c  handle the special case when z=1

      if (z.eq.one) then
        re=gamm(c)*gamm(c-a-b)/gamm(c-a)/gamm(c-b)
        return
      endif

c  transform to a new variable w which lies between 0 and 1/2

      if(z .lt. -one) then
        a1=a
        b1=c-b
        c1=a-b+1
        a2=b
        b2=c-a
        c2=b-a+1
        w=one/(one-z)
        flag=1
      elseif( (-one.le.z) .and. (z.lt.zero) )then
        a1=a
        b1=c-b
        c1=c
        a2=a1
        b2=b1
        c2=c1
        w=z/(z-one)
        flag=2
      elseif( (zero.le.z) .and. (z.le.half) )then
        a1=a
        b1=b
        c1=c
        a2=a1
        b2=b1
        c2=c1
        w=z
        flag=3
      elseif( (half.lt.z) .and. (z.le.one) )then
        a1=a
        b1=b
        c1=a+b-c+1
        a2=c-a
        b2=c-b
        c2=c-a-b+1
        w=one-z
        flag=4
      elseif( (one.lt.z) .and. (z.le.two) )then
        a1=a
        b1=a-c+1
        c1=a+b-c+1
        a2=c-a
        b2=1-a
        c2=c-a-b+1
        w=one-(one/z)
        flag=5
      elseif(two.lt.z)then
        a1=a
        b1=a-c+1
        c1=a-b+1
        a2=b-c+1
        b2=b
        c2=b-a+1
        w=one/z
        flag=6
      endif

c  compute the hypergeometric function of z via the transformation
c  theory

      if (flag .eq. 1)then
        k=nint(dble(a-b))
        test=a-b-dble(k)
        if (dabs(test).lt.tol) then
          fix=.true.
          flag2=0
          if (a.lt.b) then
            temp1=a
            temp2=b
            b=temp1
            a=temp2
            flag2=1
          endif
          k=nint(dble(a-b))
          eps=a-b-dble(k)
          call fix1(a,b,c,n,k,f1,w,machep,eps)
          do m=n+5,nmax,5
            call fix1(a,b,c,m,k,f2,w,machep,eps)
            test=dabs(f1-f2)
            if(test.le.machep)go to 30
            f1=f2
          end do
          write(*,*)'fix1 warning: not completely converged'
  30      re=f2    
          if (flag2.eq.1) then
            a=temp1
            b=temp2
          endif
        else
          call hyper(w,a1,b1,c1,f1,machep)
          call hyper(w,a2,b2,c2,f2,machep)

            x1=b
            coeff1=one
  1         if (x1.lt.one) then
              coeff1=coeff1*x1
              x1=x1+one
              go to 1
            endif

            x2=c-a
            coeff2=one
  2         if (x2.lt.one) then
              coeff2=coeff2*x2
              x2=x2+one
              go to 2
            endif

            x3=a
            coeff3=one
  3         if (x3.lt.one) then
              coeff3=coeff3*x3
              x3=x3+one
              go to 3
            endif

            x4=c-b
            coeff4=one
  4         if (x4.lt.one) then
              coeff4=coeff4*x4
              x4=x4+one
              go to 4
            endif

          re=(w**a)*gamm(c)*gamm(b-a)*coeff1*coeff2/gamm(x1)/gamm(x2)*f1
     #     +(w**b)*gamm(c)*gamm(a-b)*coeff3*coeff4/gamm(x3)/gamm(x4)*f2

        endif
      elseif (flag .eq. 2)then
        call hyper(w,a1,b1,c1,f1,machep)
        re=((one-w)**a)*f1
      elseif (flag .eq. 3)then
        call hyper(w,a1,b1,c1,f1,machep)
        re=f1
      elseif (flag .eq. 4)then
        k=nint(dble(c-a-b))
        test=c-a-b-dble(k)
        if (dabs(test).lt.tol) then
          fix=.true.
          if (k.ge.zero) then
            eps=c-a-b-dble(k)
            call fix4a(a,b,c,n,k,f1,w,machep,eps)
            do m=n+5,nmax,5
            call fix4a(a,b,c,m,k,f2,w,machep,eps)
            test=dabs(f1-f2)
            if(test.le.machep)go to 31
            f1=f2
            end do
            write(*,*)'fix4a warning: not completely converged'
  31        re=f2   
          else
            k=-k
            eps=c-a-b+dble(k)
            call fix4b(a,b,c,n,k,f1,w,machep,eps)
            do m=n+5,nmax,5
            call fix4b(a,b,c,m,k,f2,w,machep,eps)
            test=dabs(f1-f2)
            if(test.le.machep)go to 32
            f1=f2
          end do
          write(*,*)'fix4b warning: not completely converged'
  32      re=f2   
          endif
        else
          call hyper(w,a1,b1,c1,f1,machep)
          call hyper(w,a2,b2,c2,f2,machep)

            x1=c-a
            coeff1=one
  5         if (x1.lt.one) then
              coeff1=coeff1*x1
              x1=x1+one
              go to 5
            endif

            x2=c-b
            coeff2=one
  6         if (x2.lt.one) then
              coeff2=coeff2*x2
              x2=x2+one
              go to 6
            endif

            x3=a
            coeff3=one
  7         if (x3.lt.one) then
              coeff3=coeff3*x3
              x3=x3+one
              go to 7
            endif

            x4=b
            coeff4=one
  8         if (x4.lt.one) then
              coeff4=coeff4*x4
              x4=x4+one
              go to 8
            endif

          re=gamm(c)*gamm(c-a-b)*coeff1*coeff2/gamm(x1)/gamm(x2)*f1
     #       +w**(c-a-b)*gamm(c)*gamm(a+b-c)*coeff3*coeff4/gamm(x3)
     #       /gamm(x4)*f2

        endif
      elseif (flag .eq. 5)then
        k=nint(dble(c-a-b))
        test=c-a-b-dble(k)
        if (dabs(test).lt.tol) then
          fix=.true.
          if (k.ge.zero) then
            eps=c-a-b-dble(k)
            call fix5a(a,b,c,n,k,f1,im,w,machep,eps,pi)
            do m=n+5,nmax,5
            call fix5a(a,b,c,m,k,f2,im,w,machep,eps,pi)
            test=dabs(f1-f2)
            if(test.le.machep)go to 33
            f1=f2
            end do
            write(*,*)'fix5a warning: not completely converged'
  33        re=f2   
          else
            k=-k
            eps=c-a-b+dble(k)
            call fix5b(a,b,c,n,k,f1,im,w,machep,eps,pi)
            do m=n+5,nmax,5
            call fix5b(a,b,c,m,k,f2,im,w,machep,eps,pi)
            test=dabs(f1-f2)
            if(test.le.machep)go to 34
            f1=f2
            end do
            write(*,*)'fix5b warning: not completely converged'
  34        re=f2   
          endif
        else
          call hyper(w,a1,b1,c1,f1,machep)
          call hyper(w,a2,b2,c2,f2,machep)

            x1=c-a
            coeff1=one
  11        if (x1.lt.one) then
              coeff1=coeff1*x1
              x1=x1+one
              go to 11
            endif

            x2=c-b
            coeff2=one
  12        if (x2.lt.one) then
              coeff2=coeff2*x2
              x2=x2+one
              go to 12
            endif

            x3=a
            coeff3=one
  13        if (x3.lt.one) then
              coeff3=coeff3*x3
              x3=x3+one
              go to 13
            endif

            x4=b
            coeff4=one
  14        if (x4.lt.one) then
              coeff4=coeff4*x4
              x4=x4+one
              go to 14
            endif

          re=(one-w)**a*gamm(c)*gamm(c-a-b)*coeff1*coeff2/gamm(x1)
     #       /gamm(x2)*f1+w**(c-a-b)*(one-w)**b*dcos(pi*(c-a-b))
     #       *gamm(c)*gamm(a+b-c)*coeff3*coeff3/gamm(x3)/gamm(x4)*f2

        endif
      elseif (flag .eq. 6)then
        k=nint(dble(a-b))
        test=a-b-dble(k)
        if (dabs(test).lt.tol) then
          fix=.true.
          flag2=0
          if (a.lt.b) then
            temp1=a
            temp2=b
            b=temp1
            a=temp2
            flag2=1
          endif
          k=nint(dble(a-b))
          eps=a-b-dble(k)
          call fix6(a,b,c,n,k,f1,im,w,machep,eps,pi)
          do m=n+5,nmax,5
          call fix6(a,b,c,m,k,f2,im,w,machep,eps,pi)
          test=dabs(f1-f2)
          if(test.le.machep)go to 35
          f1=f2
          end do
          write(*,*)'fix6 warning: not completely converged'
  35      re=f2   
          if (flag2.eq.1) then
            a=temp1
            b=temp2
          endif
        else
          call hyper(w,a1,b1,c1,f1,machep)
          call hyper(w,a2,b2,c2,f2,machep)

            x1=b
            coeff1=one
  15        if (x1.lt.one) then
              coeff1=coeff1*x1
              x1=x1+one
              go to 15
            endif

            x2=c-a
            coeff2=one
  16        if (x2.lt.one) then
              coeff2=coeff2*x2
              x2=x2+one
              go to 16
            endif

            x3=a
            coeff3=one
  17        if (x3.lt.one) then
              coeff3=coeff3*x3
              x3=x3+one
              go to 17
            endif

            x4=c-b
            coeff4=one
  18        if (x4.lt.one) then
              coeff4=coeff4*x4
              x4=x4+one
              go to 18
            endif

          re=w**a*dcos(pi*a)*gamm(c)*gamm(b-a)*coeff1*coeff2/gamm(x1)
     #       /gamm(x2)*f1+w**b*dcos(pi*b)*gamm(c)*gamm(a-b)*coeff3
     #       *coeff4/gamm(x3)/gamm(x4)*f2

        endif
      endif

c     if(fix)then
c       write(20,'(2(a6,2x,i3,2x))')'case=',flag,'m=',m
c     endif

      return
      end


c***********************************************************************
c
c  subroutine name     - hyper
c
c  computation
c  performed           - calculates the hypergeometric function,
c                        f(a,b;c;w), from its power series for 0<w<.5.
c
c  usage               - call hyper(w,a,b,c,f)
c
c  arguments        
c                   w  - the transformed independent variable.
c
c                a,b,c - the parameters of the hypergeometric function.
c
c                   f  - the computed value of the hypergeometric
c                        function which is returned to the caller.
c
c  precision           - double
c
c  language            - fortran
c
c***********************************************************************

      subroutine hyper(w,a,b,c,f,machep)
      implicit none
      integer i,m,n,nmax,k,k0,k1
      parameter (nmax=100)
      real*8  a,b,c,w,f,alpha0,alpha1,rn,gamm,term,machep,binom
      common /bcoeff/binom(5151)

c  compute the number of sums needed to get good convergence

      alpha1=a+b-c
      k1=nint(dble(alpha1))

      do 10 n=1,nmax

        rn=0.d0
        alpha0=(a+n+1)*(b+n+1)/(c+n+1)-(n+1)
        k0=nint(dble(alpha0))
        k=max(k0,k1)
        if (k.le.1) k=1
        if (n+k.ge.100) then
           write(*,*)'error in hyp:  binomial coefficient routine
     #no longer valid'

          return
        endif

        do m=0,k
          rn=rn+binom((n+k+1)*(n+k+2)/2+m+1)
        end do   

        term=1.d0
        do i=1,n+1
          term=(a+i-1)*(b+i-1)/(c+i-1)/(k+i)*term
        end do

        rn=rn*term*(w**(n+1))/(1-w)
        if (dabs(rn).lt.machep) go to 100

 10   continue

      write(*,*)'error in hyp:  nmax not large enough'
      return

 100  continue
c     write(20,'(2(a6,i3,5x))')'n=',n

c  evaluate the hypergeometric function of w from its power series

      term=1.d0
      f=1.d0
      do 20 k=1,n
        term=term*(a+k-1)*(b+k-1)/(c+k-1)*w/k
        f=f+term
  20  continue

      return
      end

c***********************************************************************
c
c  function name      - gamm
c
c  computation
c  performed          - calculates the gamma function
c
c  usage              - gamm(x)
c
c  argument         x - any real number (excluding negative integers).
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      function gamm(x)

      real*8  zero,one
      parameter(zero=0.d0,one=1.d0)
      real*8  x,xx,coeff,gamm,g

c  scale input variable and change it's name
c  so that it does not get destroyed

      xx=x-one
      coeff=one

  100 if ( (zero.le.xx) .and. (xx.le.one) ) go to 200

        if (xx.lt.zero) then
          xx=xx+one
          coeff=coeff/xx
          go to 100
        else
          coeff=xx*coeff
          xx=xx-one
          go to 100
        endif

  200 gamm=coeff*g(xx)

      return
      end

c***********************************************************************
c
c  function name     - g
c
c  computation
c  performed         - calculates gamma(xx+1) for xx in the interval
c                      [0,1] using clenshaw's recurrence formula with
c                      tchebychev polynomials and the tabulated
c                      expansion coefficients.
c
c  usage             - g(xx)
c
c  argument       xx - scaled value for 'x' in 'gamm(x)'.
c
c  precision         - double
c
c  language          - fortran
c
c***********************************************************************

      function g(xx)

      real*8  zero,one,two
      parameter (zero=0.d0,one=1.d0,two=2.d0)
      real*8  c(0:41),xx,y,y1,y2,g

c  use clenshaw recurrence scheme with tchebychev polynomials
c  and the expansion coefficients tabulated in 'cheb' for 0<xx<1 .

      call cheb(c,41,1)

      y1=zero
      y2=zero

      do 10 k=41,1,-1
        y=two*(two*xx-one)*y1-y2+c(k)
        y2=y1
        y1=y
   10 continue

      g=-y2+(two*xx-one)*y1+c(0)

      return
      end

c**********************************************************************
c
c  subroutine name    - fix1
c
c  computation
c  performed          - calculates the hypergeometric function for
c                       z less than -1 when a-b is near an integer.
c
c  usage              - call fix1(a,b,c,n,k,f,w,machep,eps)
c
c  arguments    a,b,c - parameters of the hypergeometric function.
c
c                  n  - the upper limit of the finite series expansion
c                       of the hypergeometric function.
c
c                  k  - equals the nearest integer of a-b.
c
c                  f  - computed value of the hypergeometric function.
c
c                  w  - transformed independent variable.
c
c              machep - equals machine epsilon.
c
c                eps  - equals a-b-k.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine fix1(a,b,c,n,k,f,w,machep,eps)

      real*8  zero,one,two,four,eighth,seven,eight,sxteen
      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,eighth=1.d0/8.d0,
     #           seven=7.d0,eight=8.d0,sxteen=16.d0,nmax=100)
      real*8   a,b,c,w,f,eps,gamm,machep,test,arg,rn,sum,et1,et2,
     #         term1,term2,term3,term4,term5,term6,term7,term8,
     #         temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,
     #         coeff,coeff1,coeff2,coeff3,coeff4,x,x1,x2,x3,x4,
     #         t1(0:80),t2(0:80),t3(0:80),t4(0:80),c1(0:80),c2(0:80),
     #         c3(0:80),c4(0:80),f1(0:80),f2(0:80),f3(0:80),f4(0:80),
     #         g1(0:nmax),g2(0:nmax),g3(0:nmax),g4(0:nmax),g5(0:nmax),
     #         fff1(0:nmax),ff1(0:nmax),fff2(0:nmax),ff2(0:nmax),
     #         ff3(0:nmax),ff4(0:nmax),poch1(0:nmax),poch2(0:nmax),
     #         e1(0:nmax),e2(0:nmax),e3(0:nmax),e4(0:nmax)

      integer  flag

c  calculate the extra terms

      x=b-one
      sum=zero
      coeff=one
      flag=0
  1   if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 1
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
  2     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 2
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 1
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 3 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 3      continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 4 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 4      continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 5 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 5      continue
      endif

      if (flag.eq.0) then
        x1=b
        coeff1=one
 6      if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 6
        endif
        x2=b+eps
        coeff2=one
 7      if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 7
        endif
        temp=sum+coeff*temp
        et1=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 8      if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 8
        endif
        x2=x+one+eps
        coeff2=one
 9      if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 9
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et1=sum+coeff*temp
      endif
      et1=-et1
c     write(10,*)et1,(one/gamm(a-dble(k)-eps)-one/gamm(a-dble(k)))
c    #                                           /eps
      x=c-a+dble(k)-one
      sum=zero
      coeff=one
      flag=0
  10  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 10
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 11     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 11
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 10
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c2,41,1)
        t2(0)=one
        t2(1)=two*(x+eps)-one
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 12 i=2,41
          t2(i)=(four*(x+eps)-two)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-two)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
 12     continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c2,55,2)
        t2(0)=one
        t2(1)=two*(x+eps)
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 13 i=2,55
          t2(i)=four*(x+eps)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+four*x*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
 13     continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c2,34,3)
        t2(0)=one
        t2(1)=two*(x+eps)-two
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 14 i=2,34
          t2(i)=(four*(x+eps)-four)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-four)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
 14     continue
      endif

      if (flag.eq.0) then
        x1=c-a+dble(k)
        coeff1=one
 15     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 15
        endif
        x2=c-a+dble(k)+eps
        coeff2=one
 16     if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 16
        endif
        temp2=sum+coeff*temp2
        et2=-temp2*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 17     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 17
        endif
        x2=x+one+eps
        coeff2=one
 18     if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 18
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et2=sum+coeff*temp2
      endif

c     write(10,*)et2,(one/gamm(c-a+dble(k)+eps)-one/gamm(c-a+dble(k)))
c    #                                                         /eps

c  calculate the f-functions

      x1=a
      coeff1=one
 20   if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 20
      endif

      x2=a-dble(k)
      coeff2=one
 21   if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 21
      endif

      x3=c-a
      coeff3=one
 22   if (x3.lt.one) then
        coeff3=x3*coeff3
        x3=x3+one
        go to 22
      endif

      x4=c-a+dble(k)
      coeff4=one
 23   if (x4.lt.one) then
        coeff4=x4*coeff4
        x4=x4+one
        go to 23
      endif

      coeff=one
      arg=-eps-dble(k)
 24   if (arg.lt.-eps) then
        coeff=coeff/arg
        arg=arg+one
        go to 24
      endif

      fff1(0)=one
      fff2(0)=one
      ff1(0)=one
      ff2(0)=one
      do 25 i=1,k
        fff1(0)=(b+dble(i-1))*fff1(0)
        ff2(0)=(c-a+dble(i-1))*ff2(0)
 25   continue

      fff1(0)=fff1(0)*coeff1/gamm(x1)
      fff2(0)=coeff3/gamm(x3)
      ff1(0)=coeff2/gamm(x2)
      ff2(0)=ff2(0)*coeff4/gamm(x4)
      ff3(0)=(-1)**(k+1)*gamm(one-eps)*coeff
      ff4(0)=gamm(one+eps)

c     do 26 i=1,n
c       fff1(i)=(b+dble(k+i-1))*fff1(i-1)
c       fff2(i)=(c-b+dble(i-1))*fff2(i-1)
c       ff1(i)=(a+dble(i-1))*ff1(i-1)
c       ff2(i)=(c-a+dble(k+i-1))*ff2(i-1)
c       ff3(i)=ff3(i-1)/(eps+dble(i)+dble(k))
c       ff4(i)=ff4(i-1)/(dble(i)-eps)
c26   continue

c     do 27 i=0,n
c       write(10,*)'fff1=',fff1(i),gamm(b+i+k)/gamm(a)/gamm(b)
c       write(10,*)'ff1=',ff1(i),gamm(a+i)/gamm(a)/gamm(a-k)
c       write(10,*)'fff2=',fff2(i),gamm(c-b+i)/gamm(c-a)/gamm(c-b)
c       write(10,*)'ff2=',ff2(i),gamm(c-a+i+k)/gamm(c-a)/gamm(c-a+k)
c       write(10,*)'ff3=',ff3(i),(-1)**(k+i)*eps*gamm(-k-i-eps)
c       write(10,*)'ff4=',ff4(i),(-1)**i*eps*gamm(eps-i)
c27   continue

c   calculate  g1,g2

      x1=a
      coeff1=one
 100  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 100
      endif

      x2=c-a
      coeff2=one
 101  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 101
      endif

      g1(0)=zero
      g2(0)=zero
      poch1(0)=one
      poch2(0)=one
      do 102 i=1,k
        g1(0)=g1(0)*(a-eps+dble(i-k-1))-poch1(0)
        poch1(0)=poch1(0)*(a+dble(i-k-1))
 102  continue

      g1(0)=g1(0)*coeff1/gamm(x1)
      g2(0)=g2(0)*coeff2/gamm(x2)
      poch1(0)=poch1(0)*coeff1/gamm(x1)
      poch2(0)=poch2(0)*coeff2/gamm(x2)
      do 103 i=1,n
        poch1(i)=(a+dble(i-1))*poch1(i-1)
        poch2(i)=(c-a+dble(k+i-1))*poch2(i-1)
        g1(i)=g1(i-1)*(a-eps+dble(i-1))-poch1(i-1)
        g2(i)=g2(i-1)*(c-a+eps+dble(k+i-1))+poch2(i-1)
 103  continue

c     do 104 i=0,n
c       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
c       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
c104  continue

c  calculate  g3,g4,g5

      x3=zero
      call cheb(c3,55,2)
      t3(0)=one
      t3(1)=two*(x3-eps)
      f3(0)=zero
      f3(1)=-two
      g3(0)=c3(1)*f3(1)

      x4=zero
      call cheb(c4,55,2)
      t4(0)=one
      t4(1)=two*(x4+eps)
      f4(0)=zero
      f4(1)=two
      g4(0)=c4(1)*f4(1)

      do 105 i=2,55
        t3(i)=four*(x3-eps)*t3(i-1)-t3(i-2)
        t4(i)=four*(x4+eps)*t4(i-1)-t4(i-2)
        f3(i)=-four*t3(i-1)+four*x3*f3(i-1)-f3(i-2)
        f4(i)=four*t4(i-1)+four*x4*f4(i-1)-f4(i-2)
        g3(0)=g3(0)+c3(i)*f3(i)
        g4(0)=g4(0)+c4(i)*f4(i)
 105  continue

      g3(0)=-g3(0)
      do 106 i=-k,-1
        g3(0)=(g3(0)+one/gamm(dble(k+i+2)))/(dble(k+i+1)+eps)
 106  continue

      test=dabs(eps*dlog(w))
      temp=-dlog(w)
      if (eps.le.zero) then
         if (test.ge.eighth) then
           temp=(one-exp(test))/eps
         else
           i=1
 107       rn=(eps**(i)*(dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 108
           temp=temp-rn
           i=i+1
           go to 107
         endif
 108     g5(0)=w**(a-eps)*temp
      else
         if (test.ge.eighth) then
           temp=(exp(test)-one)/eps
         else
           i=1
 109       rn=(eps**(i)*(-dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 110
           temp=temp+rn
           i=i+1
           go to 109
         endif
 110     g5(0)=(w**a)*temp
      endif

c     write(10,*)g3(0),(-1)**k*gamm(-k-eps)+one/eps/gamm(dble(k+1))
c     write(10,*)g4(0),gamm(eps)-one/eps
c     write(10,*)g5(0),w**(a-eps)/eps-w**a/eps

      e1(0)=one/gamm(dble(k+1))
      e2(0)=one
      do 120 i=1,n
        e1(i)=e1(i-1)/dble(k+i)
        e2(i)=e2(i-1)/dble(i)
        g3(i)=(g3(i-1)+e1(i))/(dble(k+i)+eps)
        g4(i)=(g4(i-1)+e2(i))/(dble(i)-eps)
        g5(i)=w*g5(i-1)
 120  continue

      e1(0)=one
      e2(0)=one
      e3(0)=one
      e4(0)=one
      do 130 i=1,k
        e2(0)=(c-a+dble(i-1))/dble(i)*e2(0)
        e4(0)=e4(0)/dble(i)
 130  continue

c     do 140 i=1,n
c       e1(i)=(a+dble(i-1))/dble(i)*e1(i-1)
c       e2(i)=(c-a+dble(k+i-1))/dble(k+i)*e2(i-1)
c       e3(i)=e3(i-1)/dble(i)
c       e4(i)=e4(i-1)/dble(k+i)
c140  continue

c   put everything back together again

      term1=-gamm(c)*w**a*e3(0)*fff2(0)*ff3(0)*(-1)**k
      term2=gamm(c)*w**a*e3(0)*fff1(0)*ff3(0)*(-1)**k
      term3=gamm(c)*w**a*e3(0)*fff1(0)*ff2(0)*(-1)**k
      term4=gamm(c)*w**a*e4(0)*fff1(0)*ff2(0)*(-1)**k
      term5=gamm(c)*e4(0)*fff1(0)*ff2(0)*ff4(0)*(-1)**k
      term6=gamm(c)*w**a*e1(0)*fff2(0)*ff3(0)*(-1)**k
      term7=gamm(c)*w**(a-eps)*e2(0)*fff1(0)*ff4(0)*(-1)**k

      temp=g1(0)*term1+g2(0)*term2+g3(0)*term3+g4(0)*term4
     #                                       +g5(0)*term5
      temp1=term6
      temp2=term7
      do 150 i=1,n
        term1=term1*w*(c-b+dble(i-1))/(eps+dble(i+k))/dble(i)
        term2=term2*w*(b+dble(k+i-1))/(eps+dble(i+k))/dble(i)
        term3=term3*w*(b+dble(k+i-1))*(c-a+dble(k+i-1))/dble(i)
        term4=term4*w*(b+dble(k+i-1))*(c-a+dble(k+i-1))/dble(k+i)
        term5=term5*(b+dble(k+i-1))*(c-a+dble(k+i-1))/dble(k+i)
     #                                       /(dble(i)-eps)
        term6=term6*w*(a+dble(i-1))/dble(i)*(c-b+dble(i-1))
     #                                       /(eps+dble(i+k))
        term7=term7*w*(c-a+dble(k+i-1))/dble(k+i)*(b+dble(k+i-1))
     #                                       /(dble(i)-eps)
        temp=temp+g1(i)*term1+g2(i)*term2+g3(i)*term3+g4(i)*term4
     #                                       +g5(i)*term5
        temp1=temp1+term6
        temp2=temp2+term7
 150  continue

c  calculate the finite series term

      poch1(0)=one
      poch2(0)=one
      do 160 i=1,k-1
        poch1(i)=(b+dble(i-1))*poch1(i-1)
        poch2(i)=(c-a+dble(i-1))*poch2(i-1)
 160  continue

      temp6=zero
      do 170 i=0,k-1
        temp6=temp6+poch1(i)*poch2(i)*gamm(dble(k-i)+eps)
     #                  /gamm(dble(i+1))*(-w)**i
 170  continue

      x1=a
      coeff1=one
 180  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 180
      endif

      x2=c-b
      coeff2=one
 190  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 190
      endif

      f=temp+et1*temp1+et2*temp2+coeff1*coeff2/gamm(x1)/gamm(x2)
     #                             *gamm(c)*w**(a-eps-dble(k))*temp6

c  alternative method  (must also turn on the individual functions)

c     temp3=zero
c     temp4=zero
c     temp5=zero
c     do 200 i=0,n
c       term1=-gamm(c)*w**(a+dble(i))*e3(i)*g1(i)*fff2(i)*ff3(i)*(-1)**k
c       term2=gamm(c)*w**(a+dble(i))*e3(i)*fff1(i)*g2(i)*ff3(i)*(-1)**k
c       term3=gamm(c)*w**(a+dble(i))*e3(i)*fff1(i)*ff2(i)*g3(i)*(-1)**k
c       term4=gamm(c)*w**(a+dble(i))*e4(i)*fff1(i)*ff2(i)*g4(i)*(-1)**k
c       term5=gamm(c)*fff1(i)*e4(i)*ff2(i)*ff4(i)*g5(i)*(-1)**k
c       temp3=temp3+term1+term2+term3+term4+term5
c       temp4=temp4+gamm(c)*w**(a+dble(i))*e1(i)*fff2(i)*ff3(i)*(-1)**k
c       temp5=temp5+gamm(c)*w**(a+dble(i)-eps)*e2(i)*fff1(i)*ff4(i)
c    #                                                       *(-1)**k
c200  continue
c     write(10,*)'temp=',temp,temp3
c     write(10,*)'temp1=',temp1,temp4
c     write(10,*)'temp2=',temp2,temp5

c     x=temp3+et1*temp4+et2*temp5+coeff1*coeff2/gamm(x1)/gamm(x2)
c    #                             *gamm(c)*w**(a-eps-dble(k))*temp6
c     write(10,*)'f=',f,x

      return
      end

c**********************************************************************
c
c  subroutine name    - fix4a
c
c  computation
c  performed          - calculates the hypergeometric function for z
c                       in the interval (.5,1) when c-a-b is near a
c                       positive integer.
c
c  usage              - call fix4a(a,b,c,n,k,f,w,machep,eps)
c
c  arguments    a,b,c - parameters of the hypergeometric function.
c
c                  n  - the upper limit of the finite series expansion
c                       of the hypergeometric function.
c
c                  k  - equals the nearest integer of c-a-b.
c
c                  f  - computed value of the hypergeometric function.
c
c                  w  - transformed independent variable.
c
c              machep - equals machine epsilon.
c
c                eps  - equals c-a-b-k.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine fix4a(a,b,c,n,k,f,w,machep,eps)

      real*8  zero,one,two,four,eighth,seven,eight,sxteen
      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,eighth=1.d0/8.d0,
     #           seven=7.d0,eight=8.d0,sxteen=16.d0,nmax=100)
      real*8   a,b,c,w,f,gamm,eps,machep,test,arg,rn,sum,et1,et2,
     #         term1,term2,term3,term4,term5,term6,term7,term8,
     #         temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,
     #         coeff,coeff1,coeff2,coeff3,coeff4,x,x1,x2,x3,x4,
     #         t1(0:80),t2(0:80),t3(0:80),t4(0:80),c1(0:80),c2(0:80),
     #         c3(0:80),c4(0:80),f1(0:80),f2(0:80),f3(0:80),f4(0:80),
     #         g1(0:nmax),g2(0:nmax),g3(0:nmax),g4(0:nmax),g5(0:nmax),
     #         fff1(0:nmax),ff1(0:nmax),fff2(0:nmax),ff2(0:nmax),
     #         ff3(0:nmax),ff4(0:nmax),poch1(0:nmax),poch2(0:nmax),
     #         e1(0:nmax),e2(0:nmax),e3(0:nmax),e4(0:nmax)

      integer  flag

c  calculate the extra terms

      x=a+dble(k)-one
      sum=zero
      coeff=one
      flag=0
  1   if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 1
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
  2     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 2
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 1
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 3 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
  3     continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 4 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
  4     continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 5 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
  5     continue
      endif

      if (flag.eq.0) then
        x1=a+dble(k)
        coeff1=one
  6     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 6
        endif
        x2=a+dble(k)+eps
        coeff2=one
  7     if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 7
        endif
        temp=sum+coeff*temp
        et1=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
  8     if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 8
        endif
        x2=x+one+eps
        coeff2=one
  9     if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 9
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et1=sum+coeff*temp
      endif

c     write(10,*)et1,(one/gamm(a+k+eps)-one/gamm(a+k))/eps

      x=b+dble(k)-one
      sum=zero
      coeff=one
      flag=0
  10  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 10
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
  11    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 11
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 10
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c2,41,1)
        t2(0)=one
        t2(1)=two*(x+eps)-one
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 12 i=2,41
          t2(i)=(four*(x+eps)-two)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-two)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
  12    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c2,55,2)
        t2(0)=one
        t2(1)=two*(x+eps)
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 13 i=2,55
          t2(i)=four*(x+eps)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+four*x*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
  13    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c2,34,3)
        t2(0)=one
        t2(1)=two*(x+eps)-two
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 14 i=2,34
          t2(i)=(four*(x+eps)-four)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-four)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
  14    continue
      endif

      if (flag.eq.0) then
        x1=b+dble(k)
        coeff1=one
  15    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 15
        endif
        x2=b+dble(k)+eps
        coeff2=one
  16    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 16
        endif
        temp2=sum+coeff*temp2
        et2=-temp2*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
  17    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 17
        endif
        x2=x+one+eps
        coeff2=one
  18    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 18
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et2=sum+coeff*temp2
      endif

c     write(10,*)et2,(one/gamm(b+k+eps)-one/gamm(b+k))/eps

c  calculate the f-functions

      x1=a
      coeff1=one
  20  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 20
      endif

      x2=b
      coeff2=one
  21  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 21
      endif

      coeff=one
      arg=-eps-dble(k)
  22  if (arg.lt.-eps) then
        coeff=coeff/arg
        arg=arg+one
        go to 22
      endif

      ff1(0)=coeff1/gamm(x1)
      ff2(0)=coeff2/gamm(x2)
      fff1(0)=coeff1/gamm(x1)
      fff2(0)=coeff2/gamm(x2)
      ff3(0)=(-1)**(k+1)*gamm(one-eps)*coeff
      ff4(0)=gamm(one+eps)

c     do 23 i=1,n
c       ff1(i)=(a+dble(k+i-1))*ff1(i-1)
c       ff2(i)=(b+dble(k+i-1))*ff2(i-1)
c       fff1(i)=(c-b+dble(i-1))*fff1(i-1)
c       fff2(i)=(c-a+dble(i-1))*fff2(i-1)
c       ff3(i)=ff3(i-1)/(eps+dble(i)+dble(k))
c       ff4(i)=ff4(i-1)/(dble(i)-eps)
c 23  continue

c     do 24 i=0,n
c       write(10,*)'fff1=',fff1(i),gamm(c-b+i)/gamm(a)/gamm(c-b)
c       write(10,*)'ff1=',ff1(i),gamm(a+k+i)/gamm(a)/gamm(a+k)
c       write(10,*)'fff2=',fff2(i),gamm(c-a+i)/gamm(b)/gamm(c-a)
c       write(10,*)'ff2=',ff2(i),gamm(b+k+i)/gamm(b)/gamm(b+k)
c       write(10,*)'ff3=',ff3(i),(-1)**(k+i)*eps*gamm(-k-i-eps)
c       write(10,*)'ff4=',ff4(i),(-1)**i*eps*gamm(eps-i)
c 24  continue

c   calculate  g1,g2

      g1(0)=zero
      g2(0)=zero
      poch1(0)=coeff1/gamm(x1)
      poch2(0)=coeff2/gamm(x2)
      do 100 i=1,n
        poch1(i)=(a+dble(k+i-1))*poch1(i-1)
        poch2(i)=(b+dble(k+i-1))*poch2(i-1)
        g1(i)=g1(i-1)*(a+eps+dble(k+i-1))+poch1(i-1)
        g2(i)=g2(i-1)*(b+eps+dble(k+i-1))+poch2(i-1)
 100  continue

c     do 101 i=0,n
c       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
c       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
c101  continue

c  calculate  g3,g4,g5

      x3=zero
      call cheb(c3,55,2)
      t3(0)=one
      t3(1)=two*(x3-eps)
      f3(0)=zero
      f3(1)=-two
      g3(0)=c3(1)*f3(1)

      x4=zero
      call cheb(c4,55,2)
      t4(0)=one
      t4(1)=two*(x4+eps)
      f4(0)=zero
      f4(1)=two
      g4(0)=c4(1)*f4(1)

      do 105 i=2,55
        t3(i)=four*(x3-eps)*t3(i-1)-t3(i-2)
        t4(i)=four*(x4+eps)*t4(i-1)-t4(i-2)
        f3(i)=-four*t3(i-1)+four*x3*f3(i-1)-f3(i-2)
        f4(i)=four*t4(i-1)+four*x4*f4(i-1)-f4(i-2)
        g3(0)=g3(0)+c3(i)*f3(i)
        g4(0)=g4(0)+c4(i)*f4(i)
 105  continue

      g3(0)=-g3(0)
      do 106 i=-k,-1
        g3(0)=(g3(0)+one/gamm(dble(k+i+2)))/(dble(k+i+1)+eps)
 106  continue

      test=eps*dlog(w)
      temp=dlog(w)
         if (dabs(test).ge.eighth) then
           temp=(exp(test)-one)/eps
         else
           i=1
 107       rn=(eps**(i)*(dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 108
           temp=temp+rn
           i=i+1
           go to 107
         endif
 108     g5(0)=(w**k)*temp

c     write(10,*)g3(0),(-1)**k*gamm(-k-eps)+one/eps/gamm(dble(k+1))
c     write(10,*)g4(0),gamm(eps)-one/eps
c     write(10,*)g5(0),w**(k+eps)/eps-w**k/eps

      do 120 i=1,n
        g3(i)=(g3(i-1)+one/gamm(dble(k+i+1)))/(dble(k+i)+eps)
        g4(i)=(g4(i-1)+one/gamm(dble(i+1)))/(dble(i)-eps)
        g5(i)=w*g5(i-1)
 120  continue

      e1(0)=one
      e2(0)=one
      e3(0)=-one
      e4(0)=one
      do 130 i=1,k
        e1(0)=(a+dble(i-1))*e1(0)
        e2(0)=(b+dble(i-1))*e2(0)
        e3(0)=e3(0)/dble(i)
 130  continue

c     do 140 i=1,n
c       e1(i)=(a+dble(k+i-1))*e1(i-1)
c       e2(i)=(b+dble(k+i-1))*e2(i-1)
c       e3(i)=e3(i-1)/dble(k+i)
c       e4(i)=e4(i-1)/dble(i)
c140  continue

c  put everything back together again

      term1=gamm(c)*(-1)**k*fff2(0)*ff3(0)*e4(0)*w**(c-a-b)
      term2=gamm(c)*(-1)**k*ff1(0)*ff3(0)*e4(0)*w**(c-a-b)
      term3=gamm(c)*(-1)**k*ff1(0)*ff2(0)*e4(0)*w**(c-a-b)
      term4=-gamm(c)*(-1)**k*ff1(0)*ff2(0)*e3(0)*w**(c-a-b)
      term5=gamm(c)*(-1)**k*ff1(0)*ff2(0)*e3(0)*ff4(0)
      term6=-gamm(c)*(-w)**k*et1*e1(0)*ff2(0)*e3(0)*ff4(0)
      term7=-gamm(c)*(-w)**k*et2*ff1(0)*e2(0)*e3(0)*ff4(0)
      term8=-gamm(c)*(-w)**k*eps*et1*et2*e1(0)*e2(0)*e3(0)*ff4(0)

      temp=g1(0)*term1+g2(0)*term2+g3(0)*term3+g4(0)*term4
     #                                       +g5(0)*term5
      temp1=term6
      temp2=term7
      temp3=term8

      do 150 i=1,n
        term1=term1*w*(b+eps+dble(k+i-1))/(eps+dble(i+k))/dble(i)
        term2=term2*w*(a+dble(k+i-1))/(eps+dble(i+k))/dble(i)
        term3=term3*w*(a+dble(k+i-1))*(b+dble(k+i-1))/dble(i)
        term4=term4*w*(a+dble(k+i-1))*(b+dble(k+i-1))/dble(k+i)
        term5=term5*(a+dble(k+i-1))*(b+dble(k+i-1))/dble(k+i)
     #                                       /(dble(i)-eps)
        term6=term6*w*(a+dble(k+i-1))*(b+dble(k+i-1))/dble(k+i)
     #                                       /(dble(i)-eps)
        term7=term7*w*(a+dble(k+i-1))*(b+dble(k+i-1))/dble(k+i)
     #                                       /(dble(i)-eps)
        term8=term8*w*(a+dble(k+i-1))*(b+dble(k+i-1))/dble(k+i)
     #                                       /(dble(i)-eps)
        temp=temp+g1(i)*term1+g2(i)*term2+g3(i)*term3+g4(i)*term4
     #                                       +g5(i)*term5
        temp1=temp1+term6
        temp2=temp2+term7
        temp3=temp3+term8
 150  continue

c  calculate the finite series term

      poch1(0)=one
      poch2(0)=one
      do 160 i=1,k-1
        poch1(i)=(a+dble(i-1))*poch1(i-1)
        poch2(i)=(b+dble(i-1))*poch2(i-1)
 160  continue

      temp4=zero
      do 170 i=0,k-1
        temp4=temp4+poch1(i)*poch2(i)*gamm(eps+dble(k-i))*(-w)**i
     #                                      /gamm(dble(i+1))
 170  continue

      x1=c-a
      coeff1=one
 180  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 180
      endif

      x2=c-b
      coeff2=one
 190  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 190
      endif

      temp4=temp4*gamm(c)*coeff1*coeff2/gamm(x1)/gamm(x2)

      f=temp+temp1+temp2+temp3+temp4

c  alternative method  (must also turn on the individual functions)

c     temp5=zero
c     temp6=zero
c     temp7=zero
c     temp8=zero
c     do 200 i=0,n
c       term1=gamm(c)*(-1)**k*e4(i)*g1(i)*fff2(i)*ff3(i)
c    #                                     *w**(dble(k+i)+eps)
c       term2=gamm(c)*(-1)**k*e4(i)*g2(i)*ff1(i)*ff3(i)
c    #                                     *w**(dble(k+i)+eps)
c       term3=gamm(c)*(-1)**k*e4(i)*g3(i)*ff1(i)*ff2(i)
c    #                                     *w**(dble(k+i)+eps)
c       term4=-gamm(c)*(-1)**k*e3(i)*g4(i)*ff1(i)*ff2(i)
c    #                                     *w**(dble(k+i)+eps)
c       term5=gamm(c)*(-1)**k*e3(i)*g5(i)*ff1(i)*ff2(i)*ff4(i)

c       temp5=temp5+term1+term2+term3+term4+term5
c       temp6=temp6-gamm(c)*e3(i)*(-1)**k*et1*e1(i)*ff2(i)*ff4(i)
c    #                                            *w**(k+i)
c       temp7=temp7-gamm(c)*e3(i)*(-1)**k*et2*e2(i)*ff1(i)*ff4(i)
c    #                                            *w**(k+i)
c       temp8=temp8-gamm(c)*e3(i)*(-1)**k*eps*et1*et2*e1(i)*e2(i)*ff4(i)
c    #                                            *w**(k+i)
c200  continue
c     write(10,*)'temp=',temp,temp5
c     write(10,*)'temp1=',temp1,temp6
c     write(10,*)'temp2=',temp2,temp7
c     write(10,*)'temp3=',temp3,temp8

c     x=temp5+temp6+temp7+temp8+temp4
c     write(10,*)'f=',f,x

      return
      end

c**********************************************************************
c
c  subroutine name    - fix4b
c
c  computation
c  performed          - calculates the hypergeometric function for z
c                       in the interval (.5,1) when c-a-b is near a
c                       negative integer.
c
c  usage              - call fix4b(a,b,c,n,k,f,w,machep,eps)
c
c  arguments    a,b,c - parameters of the hypergeometric function.
c
c                  n  - the upper limit of the finite series expansion
c                       of the hypergeometric function.
c
c                  k  - equals the nearest integer of a+b-c.
c
c                  f  - computed value of the hypergeometric function.
c
c                  w  - transformed independent variable.
c
c              machep - equals machine epsilon.
c
c                eps  - equals c-a-b+k.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine fix4b(a,b,c,n,k,f,w,machep,eps)

      real*8  zero,one,two,four,eighth,seven,eight,sxteen
      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,eighth=1.d0/8.d0,
     #           seven=7.d0,eight=8.d0,sxteen=16.d0,nmax=100)
      real*8   a,b,c,w,f,gamm,eps,machep,test,arg,rn,sum,et1,et2,
     #         term1,term2,term3,term4,term5,term6,term7,term8,
     #         temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,
     #         coeff,coeff1,coeff2,coeff3,coeff4,x,x1,x2,x3,x4,
     #         t1(0:80),t2(0:80),t3(0:80),t4(0:80),c1(0:80),c2(0:80),
     #         c3(0:80),c4(0:80),f1(0:80),f2(0:80),f3(0:80),f4(0:80),
     #         g1(0:nmax),g2(0:nmax),g3(0:nmax),g4(0:nmax),g5(0:nmax),
     #         fff1(0:nmax),ff1(0:nmax),fff2(0:nmax),ff2(0:nmax),
     #         ff3(0:nmax),ff4(0:nmax),poch1(0:nmax),poch2(0:nmax),
     #         e1(0:nmax),e2(0:nmax),e3(0:nmax),e4(0:nmax)

      integer  flag

c  calculate the extra terms

      x=a-dble(k)-one
      sum=zero
      coeff=one
      flag=0
   1  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 1
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
   2    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 2
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 1
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 3 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
   3    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 4 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
   4    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 5 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
   5    continue
      endif

      if (flag.eq.0) then
        x1=a-dble(k)
        coeff1=one
   6    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 6
        endif
        x2=a-dble(k)+eps
        coeff2=one
   7    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 7
        endif
        temp=sum+coeff*temp
        et1=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
   8    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 8
        endif
        x2=x+one+eps
        coeff2=one
   9    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 9
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et1=sum+coeff*temp
      endif

c     write(10,*)et1,(one/gamm(a-k+eps)-one/gamm(a-k))/eps

      x=b-dble(k)-one
      sum=zero
      coeff=one
      flag=0
  10  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 10
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
  11    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 11
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 10
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c2,41,1)
        t2(0)=one
        t2(1)=two*(x+eps)-one
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 12 i=2,41
          t2(i)=(four*(x+eps)-two)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-two)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
  12    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c2,55,2)
        t2(0)=one
        t2(1)=two*(x+eps)
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 13 i=2,55
          t2(i)=four*(x+eps)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+four*x*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
  13    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c2,34,3)
        t2(0)=one
        t2(1)=two*(x+eps)-two
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 14 i=2,34
          t2(i)=(four*(x+eps)-four)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-four)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
  14    continue
      endif

      if (flag.eq.0) then
        x1=b-dble(k)
        coeff1=one
  15    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 15
        endif
        x2=b-dble(k)+eps
        coeff2=one
  16    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 16
        endif
        temp2=sum+coeff*temp2
        et2=-temp2*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
  17    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 17
        endif
        x2=x+one+eps
        coeff2=one
  18    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 18
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et2=sum+coeff*temp2
      endif

c     write(10,*)et2,(one/gamm(b-k+eps)-one/gamm(b-k))/eps

c  calculate the f-functions

      x1=a
      coeff1=one
  20  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 20
      endif

      x2=b
      coeff2=one
  21  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 21
      endif

      x3=a-dble(k)
      coeff3=one
  22  if (x3.lt.one) then
        coeff3=x3*coeff3
        x3=x3+one
        go to 22
      endif

      x4=b-dble(k)
      coeff4=one
  23  if (x4.lt.one) then
        coeff4=x4*coeff4
        x4=x4+one
        go to 23
      endif

      coeff=one
      arg=eps-dble(k)
  24  if (arg.lt.eps) then
        coeff=coeff/arg
        arg=arg+one
        go to 24
      endif

      fff1(0)=one
      fff2(0)=one
      ff1(0)=one
      ff2(0)=one
      do 25 i=1,k
        fff1(0)=(c-b+dble(i-1))*fff1(0)
        fff2(0)=(c-a+dble(i-1))*fff2(0)
  25  continue

      fff1(0)=fff1(0)*coeff1/gamm(x1)
      fff2(0)=fff2(0)*coeff2/gamm(x2)
      ff1(0)=ff1(0)*coeff3/gamm(x3)
      ff2(0)=ff2(0)*coeff4/gamm(x4)
      ff3(0)=-gamm(one-eps)
      ff4(0)=(-1)**k*gamm(one+eps)*coeff

c     do 26 i=1,n
c       fff1(i)=(c-b+dble(k+i-1))*fff1(i-1)
c       fff2(i)=(c-a+dble(k+i-1))*fff2(i-1)
c       ff1(i)=(a+dble(i-1))*ff1(i-1)
c       ff2(i)=(b+dble(i-1))*ff2(i-1)
c       ff3(i)=ff3(i-1)/(eps+dble(i))
c       ff4(i)=ff4(i-1)/(dble(k+i)-eps)
c 26  continue

c     do 27 i=0,n
c       write(10,*)'fff1=',fff1(i),gamm(a+eps+i)/gamm(a)/gamm(c-b)
c       write(10,*)'fff2=',fff2(i),gamm(b+eps+i)/gamm(b)/gamm(c-a)
c       write(10,*)'ff1=',ff1(i),gamm(a+i)/gamm(a)/gamm(a-k)
c       write(10,*)'ff2=',ff2(i),gamm(b+i)/gamm(b)/gamm(b-k)
c       write(10,*)'ff3=',ff3(i),(-1)**i*eps*gamm(-i-eps)
c       write(10,*)'ff4=',ff4(i),(-1)**(k+i)*eps*gamm(eps-k-i)
c 27  continue

c   calculate  g1,g2

      x1=a
      coeff1=one
 100  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 100
      endif

      x2=b
      coeff2=one
 101  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 101
      endif

      g1(0)=zero
      g2(0)=zero
      poch1(0)=one
      poch2(0)=one
      do 102 i=1,k
        g1(0)=g1(0)*(a+eps+dble(i-k-1))+poch1(0)
        g2(0)=g2(0)*(b+eps+dble(i-k-1))+poch2(0)
        poch1(0)=poch1(0)*(a+dble(i-k-1))
        poch2(0)=poch2(0)*(b+dble(i-k-1))
 102  continue

      g1(0)=g1(0)*coeff1/gamm(x1)
      g2(0)=g2(0)*coeff2/gamm(x2)
      poch1(0)=poch1(0)*coeff1/gamm(x1)
      poch2(0)=poch2(0)*coeff2/gamm(x2)
      do 103 i=1,n
        poch1(i)=(a+i-1)*poch1(i-1)/i
        poch2(i)=(b+i-1)*poch2(i-1)/i
        g1(i)=(g1(i-1)*(a+eps+i-1)+poch1(i-1))/i
        g2(i)=(g2(i-1)*(b+eps+i-1)+poch2(i-1))/i
 103  continue

c     do 104 i=0,n
c       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps/gamma(i+1.0)
c       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps/gamma(i+1.0)
c104  continue

c  calculate  g3,g4,g5

      x3=zero
      call cheb(c3,55,2)
      t3(0)=one
      t3(1)=two*(x3-eps)
      f3(0)=zero
      f3(1)=-two
      g3(0)=c3(1)*f3(1)

      x4=zero
      call cheb(c4,55,2)
      t4(0)=one
      t4(1)=two*(x4+eps)
      f4(0)=zero
      f4(1)=two
      g4(0)=c4(1)*f4(1)

      do 105 i=2,55
        t3(i)=four*(x3-eps)*t3(i-1)-t3(i-2)
        t4(i)=four*(x4+eps)*t4(i-1)-t4(i-2)
        f3(i)=-four*t3(i-1)+four*x3*f3(i-1)-f3(i-2)
        f4(i)=four*t4(i-1)+four*x4*f4(i-1)-f4(i-2)
        g3(0)=g3(0)+c3(i)*f3(i)
        g4(0)=g4(0)+c4(i)*f4(i)
 105  continue

      g3(0)=-g3(0)
      do 106 i=-k,-1
        g4(0)=(g4(0)+one/gamm(dble(k+i+2)))/(dble(k+i+1)-eps)
 106  continue

      test=eps*dlog(w)
      temp=dlog(w)
         if (dabs(test).ge.eighth) then
           temp=(exp(test)-one)/eps
         else
           i=1
 107       rn=(eps**(i)*(dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 108
           temp=temp+rn
           i=i+1
           go to 107
         endif
 108     g5(0)=temp

c     write(10,*)g3(0),gamm(-eps)+one/eps
c     write(10,*)g4(0),(-1)**k*gamm(eps-dble(k))-one/eps/gamm(dble(k+1))
c     write(10,*)g5(0),w**eps/eps-one/eps

      do 120 i=1,n
        temp=one/gamm(dble(k+1))
        do 121 j=1,i
          temp=temp*dble(j)/dble(k+j)
 121    continue
        g3(i)=(g3(i-1)*dble(i)+one)/(dble(i)+eps)
        g4(i)=(g4(i-1)*dble(i)+temp)/(dble(k+i)-eps)
        g5(i)=w*g5(i-1)
 120  continue

      e1(0)=one
      e2(0)=one
      e3(0)=-one
      e4(0)=one
      do 130 i=1,k
        e4(0)=e4(0)/dble(i)
 130  continue

c     do 140 i=1,n
c       e1(i)=(a+dble(i-1))*e1(i-1)
c       e2(i)=(b+dble(i-1))*e2(i-1)
c       e3(i)=e3(i-1)/dble(i)
c       e4(i)=e4(i-1)/dble(k+i)
c140  continue

c  put everything back together again

      term1=gamm(c)*(-1)**k*fff2(0)*ff3(0)*e4(0)*w**eps
      term2=gamm(c)*(-1)**k*ff1(0)*ff3(0)*e4(0)*w**eps
      term3=gamm(c)*(-1)**k*ff1(0)*ff2(0)*e4(0)*w**eps
      term4=-gamm(c)*(-1)**k*ff1(0)*ff2(0)*e3(0)*w**eps
      term5=gamm(c)*(-1)**k*ff1(0)*ff2(0)*e3(0)*ff4(0)
      term6=-gamm(c)*(-1)**k*et1*e1(0)*ff2(0)*e3(0)*ff4(0)
      term7=-gamm(c)*(-1)**k*et2*ff1(0)*e2(0)*e3(0)*ff4(0)
      term8=-gamm(c)*(-1)**k*eps*et1*et2*e1(0)*e2(0)*e3(0)*ff4(0)

      temp=g1(0)*term1+g2(0)*term2+g3(0)*term3+g4(0)*term4
     #                                       +g5(0)*term5
      temp1=term6
      temp2=term7
      temp3=term8

      do 150 i=1,n
        term1=term1*w*(b+eps+dble(i-1))/(eps+dble(i))*dble(i)/dble(k+i)
        term2=term2*w*(a+dble(i-1))/(eps+dble(i))*dble(i)/dble(k+i)
        term3=term3*w*(a+dble(i-1))/dble(i)*(b+dble(i-1))/dble(k+i)
        term4=term4*w*(a+dble(i-1))/dble(i)*(b+dble(i-1))/dble(i)
        term5=term5*(a+dble(i-1))/dble(i)*(b+dble(i-1))/(dble(k+i)-eps)
        term6=term6*w*(a+dble(i-1))/dble(i)*(b+dble(i-1))
     #                                                 /(dble(k+i)-eps)
        term7=term7*w*(a+dble(i-1))/dble(i)*(b+dble(i-1))
     #                                                 /(dble(k+i)-eps)
        term8=term8*w*(a+dble(i-1))/dble(i)*(b+dble(i-1))
     #                                                 /(dble(k+i)-eps)
        temp=temp+g1(i)*term1+g2(i)*term2
     #        +g3(i)*term3+g4(i)*term4+g5(i)*term5
        temp1=temp1+term6
        temp2=temp2+term7
        temp3=temp3+term8
 150  continue

c  calculate the finite series term

      poch1(0)=one
      poch2(0)=one
      do 160 i=1,k-1
        poch1(i)=(c-a+dble(i-1))*poch1(i-1)
        poch2(i)=(c-b+dble(i-1))*poch2(i-1)
 160  continue

      temp4=zero
      do 170 i=0,k-1
        temp4=temp4+poch1(i)*poch2(i)*gamm(-eps+dble(k-i))*(-1)**i
     #                  *w**(eps+dble(i-k))/gamm(dble(i+1))
 170  continue

      x1=a
      coeff1=one
 180  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 180
      endif

      x2=b
      coeff2=one
 190  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 190
      endif

      temp4=temp4*gamm(c)*coeff1*coeff2/gamm(x1)/gamm(x2)

      f=temp+temp1+temp2+temp3+temp4

c  alternative method (must also turn on the individual functions)

c     temp5=zero
c     temp6=zero
c     temp7=zero
c     temp8=zero
c     do 200 i=0,n
c       term1=w**(dble(i)+eps)/gamm(dble(k+i+1))
c    #                      *g1(i)*fff2(i)*ff3(i)
c       term2=w**(dble(i)+eps)/gamm(dble(k+i+1))
c    #                      *g2(i)*ff1(i)*ff3(i)
c       term3=w**(dble(i)+eps)/gamm(dble(k+i+1))
c    #                      *g3(i)*ff1(i)*ff2(i)
c       term4=w**(dble(i)+eps)/gamm(dble(i+1))
c    #                      *g4(i)*ff1(i)*ff2(i)
c       term5=-ff1(i)/gamm(dble(i+1))*ff2(i)*ff4(i)*g5(i)
c
c       temp5=temp5+term1+term2+term3+term4+term5
c
c       temp6=temp6+ff1(i)*et2/gamm(dble(i+1))*ff4(i)*e2(i)
c    #                                 *w**(dble(i))
c       temp7=temp7+ff2(i)*et1/gamm(dble(i+1))*ff4(i)*e1(i)
c    #                                 *w**(dble(i))
c       temp8=temp8+e1(i)*et1*et2*eps/gamm(dble(i+1))*ff4(i)
c    #                                 *w**(dble(i))*e2(i)
c200  continue
c     write(10,*)'temp=',temp,temp5*gamm(c)*(-1)**k
c     write(10,*)'temp1=',temp1,temp7*gamm(c)*(-1)**k
c     write(10,*)'temp2=',temp2,temp6*gamm(c)*(-1)**k
c     write(10,*)'temp3=',temp3,temp8*gamm(c)*(-1)**k
c
c     x=(-1)**k*gamm(c)*(temp5+temp6+temp7+temp8)+temp4
c
c     write(10,*)'f=',f,x

      return
      end

c**********************************************************************
c
c  subroutine name    - fix5a
c
c  computation
c  performed          - calculates the hypergeometric function for z
c                       in the interval (1,2) when c-a-b is near a
c                       positive integer.
c
c  usage              - call fix5a(a,b,c,n,k,re,im,w,machep,eps,pi)
c
c  arguments    a,b,c - parameters of the hypergeometric function.
c
c                  n  - the upper limit of the finite series expansion
c                       of the hypergeometric function.
c
c                  k  - equals the nearest integer of c-a-b.
c
c               re,im - computed values for the real and imaginary parts
c                       of the hypergeometric function.
c
c                  w  - transformed independent variable.
c
c              machep - equals machine epsilon.
c
c                eps  - equals c-a-b-k.
c
c                 pi  - equals 3.1415... to machine accuracy.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine fix5a(a,b,c,n,k,re,im,w,machep,eps,pi)

      real*8  zero,one,two,four,eighth,seven,eight,sxteen
      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,eighth=1.d0/8.d0,
     #           seven=7.d0,eight=8.d0,sxteen=16.d0,nmax=100)
      real*8   a,b,c,w,re,im,gamm,temp,temp2,g1(0:nmax),g2,
     #         g3(0:nmax),g4(0:nmax),g5(0:nmax),x,x1,x2,x3,x4,psi,rn,
     #         t1(0:80),t2(0:80),t3(0:80),t4(0:80),test,machep,pi,
     #         f1(0:80),f2(0:80),f3(0:80),f4(0:80),ff3(0:nmax),eps,
     #         ff4(0:nmax),coeff1,coeff2,c1(0:80),c2(0:80),c3(0:80),
     #         c4(0:80),sum,term1,term2,term3,term4,term5,poch1(0:nmax),
     #         coeff,temp1,et1,et2,e1(0:nmax),e2(0:nmax),e3(0:nmax),
     #         ff1(0:nmax),fff1(0:nmax),coeff3,coeff4,f(0:nmax),error,
     #         poch2(0:nmax)

      x3=zero
      call cheb(c3,55,2)
      t3(0)=one
      t3(1)=two*(x3+eps)
      f3(0)=zero
      f3(1)=two
      g3(0)=c3(1)*f3(1)

      x4=zero
      call cheb(c4,55,2)
      t4(0)=one
      t4(1)=two*(x4-eps)
      f4(0)=zero
      f4(1)=-two
      g4(0)=c4(1)*f4(1)

      do 7 i=2,55
        t3(i)=four*(x3+eps)*t3(i-1)-t3(i-2)
        t4(i)=four*(x4-eps)*t4(i-1)-t4(i-2)
        f3(i)=four*t3(i-1)+four*x3*f3(i-1)-f3(i-2)
        f4(i)=-four*t4(i-1)+four*x4*f4(i-1)-f4(i-2)
        g3(0)=g3(0)+c3(i)*f3(i)
        g4(0)=g4(0)+c4(i)*f4(i)
  7   continue

      g4(0)=-g4(0)
      do 10 i=-k,-1
        g4(0)=(g4(0)+one/gamm(dble(k+i+2)))/(dble(k+i+1)+eps)
  10  continue

      test=eps*dlog(w)
      temp=dlog(w)
         if (dabs(test).ge.eighth) then
           temp=(exp(test)-one)/eps
         else
           i=1
  20       rn=(eps**(i)*(dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 30
           temp=temp+rn
           i=i+1
           go to 20
         endif
  30     g5(0)=temp*w**k

c     write(10,*)g3(0),gamm(-eps)+one/eps
c     write(10,*)g4(0),(-1)**k*gamm(eps-dble(k))-one/eps/gamm(dble(k+1))
c     write(10,*)g5(0),w**eps/eps-one/eps

      do 60 i=1,n
        g3(i)=(g3(i-1)+one/gamm(dble(i+1)))/(dble(i)-eps)
        g4(i)=(g4(i-1)+one/gamm(dble(k+i+1)))/(dble(k+i)+eps)
        g5(i)=w*g5(i-1)
  60  continue

      do 65 i=0,n
        ff3(i)=eps*g3(i)+one/gamm(dble(i+1))
        ff4(i)=eps*g4(i)-one/gamm(dble(k+i+1))
  65  continue

c  calculate the extra terms

      x=a-one
      sum=zero
      coeff=one
      flag=0
  61  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 61
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 610    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 610
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 61
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 611 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 611    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 612 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 612    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 613 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 613    continue
      endif

      if (flag.eq.0) then
        x1=a
        coeff1=one
 614    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 614
        endif
        x2=a+eps
        coeff2=one
 615    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 615
        endif
        temp=sum+coeff*temp
        et1=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 616    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 616
        endif
        x2=x+one+eps
        coeff2=one
 617    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 617
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et1=sum+coeff*temp
      endif
      et1=-et1
c     write(10,*)et1,(one/gamm(c-b-dble(k)-eps)-one/gamm(c-b-dble(k)))
c    #                                                  /eps

      x=b-one
      sum=zero
      coeff=one
      flag=0
  62  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 62
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 620    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 620
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 62
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c2,41,1)
        t2(0)=one
        t2(1)=two*(x+eps)-one
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 621 i=2,41
          t2(i)=(four*(x+eps)-two)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-two)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
 621    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c2,55,2)
        t2(0)=one
        t2(1)=two*(x+eps)
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 622 i=2,55
          t2(i)=four*(x+eps)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+four*x*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
 622    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c2,34,3)
        t2(0)=one
        t2(1)=two*(x+eps)-two
        f2(0)=zero
        f2(1)=two
        temp2=c2(1)*f2(1)
        do 623 i=2,34
          t2(i)=(four*(x+eps)-four)*t2(i-1)-t2(i-2)
          f2(i)=four*t2(i-1)+(four*x-four)*f2(i-1)-f2(i-2)
          temp2=temp2+c2(i)*f2(i)
 623    continue
      endif

      if (flag.eq.0) then
        x1=b
        coeff1=one
 624    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 624
        endif
        x2=b+eps
        coeff2=one
 625    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 625
        endif
        temp2=sum+coeff*temp2
        et2=-temp2*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 626    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 626
        endif
        x2=x+one+eps
        coeff2=one
 627    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 627
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et2=sum+coeff*temp2
      endif

c     write(10,*)et2,(one/gamm(b+eps)-one/gamm(b))/eps

      fff1(0)=one
      do 685 i=1,k
        fff1(0)=(a+dble(i-1))*fff1(0)
 685  continue

      ff1(0)=one
      do 686 i=1,n
        fff1(i)=(a+dble(k+i-1))*fff1(i-1)
        ff1(i)=(c-b+dble(i-1))*ff1(i-1)
 686  continue

      x1=c-b
      coeff1=one
 687  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 687
      endif

      x2=c-b-dble(k)
      coeff2=one
 688  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 688
      endif

      do 691 i=0,n
        x3=b+eps-dble(i)
        coeff3=one
 689    if (x3.lt.one) then
          coeff3=x3*coeff3
          x3=x3+one
          go to 689
        endif

        x4=b-dble(i)
        coeff4=one
 690    if (x4.lt.one) then
          coeff4=x4*coeff4
          x4=x4+one
          go to 690
        endif
        f(i)=ff1(i)*coeff4/gamm(x4)
        fff1(i)=fff1(i)*coeff1*coeff3/gamm(x1)/gamm(x3)
        ff1(i)=ff1(i)*coeff2*coeff4/gamm(x2)/gamm(x4)
c       write(10,*)'fff1=',fff1(i),gamm(c-b-eps+dble(i))
c    #            /gamm(a)/gamm(b+eps-dble(i))/gamm(c-b)
c       write(10,*)'ff1=',ff1(i),gamm(c-b+dble(i))
c    #            /gamm(a+eps)/gamm(b-dble(i))/gamm(c-b)
 691  continue

c   calculate  g1

      e1(0)=zero
      poch1(0)=one
      do 697 i=1,k
        e1(0)=e1(0)*(c-b-eps+dble(i-k-1))-poch1(0)
        poch1(0)=poch1(0)*(c-b+dble(i-k-1))
 697  continue
      do 698 i=1,n
        poch1(i)=(c-b+dble(i-1))*poch1(i-1)
        e1(i)=e1(i-1)*(c-b-eps+dble(i-1))-poch1(i-1)
 698  continue

      do 700 i=0,n
        x1=b-dble(i)
        coeff1=one
 699    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 699
        endif
        e1(i)=e1(i)*coeff1/gamm(x1)
 700  continue

      e2(0)=et2
      do 702 i=1,n
        x1=b-dble(i-1)
        coeff1=one
 701    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 701
        endif
        e2(i)=e2(i-1)*(b+eps-dble(i))+coeff1/gamm(x1)
 702  continue

      e3(0)=one
      do 703 i=1,k
        e3(0)=(a+dble(i-1))*e3(0)
 703  continue

      do 704 i=1,n
        e3(i)=(a+dble(k+i-1))*e3(i-1)
 704  continue

      x1=c-b
      coeff1=one
 705  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 705
      endif

      do 706 i=0,n
        g1(i)=(e2(i)*e3(i)+e1(i))*coeff1/gamm(x1)
c       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
 706  continue

c  calculate g2

      g2=zero
      if (dabs(eps).lt..1d0) then
        i=1
 707    rn=(-1)**i*pi**(i+i)*eps**(i+i-1)/gamm(dble(i+i+1))
        if (dabs(rn).lt.machep) go to 708
        g2=g2+rn
        i=i+1
        go to 707
      else
        g2=(cos(pi*eps)-one)/eps
      endif
 708  continue
c     write(10,*)'g2=',g2,(cos(pi*eps)-one)/eps

      temp=zero
      temp1=zero
      do 70 i=0,n
        term1=-g1(i)*cos(pi*eps)/gamm(dble(i+1))
     #                      *ff4(i)*w**(eps+dble(i+k))
        term2=fff1(i)*g2/gamm(dble(i+1))
     #                      *ff4(i)*w**(eps+dble(i+k))
        term3=-fff1(i)*g3(i)*ff4(i)*w**(eps+dble(i+k))
        term4=fff1(i)*ff3(i)*g4(i)*w**(eps+dble(i+k))
        term5=-fff1(i)*ff3(i)/gamm(dble(k+i+1))*g5(i)
        temp=temp+(term1+term2+term3+term4+term5)*(-1)**i
        temp1=temp1+(et1*f(i)*cos(pi*eps)/gamm(dble(i+1))*ff4(i)
     #                              *w**(dble(i+k)+eps))*(-1)**i
  70  continue

      poch1(0)=one
      poch2(0)=one
      do 71 i=1,k-1
        poch1(i)=(a+dble(i-1))*poch1(i-1)
        poch2(i)=(one-c+a+dble(i-1))*poch2(i-1)
  71  continue

      x1=c-a
      coeff1=one
  72  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 72
      endif

      x2=c-b
      coeff2=one
  73  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 73
      endif

      temp2=zero
      do 80 i=0,k-1
        temp2=temp2+poch1(i)*poch2(i)*coeff1*coeff2/gamm(x1)/gamm(x2)
     #        *gamm(dble(k-i)+eps)/gamm(dble(i+1))*(-w)**i
  80  continue

c     term1=zero
c     do 81 i=0,k-1
c       term1=term1+gamm(a+dble(i))/gamm(a)*gamm(a-c+dble(1+i))
c    #        /gamm(a-c+one)*gamm(eps+dble(k-i))*(-w)**i/gamm(dble(i+1))
c    #        /gamm(c-a)/gamm(c-b)
c 81  continue
c     write(10,*)temp2,term1

      re=(one-w)**a*gamm(c)*(temp+temp1+temp2)

c  calculate the imaginary part

      x1=a
      coeff1=one
  90  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 90
      endif

      temp=zero
      do 91 i=0,n
        temp=temp+(-1)**i*f(i)/gamm(dble(i+1))*ff4(i)*w**(eps+dble(i+k))
     #       *coeff1/gamm(x1)
  91  continue

      if (dabs(eps).lt..1d0) then
        temp1=one
        i=1
  92    temp2=temp1+(-1)**i*(pi*eps)**(i+i)/gamm(dble(i+i+2))
        error=(temp2-temp1)/temp2
        if (dabs(error).lt.machep) go to 93
        i=i+1
        temp1=temp2
        go to 92
      else
        temp2=sin(pi*eps)/pi/eps
      endif
  93  continue
c     write(10,*)temp2,sin(pi*eps)/pi/eps

      im=-pi*temp2*temp

      return
      end

c**********************************************************************
c
c  subroutine name    - fix5b
c
c  computation
c  performed          - calculates the hypergeometric function for z
c                       in the interval (1,2) when c-a-b is near a
c                       negative integer.
c
c  usage              - call fix5b(a,b,c,n,k,re,im,w,machep,eps,pi)
c
c  arguments    a,b,c - parameters of the hypergeometric function.
c
c                  n  - the upper limit of the finite series expansion
c                       of the hypergeometric function.
c
c                  k  - equals the nearest integer of a+b-c.
c
c               re,im - computed values for the real and imaginary parts
c                       of the hypergeometric function.
c
c                  w  - transformed independent variable.
c
c              machep - equals machine epsilon.
c
c                eps  - equals c-a-b+k.
c
c                 pi  - equals 3.1415... to machine accuracy.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine fix5b(a,b,c,n,k,re,im,w,machep,eps,pi)

      real*8  zero,one,two,four,eighth,seven,eight,sxteen
      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,eighth=1.d0/8.d0,
     #           seven=7.d0,eight=8.d0,sxteen=16.d0,nmax=100)
      real*8   a,b,c,w,re,im,gamm,temp,temp2,g1(0:nmax),g2(0:nmax),
     #         g3(0:nmax),g4(0:nmax),g5(0:nmax),x,x1,x2,x3,x4,psi,rn,
     #         t1(0:80),t2(0:80),t3(0:80),t4(0:80),test,machep,pi,
     #         f1(0:80),f2(0:80),f3(0:80),f4(0:80),ff3(0:nmax),eps,
     #         ff4(0:nmax),coeff1,coeff2,c1(0:80),c2(0:80),c3(0:80),
     #         c4(0:80),sum,term1,term2,term3,term4,term5,term6,
     #         coeff,temp1,et1,et2,e1,e2(0:nmax),coeff3,coeff4,
     #         fff1(0:nmax),fff2(0:nmax),ff1(0:nmax),ff2(0:nmax),
     #         poch1(0:nmax),poch2(0:nmax),ttest,error

      integer  flag

      x3=zero
      call cheb(c3,55,2)
      t3(0)=one
      t3(1)=two*(x3-eps)
      f3(0)=zero
      f3(1)=-two
      g3(0)=c3(1)*f3(1)

      x4=zero
      call cheb(c4,55,2)
      t4(0)=one
      t4(1)=two*(x4+eps)
      f4(0)=zero
      f4(1)=two
      g4(0)=c4(1)*f4(1)

      do 7 i=2,55
        t3(i)=four*(x3-eps)*t3(i-1)-t3(i-2)
        t4(i)=four*(x4+eps)*t4(i-1)-t4(i-2)
        f3(i)=-four*t3(i-1)+four*x3*f3(i-1)-f3(i-2)
        f4(i)=four*t4(i-1)+four*x4*f4(i-1)-f4(i-2)
        g3(0)=g3(0)+c3(i)*f3(i)
        g4(0)=g4(0)+c4(i)*f4(i)
  7   continue

      g3(0)=-g3(0)
      do 10 i=-k,-1
        g4(0)=(g4(0)+one/gamm(dble(k+i+2)))/(dble(k+i+1)-eps)
  10  continue

      test=eps*dlog(w)
      temp=dlog(w)
         if (dabs(test).ge.eighth) then
           temp=(exp(test)-one)/eps
         else
           i=1
  20       rn=(eps**(i)*(dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 30
           temp=temp+rn
           i=i+1
           go to 20
         endif
  30     g5(0)=temp

c     write(10,*)g3(0),gamm(-eps)+one/eps
c     write(10,*)g4(0),(-1)**k*gamm(eps-dble(k))-one/eps/gamm(dble(k+1))
c     write(10,*)g5(0),w**eps/eps-one/eps

      do 60 i=1,n
        g3(i)=(g3(i-1)+one/gamm(dble(i+1)))/(dble(i)+eps)
        g4(i)=(g4(i-1)+one/gamm(dble(k+i+1)))/(dble(k+i)-eps)
        g5(i)=w*g5(i-1)
  60  continue

      do 65 i=0,n
        ff3(i)=eps*g3(i)-one/gamm(dble(i+1))
        ff4(i)=eps*g4(i)+one/gamm(dble(k+i+1))
  65  continue

c  calculate the extra terms

      x=a-dble(k)-one
      sum=zero
      coeff=one
      flag=0
  61  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 61
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 610    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 610
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 61
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 611 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 611    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 612 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 612    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 613 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 613    continue
      endif

      if (flag.eq.0) then
        x1=a-dble(k)
        coeff1=one
 614    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 614
        endif
        x2=a-dble(k)+eps
        coeff2=one
 615    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 615
        endif
        temp=sum+coeff*temp
        et1=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 616    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 616
        endif
        x2=x+one+eps
        coeff2=one
 617    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 617
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et1=sum+coeff*temp
      endif

c     write(10,*)et1,(one/gamm(a-dble(k)+eps)-one/gamm(a-dble(k)))
c    #                                                  /eps

      x=b-dble(k)-one
      sum=zero
      coeff=one
      flag=0
  62  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 62
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 620    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 620
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 62
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 621 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 621    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 622 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 622    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 623 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 623    continue
      endif

      if (flag.eq.0) then
        x1=b-dble(k)
        coeff1=one
 624    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 624
        endif
        x2=b-dble(k)+eps
        coeff2=one
 625    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 625
        endif
        temp=sum+coeff*temp
        et2=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 626    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 626
        endif
        x2=x+one+eps
        coeff2=one
 627    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 627
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et2=sum+coeff*temp
      endif

c     write(10,*)et2,(one/gamm(b-dble(k)+eps)-one/gamm(b-dble(k)))
c    #                                                  /eps

      fff1(0)=one
      do 685 i=1,k
        fff1(0)=(c-b+dble(i-1))*fff1(0)
 685  continue

      ff1(0)=one
      e2(0)=one
      do 686 i=1,n
        fff1(i)=(c-b+dble(k+i-1))*fff1(i-1)
        ff1(i)=(a+dble(i-1))*ff1(i-1)
        e2(i)=ff1(i)
 686  continue

      x1=a
      coeff1=one
 687  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 687
      endif

      x2=a-dble(k)
      coeff2=one
 688  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 688
      endif

      do 691 i=0,n
        x3=b+eps-dble(i+k)
        coeff3=one
 689    if (x3.lt.one) then
          coeff3=x3*coeff3
          x3=x3+one
          go to 689
        endif

        x4=b-dble(i+k)
        coeff4=one
 690    if (x4.lt.one) then
          coeff4=x4*coeff4
          x4=x4+one
          go to 690
        endif
        fff1(i)=fff1(i)*coeff1/gamm(x1)
        ff1(i)=ff1(i)*coeff2/gamm(x2)
        fff2(i)=coeff3/gamm(x3)
        ff2(i)=coeff4/gamm(x4)
c       write(10,*)'fff1=',fff1(i),gamm(c-b+dble(i+k))/gamm(a)/gamm(c-b)
c       write(10,*)'ff1=',ff1(i),gamm(a+dble(i))/gamm(a)/gamm(a-dble(k))
c       write(10,*)'fff2=',fff2(i),one/gamm(b+eps-dble(k+i))
c       write(10,*)'ff2=',ff2(i),one/gamm(b-dble(k+i))
 691  continue

c   calculate  g1

      g1(0)=zero
      poch1(0)=one
      do 697 i=1,k
        g1(0)=g1(0)*(a+eps+dble(i-k-1))+poch1(0)
        poch1(0)=poch1(0)*(a+dble(i-k-1))
 697  continue
      do 698 i=1,n
        poch1(i)=(a+dble(i-1))*poch1(i-1)
        g1(i)=g1(i-1)*(a+eps+dble(i-1))+poch1(i-1)
 698  continue

      x1=a
      coeff1=one
 699  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 699
      endif
      do 700 i=0,n
        g1(i)=g1(i)*coeff1/gamm(x1)
c       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
 700  continue

c   calculate  g2

      g2(0)=et2
      do 702 i=1,n
        x1=b-dble(k+i-1)
        coeff1=one
 701    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 701
        endif
        g2(i)=g2(i-1)*(b+eps-dble(i+k))+coeff1/gamm(x1)
c       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
 702  continue

c  calculate  e1

      e1=zero
      if (dabs(eps).lt..1d0) then
        i=1
 703    rn=(-1)**i*pi**(i+i)*eps**(i+i-1)/gamm(dble(i+i+1))
        if (dabs(rn).lt.machep) go to 704
        e1=e1+rn
        i=i+1
        go to 703
      else
        e1=(cos(pi*eps)-one)/eps
      endif
 704  continue
c     write(10,*)'e1=',e1,(cos(pi*eps)-one)/eps

c  put everything back together again

      ttest=zero
      temp=zero
      temp1=zero
      do 70 i=0,n
        term1=g1(i)/gamm(dble(k+i+1))*ff2(i)*ff3(i)
     #                   *cos(pi*eps)*w**(dble(i)+eps)
        term2=-ff1(i)/gamm(dble(k+i+1))*g2(i)*ff3(i)
     #                   *cos(pi*eps)*w**(dble(i)+eps)
        term3=ff1(i)/gamm(dble(k+i+1))*fff2(i)*g3(i)
     #                   *cos(pi*eps)*w**(dble(i)+eps)
        term4=ff1(i)/gamm(dble(i+1))*fff2(i)*g4(i)
     #                   *cos(pi*eps)*w**(dble(i)+eps)
        term5=-ff1(i)/gamm(dble(i+1))*ff4(i)*fff2(i)
     #                   *cos(pi*eps)*g5(i)
        term6=-ff1(i)/gamm(dble(i+1))*ff4(i)*fff2(i)
     #                   *w**(dble(i))*e1
        temp=temp+(term1+term2+term3+term4+term5+term6)*(-1)**(k+i)
        temp1=temp1+e2(i)/gamm(dble(i+1))*et1*fff2(i)
     #                       *ff4(i)*w**(dble(i))*(-1)**(k+i)
c       ttest=ttest+(-1)**(k+i)*(cos(pi*eps)*fff1(i)*ff2(i)*ff3(i)
c    #    /gamm(dble(k+i+1))*w**(dble(i)+eps)+ff1(i)*fff2(i)
c    #    /gamm(dble(i+1))*ff4(i)*w**(dble(i)))/eps
c       write(10,*)temp,ttest
c       ttest=ttest+(-1)**(k+i)*gamm(a+dble(i))/gamm(a)*fff2(i)
c    #    /gamm(dble(i+1))*ff4(i)*w**(dble(i))*et1
c       write(10,*)temp1,ttest
c       ttest=ttest+gamm(dble(i+k+1)-b)/gamm(one-b)*gamm(c-b+dble(i+k))
c    #    /gamm(c-b)/gamm(a)/gamm(b)*gamm(-eps-dble(i))*(-1)**(i+k)
c    #    *w**(eps+dble(i))/gamm(dble(i+k+1))*cos(pi*(eps-dble(k)))
c    #    +gamm(a+dble(i))/gamm(a)*gamm(a-c+dble(i+1))/gamm(a-c+one)
c    #    /gamm(c-a)/gamm(c-b)*gamm(eps-dble(k+i))*(-w)**i
c    #    /gamm(dble(i+1))
c       write(10,*)temp+temp1,ttest
  70  continue

      poch1(0)=one
      poch2(0)=one
      do 71 i=1,k-1
        poch1(i)=(c-b+dble(i-1))*poch1(i-1)
        poch2(i)=(dble(i)-b)*poch2(i-1)
  71  continue

      x1=a
      coeff1=one
  72  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 72
      endif

      x2=b
      coeff2=one
  73  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 73
      endif

      temp2=zero
      do 80 i=0,k-1
        temp2=temp2+coeff1*coeff2/gamm(x1)/gamm(x2)*poch1(i)*poch2(i)
     #  *gamm(dble(k-i)-eps)/gamm(dble(i+1))*w**(eps+dble(i-k))*(-1)**i
  80  continue

c     term1=zero
c     do 81 i=0,k-1
c       term1=term1+gamm(dble(i+1)-b)/gamm(a)*gamm(c-b+dble(i))/gamm(b)
c    #        /gamm(one-b)/gamm(c-b)*(-1)**i/gamm(dble(i+1))
c    #        *gamm(dble(k-i)-eps)*w**(eps+dble(i-k))
c 81  continue
c     write(10,*)temp2,term1

      re=gamm(c)*(one-w)**a*(temp+temp1+temp2*cos(pi*(eps-dble(k))))
c     write(10,*)re,(ttest+term1*cos(pi*(eps-dble(k))))
c    #             *gamm(c)*(one-w)**a

c  calculate the imaginary part

      im=temp2*sin(pi*(eps-dble(k)))

      poch1(0)=one
      poch2(0)=one
      do 90 i=1,k
        poch1(0)=(c-b+dble(i-1))*poch1(0)
        poch2(0)=(dble(i)-b)*poch2(0)
  90  continue
      do 91 i=1,n
        poch1(i)=(c-b+dble(k+i-1))*poch1(i-1)
        poch2(i)=(dble(k+i)-b)*poch2(i-1)
  91  continue

      temp=zero
      do 92 i=0,n
        temp=temp+poch1(i)/gamm(dble(k+i+1))*ff3(i)
     #       *w**(eps+dble(i))*poch2(i)
  92  continue

      x1=a
      coeff1=one
  93  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 93
      endif

      x2=b
      coeff2=one
  94  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 94
      endif

      if (dabs(eps).lt..1d0) then
        temp1=one
        i=1
  95    temp2=temp1+(-1)**i*(pi*eps)**(i+i)/gamm(dble(i+i+2))
        error=(temp2-temp1)/temp2
        if (dabs(error).lt.machep) go to 96
        i=i+1
        temp1=temp2
        go to 95
      else
        temp2=sin(pi*eps)/pi/eps
      endif
  96  continue
c     write(10,*)temp2,sin(pi*eps)/pi/eps

      im=im+temp*coeff1*coeff2/gamm(x1)/gamm(x2)*pi*temp2
      im=-im

      return
      end

c**********************************************************************
c
c  subroutine name    - fix6
c
c  computation
c  performed          - calculates the hypergeometric function for z
c                       greater than 2 when a-b is near an integer.
c
c  usage              - call fix6(a,b,c,n,k,re,im,w,machep,eps,pi)
c
c  arguments    a,b,c - parameters of the hypergeometric function.
c
c                  n  - the upper limit of the finite series expansion
c                       of the hypergeometric function.
c
c                  k  - equals the nearest integer of a-b.
c
c               re,im - computed values for the real and imaginary parts
c                       of the hypergeometric function.
c
c                  w  - transformed independent variable.
c
c              machep - equals machine epsilon.
c
c                eps  - equals a-b-k.
c
c                 pi  - equals 3.1415... to machine accuracy.
c
c  precision          - double
c
c  language           - fortran
c
c***********************************************************************

      subroutine fix6(a,b,c,n,k,re,im,w,machep,eps,pi)

      real*8  zero,one,two,four,eighth,seven,eight,sxteen
      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,eighth=1.d0/8.d0,
     #           seven=7.d0,eight=8.d0,sxteen=16.d0,nmax=100)
      real*8   a,b,c,w,re,im,gamm,temp,temp2,g1(0:nmax),g2(0:nmax),
     #         g3(0:nmax),g4(0:nmax),g5(0:nmax),x,x1,x2,x3,x4,psi,rn,
     #         t1(0:80),t2(0:80),t3(0:80),t4(0:80),test,machep,pi,
     #         f1(0:80),f2(0:80),f3(0:80),f4(0:80),ff3(0:nmax),eps,
     #         ff4(0:nmax),coeff1,coeff2,c1(0:80),c2(0:80),c3(0:80),
     #         c4(0:80),sum,term1,term2,term3,term4,term5,et1,et2,error,
     #         term6,temp1,coeff,coeff3,coeff4,fff1(0:nmax),ff1(0:nmax),
     #         fff2(0:nmax),ff2(0:nmax),poch1(0:nmax),poch2(0:nmax),e1

      integer  flag

      x3=zero
      call cheb(c3,55,2)
      t3(0)=one
      t3(1)=two*(x3+eps)
      f3(0)=zero
      f3(1)=two
      g3(0)=c3(1)*f3(1)

      x4=zero
      call cheb(c4,55,2)
      t4(0)=one
      t4(1)=two*(x4-eps)
      f4(0)=zero
      f4(1)=-two
      g4(0)=c4(1)*f4(1)

      do 7 i=2,55
        t3(i)=four*(x3+eps)*t3(i-1)-t3(i-2)
        t4(i)=four*(x4-eps)*t4(i-1)-t4(i-2)
        f3(i)=four*t3(i-1)+four*x3*f3(i-1)-f3(i-2)
        f4(i)=-four*t4(i-1)+four*x4*f4(i-1)-f4(i-2)
        g3(0)=g3(0)+c3(i)*f3(i)
        g4(0)=g4(0)+c4(i)*f4(i)
  7   continue

      g4(0)=-g4(0)
      do 10 i=-k,-1
        g4(0)=(g4(0)+one/gamm(dble(k+i+2)))/(dble(k+i+1)+eps)
  10  continue

      test=-eps*dlog(w)
      temp=-dlog(w)
         if (dabs(test).ge.eighth) then
           temp=(exp(test)-one)/eps
         else
           i=1
  20       rn=(eps**(i)*(-dlog(w))**(i+1))/gamm(dble(i+2))
           if (dabs(rn).lt.machep) go to 30
           temp=temp+rn
           i=i+1
           go to 20
         endif
  30     g5(0)=temp*w**a

c     write(10,*)g3(0),gamm(eps)-one/eps
c     write(10,*)g4(0),(-1)**k*gamm(-eps-k)+one/eps/gamm(dble(k+1))
c     write(10,*)g5(0),w**(a-eps)/eps-w**a/eps

      do 60 i=1,n
        g3(i)=(g3(i-1)+one/gamm(dble(i+1)))/(dble(i)-eps)
        g4(i)=(g4(i-1)+one/gamm(dble(k+i+1)))/(dble(k+i)+eps)
        g5(i)=w*g5(i-1)
  60  continue

      do 65 i=0,n
        ff3(i)=eps*g3(i)+one/gamm(dble(i+1))
        ff4(i)=eps*g4(i)-one/gamm(dble(k+i+1))
  65  continue

c  calculate the extra terms

      x=b-one
      sum=zero
      coeff=one
      flag=0
  61  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 61
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 610    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 610
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 61
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 611 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 611    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 612 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 612    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 613 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 613    continue
      endif

      if (flag.eq.0) then
        x1=b
        coeff1=one
 614    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 614
        endif
        x2=b+eps
        coeff2=one
 615    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 615
        endif
        temp=sum+coeff*temp
        et1=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 616    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 616
        endif
        x2=x+one+eps
        coeff2=one
 617    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 617
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et1=sum+coeff*temp
      endif
      et1=-et1
c     write(10,*)et1,(one/gamm(a-dble(k)-eps)-one/gamm(a-dble(k)))
c    #                                                  /eps

      x=c-a+k-one
      sum=zero
      coeff=one
      flag=0
  62  if (x.gt.one) then
        sum=sum+coeff*gamm(x+eps)
        coeff=coeff*x
        x=x-one
        go to 62
      elseif (x.lt.zero) then
        x1=x+eps+two
        coeff1=one
 620    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 620
        endif
        sum=sum+coeff*coeff1/gamm(x1)
        coeff=coeff*(x+one)
        x=x+one
        flag=1
        go to 62
      endif

      if ((x .ge. .25d0).and.(x .le. .75d0)) then
        call cheb(c1,41,1)
        t1(0)=one
        t1(1)=two*(x+eps)-one
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 621 i=2,41
          t1(i)=(four*(x+eps)-two)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-two)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 621    continue
      elseif ((x .ge. 0.d0).and.(x .lt. .25d0)) then
        call cheb(c1,55,2)
        t1(0)=one
        t1(1)=two*(x+eps)
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 622 i=2,55
          t1(i)=four*(x+eps)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+four*x*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 622    continue
      elseif ((x .gt. .75d0).and.(x .le. 1.d0)) then
        call cheb(c1,34,3)
        t1(0)=one
        t1(1)=two*(x+eps)-two
        f1(0)=zero
        f1(1)=two
        temp=c1(1)*f1(1)
        do 623 i=2,34
          t1(i)=(four*(x+eps)-four)*t1(i-1)-t1(i-2)
          f1(i)=four*t1(i-1)+(four*x-four)*f1(i-1)-f1(i-2)
          temp=temp+c1(i)*f1(i)
 623    continue
      endif

      if (flag.eq.0) then
        x1=c-a+dble(k)
        coeff1=one
 624    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 624
        endif
        x2=c-a+dble(k)+eps
        coeff2=one
 625    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 625
        endif
        temp=sum+coeff*temp
        et2=-temp*coeff1*coeff2/gamm(x1)/gamm(x2)
      elseif (flag.eq.one) then
        x1=x+one
        coeff1=one
 626    if (x1.lt.one) then
          coeff1=x1*coeff1
          x1=x1+one
          go to 626
        endif
        x2=x+one+eps
        coeff2=one
 627    if (x2.lt.one) then
          coeff2=x2*coeff2
          x2=x2+one
          go to 627
        endif
        coeff=-coeff*coeff1*coeff2/gamm(x1)/gamm(x2)
        et2=sum+coeff*temp
      endif
      et2=-et2
c     write(10,*)et2,(one/gamm(c-b-eps)-one/gamm(c-b))/eps
c
      fff1(0)=one
      fff2(0)=one
      ff2(0)=one
      do 685 i=1,k
        fff1(0)=(b+dble(i-1))*fff1(0)
        fff2(0)=(b-c+eps+dble(i))*fff2(0)
        ff2(0)=(b-c+dble(i))*ff2(0)
 685  continue

      ff1(0)=one
      do 686 i=1,n
        fff1(i)=(b+dble(k+i-1))*fff1(i-1)
        fff2(i)=(b-c+eps+dble(k+i))*fff2(i-1)
        ff1(i)=(a+dble(i-1))*ff1(i-1)
        ff2(i)=(b-c+dble(k+i))*ff2(i-1)
 686  continue

      x1=a
      coeff1=one
 687  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 687
      endif

      x2=a-dble(k)
      coeff2=one
 688  if (x2.lt.one) then
        coeff2=x2*coeff2
        x2=x2+one
        go to 688
      endif

      x3=c-b-eps
      coeff3=one
 689  if (x3.lt.one) then
        coeff3=x3*coeff3
        x3=x3+one
        go to 689
      endif

      x4=c-b
      coeff4=one
 690  if (x4.lt.one) then
        coeff4=x4*coeff4
        x4=x4+one
        go to 690
      endif

      do 691 i=0,n
        fff1(i)=fff1(i)*coeff1/gamm(x1)
        ff1(i)=ff1(i)*coeff2/gamm(x2)
        fff2(i)=fff2(i)*coeff3/gamm(x3)
        ff2(i)=ff2(i)*coeff4/gamm(x4)
c       write(10,*)'fff1=',fff1(i),gamm(b+dble(i+k))/gamm(a)/gamm(b)
c       write(10,*)'ff1=',ff1(i),gamm(a+dble(i))/gamm(a)/gamm(a-dble(k))
c       write(10,*)'fff2=',fff2(i),gamm(a-c+dble(i+1))/gamm(c-b-eps)
c    #                               /gamm(one-c+b+eps)
c       write(10,*)'ff2=',ff2(i),gamm(b-c+dble(i+k+1))/gamm(c-b)
c    #                               /gamm(one-c+b)
 691  continue

c   calculate  g1

      g1(0)=zero
      poch1(0)=one
      do 697 i=1,k
        g1(0)=g1(0)*(a-eps+dble(i-k-1))-poch1(0)
        poch1(0)=poch1(0)*(a+dble(i-k-1))
 697  continue
      do 698 i=1,n
        poch1(i)=(a+dble(i-1))*poch1(i-1)
        g1(i)=g1(i-1)*(a-eps+dble(i-1))-poch1(i-1)
 698  continue

      x1=a
      coeff1=one
 699  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 699
      endif
      do 700 i=0,n
        g1(i)=g1(i)*coeff1/gamm(x1)
c       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
 700  continue

c   calculate  g2

      g2(0)=zero
      poch2(0)=one
      do 701 i=1,k
        g2(0)=g2(0)*(b-c+eps+dble(i))+poch2(0)
        poch2(0)=poch2(0)*(b-c+dble(i))
 701  continue
      do 702 i=1,n
        poch2(i)=(b-c+dble(i+k))*poch2(i-1)
        g2(i)=g2(i-1)*(b-c+eps+dble(i+k))+poch2(i-1)
 702  continue

      x1=c-b
      coeff1=one
 703  if (x1.lt.one) then
        coeff1=x1*coeff1
        x1=x1+one
        go to 703
      endif

      poch2(0)=one
      do 704 i=1,k
        poch2(0)=(b-c+eps+dble(i))*poch2(0)
 704  continue
      do 705 i=1,n
        poch2(i)=(b-c+eps+dble(i+k))*poch2(i-1)
 705  continue

      do 706 i=0,n
        g2(i)=et2*poch2(i)+g2(i)*coeff1/gamm(x1)
c       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
 706  continue

c  calculate  e1

      e1=zero
      if (dabs(eps).lt..1d0) then
        i=1
 707    rn=(-1)**i*pi**(i+i)*eps**(i+i-1)/gamm(dble(i+i+1))
        if (dabs(rn).lt.machep) go to 708
        e1=e1+rn
        i=i+1
        go to 707
      else
        e1=(cos(pi*eps)-one)/eps
      endif
 708  continue
c     write(10,*)'e1=',e1,(cos(pi*eps)-one)/eps

      poch1(0)=one
      poch2(0)=one
      do 709 i=1,n
        poch1(i)=(a+dble(i-1))*poch1(i-1)
        poch2(i)=(a-c+dble(i))*poch2(i-1)
 709  continue

c  put everything back together again

      temp=zero
      temp1=zero
      temp2=zero
      do 70 i=0,n
        term1=-g1(i)/gamm(dble(i+1))*fff2(i)*ff4(i)
     #                     *w**(a+dble(i))*cos(pi*eps)
        term2=fff1(i)/gamm(dble(i+1))*g2(i)*ff4(i)
     #                     *w**(a+dble(i))*cos(pi*eps)
        term3=-fff1(i)*ff2(i)*g3(i)*ff4(i)
     #                     *w**(a+dble(i))*cos(pi*eps)
        term4=fff1(i)*ff2(i)*ff3(i)*g4(i)
     #                     *w**(a+dble(i))*cos(pi*eps)
        term5=fff1(i)/gamm(dble(k+i+1))*ff2(i)*ff3(i)
     #                     *g5(i)*cos(pi*eps)
        term6=-fff1(i)/gamm(dble(k+i+1))*ff2(i)*ff3(i)
     #                     *w**(a-eps+dble(i))*e1
        temp=temp+term1+term2+term3+term4+term5+term6
        temp1=temp1+poch1(i)/gamm(dble(i+1))*poch2(i)*ff4(i)
     #                                      *w**(a+dble(i))
        temp2=temp2+poch1(i)/gamm(dble(i+1))*fff2(i)*ff4(i)
     #                                      *w**(a+dble(i))
  70  continue

      x1=b
      coeff1=one
  71  if (x1.lt.one) then
        coeff1=coeff1*x1
        x1=x1+one
        go to 71
      endif

      x2=c-a
      coeff2=one
  72  if (x2.lt.one) then
        coeff2=coeff2*x2
        x2=x2+one
        go to 72
      endif

      term1=temp*gamm(c)*cos(pi*b)*(-1)**k
      term2=-temp1*gamm(c)*sin(pi*b)*coeff1/gamm(x1)*coeff2/gamm(x2)
      term3=temp2*gamm(c)*cos(pi*b)*cos(pi*eps)*(-1)**k*et1
      term4=temp*gamm(c)*sin(pi*b)*(-1)**k
      term5=temp1*gamm(c)*cos(pi*b)*coeff1/gamm(x1)*coeff2/gamm(x2)
      term6=term6*gamm(c)*sin(pi*b)*cos(pi*eps)*(-1)**k*et1

      if (dabs(eps).lt..1d0) then
        temp1=one
        i=1
  80    temp2=temp1+(-1)**i*(pi*eps)**(i+i)/gamm(dble(i+i+2))
        error=(temp2-temp1)/temp2
        if (dabs(error).lt.machep) go to 81
        i=i+1
        temp1=temp2
        go to 80
      else
        temp2=sin(pi*eps)/pi/eps
      endif
  81  continue
c     write(10,*)temp2,sin(pi*eps)/(pi*eps)

      term2=term2*pi*temp2
      term5=term5*pi*temp2

      re=term1+term2+term3
      im=term4+term5+term6

c  calculate the finite series contribution

      poch1(0)=one
      poch2(0)=one
      do 82 i=1,n
        poch1(i)=(b+dble(i-1))*poch1(i-1)
        poch2(i)=(b-c+dble(i))*poch2(i-1)
  82  continue

      temp=zero
      do 83 i=0,k-1
        temp=temp+poch1(i)*poch2(i)/gamm(dble(i+1))*gamm(eps+dble(k-i))
     #                              *(-1)**i*w**(b+dble(i))
  83  continue

      x1=a
      coeff1=one
  84  if (x1.lt.one) then
        coeff1=coeff1*x1
        x1=x1+one
        go to 84
      endif

      x2=c-b
      coeff2=one
  85  if (x2.lt.one) then
        coeff2=coeff2*x2
        x2=x2+one
        go to 85
      endif

      temp=temp*gamm(c)*coeff1/gamm(x1)*coeff2/gamm(x2)

      re=re+temp*cos(pi*b)
      im=im+temp*sin(pi*b)

      return
      end

c**********************************************************************
*
*   subroutine name     - geteps
*
*   computation
*   performed           - compute the smallest number machep such that
*                           machep+1 is not equal to 1 in the finite
*                           precision arithmetic used by the computer.
*
*   usage               - call geteps(machep,neps)
*
*   argument     machep - double precision (output).  the smallest
*                           number such that machep+1 is not equal to 1
*                           in the finite precision arithmetic used by
*                           the computer.
*                  neps - integer (output).  machine epsilon is machep =
*                           (1/2)**neps
*
*   precision           - double
*
*   language            - fortran 77
*
************************************************************************
*
      subroutine geteps(machep,neps)
*
      real*8            machep,one,two,temp
      integer            neps
      parameter (one = 1.0d0, two = 2.0d0)
      machep = one
      neps = 0
100   continue
      machep = machep/two
      neps = neps+1
      temp = machep+one
      if (temp .ne. one) go to 100
*
      machep = two*machep
      neps = neps-1
*
      return
      end
*
c***********************************************************************

      subroutine binomc

c       a
c     (   ) = binom(a*(a+1)/2+b+1)
c       b

      double precision   binom,one
      common /bcoeff/binom(5151)

      maxnll=100

      if (maxnll .lt. 0) go to 300
      if (maxnll .gt. 100) go to 300
      one = 1.0d0
      binom(1) = one
      if (maxnll .eq. 0) go to 300
      binom(2) = one
      binom(3) = one
      if (maxnll .eq. 1) go to 300
      ij = 4
      imax = maxnll+1
      do 200 i = 3,imax
         ii = ((i-1)*(i-2))/2
         binom(ij) = one
         ij = ij+1
         jmax = i-1
         do 100 j = 2,jmax
            binom(ij) = binom(ii+j-1)+binom(ii+j)
            ij = ij+1
100      continue
         binom(ij) = one
         ij = ij+1
200   continue
c
300   continue
c
      return
      end

c***********************************************************************
c
c  subroutine name    - cheb
c
c  computation
c  performed          - tabulates the tchebychev coefficients which
c                       were computed by the program 'tcheb2'.  the
c                       three sets of coefficients correspond to
c                       the three gamma function expansions shown in
c                       equations (4.35),(4.36), and (4.37). see
c                       'tcheb2' for additional documentation.
c
c  usage              - call cheb(c,n,flag)
c
c  arguments       c  - the array (output) which contains the
c                       tchebychev coefficients.
c
c                  n  - the dimension (input) of the array 'c'.
c
c                flag - the parameter (input) which tells the sub-
c                       routine which tchebychev coefficients to
c                       return to the caller.
c
c  precision          - double (although the coefficients are
c                               accurate to quadruple)
c
c  language           - fortran 77
c
c***********************************************************************

      subroutine cheb(c,n,flag)

      real*8  c(0:n)
      integer  flag

      if (flag.eq.1) go to 100
      if (flag.eq.2) go to 200
      if (flag.eq.3) go to 300

c  tchebychev expansion coefficients for the range, 0<x<1

 100  c(0) =  0.94178559779549466571096003120435196d+00
      c(1) =  0.44153813248410067571913157711414607d-02
      c(2) =  0.56850436815993633786326645888162378d-01
      c(3) = -0.42198353964185605010125001866024699d-02
      c(4) =  0.13268081812124602205840067963889683d-02
      c(5) = -0.18930245297988804325239470239464680d-03
      c(6) =  0.36069253274412452565780822094442805d-04
      c(7) = -0.60567619044608642184855483216922771d-05
      c(8) =  0.10558295463022833447318234541645507d-05
      c(9) = -0.18119673655423840482918555144273961d-06
      c(10)=  0.31177249647153222777902517006137963d-07
      c(11)= -0.53542196390196871408740949118221475d-08
      c(12)=  0.91932755198595889468877475468573503d-09
      c(13)= -0.15779412802883397617671187106425584d-09
      c(14)=  0.27079806229349545432695717700017206d-10
      c(15)= -0.46468186538257301439531283506784063d-11
      c(16)=  0.79733501920074196555512936759234830d-12
      c(17)= -0.13680782098309160264738694164685656d-12
      c(18)=  0.23473194865638006534799539031857605d-13
      c(19)= -0.40274326149490669507857892267787757d-14
      c(20)=  0.69100517473721009958174457696435176d-15
      c(21)= -0.11855845002219929396593062972684083d-15
      c(22)=  0.20341485424963760969383490105975402d-16
      c(23)= -0.34900543417173691101844936408331408d-17
      c(24)=  0.59879938564842634972645168624438135d-18
      c(25)= -0.10273780578716378747008169519685451d-18
      c(26)=  0.17627028160574041125936108594612916d-19
      c(27)= -0.30243206536626379817809691872233988d-20
      c(28)=  0.51889146600668142375785699199940389d-21
      c(29)= -0.89027708392150216484577040964212789d-22
      c(30)=  0.15274740724470977041487116294681806d-22
      c(31)= -0.26207312865170684216151526387496724d-23
      c(32)=  0.44964644619824783627762340991300087d-24
      c(33)= -0.77147147879836211531329396406348717d-25
      c(34)=  0.13236365808260955301316348853544449d-25
      c(35)= -0.22709797413377406198008958539204735d-26
      c(36)=  0.38966913277073699893252807432563276d-27
      c(37)= -0.66795989154793901466615113245736539d-28
      c(38)=  0.11456694360946249087722449327564468d-28
      c(39)= -0.20956088513945987438866120550893160d-29
      c(40)=  0.34345153487326051089311279207743562d-30
      c(41)= -0.74448389617685196161619686887550341d-31
      return

c  tchebychev expansion coefficients for the range,  -.5<x<.5

 200  c(0) =  0.11528686913857579339872890819003657d+01
      c(1) = -0.39836641427188668813550502856567435d+00
      c(2) =  0.16381491849746834445969671065563396d+00
      c(3) = -0.41349972584595838242416447164595642d-01
      c(4) =  0.11739888104509743948748485834561229d-01
      c(5) = -0.31509159742825717845846783104528302d-02
      c(6) =  0.85084809366682540330028115184077086d-03
      c(7) = -0.22845443192182297253614554810213881d-03
      c(8) =  0.61296656896858907270916323759970391d-04
      c(9) = -0.16433766723011959082591541534833589d-04
      c(10)=  0.44046701847148520660258125028242579d-05
      c(11)= -0.11803851479587223345492859134791582d-05
      c(12)=  0.31630339312403588488305625683201151d-06
      c(13)= -0.84755796666686117564957022251013564d-07
      c(14)=  0.22710572677209079780536954678987573d-07
      c(15)= -0.60853209609268373214751556259951644d-08
      c(16)=  0.16305620921375867864482570008163625d-08
      c(17)= -0.43690846345047718022878883179027790d-09
      c(18)=  0.11706935476739890379554689241357534d-09
      c(19)= -0.31368649843198552351255033209421610d-10
      c(20)=  0.84052057618382692960217222664957228d-11
      c(21)= -0.22521682699590609081199019088965996d-11
      c(22)=  0.60346669123807723976181127096882828d-12
      c(23)= -0.16169841538137032176079290114309245d-12
      c(24)=  0.43326960175123609635570088625382667d-13
      c(25)= -0.11609424034675431553315176322024985d-13
      c(26)=  0.31107358004300087572452155428660087d-14
      c(27)= -0.83351914632193111475558815401948979d-15
      c(28)=  0.22334078222557889355389486422061460d-15
      c(29)= -0.59843982246058550382747881611851515d-16
      c(30)=  0.16035146716190080240936859943115090d-16
      c(31)= -0.42966046133076898235808019603294715d-17
      c(32)=  0.11512717363557431988678458870224873d-17
      c(33)= -0.30848233202835882015258583966299712d-18
      c(34)=  0.82657591746540727258216017499064442d-19
      c(35)= -0.22148034956862123422799663231945171d-19
      c(36)=  0.59345480806145642339133686333296721d-20
      c(37)= -0.15901573656881585725893714030807897d-20
      c(38)=  0.42608138203898096080539369435375448d-21
      c(39)= -0.11416816226321087557458906349840213d-21
      c(40)=  0.30591266842950015571055286508657438d-22
      c(41)= -0.81969053674548061989664444282339330d-23
      c(42)=  0.21963543471485197662543467891802004d-23
      c(43)= -0.58851140572211577956963471197095354d-24
      c(44)=  0.15769121438531798083082131134888596d-24
      c(45)= -0.42253211944581570323425035302537635d-25
      c(46)=  0.11321706791574145306428072576766804d-25
      c(47)= -0.30335842761477973373797446515125892d-26
      c(48)=  0.81281383350578045680446098123885346d-27
      c(49)= -0.21782407988772728568103833180457024d-27
      c(50)=  0.58395544064782062129754390403734767d-28
      c(51)= -0.15729062977489325257494410942884130d-28
      c(52)=  0.42390612257722955199550993363196147d-29
      c(53)= -0.11242203351086692027388616387423238d-29
      c(54)=  0.27892280419588143241883200553486195d-30
      c(55)= -0.75766427928255356179910217971637866d-31
      return

c  tchebychev expansion coefficients for the range,  .5<x<1.5

 300  c(0) =  0.10532770878177862619534128247576828d+01
      c(1) =  0.21902166104535936497306369004840667d+00
      c(2) =  0.53885821783347712865216341722976574d-01
      c(3) =  0.25387290658986838596948519579519148d-02
      c(4) =  0.61466596479014144199820446583715941d-03
      c(5) = -0.32319247384294465724865638122474435d-05
      c(6) =  0.60054921157267140200751871810266970d-05
      c(7) = -0.41824428090189489334617924547407754d-06
      c(8) =  0.74607235650174366232051332482639985d-07
      c(9) = -0.84349526185192483560074198183789434d-08
      c(10)=  0.11322169721817117406057072389666464d-08
      c(11)= -0.14175349900034682206860980369914924d-09
      c(12)=  0.18156967683771854495445069753509525d-10
      c(13)= -0.23052163748763990586386231147733255d-11
      c(14)=  0.29327030584105892891631030300077869d-12
      c(15)= -0.37268590170679729030689484336505900d-13
      c(16)=  0.47360432581610222494078892575939043d-14
      c(17)= -0.60172423075766780010690060490450222d-15
      c(18)=  0.76443979970650480527157880770622904d-16
      c(19)= -0.97108892590783757664936380167684001d-17
      c(20)=  0.12335488659810502174628042595177563d-17
      c(21)= -0.15668997427423797214874298423999374d-18
      c(22)=  0.19902969432180950952170748993213290d-19
      c(23)= -0.25280701093316992983208535829903356d-20
      c(24)=  0.32111217127088658654008440525466587d-21
      c(25)= -0.40787027055654288157193053732139852d-22
      c(26)=  0.51806681115442807351458062924762066d-23
      c(27)= -0.65803415226414646040514708695329147d-24
      c(28)=  0.83581632724068042390791744946381128d-25
      c(29)= -0.10616267321620223331012310816058461d-25
      c(30)=  0.13484159784261929973156667845986312d-26
      c(31)= -0.17130640476670792317750095910458264d-27
      c(32)=  0.21720215147689502411187819143753676d-28
      c(33)= -0.27633054946463729557612727034555572d-29
      c(34)=  0.26664265210535867308016959008022142d-30
      return
      end
