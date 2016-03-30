subroutine qng ( f, ncomp, a, b, epsabs, epsrel, result, abserr, neval, ier )

  !*****************************************************************************80
  !
  !! QNG estimates an integral, using non-adaptive integration.
  !
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A,B),
  !    hopefully satisfying
  !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
  !
  !    The routine is a simple non-adaptive automatic integrator, based on
  !    a sequence of rules with increasing degree of algebraic
  !    precision (Patterson, 1968).
  !
  !  Author:
  !
  !    Robert Piessens, Elise de Doncker-Kapenger, 
  !    Christian Ueberhuber, David Kahaner
  !
  !  Reference:
  !
  !    Robert Piessens, Elise de Doncker-Kapenger, 
  !    Christian Ueberhuber, David Kahaner,
  !    QUADPACK, a Subroutine Package for Automatic Integration,
  !    Springer Verlag, 1983
  !
  !  Parameters:
  !
  !    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
  !      function f ( x )
  !      real ( kind = 4 ) f
  !      real ( kind = 4 ) x
  !    which evaluates the integrand function.
  !
  !    Input, real ( kind = 4 ) A, B, the limits of integration.
  !
  !    Input, real ( kind = 4 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
  !    RESULT is obtained by applying the 21-point Gauss-Kronrod rule (RES21)
  !    obtained  by optimal addition of abscissae to the 10-point Gauss rule
  !    (RES10), or by applying the 43-point rule (RES43) obtained by optimal
  !    addition of abscissae to the 21-point Gauss-Kronrod rule, or by 
  !    applying the 87-point rule (RES87) obtained by optimal addition of
  !    abscissae to the 43-point rule.
  !
  !    Output, real ( kind = 4 ) ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer ( kind = 4 ) NEVAL, the number of times the integral was evaluated.
  !
  !           ier    - ier = 0 normal and reliable termination of the
  !                            routine. it is assumed that the requested
  !                            accuracy has been achieved.
  !                    ier > 0 abnormal termination of the routine. it is
  !                            assumed that the requested accuracy has
  !                            not been achieved.
  !                    ier = 1 the maximum number of steps has been
  !                            executed. the integral is probably too
  !                            difficult to be calculated by qng.
  !                        = 6 the input is invalid, because
  !                            epsabs < 0 and epsrel < 0,
  !                            result, abserr and neval are set to zero.
  !
  !  Local Parameters:
  !
  !           centr  - mid point of the integration interval
  !           hlgth  - half-length of the integration interval
  !           fcentr - function value at mid point
  !           absc   - abscissa
  !           fval   - function value
  !           savfun - array of function values which have already
  !                    been computed
  !           res10  - 10-point Gauss result
  !           res21  - 21-point Kronrod result
  !           res43  - 43-point result
  !           res87  - 87-point result
  !           resabs - approximation to the integral of abs(f)
  !           resasc - approximation to the integral of abs(f-i/(b-a))
  !
  implicit none

  interface
     function f(x, ncomp)
       real ( kind = 8) x
       integer ( kind = 8) ncomp
       real ( kind = 8), dimension(ncomp) :: f
     end function f
  end interface
  
  integer ( kind = 8) ncomp
  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr(ncomp)
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  !real ( kind = 8 ), external :: f
  real ( kind = 8 ) fcentr(ncomp)
  real ( kind = 8 ) fval(ncomp)
  real ( kind = 8 ) fval1(ncomp)
  real ( kind = 8 ) fval2(ncomp)
  real ( kind = 8 ) fv1(ncomp,5)
  real ( kind = 8 ) fv2(ncomp,5)
  real ( kind = 8 ) fv3(ncomp,5)
  real ( kind = 8 ) fv4(ncomp,5)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) ier(ncomp)
  integer ( kind = 8 ) ipx
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l
  integer ( kind = 8 ) neval
  real ( kind = 8 ) result(ncomp)
  real ( kind = 8 ) res10(ncomp)
  real ( kind = 8 ) res21(ncomp)
  real ( kind = 8 ) res43(ncomp)
  real ( kind = 8 ) res87(ncomp)
  real ( kind = 8 ) resabs(ncomp)
  real ( kind = 8 ) resasc(ncomp)
  real ( kind = 8 ) reskh(ncomp)
  real ( kind = 8 ) savfun(ncomp,21)
  real ( kind = 8 ) w10(5)
  real ( kind = 8 ) w21a(5)
  real ( kind = 8 ) w21b(6)
  real ( kind = 8 ) w43a(10)
  real ( kind = 8 ) w43b(12)
  real ( kind = 8 ) w87a(21)
  real ( kind = 8 ) w87b(23)
  real ( kind = 8 ) x1(5)
  real ( kind = 8 ) x2(5)
  real ( kind = 8 ) x3(11)
  real ( kind = 8 ) x4(22)
  !
  !           the following data statements contain the abscissae
  !           and weights of the integration rules used.
  !
  !           x1      abscissae common to the 10-, 21-, 43- and 87-point
  !                   rule
  !           x2      abscissae common to the 21-, 43- and 87-point rule
  !           x3      abscissae common to the 43- and 87-point rule
  !           x4      abscissae of the 87-point rule
  !           w10     weights of the 10-point formula
  !           w21a    weights of the 21-point formula for abscissae x1
  !           w21b    weights of the 21-point formula for abscissae x2
  !           w43a    weights of the 43-point formula for absissae x1, x3
  !           w43b    weights of the 43-point formula for abscissae x3
  !           w87a    weights of the 87-point formula for abscissae x1,
  !                   x2 and x3
  !           w87b    weights of the 87-point formula for abscissae x4
  !
  data x1(1),x1(2),x1(3),x1(4),x1(5)/ &
       9.739065285171717D-01,     8.650633666889845D-01, &
       6.794095682990244D-01,     4.333953941292472D-01, &
       1.488743389816312D-01/
  data x2(1),x2(2),x2(3),x2(4),x2(5)/ &
       9.956571630258081D-01,     9.301574913557082D-01, &
       7.808177265864169D-01,     5.627571346686047D-01, &
       2.943928627014602D-01/
  data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),x3(9),x3(10), &
       x3(11)/ &
       9.993333609019321D-01,     9.874334029080889D-01, &
       9.548079348142663D-01,     9.001486957483283D-01, &
       8.251983149831142D-01,     7.321483889893050D-01, &
       6.228479705377252D-01,     4.994795740710565D-01, &
       3.649016613465808D-01,     2.222549197766013D-01, &
       7.465061746138332D-02/
  data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),x4(10), &
       x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),x4(19), &
       x4(20),x4(21),x4(22)/         9.999029772627292D-01, &
       9.979898959866787D-01,     9.921754978606872D-01, &
       9.813581635727128D-01,     9.650576238583846D-01, &
       9.431676131336706D-01,     9.158064146855072D-01, &
       8.832216577713165D-01,     8.457107484624157D-01, &
       8.035576580352310D-01,     7.570057306854956D-01, &
       7.062732097873218D-01,     6.515894665011779D-01, &
       5.932233740579611D-01,     5.314936059708319D-01, &
       4.667636230420228D-01,     3.994248478592188D-01, &
       3.298748771061883D-01,     2.585035592021616D-01, &
       1.856953965683467D-01,     1.118422131799075D-01, &
       3.735212339461987D-02/
  data w10(1),w10(2),w10(3),w10(4),w10(5)/ &
       6.667134430868814D-02,     1.494513491505806D-01, &
       2.190863625159820D-01,     2.692667193099964D-01, &
       2.955242247147529D-01/
  data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/ &
       3.255816230796473D-02,     7.503967481091995D-02, &
       1.093871588022976D-01,     1.347092173114733D-01, &
       1.477391049013385D-01/
  data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/ &
       1.169463886737187D-02,     5.475589657435200D-02, &
       9.312545458369761D-02,     1.234919762620659D-01, &
       1.427759385770601D-01,     1.494455540029169D-01/
  data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7), &
       w43a(8),w43a(9),w43a(10)/     1.629673428966656D-02, &
       3.752287612086950D-02,     5.469490205825544D-02, &
       6.735541460947809D-02,     7.387019963239395D-02, &
       5.768556059769796D-03,     2.737189059324884D-02, &
       4.656082691042883D-02,     6.174499520144256D-02, &
       7.138726726869340D-02/
  data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),w43b(7), &
       w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/ &
       1.844477640212414D-03,     1.079868958589165D-02, &
       2.189536386779543D-02,     3.259746397534569D-02, &
       4.216313793519181D-02,     5.074193960018458D-02, &
       5.837939554261925D-02,     6.474640495144589D-02, &
       6.956619791235648D-02,     7.282444147183321D-02, &
       7.450775101417512D-02,     7.472214751740301D-02/
  data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),w87a(7), &
       w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),w87a(13),w87a(14), &
       w87a(15),w87a(16),w87a(17),w87a(18),w87a(19),w87a(20),w87a(21)/ &
       8.148377384149173D-03,     1.876143820156282D-02, &
       2.734745105005229D-02,     3.367770731163793D-02, &
       3.693509982042791D-02,     2.884872430211531D-03, &
       1.368594602271270D-02,     2.328041350288831D-02, &
       3.087249761171336D-02,     3.569363363941877D-02, &
       9.152833452022414D-04,     5.399280219300471D-03, &
       1.094767960111893D-02,     1.629873169678734D-02, &
       2.108156888920384D-02,     2.537096976925383D-02, &
       2.918969775647575D-02,     3.237320246720279D-02, &
       3.478309895036514D-02,     3.641222073135179D-02, &
       3.725387550304771D-02/
  data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7), &
       w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14), &
       w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),w87b(21), &
       w87b(22),w87b(23)/            2.741455637620724D-04, &
       1.807124155057943D-03,     4.096869282759165D-03, &
       6.758290051847379D-03,     9.549957672201647D-03, &
       1.232944765224485D-02,     1.501044734638895D-02, &
       1.754896798624319D-02,     1.993803778644089D-02, &
       2.219493596101229D-02,     2.433914712600081D-02, &
       2.637450541483921D-02,     2.828691078877120D-02, &
       3.005258112809270D-02,     3.164675137143993D-02, &
       3.305041341997850D-02,     3.425509970422606D-02, &
       3.526241266015668D-02,     3.607698962288870D-02, &
       3.669860449845609D-02,     3.712054926983258D-02, &
       3.733422875193504D-02,     3.736107376267902D-02/
  !
  !  Test on validity of parameters.
  !
  do k = 1,ncomp
     result(k) = 0.0D+00
     abserr(k) = 0.0D+00
  end do

  neval = 0

  if ( epsabs < 0.0D+00 .and. epsrel < 0.0D+00 ) then
     do k=1,ncomp
        ier(k) = 6
     end do
     return
  end if

  hlgth = 5.0D-01 * ( b - a )
  dhlgth = abs ( hlgth )
  centr = 5.0D-01 * ( b + a )
  fcentr = f(centr,ncomp)
  neval = 21
  do k=1,ncomp
     ier(k) = 1
  end do
  !
  !  Compute the integral using the 10- and 21-point formula.
  !

  do l = 1, 3

     if ( l == 1 ) then

        do k = 1,ncomp
           res10(k) = 0.0D+00
        end do

        res21 = w21b(6) * fcentr
        resabs = w21b(6) * abs(fcentr)


        do k = 1, 5
           absc = hlgth * x1(k)
           fval1 = f(centr+absc,ncomp)
           fval2 = f(centr-absc,ncomp)
           fval = fval1 + fval2
           res10 = res10 + w10(k)*fval
           res21 = res21 + w21a(k)*fval
           resabs = resabs + w21a(k)*(abs(fval1)+abs(fval2))
           savfun(:,k) = fval
           fv1(:,k) = fval1
           fv2(:,k) = fval2
        end do

        ipx = 5

        do k = 1, 5
           ipx = ipx + 1
           absc = hlgth * x2(k)
           fval1 = f(centr+absc,ncomp)
           fval2 = f(centr-absc,ncomp)
           fval = fval1 + fval2
           res21 = res21 + w21b(k) * fval
           resabs = resabs + w21b(k) * ( abs ( fval1 ) + abs ( fval2 ) )
           savfun(:,ipx) = fval
           fv3(:,k) = fval1
           fv4(:,k) = fval2
        end do
        !
        !  Test for convergence.
        !
        result = res21 * hlgth
        resabs = resabs * dhlgth
        reskh = 5.0D-01 * res21
        resasc = w21b(6) * abs ( fcentr - reskh )

        do k = 1, 5
           resasc = resasc+w21a(k)*(abs(fv1(:,k)-reskh)+abs(fv2(:,k)-reskh)) &
                +w21b(k)*(abs(fv3(:,k)-reskh)+abs(fv4(:,k)-reskh))
        end do

        abserr = abs ( ( res21 - res10 ) * hlgth )
        resasc = resasc * dhlgth
        !
        !  Compute the integral using the 43-point formula.
        !
     else if ( l == 2 ) then

        res43 = w43b(12)*fcentr
        neval = 43

        do k = 1, 10
           res43 = res43 + savfun(:,k) * w43a(k)
        end do

        do k = 1, 11
           ipx = ipx + 1
           absc = hlgth * x3(k)
           fval = f(absc+centr,ncomp) + f(centr-absc,ncomp)
           res43 = res43 + fval * w43b(k)
           savfun(:,ipx) = fval
        end do
        !
        !  Test for convergence.
        !
        result = res43 * hlgth
        abserr = abs((res43-res21)*hlgth)
        !
        !  Compute the integral using the 87-point formula.
        !
     else if ( l == 3 ) then
        
        res87 = w87b(23) * fcentr
        neval = 87

        do k = 1, 21
           res87 = res87 + savfun(:,k) * w87a(k)
        end do

        do k = 1, 22
           absc = hlgth * x4(k)
           res87 = res87 + w87b(k) * ( f(absc+centr,ncomp) + f(centr-absc,ncomp) )
        end do

        result = res87 * hlgth
        abserr = abs ( ( res87 - res43) * hlgth )

     end if

     do k = 1,ncomp
        if ( resasc(k) /= 0.0D+00.and.abserr(k) /= 0.0D+00 ) then
           abserr(k) = resasc(k) * min ( 1.0D+00,(2.0D+02*abserr(k)/resasc(k))**1.5D+00)
        end if

        if ( resabs(k)> tiny ( resabs(k)) / ( 5.0D+01 * epsilon ( resabs(k)) ) ) then
           abserr(k) = max (( epsilon ( resabs(k)) *5.0D+01) * resabs(k), abserr(k) )
        end if

        
     end do

     if ( abserr(k) <= max ( epsabs, epsrel*abs(result(k)))) then
        ier(k) = 0
     end if
     
  end do

  return
end subroutine qng
