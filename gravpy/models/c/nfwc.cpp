
#include nfwc.h

void nfw(double** res, double* x, double* y, double* modelargs, int len){
  double res[len][6];
  double x  [len];
  double y  [len];
  double modelargs[6];
  
  
  for (int i = 0; i < len; i++){
    single_eval(res[i], x[i], y[i], modelargs)
  }
  
}

void single_eval(double res[6], double x, double y, double modelargs[6]){
  double b,x0,y0,e,te,s,a;
  double r,sint,cost,front,xx,dx,phir_r,phirr,phix,phiy,phixx,phiyy,phixy,pot,t1,t2;
  double temparr[6];
  
  b = modelargs[0];
  e = modelargs[3];
  s = modelargs[5];

  r = sqrt(x*x+y*y);

  if (r > 0.0){
    cost = x/r;
    sint = y/r;
  }
  else{
    cost = 1.0;
    sint = 0.0;
  }

  front = (1.0-e)*b;

  if (abs(e) < 1.0e-6){
    xx = r/s;
    dx = xx-1.0;
    phir_r = front*nfw0phir(xx)/xx;

    if (abs(dx) < 1.0e-4){
      phirr = -phir_r + 4.0*front*(1.0/3.0-0.4*dx+13.0/35.0*dx*dx-20.0/63.0*dx*dx*dx);
    }
    else {
      phirr = -phir_r + 4.0*front*(1.0-nfwFfunc(xx))/(xx*xx-1.0);
    }

    phix  = phir_r*x;
    phiy  = phir_r*y;
    phixx = phir_r*sint*sint + phirr*cost*cost;
    phiyy = phir_r*cost*cost + phirr*sint*sint;
    phixy = (phirr-phir_r)*sint*cost;

    t1 = log(xx/2.0);
    if (xx < 1.0e-2)  {
      pot = -0.5*xx*xx*(t1+(1.0+3.0*t1)*xx*xx/8.0+(3.0/32.0+5.0/24.0*t1)*xx*xx*xx*xx);}
    else if (x <= 1.0){
      t2  = atanh(sqrt(1.0-xx*xx));
      pot = t1*t1-t2*t2;}
    else{
      t2  = atanh(sqrt(xx*xx-1.0));
      pot = t1*t1+t2*t2;}

    pot *= 2*front*s*s;
      
  }
  else{
    nfw_integral(temparr,1.0e-8,1.0,r,sint,cost,e,s);

    pot   = front*0.5     *temparr[1];
    phix  = front*x*temparr[2];
    phiy  = front*y*temparr[3];
    phixx = front*(temparr[2]+2.0*x*x*temparr[4]);
    phiyy = front*(temparr[3]+2.0*y*y*temparr[5]);
    phixy = front*(           2.0*x*y*temparr[6]);
      
  }    

  res = {pot,phix,phiy,phixx,phiyy,phixy};
    
}

double nfw0phir(double x){
  double tmp,nfwFfunc(double);
  if (x == 0.0) return 0.0;
  nfw_F = nfwFfunc(x);
  if (x < 1.0e-4) {
    tmp = -2.0*(x*(1.0+0.75*x*x)*log(0.5*x)+0.5*x*(1.0+0.875*x*x));
  } else {
    tmp = 4.0/x*(log(0.5*x)+nfw_F);
  }
  return tmp;
}

double nfwFfunc(double x)
{
  double l2x,dx,tmp,ans,atanh(double);

  dx = x*x - 1.0;
  if (fabs(x)<1.0e-2) {
    /* series with error O[x^6] */
    l2x  = log(2.0/x);
    ans  = l2x;
    ans += x*x*(0.5*l2x-0.25);
    ans += x*x*x*x*(0.375*l2x-0.21875);
  } else if (fabs(dx)<1.0e-2) {
    /* series with error O[dx^6] */
    ans  = 1.0;
    ans -= dx/3.0;
    ans += dx*dx/5.0;
    ans -= dx*dx*dx/7.0;
    ans += dx*dx*dx*dx/9.0;
    ans -= dx*dx*dx*dx*dx/11.0;
  } else if (x < 1.0) {
    tmp = sqrt(1.0-x*x);
    ans = atanh(tmp)/tmp;
  } else {
    tmp = sqrt(x*x-1.0);
    ans = atan(tmp)/tmp;
  }
  return ans;
}

void nfwkap(double u, double* mphiu, double* k, double* kp)
{
    double q,t0,t1,t2,x2,dx,phir,nfw0phir(double);

    /* utility */
    q  = 1.0-ep;
    t0 = cost*cost+sint*sint/(1.0-(1.0-q*q)*u);
    t1 = t0*r*r;
    t2 = t1*u;
    x2 = t2/(ss*ss);
    dx = x2-1.0;

    /* F function computed in here */
    phir = ss*nfw0phir(sqrt(x2));

    *mphiu = sqrt(t2)*phir/u;

    if (fabs(dx)<1.0e-4) {
      /* deflection integrals */
      *k  = 2.0/3.0-0.4*dx+2.0/7.0*dx*dx-2.0/9.0*dx*dx*dx;
      /* magnification integrals */
      *kp = -0.4+4.0/7.0*dx-2.0/3.0*dx*dx+8.0/11.0*dx*dx*dx;
    } else {
      /* deflection integrals */
      *k  = 2.0*(1.0-F)/dx;
      /* magnification integrals */
      *kp = (3.0*x2*F-2.0*x2-1.0)/(x2*ss*ss*dx*dx);
    }
}
