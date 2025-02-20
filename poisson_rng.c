double poidev(double xm) {
  static double sq,alxm,g,oldm=-1.0;
  double em,t,y;
  if (xm<12.0) {
    if (xm!=oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em=-1;
    t=1.0;
    do {
      em+=1.0;
      t*=UNIFORM_01;
    }
    while (t>g);
  } else {
    if (xm!=oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-lgamma(xm+1.0);
    }
    do {
        do {
          y=tan(3.1415926535897932385*UNIFORM_01);
          em=sq*y+xm;
        } while (em<0.0);
        em=floor(em);
        t=0.9*(1.0+y*y)*exp(em*alxm-lgamma(em+1.0)-g);
    } while (UNIFORM_01>t);
  } 
  return em;
}
