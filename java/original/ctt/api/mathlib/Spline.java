
/*--------------------------------------------------------------------------*
 * Spline3
 *
 * Copyright 2004 Battelle Memorial Institute. All rights reserved.
 * This computer program contains the confidential information of Battelle.
 * Use or duplication in any form is subject to conditions set forth in a
 * separate license agreement. Unauthorized use or duplication may violate
 * the license agreement, state, federal and/or international laws including
 * the Copyright Laws of the United States and of other international
 * jurisdictions.
 *
 * Author: mooneyd
 *
 * Date: 4/6/2004
 *
 * Class to compute B-spline and their derivatives including fitting B-splines to data
 * using weighted least squares
 * 
 * The core algorithms of this class are based on the algorithms and FORTRAN implementations in 
 * de Boor's A Practical Guide to Splines
 * 
 * References:  math
 * 
 * Revision: 
 *
 * History: 
 * 
 *
\*--------------------------------------------------------------------------*/
package ctt.api.mathlib;

public class Spline {
//	  subroutine l2err ( prfun , ftau , error )
//c  from  * a practical guide to splines *  by c. de boor    
//c  this routine is to be called in the main program  l 2 m a i n .
//calls subprogram  ppvalu(interv)
//c  this subroutine computes various errors of the current l2-approxi-
//c  mation , whose pp-repr. is contained in common block  approx  ,
//c  to the given data contained in common block  data . it prints out
//c  the average error  e r r l 1 , the l2-error  e r r l 2,  and the
//c  maximum error  e r r m a x .
//c
//c******  i n p u t  ******
//c  prfun  a hollerith string.  if prfun = 'ON', the routine prints out
//c          the value of the approximation as well as its error at
//c          every data point.
//c
//c******  o u t p u t  ******
//c  ftau(1), ..., ftau(ntau),  with  ftau(i)  the approximation  f at
//c          tau(i), all i.
//c  error(1), ..., error(ntau),  with  error(i) = scale*(g - f)
//c          at tau(i), all i. here,  s c a l e  equals  1. in case
//c          prfun .ne. 'ON' , or the abs.error is greater than 100 some-
//c          where. otherwise, s c a l e  is such that the maximum of
//c          abs(error))  over all  i  lies between  10  and  100. this
//c          makes the printed output more illustrative.
//c
//	  integer prfun,   ie,k,l,ll,lpkmax,ltkmax,ntau,ntmax,on
//	  real ftau(1),error(1),  break,coef,err,errmax,errl1,errl2
//	 *                       ,gtau,scale,tau,totalw,weight
//	  dimension ftau(ntau),error(ntau)
//	  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
//	  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
//	  common /approx/ break(lpkmax),coef(ltkmax),l,k
//C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
//C     common /approx/ break(100),coef(2000),l,k
//	  data on /'ON'/
//	  errl1 = 0.
//	  errl2 = 0.
//	  errmax = 0.
//	  do 10 ll=1,ntau
//		 ftau(ll) = ppvalu (break, coef, l, k, tau(ll), 0 )
//		 error(ll) = gtau(ll) - ftau(ll)
//		 err = abs(error(ll))
//		 if (errmax .lt. err)   errmax = err
//		 errl1 = errl1 + err*weight(ll)
//   10    errl2 = errl2 + err**2*weight(ll)
//	  errl1 = errl1/totalw
//	  errl2 = sqrt(errl2/totalw)
//	  print 615,errl2,errl1,errmax
//  615 format(///' least square error =',e20.6/
//	 1          ' average error      =',e20.6/
//	 2          ' maximum error      =',e20.6//)
//	  if (prfun .ne. on)                return
//c     **  scale error curve and print  **
//	  ie = 0
//	  scale = 1.
//	  if (errmax .ge. 10.)              go to 18
//	  do 17 ie=1,9
//		 scale = scale*10.
//		 if (errmax*scale .ge. 10.)     go to 18
//   17    continue
//   18 do 19 ll=1,ntau
//   19    error(ll) = error(ll)*scale
//	  print 620,ie,(ll,tau(ll),ftau(ll),error(ll),ll=1,ntau)
//  620 format(///14x,'approximation and scaled error curve'/7x,
//	 1'data point',7x,'approximation',3x,'deviation x 10**',i1/
//	 2(i4,f16.8,f16.8,f17.6))
//										return
//	  end

public static void	plotdata ( double minRange, double maxRange, int numPoints,double[] t,  int n, int k, double[] bcoef){
	int i,j,jj,left,leftmk, ll, mm;
		final int kmax=20;
		final int ntmax=200;
		int jmax=20;
		double dw;
		double biatx[] = new double[k+1];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
		int j_it[]= new int[1];
		j_it[0]=1;
		
		left=k-1;
		leftmk=0;
		double step;
		double dataPoint;
		double plotData[]=new double[20];

		step=(maxRange-minRange)/numPoints;
		
	
		for(i=0;i<=numPoints;i++){
		//left=0;
			dataPoint= minRange+i*step;
		while(dataPoint>=t[left+1] && left<n){
			left++;
			leftmk++;
		}
		bsplvb(t,k,1,dataPoint,left, biatx,deltar, deltal, j_it);
		System.out.print("@@\t"+dataPoint);
		for(int ix=0;ix<20;ix++){
			plotData[ix]=0;
		}
		for(int ix=0;ix<k;ix++){
			plotData[leftmk+ix]=biatx[ix];
		}
		double sum;
		sum=0;
		for(int ix=0;ix<n;ix++){
			sum=sum+plotData[ix]*bcoef[ix];
		}
		//System.out.print("3@\t" );
	System.out.println("\t"+ sum);
//		for(int ix=0;ix<20;ix++){
//			System.out.print("\t"+plotData[ix]);
//		}
//		System.out.print("\n" );
	//	 System.out.println("@@\t"+data[0][ll]+"\t"+left+"\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]+"\t"+biatx[3]);
		//System.out.println("\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]);
	
	}	
}

public static void	plotDerivativeData ( double minRange, double maxRange, int numPoints,double[] t,  int n, int k, double[] bcoef){
	int i,j,jj,left,leftmk, ll, mm;
		final int kmax=20;
		final int ntmax=200;
		int jmax=20;
		double dw;
		double biatx[] = new double[k+1];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
		double dbiatx[][] = new double[k][k];
		int j_it[]= new int[1];
		j_it[0]=1;
		
		left=k-1;
		leftmk=0;
		double step;
		double dataPoint;
		double plotData[]=new double[20];
		double a[][] =new double[k][k];
		step=(maxRange-minRange)/numPoints;
		
	    int nderiv=2;
		for(i=0;i<=numPoints;i++){
		//left=0;
			dataPoint= minRange+i*step;
		while(dataPoint>=t[left+1] && left<n){
			left++;
			leftmk++;
		}
		
		bsplvd(t, k,  dataPoint, left, a, dbiatx, nderiv) ;
		//bsplvb(t,k,1,dataPoint,left, biatx,deltar, deltal, j_it);
		//System.out.print("@Deriv@\t"+dataPoint);
		for(int ix=0;ix<20;ix++){
			plotData[ix]=0;
		}
		for(int ix=0;ix<k;ix++){
			plotData[leftmk+ix]=dbiatx[ix][1];
		}
		double sum;
		sum=0;
		for(int ix=0;ix<n;ix++){
			sum=sum+plotData[ix]*bcoef[ix];
		}
		//System.out.print("3@\t" );
		//System.out.println("\t"+ sum);
//		for(int ix=0;ix<20;ix++){
//			System.out.print("\t"+plotData[ix]);
//		}
//		System.out.print("\n" );
	//	 System.out.println("@@\t"+data[0][ll]+"\t"+left+"\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]+"\t"+biatx[3]);
		//System.out.println("\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]);
	
	}	
}

public static void	plotDerivativeData2 ( double minRange, double maxRange, int numPoints,double[] t,  int n, int k, double[] bcoef){
		int i,j,jj,left,leftmk, ll, mm;
		int nderiv=2;
		final int kmax=20;
		final int ntmax=200;
		int jmax=20;
		double dw;
		double biatx[] = new double[k+1];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
		double dbiatx[][] = new double[k][nderiv];
		int j_it[]= new int[1];
		j_it[0]=1;
		
		left=k-1;
		leftmk=0;
		double step;
		double dataPoint;
		double plotData[]=new double[20];
		double a[][] =new double[k][k];
		step=(maxRange-minRange)/numPoints;
		
	
		for(i=0;i<=numPoints;i++){
		//left=0;
			dataPoint= minRange+i*step;
		while(dataPoint>=t[left+1] && left<n){
			left++;
			leftmk++;
		}
		
		bsplvd(t, k,  dataPoint, left, a, dbiatx, nderiv) ;
		//bsplvb(t,k,1,dataPoint,left, biatx,deltar, deltal, j_it);
		System.out.print("@Deriv@\t"+dataPoint);
		for(int ix=0;ix<20;ix++){
			plotData[ix]=0;
		}
		for(int ix=0;ix<k;ix++){
			plotData[leftmk+ix]=dbiatx[ix][1];
		 	System.out.print("\t" + dbiatx[ix][1]);
		}
		// System.out.print("\n" );
		double sum;
		sum=0;
		for(int ix=0;ix<n;ix++){
			sum=sum+plotData[ix]*bcoef[ix];
		}
		//System.out.print("3@\t" );
		System.out.println("\t"+ sum);
//		for(int ix=0;ix<20;ix++){
//			System.out.print("\t"+plotData[ix]);
//		}
//		System.out.print("\n" );
	//	 System.out.println("@@\t"+data[0][ll]+"\t"+left+"\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]+"\t"+biatx[3]);
		//System.out.println("\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]);
	
	}	
}

//	  subroutine l2err ( prfun , ftau , error )
//c  from  * a practical guide to splines *  by c. de boor    
//c  this routine is to be called in the main program  l 2 m a i n .
//calls subprogram  ppvalu(interv)
//c  this subroutine computes various errors of the current l2-approxi-
//c  mation , whose pp-repr. is contained in common block  approx  ,
//c  to the given data contained in common block  data . it prints out
//c  the average error  e r r l 1 , the l2-error  e r r l 2,  and the
//c  maximum error  e r r m a x .
//c
//c******  i n p u t  ******
//c  prfun  a hollerith string.  if prfun = 'ON', the routine prints out
//c          the value of the approximation as well as its error at
//c          every data point.
//c
//c******  o u t p u t  ******
//c  ftau(1), ..., ftau(ntau),  with  ftau(i)  the approximation  f at
//c          tau(i), all i.
//c  error(1), ..., error(ntau),  with  error(i) = scale*(g - f)
//c          at tau(i), all i. here,  s c a l e  equals  1. in case
//c          prfun .ne. 'ON' , or the abs.error is greater than 100 some-
//c          where. otherwise, s c a l e  is such that the maximum of
//c          abs(error))  over all  i  lies between  10  and  100. this
//c          makes the printed output more illustrative.
//c
//	  integer prfun,   ie,k,l,ll,lpkmax,ltkmax,ntau,ntmax,on
//	  real ftau(1),error(1),  break,coef,err,errmax,errl1,errl2
//	 *                       ,gtau,scale,tau,totalw,weight
//	  dimension ftau(ntau),error(ntau)
//	  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
//	  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
//	  common /approx/ break(lpkmax),coef(ltkmax),l,k
//C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
//C     common /approx/ break(100),coef(2000),l,k
//	  data on /'ON'/
//	  errl1 = 0.
//	  errl2 = 0.
//	  errmax = 0.
//	  do 10 ll=1,ntau
//		 ftau(ll) = ppvalu (break, coef, l, k, tau(ll), 0 )
//		 error(ll) = gtau(ll) - ftau(ll)
//		 err = abs(error(ll))
//		 if (errmax .lt. err)   errmax = err
//		 errl1 = errl1 + err*weight(ll)
//   10    errl2 = errl2 + err**2*weight(ll)
//	  errl1 = errl1/totalw
//	  errl2 = sqrt(errl2/totalw)
//	  print 615,errl2,errl1,errmax
//  615 format(///' least square error =',e20.6/
//	 1          ' average error      =',e20.6/
//	 2          ' maximum error      =',e20.6//)
//	  if (prfun .ne. on)                return
//c     **  scale error curve and print  **
//	  ie = 0
//	  scale = 1.
//	  if (errmax .ge. 10.)              go to 18
//	  do 17 ie=1,9
//		 scale = scale*10.
//		 if (errmax*scale .ge. 10.)     go to 18
//   17    continue
//   18 do 19 ll=1,ntau
//   19    error(ll) = error(ll)*scale
//	  print 620,ie,(ll,tau(ll),ftau(ll),error(ll),ll=1,ntau)
//  620 format(///14x,'approximation and scaled error curve'/7x,
//	 1'data point',7x,'approximation',3x,'deviation x 10**',i1/
//	 2(i4,f16.8,f16.8,f17.6))
//										return
//	  end

public static double splineEval( double x,  double[] t, int n, int k, double[] bcoef){
		/*
		 * Inputs:
		 * x is array of points at which the spline is to be evaluated
		 
		 * t is the knot sequence
		 * n is the dimension of the spline space...returned from l2knts
		 * k is the order of the splines plus 1
		 * bcoef is the coeffecients defining the spline 
		 * 
		 * Ouputs:
		 * y is the 
		 * 
		 * Return code
		 * 0 or positive indicates successful run
		 * negative values indicate an error condition
		 * 
		 */
		int i,j,jj,left,leftmk, ll, mm;
		final int kmax=20;
		final int ntmax=200;
		int jmax=20;
		double dw;
		double biatx[] = new double[k+1];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
		int j_it[]= new int[1];
		j_it[0]=1;
		
		left=k-1;
		leftmk=0;
		double step;
		double dataPoint;
		double plotData[]=new double[20];
		double sum;

			dataPoint= x ;
			while(dataPoint>=t[left+1] && left<n){
				left++;
				leftmk++;
			}
			bsplvb(t,k,1,dataPoint,left, biatx,deltar, deltal, j_it);
			//System.out.print("@@\t"+dataPoint);
			
			for(int ix=0;ix<20;ix++){
				plotData[ix]=0;
			}
			
			for(int ix=0;ix<k;ix++){
				plotData[leftmk+ix]=biatx[ix];
			}
			
	
			sum=0;
			for(int ix=0;ix<n;ix++){
				sum=sum+plotData[ix]*bcoef[ix];
			}
			//System.out.println("\t"+ sum);
		
//		}	
		return sum;
}

//	 subroutine l2err ( prfun , ftau , error )
//c  from  * a practical guide to splines *  by c. de boor    
//c  this routine is to be called in the main program  l 2 m a i n .
//calls subprogram  ppvalu(interv)
//c  this subroutine computes various errors of the current l2-approxi-
//c  mation , whose pp-repr. is contained in common block  approx  ,
//c  to the given data contained in common block  data . it prints out
//c  the average error  e r r l 1 , the l2-error  e r r l 2,  and the
//c  maximum error  e r r m a x .
//c
//c******  i n p u t  ******
//c  prfun  a hollerith string.  if prfun = 'ON', the routine prints out
//c          the value of the approximation as well as its error at
//c          every data point.
//c
//c******  o u t p u t  ******
//c  ftau(1), ..., ftau(ntau),  with  ftau(i)  the approximation  f at
//c          tau(i), all i.
//c  error(1), ..., error(ntau),  with  error(i) = scale*(g - f)
//c          at tau(i), all i. here,  s c a l e  equals  1. in case
//c          prfun .ne. 'ON' , or the abs.error is greater than 100 some-
//c          where. otherwise, s c a l e  is such that the maximum of
//c          abs(error))  over all  i  lies between  10  and  100. this
//c          makes the printed output more illustrative.
//c
//	  integer prfun,   ie,k,l,ll,lpkmax,ltkmax,ntau,ntmax,on
//	  real ftau(1),error(1),  break,coef,err,errmax,errl1,errl2
//	 *                       ,gtau,scale,tau,totalw,weight
//	  dimension ftau(ntau),error(ntau)
//	  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
//	  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
//	  common /approx/ break(lpkmax),coef(ltkmax),l,k
//C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
//C     common /approx/ break(100),coef(2000),l,k
//	  data on /'ON'/
//	  errl1 = 0.
//	  errl2 = 0.
//	  errmax = 0.
//	  do 10 ll=1,ntau
//		 ftau(ll) = ppvalu (break, coef, l, k, tau(ll), 0 )
//		 error(ll) = gtau(ll) - ftau(ll)
//		 err = abs(error(ll))
//		 if (errmax .lt. err)   errmax = err
//		 errl1 = errl1 + err*weight(ll)
//   10    errl2 = errl2 + err**2*weight(ll)
//	  errl1 = errl1/totalw
//	  errl2 = sqrt(errl2/totalw)
//	  print 615,errl2,errl1,errmax
//  615 format(///' least square error =',e20.6/
//	 1          ' average error      =',e20.6/
//	 2          ' maximum error      =',e20.6//)
//	  if (prfun .ne. on)                return
//c     **  scale error curve and print  **
//	  ie = 0
//	  scale = 1.
//	  if (errmax .ge. 10.)              go to 18
//	  do 17 ie=1,9
//		 scale = scale*10.
//		 if (errmax*scale .ge. 10.)     go to 18
//   17    continue
//   18 do 19 ll=1,ntau
//   19    error(ll) = error(ll)*scale
//	  print 620,ie,(ll,tau(ll),ftau(ll),error(ll),ll=1,ntau)
//  620 format(///14x,'approximation and scaled error curve'/7x,
//	 1'data point',7x,'approximation',3x,'deviation x 10**',i1/
//	 2(i4,f16.8,f16.8,f17.6))
//										return
//	  end

public static double splineDerivativeEval(int nderiv, double x,  int numPoints, double[] t, int n, int k, double[] bcoef){
		/*
		 * Inputs:
		 * nderiv is the order of the derviative to be returned
		 * x is point which the spline is to be evaluated
		 * numPoints is the number of points in x
		 * t is the knot sequence
		 * n is the dimension of the spline space...returned from l2knts
		 * k is the order of the splines plus 1
		 * bcoef is the coeffecients defining the spline 
		 * 
		 * Ouputs:
		 * method returns value of indicated derivative (nderiv) at x
		 * 
		 */
		 
		int i,j,jj,left,leftmk, ll, mm;
		final int kmax=20;
		final int ntmax=200;
		int jmax=20;
		double dw;
		double biatx[] = new double[k+1];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
	    double dbiatx[][] = new double[k][nderiv];
		double a[][] =new double[k][k];
		
		int j_it[]= new int[1];
		j_it[0]=1;
		
		left=k-1;
		leftmk=0;
		double step;
		double dataPoint;
		double plotData[]=new double[20];
		double sum;
		
		dataPoint= x ;
		while(dataPoint>=t[left+1] && left<n){
			left++;
			leftmk++;
		}			 
		
		bsplvd(t, k,  dataPoint, left, a, dbiatx, nderiv) ;
		//System.out.print("@Deriv@\t"+dataPoint);
		for(int ix=0;ix<20;ix++){
			plotData[ix]=0;
		}
		for(int ix=0;ix<k;ix++){
			plotData[leftmk+ix]=dbiatx[ix][nderiv-1];
		//	System.out.print("\t" + dbiatx[ix][nderiv-1]);
		}
		
		sum=0;
		for(int ix=0;ix<n;ix++){
			sum=sum + plotData[ix]*bcoef[ix];
		}

		//System.out.println("\t"+ sum);
	
		for(int ix=0;ix<k;ix++){
			plotData[leftmk+ix] = biatx[ix];
		}
		
		return sum;
}
	public static void	l2knts (double[] brake, int l, int k, double[] t, int n[] ){
	
//	  subroutine l2knts ( break, l, k, t, n )
//c  from  * a practical guide to splines *  by c. de boor    
//c  to be called in main program  l 2 m a i n .
//converts the breakpoint sequence  b r e a k   into a corresponding knot
//c  sequence  t  to allow the repr. of a pp function of order  k  with
//c  k-2 continuous derivatives as a spline of order  k  with knot
//c  sequence  t . this means that
//c  t(1), ..., t(n+k) =  break(1) k times, then break(i), i=2,...,l, each
//c                       once, then break(l+1) k times .
//c  therefore,  n = k-1 + l.
//c
//c******  i n p u t  ******
//c  k     the order
//c  l     the number of polynomial pieces
//c  break(1), ...,break(l+1)  the breakpoint sequence
//c
//c******  o u t p u t  ******
//c  t(1),...,t(n+k)   the knot sequence
//c  n     the dimension of the corresp. spline space of order  k .
//c
//	  integer k,l,n,   i,km1
		int i, kml;
//	  real break(1),t(1)
//c     dimension break(l+1),t(n+k)
//	  km1 = k - 1
	  kml=k - 1;
//	  do 5 i=1,km1
//	5    t(i) = break(1)
//	  do 6 i=1,l
//	6    t(km1+i) = break(i)
//	  n = km1 + l
		n[0]=(kml + l);
//	  do 7 i=1,k
//	7    t(n+i) = break(l+1)
//										return
//	  end
		for(i=0;i<k;i++){
			t[i]=brake[0];
		}
		for(i=0;i<l;i++){
			t[kml+i+1]=brake[i];
		}
		for(i=0;i<k;i++){
			t[n[0] +i+1]=brake[l-1];
		}
		//n=Integer.valueOf(5);
		return;
	}
	


	public static void	l5main ( double[][] data, int ntau2, double[] weight, double[] t){	
		
//	c  main program for least-squares approximation by splines
//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
//	c
//	c  the program, though ostensibly written for l2-approximation, is typ-
//	c  ical for programs constructing a pp approximation to a function gi-
//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
//	c  easily be replaced by one carrying out interpolation or some other
//	c  form of approximation.
//	c
//	c******  i n p u t  ******
//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
//	c  be approximated and the order and breakpoint sequence of the pp ap-
//	c  proximating function to be used. further,  s e t d a t  is expected
//	c  to  t e r m i n a t e  the run (for lack of further input or because
//	c   i c o u n t  has reached a critical value).
//	c     the number  n t i m e s  is read in in the main program. it speci
//	c  fies the number of passes through the knot improvement algorithm in
//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
//	c  that, on the average, addbrk knots are to be added per pass through
//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
//	c  every third pass (as long as  ntimes .lt. 50).
//	c
//	c******  p r i n t e d  o u t p u t  ******
//	c  is governed by the three print control hollerith strings
//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
//	c  p r f u n  = 'on'  gives printout of approximation and error at
//	c                     every data point.
//	c  the order  k , the number of pieces  l, and the interior breakpoints
//	c  are always printed out as are (in l2err) the mean, mean square, and
//	c  maximum errors in the approximation.
//	c
//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
		int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
		int n[]= new int[1];
//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
		final int lpkmax=100;
		final int ntmax=200;
		final int ltkmax=2000;
//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
//		 *    ,t(ntmax),tau,totalw,weight
		double addbrk,  gtau;
		double bcoef[] = new double[lpkmax];
		double q[][]	= new double[10][10];
		double scrtch[] = new double[ntmax];
		
		
//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
//	C    *    ,t(200),tau,totalw,weight
//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
//	c     ioned here only because it might otherwise become undefined be-
//	c     tween calls to those subroutines.
//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
		double brake[]=new double[lpkmax];
		double coef[]=new double[ltkmax];
//	C     common /approx/ break(100),coef(2000),l,k
//	c     common /approx/ also occurs in setdat and l2err.
//		  data on /'ON'/
//	c
//		  icount = 0
		icount=0;
		
//	c        i c o u n t  provides communication with the data-input-and-
//	c     termination routine  s e t d a t . it is initialized to  0  to
//	c     signal to setdat when it is being called for the first time. after
//	c     that, it is up to setdat to use icount for keeping track of the
//	c     passes through setdat .
//	c
//	c     information about the function to be approximated and order and
//	c     breakpoint sequence of the approximating pp functions is gathered
//	c     by a
//		1 call setdat(icount)

		ntau = 9;
		int ntaum1 = ntau-1;
		data[0][0] =0.0;
		data[0][1] =30.0;
		data[0][2] =45.0;
		data[0][3] =60.0;
		data[0][4] =90.0;
		data[0][5] =120.0;
		data[0][6] =135.0;
		data[0][7] =150.0;
		data[0][8] =180.0;

		data[1][0] =0.0;
		data[1][1] =0.5;
		data[1][2] =.707;
		data[1][3] =.866;
		data[1][4] =1.0;
		data[1][5] =.866;
		data[1][6] =.707;
		data[1][7] =.5;
		data[1][8] =0.0;
		
		weight[0] = 2.;
		weight[1] = 2.;
		weight[2] = 1.;
		weight[3] = 1.;
		weight[4] = 2.;
		weight[5] = 1.;
		weight[6] = 1.;
		weight[7] = 2.;
		weight[8] = 2.;
				
//		for (i=0;i< ntau;i++){
////			data[0][i] = 1. - Math.pow(.5,(i));
////			data[1][i] = Math.pow(data[0][i],2)  + 1.;
//			weight[i] = 1.;
			
//			System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
//		}
	     

		int totalw = ntau;
		l = 5;
		int lp1 = l+1;
		double step = 1./((double)l);
		k = 4;
//		for(i=0;i<lp1;i++){
//			brake[i]=i*step;
//			System.out.println("\t"+brake[i]);
//		}
		
		brake[0]=-1;
		brake[1]=50.0;
		brake[2]=75.0;
		brake[3]=145.0;
		brake[4]=181.;





//	c
//	c     breakpoints are translated into knots, and the number  n  of
//	c     b-splines to be used is obtained by a
//		  call l2knts ( break, l, k, t, n )
			n[0]=0;
			l2knts(brake,l,k-1,t,n);
			System.out.println("n\t"+n[0]);
//	c
//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
//	c     as well as the print controls  p r b c o ,  p r p c o  and
//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
//	c     not, with an increase of  addbrk  knots for every pass .
//		  print 600
//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
//		  read 500,ntimes,addbrk,prbco,prpco,prfun
//	  500 format(i3,f10.5/3a2)
//	c
//		  lbegin = l
//		  nt = 0
//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
//	c        ed by a
//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
//	  609    format(//' b-spline coefficients'/(4e20.10))
lbegin=1;
nt=0;
//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){

l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
System.out.println(ntau);
for(int ix=0; ix<ntau; ix++){
	System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
}
for(int ix=0; ix < n[0]; ix++){
	System.out.println(t[ix]);
}

for(int ix=0; ix < n[0]; ix++){
	System.out.println(t[ix]);
}
//System.out.println(n[0]);
//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]);
//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]);
  plotdata ( 0., 180.0, 18,t,n[0], k,bcoef);
  plotDerivativeData(0., 180.0, 18,t,n[0], k, bcoef);

//	c
//	c        convert the b-repr. of the approximation to pp repr.
//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
//			 print 610, k, l, (break(ll),ll=2,l)
//	  610    format(//' approximation by splines of order',i3,' on ',
//		 *         i3,' intervals. breakpoints -'/(4e20.10))
//			 if (prpco .ne. on)             go to 15
//			 print 611
//	  611    format(/' pp-representation for approximation')
//			 do 12 i=1,l
//				ii = (i-1)*k
//	   12       print 613,break(i),(coef(ii+j),j=1,k)
//	  613    format(f9.3,4e20.10/(11x,4e20.10))
//	c
//	c        compute and print out various error norms by a
//	   15    call l2err ( prfun, scrtch, q )
//	c
//	c        if newnot has been applied less than  n t i m e s  times, try
//	c        it again to obtain, from the current approx. a possibly improv-
//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
//	c        the average) than the current approximation has.
//	c           if only an increase in breakpoints is wanted, without the
//	c        adjustment that newnot provides, a fake newnot routine could be
//	c        used here which merely returns the breakpoints for  l n e w
//	c        equal intervals .
//			 if (nt .ge. ntimes)            go to 1
//			 lnew = lbegin + float(nt)*addbrk
//			 call newnot (break, coef, l, k, scrtch, lnew, t )
//			 call l2knts ( scrtch, lnew, k, t, n )
//			 nt = nt + 1
//											go to 10
//		  end
	}

		public static void	l6main ( double[][] data, int ntau2, double[] weight, double[] t){	
			
	//	c  main program for least-squares approximation by splines
	//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
	//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
	//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
	//	c
	//	c  the program, though ostensibly written for l2-approximation, is typ-
	//	c  ical for programs constructing a pp approximation to a function gi-
	//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
	//	c  easily be replaced by one carrying out interpolation or some other
	//	c  form of approximation.
	//	c
	//	c******  i n p u t  ******
	//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
	//	c  be approximated and the order and breakpoint sequence of the pp ap-
	//	c  proximating function to be used. further,  s e t d a t  is expected
	//	c  to  t e r m i n a t e  the run (for lack of further input or because
	//	c   i c o u n t  has reached a critical value).
	//	c     the number  n t i m e s  is read in in the main program. it speci
	//	c  fies the number of passes through the knot improvement algorithm in
	//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
	//	c  that, on the average, addbrk knots are to be added per pass through
	//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
	//	c  every third pass (as long as  ntimes .lt. 50).
	//	c
	//	c******  p r i n t e d  o u t p u t  ******
	//	c  is governed by the three print control hollerith strings
	//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
	//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
	//	c  p r f u n  = 'on'  gives printout of approximation and error at
	//	c                     every data point.
	//	c  the order  k , the number of pieces  l, and the interior breakpoints
	//	c  are always printed out as are (in l2err) the mean, mean square, and
	//	c  maximum errors in the approximation.
	//	c
	//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
	//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
			int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
			int n[]= new int[1];
	//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
			final int lpkmax=100;
			final int ntmax=200;
			final int ltkmax=2000;
	//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
	//		 *    ,t(ntmax),tau,totalw,weight
			double addbrk,  gtau;
			double bcoef[] = new double[lpkmax];
			double q[][]	= new double[10][10];
			double scrtch[] = new double[ntmax];
			
			
	//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
	//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
	//	C    *    ,t(200),tau,totalw,weight
	//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
	//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
	//	c     ioned here only because it might otherwise become undefined be-
	//	c     tween calls to those subroutines.
	//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
			double brake[]=new double[lpkmax];
			double coef[]=new double[ltkmax];
	//	C     common /approx/ break(100),coef(2000),l,k
	//	c     common /approx/ also occurs in setdat and l2err.
	//		  data on /'ON'/
	//	c
	//		  icount = 0
			icount=0;
			
	//	c        i c o u n t  provides communication with the data-input-and-
	//	c     termination routine  s e t d a t . it is initialized to  0  to
	//	c     signal to setdat when it is being called for the first time. after
	//	c     that, it is up to setdat to use icount for keeping track of the
	//	c     passes through setdat .
	//	c
	//	c     information about the function to be approximated and order and
	//	c     breakpoint sequence of the approximating pp functions is gathered
	//	c     by a
	//		1 call setdat(icount)
	
			ntau = 9;
			int ntaum1 = ntau-1;
			data[0][0] =0.0;
			data[0][1] =30.0;
			data[0][2] =45.0;
			data[0][3] =60.0;
			data[0][4] =90.0;
			data[0][5] =120.0;
			data[0][6] =135.0;
			data[0][7] =150.0;
			data[0][8] =180.0;
	
			data[1][0] =0.0;
			data[1][1] =0.5;
			data[1][2] =.707;
			data[1][3] =.866;
			data[1][4] =1.0;
			data[1][5] =.866;
			data[1][6] =.707;
			data[1][7] =.5;
			data[1][8] =0.0;
			
			weight[0] = 2.;
			weight[1] = 2.;
			weight[2] = 1.;
			weight[3] = 1.;
			weight[4] = 2.;
			weight[5] = 1.;
			weight[6] = 1.;
			weight[7] = 2.;
			weight[8] = 2.;
					
	//		for (i=0;i< ntau;i++){
	////			data[0][i] = 1. - Math.pow(.5,(i));
	////			data[1][i] = Math.pow(data[0][i],2)  + 1.;
	//			weight[i] = 1.;
				
	//			System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
	//		}
		     
	
			int totalw = ntau;
			l = 5;
			int lp1 = l+1;
			double step = 1./((double)l);
			k = 4;
	//		for(i=0;i<lp1;i++){
	//			brake[i]=i*step;
	//			System.out.println("\t"+brake[i]);
	//		}
			
			brake[0]=-1;
			brake[1]=50.0;
			brake[2]=75.0;
			brake[3]=145.0;
			brake[4]=181.;
	
	
	
	
	
	//	c
	//	c     breakpoints are translated into knots, and the number  n  of
	//	c     b-splines to be used is obtained by a
	//		  call l2knts ( break, l, k, t, n )
				n[0]=0;
				l2knts(brake,l,k-1,t,n);
				System.out.println("n\t"+n[0]);
	//	c
	//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
	//	c     as well as the print controls  p r b c o ,  p r p c o  and
	//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
	//	c     not, with an increase of  addbrk  knots for every pass .
	//		  print 600
	//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
	//		  read 500,ntimes,addbrk,prbco,prpco,prfun
	//	  500 format(i3,f10.5/3a2)
	//	c
	//		  lbegin = l
	//		  nt = 0
	//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
	//	c        ed by a
	//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
	//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
	//	  609    format(//' b-spline coefficients'/(4e20.10))
	lbegin=1;
	nt=0;
	//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){
	
	l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
	System.out.println(ntau);
	for(int ix=0; ix<ntau; ix++){
		System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
	}
	for(int ix=0; ix < n[0]; ix++){
		System.out.println(t[ix]);
	}
	
	for(int ix=0; ix < n[0]; ix++){
		System.out.println(t[ix]);
	}
	//System.out.println(n[0]);
	//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]);
	//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]);
	  plotdata ( 0., 180.0, 100,t,n[0], k,bcoef);
	  plotDerivativeData2(0., 180.0, 100,t,n[0], k, bcoef);
	
	//	c
	//	c        convert the b-repr. of the approximation to pp repr.
	//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
	//			 print 610, k, l, (break(ll),ll=2,l)
	//	  610    format(//' approximation by splines of order',i3,' on ',
	//		 *         i3,' intervals. breakpoints -'/(4e20.10))
	//			 if (prpco .ne. on)             go to 15
	//			 print 611
	//	  611    format(/' pp-representation for approximation')
	//			 do 12 i=1,l
	//				ii = (i-1)*k
	//	   12       print 613,break(i),(coef(ii+j),j=1,k)
	//	  613    format(f9.3,4e20.10/(11x,4e20.10))
	//	c
	//	c        compute and print out various error norms by a
	//	   15    call l2err ( prfun, scrtch, q )
	//	c
	//	c        if newnot has been applied less than  n t i m e s  times, try
	//	c        it again to obtain, from the current approx. a possibly improv-
	//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
	//	c        the average) than the current approximation has.
	//	c           if only an increase in breakpoints is wanted, without the
	//	c        adjustment that newnot provides, a fake newnot routine could be
	//	c        used here which merely returns the breakpoints for  l n e w
	//	c        equal intervals .
	//			 if (nt .ge. ntimes)            go to 1
	//			 lnew = lbegin + float(nt)*addbrk
	//			 call newnot (break, coef, l, k, scrtch, lnew, t )
	//			 call l2knts ( scrtch, lnew, k, t, n )
	//			 nt = nt + 1
	//											go to 10
	//		  end
		}

			public static void	l7main ( double[][] data, int ntau2, double[] weight, double[] t){	
				
		//	c  main program for least-squares approximation by splines
		//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
		//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
		//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
		//	c
		//	c  the program, though ostensibly written for l2-approximation, is typ-
		//	c  ical for programs constructing a pp approximation to a function gi-
		//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
		//	c  easily be replaced by one carrying out interpolation or some other
		//	c  form of approximation.
		//	c
		//	c******  i n p u t  ******
		//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
		//	c  be approximated and the order and breakpoint sequence of the pp ap-
		//	c  proximating function to be used. further,  s e t d a t  is expected
		//	c  to  t e r m i n a t e  the run (for lack of further input or because
		//	c   i c o u n t  has reached a critical value).
		//	c     the number  n t i m e s  is read in in the main program. it speci
		//	c  fies the number of passes through the knot improvement algorithm in
		//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
		//	c  that, on the average, addbrk knots are to be added per pass through
		//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
		//	c  every third pass (as long as  ntimes .lt. 50).
		//	c
		//	c******  p r i n t e d  o u t p u t  ******
		//	c  is governed by the three print control hollerith strings
		//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
		//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
		//	c  p r f u n  = 'on'  gives printout of approximation and error at
		//	c                     every data point.
		//	c  the order  k , the number of pieces  l, and the interior breakpoints
		//	c  are always printed out as are (in l2err) the mean, mean square, and
		//	c  maximum errors in the approximation.
		//	c
		//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
		//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
				int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
				int n[]= new int[1];
		//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
				final int lpkmax=100;
				final int ntmax=200;
				final int ltkmax=2000;
		//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
		//		 *    ,t(ntmax),tau,totalw,weight
				double addbrk,  gtau;
				double bcoef[] = new double[lpkmax];
				double q[][]	= new double[10][10];
				double scrtch[] = new double[ntmax];
				
				
		//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
		//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
		//	C    *    ,t(200),tau,totalw,weight
		//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
		//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
		//	c     ioned here only because it might otherwise become undefined be-
		//	c     tween calls to those subroutines.
		//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
				double brake[]=new double[lpkmax];
				double coef[]=new double[ltkmax];
		//	C     common /approx/ break(100),coef(2000),l,k
		//	c     common /approx/ also occurs in setdat and l2err.
		//		  data on /'ON'/
		//	c
		//		  icount = 0
				icount=0;
				
		//	c        i c o u n t  provides communication with the data-input-and-
		//	c     termination routine  s e t d a t . it is initialized to  0  to
		//	c     signal to setdat when it is being called for the first time. after
		//	c     that, it is up to setdat to use icount for keeping track of the
		//	c     passes through setdat .
		//	c
		//	c     information about the function to be approximated and order and
		//	c     breakpoint sequence of the approximating pp functions is gathered
		//	c     by a
		//		1 call setdat(icount)
		
				ntau = 9;
				int ntaum1 = ntau-1;
				data[0][0] =0.0;
				data[0][1] =30.0;
				data[0][2] =45.0;
				data[0][3] =60.0;
				data[0][4] =90.0;
				data[0][5] =120.0;
				data[0][6] =135.0;
				data[0][7] =150.0;
				data[0][8] =180.0;
		
				data[1][0] =0.0;
				data[1][1] =0.5;
				data[1][2] =.707;
				data[1][3] =.866;
				data[1][4] =1.0;
				data[1][5] =.866;
				data[1][6] =.707;
				data[1][7] =.5;
				data[1][8] =0.0;
				
				weight[0] = 2.;
				weight[1] = 2.;
				weight[2] = 1.;
				weight[3] = 1.;
				weight[4] = 2.;
				weight[5] = 1.;
				weight[6] = 1.;
				weight[7] = 2.;
				weight[8] = 2.;
						
		//		for (i=0;i< ntau;i++){
		////			data[0][i] = 1. - Math.pow(.5,(i));
		////			data[1][i] = Math.pow(data[0][i],2)  + 1.;
		//			weight[i] = 1.;
					
		//			System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
		//		}
			     
		
				int totalw = ntau;
				l = 2;
				int lp1 = l+1;
				double step = 1./((double)l);
				k = 4;
		//		for(i=0;i<lp1;i++){
		//			brake[i]=i*step;
		//			System.out.println("\t"+brake[i]);
		//		}
				
				brake[0]=-1;
				brake[1]=181.0;

		
		
		
		
		
		//	c
		//	c     breakpoints are translated into knots, and the number  n  of
		//	c     b-splines to be used is obtained by a
		//		  call l2knts ( break, l, k, t, n )
					n[0]=0;
					l2knts(brake,l,k-1,t,n);
					System.out.println("n\t"+n[0]);
		//	c
		//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
		//	c     as well as the print controls  p r b c o ,  p r p c o  and
		//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
		//	c     not, with an increase of  addbrk  knots for every pass .
		//		  print 600
		//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
		//		  read 500,ntimes,addbrk,prbco,prpco,prfun
		//	  500 format(i3,f10.5/3a2)
		//	c
		//		  lbegin = l
		//		  nt = 0
		//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
		//	c        ed by a
		//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
		//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
		//	  609    format(//' b-spline coefficients'/(4e20.10))
		lbegin=1;
		nt=0;
		//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){
		
		l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
		System.out.println(ntau);
		for(int ix=0; ix<ntau; ix++){
			System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
		}
		for(int ix=0; ix < n[0]; ix++){
			System.out.println(t[ix]);
		}
		
		for(int ix=0; ix < n[0]; ix++){
			System.out.println(t[ix]);
		}
		//System.out.println(n[0]);
		//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]);
		//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]);
		  plotdata ( 0., 180.0, 100,t,n[0], k,bcoef);
		  plotDerivativeData2(0., 180.0, 100,t,n[0], k, bcoef);
		
		//	c
		//	c        convert the b-repr. of the approximation to pp repr.
		//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
		//			 print 610, k, l, (break(ll),ll=2,l)
		//	  610    format(//' approximation by splines of order',i3,' on ',
		//		 *         i3,' intervals. breakpoints -'/(4e20.10))
		//			 if (prpco .ne. on)             go to 15
		//			 print 611
		//	  611    format(/' pp-representation for approximation')
		//			 do 12 i=1,l
		//				ii = (i-1)*k
		//	   12       print 613,break(i),(coef(ii+j),j=1,k)
		//	  613    format(f9.3,4e20.10/(11x,4e20.10))
		//	c
		//	c        compute and print out various error norms by a
		//	   15    call l2err ( prfun, scrtch, q )
		//	c
		//	c        if newnot has been applied less than  n t i m e s  times, try
		//	c        it again to obtain, from the current approx. a possibly improv-
		//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
		//	c        the average) than the current approximation has.
		//	c           if only an increase in breakpoints is wanted, without the
		//	c        adjustment that newnot provides, a fake newnot routine could be
		//	c        used here which merely returns the breakpoints for  l n e w
		//	c        equal intervals .
		//			 if (nt .ge. ntimes)            go to 1
		//			 lnew = lbegin + float(nt)*addbrk
		//			 call newnot (break, coef, l, k, scrtch, lnew, t )
		//			 call l2knts ( scrtch, lnew, k, t, n )
		//			 nt = nt + 1
		//											go to 10
		//		  end
			}

		public static void	l2a_main ( double[][] data, int ntau2, double[] weight, double[] t){	
			
	//	c  main program for least-squares approximation by splines
	//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
	//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
	//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
	//	c
	//	c  the program, though ostensibly written for l2-approximation, is typ-
	//	c  ical for programs constructing a pp approximation to a function gi-
	//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
	//	c  easily be replaced by one carrying out interpolation or some other
	//	c  form of approximation.
	//	c
	//	c******  i n p u t  ******
	//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
	//	c  be approximated and the order and breakpoint sequence of the pp ap-
	//	c  proximating function to be used. further,  s e t d a t  is expected
	//	c  to  t e r m i n a t e  the run (for lack of further input or because
	//	c   i c o u n t  has reached a critical value).
	//	c     the number  n t i m e s  is read in in the main program. it speci
	//	c  fies the number of passes through the knot improvement algorithm in
	//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
	//	c  that, on the average, addbrk knots are to be added per pass through
	//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
	//	c  every third pass (as long as  ntimes .lt. 50).
	//	c
	//	c******  p r i n t e d  o u t p u t  ******
	//	c  is governed by the three print control hollerith strings
	//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
	//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
	//	c  p r f u n  = 'on'  gives printout of approximation and error at
	//	c                     every data point.
	//	c  the order  k , the number of pieces  l, and the interior breakpoints
	//	c  are always printed out as are (in l2err) the mean, mean square, and
	//	c  maximum errors in the approximation.
	//	c
	//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
	//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
			int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
			int n[]= new int[1];
	//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
			final int lpkmax=100;
			final int ntmax=200;
			final int ltkmax=2000;
	//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
	//		 *    ,t(ntmax),tau,totalw,weight
			double addbrk,  gtau;
			double bcoef[] = new double[lpkmax];
			double q[][]	= new double[10][10];
			double scrtch[] = new double[ntmax];
			
			
	//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
	//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
	//	C    *    ,t(200),tau,totalw,weight
	//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
	//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
	//	c     ioned here only because it might otherwise become undefined be-
	//	c     tween calls to those subroutines.
	//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
			double brake[]=new double[lpkmax];
			double coef[]=new double[ltkmax];
	//	C     common /approx/ break(100),coef(2000),l,k
	//	c     common /approx/ also occurs in setdat and l2err.
	//		  data on /'ON'/
	//	c
	//		  icount = 0
			icount=0;
			
	//	c        i c o u n t  provides communication with the data-input-and-
	//	c     termination routine  s e t d a t . it is initialized to  0  to
	//	c     signal to setdat when it is being called for the first time. after
	//	c     that, it is up to setdat to use icount for keeping track of the
	//	c     passes through setdat .
	//	c
	//	c     information about the function to be approximated and order and
	//	c     breakpoint sequence of the approximating pp functions is gathered
	//	c     by a
	//		1 call setdat(icount)
	
			ntau = 8;
			int ntaum1 = ntau-1;
			data[0][0] =.1;
			data[0][1] =.2;
			data[0][2] =.3;
			data[0][3] =.4;
			data[0][4] =.5;
			data[0][5] =.6;
			data[0][6] =.7;
			data[0][7] =.8;
	
			data[1][0] =.1;
			data[1][1] =.35;
			data[1][2] =.45;
			data[1][3] =.5;
			data[1][4] =.4;
			data[1][5] =.35;
			data[1][6] =.25;
			data[1][7] =.1;
			
			for (i=0;i< ntau;i++){
	//			data[0][i] = 1. - Math.pow(.5,(i));
	//			data[1][i] = Math.pow(data[0][i],2)  + 1.;
				weight[i] = 1.;
				
				System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
			}
		     
	
			int totalw = ntau;
			l = 5;
			int lp1 = l+1;
			double step = 1./((double)l);
			k = 3;
	//		for(i=0;i<lp1;i++){
	//			brake[i]=i*step;
	//			System.out.println("\t"+brake[i]);
	//		}
			
			brake[0]=0;
			brake[1]=0.25;
			brake[2]=.5;
			brake[3]=.75;
			brake[4]=1;
	
	
	
	
	
	//	c
	//	c     breakpoints are translated into knots, and the number  n  of
	//	c     b-splines to be used is obtained by a
	//		  call l2knts ( break, l, k, t, n )
				n[0]=0;
				l2knts(brake,l,k-1,t,n);
				System.out.println("n\t"+n[0]);
	//	c
	//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
	//	c     as well as the print controls  p r b c o ,  p r p c o  and
	//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
	//	c     not, with an increase of  addbrk  knots for every pass .
	//		  print 600
	//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
	//		  read 500,ntimes,addbrk,prbco,prpco,prfun
	//	  500 format(i3,f10.5/3a2)
	//	c
	//		  lbegin = l
	//		  nt = 0
	//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
	//	c        ed by a
	//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
	//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
	//	  609    format(//' b-spline coefficients'/(4e20.10))
	lbegin=1;
	nt=0;
	//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){
	
	l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
	plotdata ( 0., 1.0, 18,t,n[0], k,bcoef);
	plotDerivativeData2(0., 1.0, 18,t,n[0], k, bcoef);
//	System.out.println(ntau);
//	for(int ix=0; ix<ntau; ix++){
//		System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
//	}
//	for(int ix=0; ix < n[0]; ix++){
//		System.out.println(t[ix]);
//	}
//	
//	for(int ix=0; ix < n[0]; ix++){
//		System.out.println(t[ix]);
//	}
//	System.out.println(n[0]);
//	System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]);
	//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]);
	
	
	//	c
	//	c        convert the b-repr. of the approximation to pp repr.
	//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
	//			 print 610, k, l, (break(ll),ll=2,l)
	//	  610    format(//' approximation by splines of order',i3,' on ',
	//		 *         i3,' intervals. breakpoints -'/(4e20.10))
	//			 if (prpco .ne. on)             go to 15
	//			 print 611
	//	  611    format(/' pp-representation for approximation')
	//			 do 12 i=1,l
	//				ii = (i-1)*k
	//	   12       print 613,break(i),(coef(ii+j),j=1,k)
	//	  613    format(f9.3,4e20.10/(11x,4e20.10))
	//	c
	//	c        compute and print out various error norms by a
	//	   15    call l2err ( prfun, scrtch, q )
	//	c
	//	c        if newnot has been applied less than  n t i m e s  times, try
	//	c        it again to obtain, from the current approx. a possibly improv-
	//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
	//	c        the average) than the current approximation has.
	//	c           if only an increase in breakpoints is wanted, without the
	//	c        adjustment that newnot provides, a fake newnot routine could be
	//	c        used here which merely returns the breakpoints for  l n e w
	//	c        equal intervals .
	//			 if (nt .ge. ntimes)            go to 1
	//			 lnew = lbegin + float(nt)*addbrk
	//			 call newnot (break, coef, l, k, scrtch, lnew, t )
	//			 call l2knts ( scrtch, lnew, k, t, n )
	//			 nt = nt + 1
	//											go to 10
	//		  end
		}

			public static void	l2main ( double[][] data, int ntau2, double[] weight, double[] t){	
				
		//	c  main program for least-squares approximation by splines
		//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
		//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
		//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
		//	c
		//	c  the program, though ostensibly written for l2-approximation, is typ-
		//	c  ical for programs constructing a pp approximation to a function gi-
		//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
		//	c  easily be replaced by one carrying out interpolation or some other
		//	c  form of approximation.
		//	c
		//	c******  i n p u t  ******
		//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
		//	c  be approximated and the order and breakpoint sequence of the pp ap-
		//	c  proximating function to be used. further,  s e t d a t  is expected
		//	c  to  t e r m i n a t e  the run (for lack of further input or because
		//	c   i c o u n t  has reached a critical value).
		//	c     the number  n t i m e s  is read in in the main program. it speci
		//	c  fies the number of passes through the knot improvement algorithm in
		//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
		//	c  that, on the average, addbrk knots are to be added per pass through
		//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
		//	c  every third pass (as long as  ntimes .lt. 50).
		//	c
		//	c******  p r i n t e d  o u t p u t  ******
		//	c  is governed by the three print control hollerith strings
		//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
		//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
		//	c  p r f u n  = 'on'  gives printout of approximation and error at
		//	c                     every data point.
		//	c  the order  k , the number of pieces  l, and the interior breakpoints
		//	c  are always printed out as are (in l2err) the mean, mean square, and
		//	c  maximum errors in the approximation.
		//	c
		//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
		//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
				int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
				int n[]= new int[1];
		//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
				final int lpkmax=100;
				final int ntmax=200;
				final int ltkmax=2000;
		//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
		//		 *    ,t(ntmax),tau,totalw,weight
				double addbrk,  gtau;
				double bcoef[] = new double[lpkmax];
				double q[][]	= new double[10][10];
				double scrtch[] = new double[ntmax];
				
				
		//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
		//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
		//	C    *    ,t(200),tau,totalw,weight
		//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
		//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
		//	c     ioned here only because it might otherwise become undefined be-
		//	c     tween calls to those subroutines.
		//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
				double brake[]=new double[lpkmax];
				double coef[]=new double[ltkmax];
		//	C     common /approx/ break(100),coef(2000),l,k
		//	c     common /approx/ also occurs in setdat and l2err.
		//		  data on /'ON'/
		//	c
		//		  icount = 0
				icount=0;
				
		//	c        i c o u n t  provides communication with the data-input-and-
		//	c     termination routine  s e t d a t . it is initialized to  0  to
		//	c     signal to setdat when it is being called for the first time. after
		//	c     that, it is up to setdat to use icount for keeping track of the
		//	c     passes through setdat .
		//	c
		//	c     information about the function to be approximated and order and
		//	c     breakpoint sequence of the approximating pp functions is gathered
		//	c     by a
		//		1 call setdat(icount)
		
				ntau = 8;
				int ntaum1 = ntau-1;
				data[0][0] =.1;
				data[0][1] =.2;
				data[0][2] =.3;
				data[0][3] =.4;
				data[0][4] =.5;
				data[0][5] =.6;
				data[0][6] =.7;
				data[0][7] =.8;
		
				data[1][0] =.1;
				data[1][1] =.35;
				data[1][2] =.45;
				data[1][3] =.5;
				data[1][4] =.4;
				data[1][5] =.35;
				data[1][6] =.25;
				data[1][7] =.1;
				
				for (i=0;i< ntau;i++){
		//			data[0][i] = 1. - Math.pow(.5,(i));
		//			data[1][i] = Math.pow(data[0][i],2)  + 1.;
					weight[i] = 1.;
					
					System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
				}
			     
		
				int totalw = ntau;
				l = 5;
				int lp1 = l+1;
				double step = 1./((double)l);
				k = 3;
		//		for(i=0;i<lp1;i++){
		//			brake[i]=i*step;
		//			System.out.println("\t"+brake[i]);
		//		}
				
				brake[0]=0;
				brake[1]=0.25;
				brake[2]=.5;
				brake[3]=.75;
				brake[4]=1;
		
		
		
		
		
		//	c
		//	c     breakpoints are translated into knots, and the number  n  of
		//	c     b-splines to be used is obtained by a
		//		  call l2knts ( break, l, k, t, n )
					n[0]=0;
					l2knts(brake,l,k-1,t,n);
					System.out.println("n\t"+n[0]);
		//	c
		//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
		//	c     as well as the print controls  p r b c o ,  p r p c o  and
		//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
		//	c     not, with an increase of  addbrk  knots for every pass .
		//		  print 600
		//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
		//		  read 500,ntimes,addbrk,prbco,prpco,prfun
		//	  500 format(i3,f10.5/3a2)
		//	c
		//		  lbegin = l
		//		  nt = 0
		//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
		//	c        ed by a
		//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
		//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
		//	  609    format(//' b-spline coefficients'/(4e20.10))
		lbegin=1;
		nt=0;
		//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){
		
		l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
		System.out.println(ntau);
		for(int ix=0; ix<ntau; ix++){
			System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
		}
		for(int ix=0; ix < n[0]; ix++){
			System.out.println(t[ix]);
		}
		
		for(int ix=0; ix < n[0]; ix++){
			System.out.println(t[ix]);
		}
		System.out.println(n[0]);
		System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]);
		//System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]);
		
		
		//	c
		//	c        convert the b-repr. of the approximation to pp repr.
		//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
		//			 print 610, k, l, (break(ll),ll=2,l)
		//	  610    format(//' approximation by splines of order',i3,' on ',
		//		 *         i3,' intervals. breakpoints -'/(4e20.10))
		//			 if (prpco .ne. on)             go to 15
		//			 print 611
		//	  611    format(/' pp-representation for approximation')
		//			 do 12 i=1,l
		//				ii = (i-1)*k
		//	   12       print 613,break(i),(coef(ii+j),j=1,k)
		//	  613    format(f9.3,4e20.10/(11x,4e20.10))
		//	c
		//	c        compute and print out various error norms by a
		//	   15    call l2err ( prfun, scrtch, q )
		//	c
		//	c        if newnot has been applied less than  n t i m e s  times, try
		//	c        it again to obtain, from the current approx. a possibly improv-
		//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
		//	c        the average) than the current approximation has.
		//	c           if only an increase in breakpoints is wanted, without the
		//	c        adjustment that newnot provides, a fake newnot routine could be
		//	c        used here which merely returns the breakpoints for  l n e w
		//	c        equal intervals .
		//			 if (nt .ge. ntimes)            go to 1
		//			 lnew = lbegin + float(nt)*addbrk
		//			 call newnot (break, coef, l, k, scrtch, lnew, t )
		//			 call l2knts ( scrtch, lnew, k, t, n )
		//			 nt = nt + 1
		//											go to 10
		//		  end
			}
	public static void	l3main ( double[][] data, int ntau2, double[] weight, double[] t){	
		
//	c  main program for least-squares approximation by splines
//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
//	c
//	c  the program, though ostensibly written for l2-approximation, is typ-
//	c  ical for programs constructing a pp approximation to a function gi-
//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
//	c  easily be replaced by one carrying out interpolation or some other
//	c  form of approximation.
//	c
//	c******  i n p u t  ******
//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
//	c  be approximated and the order and breakpoint sequence of the pp ap-
//	c  proximating function to be used. further,  s e t d a t  is expected
//	c  to  t e r m i n a t e  the run (for lack of further input or because
//	c   i c o u n t  has reached a critical value).
//	c     the number  n t i m e s  is read in in the main program. it speci
//	c  fies the number of passes through the knot improvement algorithm in
//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
//	c  that, on the average, addbrk knots are to be added per pass through
//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
//	c  every third pass (as long as  ntimes .lt. 50).
//	c
//	c******  p r i n t e d  o u t p u t  ******
//	c  is governed by the three print control hollerith strings
//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
//	c  p r f u n  = 'on'  gives printout of approximation and error at
//	c                     every data point.
//	c  the order  k , the number of pieces  l, and the interior breakpoints
//	c  are always printed out as are (in l2err) the mean, mean square, and
//	c  maximum errors in the approximation.
//	c
//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
		int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
		int n[]= new int[1];
//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
		final int lpkmax=100;
		final int ntmax=200;
		final int ltkmax=2000;
//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
//		 *    ,t(ntmax),tau,totalw,weight
		double addbrk,  gtau;
		double bcoef[] = new double[lpkmax];
		double q[][]	= new double[ltkmax][ltkmax];
		double scrtch[] = new double[ntmax];
		
		
//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
//	C    *    ,t(200),tau,totalw,weight
//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
//	c     ioned here only because it might otherwise become undefined be-
//	c     tween calls to those subroutines.
//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
		double brake[]=new double[lpkmax];
		double coef[]=new double[ltkmax];
//	C     common /approx/ break(100),coef(2000),l,k
//	c     common /approx/ also occurs in setdat and l2err.
//		  data on /'ON'/
//	c
//		  icount = 0
		icount=0;
		
//	c        i c o u n t  provides communication with the data-input-and-
//	c     termination routine  s e t d a t . it is initialized to  0  to
//	c     signal to setdat when it is being called for the first time. after
//	c     that, it is up to setdat to use icount for keeping track of the
//	c     passes through setdat .
//	c
//	c     information about the function to be approximated and order and
//	c     breakpoint sequence of the approximating pp functions is gathered
//	c     by a
//		1 call setdat(icount)

		ntau = 10;
		int ntaum1 = ntau-1;
		
		for (i=0;i< ntau;i++){
			data[0][i] = .1*i;
			data[1][i] = Math.pow(data[0][i],2)  + 1. +Math.random()*.1;
			weight[i] = 1.;
			
			System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
		}
//		data[0][1]=.2;
//		data[1][1]=1.1;
//	     

		int totalw = ntau;
		l = 6;
		int lp1 = l+1;
		double step = 1./((double)l);
		k = 3;
		for(i=0;i<lp1;i++){
			brake[i]=i*step;
			System.out.println("\t"+brake[i]);
		}
		




//	c
//	c     breakpoints are translated into knots, and the number  n  of
//	c     b-splines to be used is obtained by a
//		  call l2knts ( break, l, k, t, n )
			n[0]=0;
			l2knts(brake,lp1,k-1,t,n);
			System.out.println("n\t"+n[0]);
//	c
//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
//	c     as well as the print controls  p r b c o ,  p r p c o  and
//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
//	c     not, with an increase of  addbrk  knots for every pass .
//		  print 600
//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
//		  read 500,ntimes,addbrk,prbco,prpco,prfun
//	  500 format(i3,f10.5/3a2)
//	c
//		  lbegin = l
//		  nt = 0
//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
//	c        ed by a
//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
//	  609    format(//' b-spline coefficients'/(4e20.10))
lbegin=1;
nt=0;
//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){

l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
System.out.println(ntau);
for(int ix=0; ix<ntau; ix++){
	System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
}
for(int ix=0; ix < n[0]; ix++){
	System.out.println(t[ix]);
}

for(int ix=0; ix < n[0]; ix++){
	System.out.println(t[ix]);
}
System.out.println(n[0]);
System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]+"\t"+bcoef[8]+"\t"+bcoef[9]);

plotdata ( 0., .9, 100,t,n[0], k,bcoef);

//	c
//	c        convert the b-repr. of the approximation to pp repr.
//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
//			 print 610, k, l, (break(ll),ll=2,l)
//	  610    format(//' approximation by splines of order',i3,' on ',
//		 *         i3,' intervals. breakpoints -'/(4e20.10))
//			 if (prpco .ne. on)             go to 15
//			 print 611
//	  611    format(/' pp-representation for approximation')
//			 do 12 i=1,l
//				ii = (i-1)*k
//	   12       print 613,break(i),(coef(ii+j),j=1,k)
//	  613    format(f9.3,4e20.10/(11x,4e20.10))
//	c
//	c        compute and print out various error norms by a
//	   15    call l2err ( prfun, scrtch, q )
//	c
//	c        if newnot has been applied less than  n t i m e s  times, try
//	c        it again to obtain, from the current approx. a possibly improv-
//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
//	c        the average) than the current approximation has.
//	c           if only an increase in breakpoints is wanted, without the
//	c        adjustment that newnot provides, a fake newnot routine could be
//	c        used here which merely returns the breakpoints for  l n e w
//	c        equal intervals .
//			 if (nt .ge. ntimes)            go to 1
//			 lnew = lbegin + float(nt)*addbrk
//			 call newnot (break, coef, l, k, scrtch, lnew, t )
//			 call l2knts ( scrtch, lnew, k, t, n )
//			 nt = nt + 1
//											go to 10
//		  end
	}	
	public static void	l4main ( double[][] data, int ntau2, double[] weight, double[] t){	
		
//	c  main program for least-squares approximation by splines
//	c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
//	calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
//	c     ,l2err(ppvalu(interv)),ppvalu*,newnot
//	c
//	c  the program, though ostensibly written for l2-approximation, is typ-
//	c  ical for programs constructing a pp approximation to a function gi-
//	c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
//	c  easily be replaced by one carrying out interpolation or some other
//	c  form of approximation.
//	c
//	c******  i n p u t  ******
//	c  is expected in  s e t d a t  (quo vide), specifying both the data to
//	c  be approximated and the order and breakpoint sequence of the pp ap-
//	c  proximating function to be used. further,  s e t d a t  is expected
//	c  to  t e r m i n a t e  the run (for lack of further input or because
//	c   i c o u n t  has reached a critical value).
//	c     the number  n t i m e s  is read in in the main program. it speci
//	c  fies the number of passes through the knot improvement algorithm in
//	c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
//	c  that, on the average, addbrk knots are to be added per pass through
//	c  newnot. for example,  addbrk = .34  would cause a knot to be added
//	c  every third pass (as long as  ntimes .lt. 50).
//	c
//	c******  p r i n t e d  o u t p u t  ******
//	c  is governed by the three print control hollerith strings
//	c  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
//	c  p r p c o  = 'on'  gives printout of pp repr. of approximation.
//	c  p r f u n  = 'on'  gives printout of approximation and error at
//	c                     every data point.
//	c  the order  k , the number of pieces  l, and the interior breakpoints
//	c  are always printed out as are (in l2err) the mean, mean square, and
//	c  maximum errors in the approximation.
//	c
//		  integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
//		 *        ,ntimes,ntmax,on,prbco,prfun,prpco
		int i, icount, ii, j, k, l, lbegin, lnw, ll,  nt, ntau, ntimes, on, prbco, prfun, prpco;
		int n[]= new int[1];
//		  parameter (lpkmax=100,ntmax=200,ltkmax=2000)
		final int lpkmax=100;
		final int ntmax=200;
		final int ltkmax=2000;
//		  real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax)
//		 *    ,t(ntmax),tau,totalw,weight
		double addbrk,  gtau;
		double bcoef[] = new double[lpkmax];
		double q[][]	= new double[ltkmax][ltkmax];
		double scrtch[] = new double[ntmax];
		
		
//		  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
//	C     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
//	C    *    ,t(200),tau,totalw,weight
//	C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
//	c     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
//	c     ioned here only because it might otherwise become undefined be-
//	c     tween calls to those subroutines.
//		  common /approx/ break(lpkmax),coef(ltkmax),l,k
		double brake[]=new double[lpkmax];
		double coef[]=new double[ltkmax];
//	C     common /approx/ break(100),coef(2000),l,k
//	c     common /approx/ also occurs in setdat and l2err.
//		  data on /'ON'/
//	c
//		  icount = 0
		icount=0;
		
//	c        i c o u n t  provides communication with the data-input-and-
//	c     termination routine  s e t d a t . it is initialized to  0  to
//	c     signal to setdat when it is being called for the first time. after
//	c     that, it is up to setdat to use icount for keeping track of the
//	c     passes through setdat .
//	c
//	c     information about the function to be approximated and order and
//	c     breakpoint sequence of the approximating pp functions is gathered
//	c     by a
//		1 call setdat(icount)

		ntau = 10;
		int ntaum1 = ntau-1;
		
		for (i=0;i< ntau;i++){
			data[0][i] = .1*i;
			data[1][i] = Math.pow(data[0][i],2)  + 1. +Math.random()*.1;
			weight[i] = 1.;
			
			System.out.println("**"+data[0][i]+"\t"+data[1][i]+"\t"+weight[i]);
		}
//		data[0][1]=.2;
//		data[1][1]=1.1;
//	     

		int totalw = ntau;
		l = 6;
		int lp1 = l+1;
		double step = 1./((double)l);
		k = 3;
		for(i=0;i<lp1;i++){
			brake[i]=i*step;
			System.out.println("\t"+brake[i]);
		}
		




//	c
//	c     breakpoints are translated into knots, and the number  n  of
//	c     b-splines to be used is obtained by a
//		  call l2knts ( break, l, k, t, n )
			n[0]=0;
			l2knts(brake,lp1,k-1,t,n);
			System.out.println("n\t"+n[0]);
//	c
//	c     the integer  n t i m e s  and the real  a d d b r k  are requested
//	c     as well as the print controls  p r b c o ,  p r p c o  and
//	c     p r f u n .  ntimes  passes  are made through the subroutine new-
//	c     not, with an increase of  addbrk  knots for every pass .
//		  print 600
//	  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
//		  read 500,ntimes,addbrk,prbco,prpco,prfun
//	  500 format(i3,f10.5/3a2)
//	c
//		  lbegin = l
//		  nt = 0
//	c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
//	c        ed by a
//	   10    call l2appr ( t, n, k, q, scrtch, bcoef )
//			 if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
//	  609    format(//' b-spline coefficients'/(4e20.10))
lbegin=1;
nt=0;
//l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){

l2appr(data,ntau,weight,t,n[0],k,q,scrtch,bcoef);
System.out.println(ntau);
for(int ix=0; ix<ntau; ix++){
	System.out.println(data[0][ix]+"\t"+data[1][ix]+"\t"+weight[ix]);
}
for(int ix=0; ix < n[0]; ix++){
	System.out.println(t[ix]);
}

for(int ix=0; ix < n[0]; ix++){
	System.out.println(t[ix]);
}
System.out.println(n[0]);
System.out.println(bcoef[0]+"\t"+bcoef[1]+"\t"+bcoef[2]+"\t"+bcoef[3]+"\t"+bcoef[4]+"\t"+bcoef[5]+"\t"+bcoef[6]+"\t"+bcoef[7]+"\t"+bcoef[8]+"\t"+bcoef[9]);

plotdata ( 0., .9, 100,t,n[0], k,bcoef);

//	c
//	c        convert the b-repr. of the approximation to pp repr.
//			 call bsplpp ( t, bcoef, n, k, q, break, coef, l )
//			 print 610, k, l, (break(ll),ll=2,l)
//	  610    format(//' approximation by splines of order',i3,' on ',
//		 *         i3,' intervals. breakpoints -'/(4e20.10))
//			 if (prpco .ne. on)             go to 15
//			 print 611
//	  611    format(/' pp-representation for approximation')
//			 do 12 i=1,l
//				ii = (i-1)*k
//	   12       print 613,break(i),(coef(ii+j),j=1,k)
//	  613    format(f9.3,4e20.10/(11x,4e20.10))
//	c
//	c        compute and print out various error norms by a
//	   15    call l2err ( prfun, scrtch, q )
//	c
//	c        if newnot has been applied less than  n t i m e s  times, try
//	c        it again to obtain, from the current approx. a possibly improv-
//	c        ed sequence of breakpoints with  addbrk  more breakpoints (on
//	c        the average) than the current approximation has.
//	c           if only an increase in breakpoints is wanted, without the
//	c        adjustment that newnot provides, a fake newnot routine could be
//	c        used here which merely returns the breakpoints for  l n e w
//	c        equal intervals .
//			 if (nt .ge. ntimes)            go to 1
//			 lnew = lbegin + float(nt)*addbrk
//			 call newnot (break, coef, l, k, scrtch, lnew, t )
//			 call l2knts ( scrtch, lnew, k, t, n )
//			 nt = nt + 1
//											go to 10
//		  end
	}	
	public static int	l2appr ( double[][] data, int ntau, double[] weight, double[] t, int n, int k, double[][] q, double[] diag, double[] bcoef ){
	
//c  from  * a practical guide to splines *  by c. de boor    
//c  to be called in main program  l 2 m a i n .
//calls subprograms  bsplvb, bchfac/slv
//c
//constructs the (weighted discrete) l2-approximation by splines of order
//c  k  with knot sequence  t(1), ..., t(n+k)  to given data points
//c  ( tau(i), gtau(i) ), i=1,...,ntau. the b-spline coefficients
//c  b c o e f   of the approximating spline are determined from the
//c  normal equations using cholesky's method.
//c
//c******  i n p u t  ******
//c  t(1), ..., t(n+k)  the knot sequence
//c  n.....the dimension of the space of splines of order k with knots t.
//c  k.....the order
//c
//c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
//c        sed by the arbitrary dimension statement for  biatx  below, but
//c        is  n o w h e r e   c h e c k e d   for.
//c
//c******  w o r k  a r r a y s  ******
//c  q....a work array of size (at least) k*n. its first  k  rows are used
//c       for the  k  lower diagonals of the gramian matrix  c .
//c  diag.....a work array of length  n  used in bchfac .
//c
//c******  i n p u t  via  c o m m o n  /data/  ******
//c  ntau.....number of data points
//c  (tau(i),gtau(i)), i=1,...,ntau     are the  ntau  data points to be
//c        fitted .
//c  weight(i), i=1,...,ntau    are the corresponding weights .
//c
//c******  o u t p u t  ******
//c  bcoef(1), ..., bcoef(n)  the b-spline coeffs. of the l2-appr.
//c
//c******  m e t h o d  ******
//c  the b-spline coefficients of the l2-appr. are determined as the sol-
//c  ution of the normal equations
//c     sum ( (b(i),b(j))*bcoef(j) : j=1,...,n)  = (b(i),g),
//c                                               i = 1, ..., n .
//c  here,  b(i)  denotes the i-th b-spline,  g  denotes the function to
//c  be approximated, and the  i n n e r   p r o d u c t  of two funct-
//c  ions  f  and  g  is given by
//c      (f,g)  :=  sum ( f(tau(i))*g(tau(i))*weight(i) : i=1,...,ntau) .
//c  the arrays  t a u  and  w e i g h t  are given in common block
//c   d a t a , as is the array  g t a u  containing the sequence
//c  g(tau(i)), i=1,...,ntau.
//c  the relevant function values of the b-splines  b(i), i=1,...,n, are
//c  supplied by the subprogram  b s p l v b .
//c     the coeff.matrix  c , with
//c           c(i,j)  :=  (b(i), b(j)), i,j=1,...,n,
//c  of the normal equations is symmetric and (2*k-1)-banded, therefore
//c  can be specified by giving its k bands at or below the diagonal. for
//c  i=1,...,n,  we store
//c   (b(i),b(j))  =  c(i,j)  in  q(i-j+1,j), j=i,...,min0(i+k-1,n)
//c  and the right side
//c   (b(i), g )  in  bcoef(i) .
//c  since b-spline values are most efficiently generated by finding sim-
//c  ultaneously the value of  e v e r y  nonzero b-spline at one point,
//c  the entries of  c  (i.e., of  q ), are generated by computing, for
//c  each ll, all the terms involving  tau(ll)  simultaneously and adding
//c  them to all relevant entries.
//	  integer k,n,   i,j,jj,kmax,left,leftmk,ll,mm,ntau,ntmax
//	  parameter (kmax=20,ntmax=200)
//	  real bcoef(n),diag(n),q(k,n),t(1),  biatx(kmax),dw,gtau,tau,weight

	int i,j,jj,left,leftmk, ll, mm;
    final int kmax=20;
    final int ntmax=200;
	int jmax=20;
	double dw;
	double biatx[] = new double[k+1];
	double deltal[] =new double[jmax];
	double deltar[] = new double[jmax];
	int j_it[]= new int[1];
	j_it[0]=1;

//C     real bcoef(n),diag(n),q(k,n),t(1),  biatx(20),dw,gtau,tau,weight
//	  dimension t(n+k)
//c former fortran standard made it impossible to specify the exact dimen-
//c  sion of  t  without the introduction of additional otherwise super-
//c  fluous arguments.
//	  common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax)
//C     common / data / ntau, tau(200),gtau(200),weight(200)
//	  do 7 j=1,n
//		 bcoef(j) = 0.
//		 do 7 i=1,k
//	7       q(i,j) = 0.

		for(j=0;j<n;j++){
			bcoef[j]=0;
			for(i=0;i<k;i++){
				q[i][j]	= 0.0;
			}
		}
//	  left = k
//	  leftmk = 0
		left=k-1;
		leftmk=0;


//	  do 20 ll=1,ntau
//c                   locate  l e f t  s.t. tau(ll) in (t(left),t(left+1))
//   10       if (left .eq. n)            go to 15
//			if (tau(ll) .lt. t(left+1)) go to 15
//			left = left+1
//			leftmk = leftmk + 1
//										go to 10

		for(ll=0;ll<ntau;ll++){
			//left=0;
			while(data[0][ll]>=t[left+1] && left<n){
				left++;
				leftmk++;
			}
			// System.out.println("* "+t[left]+"\t"+data[0][ll]+"\t"+t[left+1]+"\t"+ left );

			
			
			
//			15    call bsplvb ( t, k, 1, tau(ll), left, biatx )
			bsplvb(t,k,1,data[0][ll],left, biatx,deltar, deltal, j_it);
		//	System.out.print("***\t"+data[0][ll]);
			for(int ix=0;ix<left;ix++){
			//	System.out.print("\t0");
			}
		//	 System.out.println("***\t"+data[0][ll]+"\t"+left+"\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]+"\t"+biatx[3]);
		//	System.out.println("\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]);

//		 c        biatx(mm) contains the value of b(left-k+mm) at tau(ll).
//		 c        hence, with  dw := biatx(mm)*weight(ll), the number dw*gtau(ll)
//		 c        is a summand in the inner product
//		 c           (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
//		 c        and the number biatx(jj)*dw is a summand in the inner product
//		 c           (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
//		 c        since  (left-k+jj) - (left-k+mm) + 1  =  jj - mm + 1 .
//				do 20 mm=1,k
//				   dw = biatx(mm)*weight(ll)
//				   j = leftmk + mm
//				   bcoef(j) = dw*gtau(ll) + bcoef(j)
//				   i = 1
//				   do 20 jj=mm,k
//					  q(i,j) = biatx(jj)*dw + q(i,j)
//			20          i = i + 1

			for(mm=0;mm<k;mm++){
				dw=biatx[mm]*weight[ll];
				j=leftmk + mm;
				bcoef[j]=dw*data[1][ll] + bcoef[j];
				i=0;
				for(jj=mm;jj<k;jj++){
					q[i][j]=biatx[jj]*dw + q[i][j];
					i++;
				}
			}

		
		}


//c
//c             construct cholesky factorization for  c  in  q , then use
//c             it to solve the normal equations
//c                    c*x  =  bcoef
//c             for  x , and store  x  in  bcoef .
//	  call bchfac ( q, k, n, diag )
//	  call bchslv ( q, k, n, bcoef )
//										return
//	  end
	bchfac(q,k,n,diag);
	bchslv(q,k,n,bcoef);
	return 0;
	}


	public static void bchfac (double[][] w, int nbands, int nRow, double[]  diag ){
//	c  from  * a practical guide to splines *  by c. de boor    
//	constructs cholesky factorization
//	c                     c  =  l * d * l-transpose
//	c  with l unit lower triangular and d diagonal, for given matrix c of
//	c  order  n r o w , in case  c  is (symmetric) positive semidefinite
//	c  and  b a n d e d , having  n b a n d s  diagonals at and below the
//	c  main diagonal.
//	c
//	c******  i n p u t  ******
//	c  nRow.....is the order of the matrix  c .
//	c  nbands.....indicates its bandwidth, i.e.,
//	c          c(i,j) = 0 for i-j .ge. nbands .
//	c  w.....workarray of size (nbands,nRow)  containing the  nbands  diago-
//	c        nals in its rows, with the main diagonal in row  1 . precisely,
//	c        w(i,j)  contains  c(i+j-1,j), i=1,...,nbands, j=1,...,nRow.
//	c          for example, the interesting entries of a seven diagonal sym-
//	c        metric matrix  c  of order  9  would be stored in  w  as
//	c
//	c
//	c
//	c
//	c
//	c
//	c        all other entries of  w  not identified in this way with an en-
//	c        try of  c  are never referenced .
//	c  diag.....is a work array of length  nRow .
//	c
//	c******  o u t p u t  ******
//	c  w.....contains the cholesky factorization  c = l*d*l-transp, with
//	c        w(1,i) containing  1/d(i,i)
//	c        and  w(i,j)  containing  l(i-1+j,j), i=2,...,nbands.
//	c
//	c******  m e t h o d  ******
//	c   gauss elimination, adapted to the symmetry and bandedness of  c , is
//	c   used .
//	c     near zero pivots are handled in a special way. the diagonal ele-
//	c  ment c(n,n) = w(1,n) is saved initially in  diag(n), all n. at the n-
//	c  th elimination step, the current pivot element, viz.  w(1,n), is com-
//	c  pared with its original value, diag(n). if, as the result of prior
//	c  elimination steps, this element has been reduced by about a word
//	c  length, (i.e., if w(1,n)+diag(n) .le. diag(n)), then the pivot is de-
//	c  clared to be zero, and the entire n-th row is declared to be linearly
//	c  dependent on the preceding rows. this has the effect of producing
//	c   x(n) = 0  when solving  c*x = b  for  x, regardless of  b. justific-
//	c  ation for this is as follows. in contemplated applications of this
//	c  program, the given equations are the normal equations for some least-
//	c  squares approximation problem, diag(n) = c(n,n) gives the norm-square
//	c  of the n-th basis function, and, at this point,  w(1,n)  contains the
//	c  norm-square of the error in the least-squares approximation to the n-
//	c  th basis function by linear combinations of the first n-1 . having
//	c  w(1,n)+diag(n) .le. diag(n) signifies that the n-th function is lin-
//	c  early dependent to machine accuracy on the first n-1 functions, there
//	c  fore can safely be left out from the basis of approximating functions
//	c     the solution of a linear system
//	c                       c*x = b
//	c   is effected by the succession of the following  t w o  calls:
//	c     call bchfac ( w, nbands, nRow, diag )       , to get factorization
//	c     call bchslv ( w, nbands, nRow, b )          , to solve for x.
//	c


//		  integer nbands,nRow,   i,imax,j,jmax,n
	int i, iMax, j, jMax, n;
//		  real w(nbands,nRow),diag(nRow),   ratio
	double ratio;
	
//	if (nRow .gt. 1)                  go to 9
//	if (w(1,1) .gt. 0.) w(1,1) = 1./w(1,1)
//									  return	
	if(nRow <= 1){
		if(w[0][0]>0) {
			w[0][0] = 1/w[0][0];
			return;
			} 
	}

//	c                                        store diagonal of  c  in  diag.
//		9 do 10 n=1,nRow
//	   10    diag(n) = w(1,n)

	for(n=0;n<nRow;n++){
		diag[n]=w[0][n];
	}

//	c                                                        factorization .
//		  do 20 n=1,nRow
//			 if (w(1,n)+diag(n) .gt. diag(n)) go to 15
//			 do 14 j=1,nbands
//	   14       w(j,n) = 0.
//											go to 20
//	   15    w(1,n) = 1./w(1,n)
//			 imax = min0(nbands-1,nRow - n)
//			 if (imax .lt. 1)               go to 20
//			 jmax = imax
//			 do 18 i=1,imax
//				ratio = w(i+1,n)*w(1,n)
//				do 17 j=1,jmax
//	   17          w(j,n+i) = w(j,n+i) - w(j+i,n)*ratio
//				jmax = jmax - 1
//	   18       w(i+1,n) = ratio
//	   20    continue
//											return
//		  end
		ratio=1;
		for(n=0;n<nRow;n++){
			if(w[0][n]+diag[n] <= diag[n]){
				for(j=0;j<nbands;j++){
					w[j][n]=0;
				}
			}
			else{
				w[0][n]=1/w[0][n];   
				iMax=Math.min(nbands-1,nRow-n+1);//check 1's 
				if(iMax >= 1){
					jMax=iMax;
					for(i=0;i<iMax;i++){
						ratio=w[i+1][n]*w[0][n];
						for(j=0;j<jMax;j++){
							w[j][n+i+1]=w[j][n+i+1]-w[j+i+1][n]*ratio;
						}
						jMax=jMax-1;
						w[i+1][n]=ratio;
					}
				}
			}
		}//20
		return;
	}
	
	public static void bchslv  (double[][] w, int nbands, int nRow, double[]  b ){

//	subroutine bchslv ( w, nbands, nRow, b )
//	c  from  * a practical guide to splines *  by c. de boor    
//	c  solves the linear system     c*x = b   of order  n r o w  for  x
//	c  provided  w  contains the cholesky factorization for the banded (sym-
//	c  metric) positive definite matrix  c  as constructed in the subroutine
//	c    b c h f a c  (quo vide).
//	c
//	c******  i n p u t  ******
//	c  nRow.....is the order of the matrix  c .
//	c  nbands.....indicates the bandwidth of  c .
//	c  w.....contains the cholesky factorization for  c , as output from
//	c        subroutine bchfac  (quo vide).
//	c  b.....the vector of length  n r o w  containing the right side.
//	c
//	c******  o u t p u t  ******
//	c  b.....the vector of length  n r o w  containing the solution.
//	c
//	c******  m e t h o d  ******
//	c  with the factorization  c = l*d*l-transpose  available, where  l  is
//	c  unit lower triangular and  d  is diagonal, the triangular system
//	c  l*y = b  is solved for  y (forward substitution), y is stored in  b,
//	c  the vector  d**(-1)*y is computed and stored in  b, then the triang-
//	c  ular system  l-transpose*x = d**(-1)*y is solved for  x (backsubstit-
//	c  ution).
//		  integer nbands,nRow,   j,jmax,n,nbndm1
		int j,n,nbndm1;
		int jMax=0;
//		  real w(nbands,nRow),b(nRow)
//		if (nRow .gt. 1)                  go to 21
//		b(1) = b(1)*w(1,1)
//										  return
		if (nRow <= 1){
			b[0]=b[0]*w[0][0];
			return;
		}

//	c
//	c                 forward substitution. solve l*y = b for y, store in b.
//	   21 nbndm1 = nbands - 1
//		  do 30 n=1,nRow
//			 jmax = min0(nbndm1,nRow-n)
//			 if (jmax .lt. 1)               go to 30
//			 do 25 j=1,jmax
//	   25       b(j+n) = b(j+n) - w(j+1,n)*b(n)
//	   30    continue
//	c
		nbndm1=nbands -1;
		for(n=0;n<nRow;n++){
			jMax=Math.min(nbndm1,nRow-n);
			if(jMax>=1){
				for(j=0;j<jMax;j++){
					b[j+n+1]=b[j+n+1]-w[j+1][n]	* b[n];
				}
			}
		}//30


//	c     backsubstitution. solve l-transp.x = d**(-1)*y  for x, store in b.
//		  n = nrow    
//	   39    b(n) = b(n)*w(1,n)   
//			 jmax = min0(nbndm1,nrow-n)
//		if (jmax .lt. 1)               go to 40
//		do 35 j=1,jmax
//35       b(n) = b(n) - w(j+1,n)*b(j+n)
//40    n = n-1  
//		if (n .gt. 0)                  go to 39
//									   return
//	 end

			n=nRow-1;
			while(n>-1){
				b[n]=b[n]*w[0][n];
				jMax=Math.min(nbndm1,nRow-n);
				if(jMax >= 1){
					for(j=0;j<jMax;j++){
						b[n]=b[n]-w[j+1][n]*b[j+n+1];
					}
					//n--;
				}
				n--;
			}			
		return;
	}
	
	
	public static void bsplvb(double[] t, int jhigh, int index, double x, int left, double[] biatx, double[]deltar, double[]deltal, int[]j_it ) {
		// This is from de Boor A Practical Guide to Splines
		// This is a java conversion from a de Boor's FORTRAN routine
		// The highlighted text is essentially de Boor's with very minor alteration
		
		//*****************************************************************
		// 
		// Calculates the value of all possible nonzero b-splines at x of order 
		//		jout = max(jhigh , (j+1)*(index-1)
		// with knot sequence t. 
	
		// Input
		// t:	a knot squence of length left + jout, assumed to be nondecreasing with
		// assumption that t(left) < t(left +1)
		// division by zero will result if t(left) = t(left + 1)
		//
		// jhigh, index: 	integers which determine the order jout=max(jhigh, (j+1)*(index-1))
		// of the b-splines whose values at x are to be returned.  index is used to avoid
		// recalculations when several columns of the triangula array of b-splines values are 
		// needed (e.g., in bsplpp or in bsplvd).  Precisely 
		//
		// if index = 1 then the calculation starts from scratch and the entire triangular
		// array of b-spline values of orders 1,2,..., jhigh is generated order by order, 
		// i.e., column by column.
		// 
		// if index = 2 only the b-spline values of order j+1, j+2, ..., jout are generated,
		// the assumption being that biatx, j, delta1, deltar are, on entry, as they were on 
		// exit at the previous call.  In particular, if jhigh = 0, the jout = j=+1, i.e., just
		// the next column of the b-spline is generated
		
		// Warning...the restriction jout <= jmax (=20) is imposed arbitrarily by the 
		// dimension statment for deltal and deltar below, but is nowhere check for.
		//
		// x:	the point at which the b-splines are to be evaluated.
		// left:	an integer chosen (usually) so that t(left) <= x <= t(left+1)
		//
		// Output
		// biatx:	array of length jout, with biatx[i] containing the value at x of the polynomial
		// of order jout which agrees with the b-spline b(left-jout+i, jout, t) on the interval
		// (t(left), t(left+1)).
		//
		// Method
		// The recurrence relation
		//
		//                   x-t(i)                       t(i+j+1) - x
		//	b(i,j+1) (x) = -------------- b(i,j) (x) + ------------------ b(i+1,j)(x)
		//	                t(i+j)-t(i)                 t(i+j+1)-t(i+1)
		//
		// is used (repeatedly) to generate the (j+1)-vector b(left-j,j+1)(x),
		// ..., b(left,j+1) (x) from the j-vector b(left-j+1, j)(x), ...., 
		// b(left, j)(x) storing the new values in biatx over the old,  the 
		// facts that 
		// 		b(i,1) =1 if t(i) <= x < t(i+1)
		// and
		//   	b(i,j)(x) = 0 unless t(i) <= x <= t(i+j)
		//
		// are used.  The particular organization of the calculations follow alogrithm (8) 
		// in chapter x of the text.
		//***************************************************************************
		
		int jmax = 20;
		int i, j, jpl;
	
		double saved, term;
		
		j=j_it[0];
		if(index==1){
			j=1;
			biatx[0]=1.0;
			for(int indx=0; indx<20; indx++){
				deltar[indx]=0.;
				deltal[indx]=0.;
			}
			if(j>=jhigh){
				j_it[0]=j;
				return;
			}
		}
		
		while(j < jhigh){
			jpl = j+1;
			deltar[j-1]=t[left+j]-x;
			deltal[j-1]= x-t[left+1-j];
			saved = 0.0;
			for(i=1; i<=j; i++){
				term = biatx[i-1]/(deltar[i-1] + deltal[jpl-i-1]);
				biatx[i-1] = saved + deltar[i-1]*term;
				saved= deltal[jpl-i-1]*term;	
//				System.out.println("@i + biatx[i] " + i + "\t"+ biatx[i]);
//				System.out.println("term + biatx[i]+ saved" + term + "\t"+ biatx[i]+ "\t"+ saved);

			}
			biatx[jpl-1]=saved;
			j=jpl;
//			System.out.println("@jpl + biatx[i] " + jpl + "\t"+ biatx[jpl-1]);

		}

		j_it[0]=j;
		return;
	}
	public static void bsplvd(double[] t, int k,  double x, int left, double[][] a, double[][] dbiatx, int nderiv) {
		// This is from de Boor, A Practical Guide to Splines
		// This is a java conversion from a de Boor's FORTRAN routine
		// The highlighted text is essentially de Boor's with very minor alteration
		
		//*****************************************************************
		//		subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )
		//  from  * a practical guide to splines *  by c. de Boor (7 may 92)    
		//  calls bsplvb
		//   calculates value and deriv.s of all b-splines which do not vanish at x
		//
		//******  i n p u t  ******
		//  t     the knot array, of length left+k (at least)
		//  k     the order of the b-splines to be evaluated
		//  x     the point at which these values are sought
		//  left  an integer indicating the left endpoint of the interval of
		//        interest. the  k  b-splines whose support contains the interval
		//               (t(left), t(left+1))
		//        are to be considered.
		//  a s s u m p t i o n  - - -  it is assumed that
		//               t(left) .lt. t(left+1)
		//        division by zero will result otherwise (in  b s p l v b ).
		//        also, the output is as advertised only if
		//               t(left) .le. x .le. t(left+1) .
		//  nderiv   an integer indicating that values of b-splines and their
		//        derivatives up to but not including the  nderiv-th  are asked
		//        for. ( nderiv  is replaced internally by the integer  m h i g h
		//        in  (1,k)closest to it.)
		//
		//******  w o r k   a r e a  ******
		//  a     an array of order (k,k), to contain b-coeff.s of the derivat-
		//        ives of a certain order of the  k  b-splines of interest.
		//
		//******  o u t p u t  ******
		//  dbiatx   an array of order (k,nderiv). its entry  (i,m)contains
		//        value of  (m-1)st  derivative of  (left-k+i)-th  b-spline of
		//        order  k  for knot sequence  t , i=1,...,k, m=1,...,nderiv.
		//
		//******  m e t h o d  ******
		//  values at  x  of all the relevant b-splines of order k,k-1,...,
		//  k+1-nderiv  are generated via  bsplvb  and stored temporarily in
		//  dbiatx .  then, the b-coeffs of the required derivatives of the b-
		//  splines of interest are generated by differencing, each from the pre-
		//	ceding one of lower order, and combined with the values of b-splines
		//  of corresponding order in  dbiatx  to produce the desired values .
		//  c
		//		integer k,left,nderiv,   i,ideriv,il,j,jlow,jp1mid,kp1,kp1mm
		//	   *                        ,ldummy,m,mhigh
		int i, ideriv, il,j ,jlow, jp1mid, kp1, kp1mm, ldummy, m, mhigh;
		//		real a(k,k),dbiatx(k,nderiv),t(1),x,   factor,fkp1mm,sum
		double factor, fkp1mm, sum;
		mhigh=Math.max(Math.min(nderiv,k),1);
		//		mhigh = max0(min0(nderiv,k),1)
		//      mhigh is usually equal to nderiv.
	 	kp1 = k+1;
//public static void bsplvb(t, jhigh, int index, double x, int left, double[] biatx, double[]deltar, double[]deltal, int[]j_it ) {
		int jhigh=kp1-mhigh;
		int jmax=20;
		double biatx[] = new double[k+1];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
		//double a[][] =new double[k][k];
		int j_it[]= new int[1];
		for(int ix =0 ; ix< k; ix++)
			for(int jx=0; jx<nderiv; jx++)
				dbiatx[ix][jx]=0.0;
				
	 	bsplvb(t,kp1-mhigh,1,x,left,biatx ,deltar, deltal,j_it);

	 	if(mhigh == 1) {
	 		return;
	 	}
		//		if (mhigh .eq. 1)                 go to 99
		//  c     the first column of  dbiatx  always contains the b-spline values
		//  c     for the current order. these are stored in column k+1-current
		//  c     order  before  bsplvb  is called to put values for the next
		//  c     higher order on top of it.
		 		ideriv = mhigh;
		//		do 15 m=2,mhigh
		//		   jp1mid = 1
		//		   do 11 j=ideriv,k
		//			  dbiatx(j,ideriv) = dbiatx(jp1mid,1)
		//	 11       jp1mid = jp1mid + 1
		//		   ideriv = ideriv - 1
		//		   call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
		//	 15    continue
			for(m=2;m<=mhigh+1;m++){
				jp1mid=0;//check
				for(j=ideriv; j<=k; j++){
					dbiatx[j-1][ideriv-1] = biatx[jp1mid];
					jp1mid++;
				}
				ideriv = ideriv - 1;
				bsplvb(t,kp1-ideriv,1,x,left,biatx,deltar, deltal,j_it);
			}
				
		//  c
		//  c     at this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j) for
		//  c     i=j,...,k and j=1,...,mhigh ('=' nderiv). in particular, the
		//  c     first column of  dbiatx  is already in final form. to obtain cor-
		//  c     responding derivatives of b-splines in subsequent columns, gene-
		//  c     rate their b-repr. by differencing, then evaluate at  x.
		//  c
		//		jlow = 1
		//		do 20 i=1,k
		//		   do 19 j=jlow,k
		//	 19       a(j,i) = 0.
		//		   jlow = i
		//	 20    a(i,i) = 1.
		jlow=0;
		for(i=0;i<k;i++){
			for(j=jlow;j<k;j++){
				a[j][i]=0.0;
			}
			jlow=i;
			a[i][i]=1.0;
		}
		
		
		
		//  c     at this point, a(.,j) contains the b-coeffs for the j-th of the
		//  c     k  b-splines of interest here.
		//  c
		//		do 40 m=2,mhigh
		//		   kp1mm = kp1 - m
		//		   fkp1mm = float(kp1mm)
		//		   il = left
		//		   i = k
			for(m=2;m<=mhigh;m++){
				kp1mm=kp1-m;
				fkp1mm=(double)kp1mm;
				il=left;
				i=k;
				

		//  c
		//  c        for j=1,...,k, construct b-coeffs of  (m-1)st  derivative of
		//  c        b-splines from those for preceding derivative by differencing
		//  c        and store again in  a(.,j) . the fact that  a(i,j) = 0  for
		//  c        i .lt. j  is used.
				for(ldummy=1; ldummy<=kp1mm; ldummy++){
				
		//		   do 25 ldummy=1,kp1mm
		//			  factor = fkp1mm/(t(il+kp1mm) - t(il))
					factor = fkp1mm/(t[il+kp1mm]-t[il]);
		//  c           the assumption that t(left).lt.t(left+1) makes denominator
		//  c           in  factor  nonzero.
		//			  do 24 j=1,i
		//	 24          a(i,j) = (a(i,j) - a(i-1,j))*factor
		//			  il = il - 1
		//	 25       i = i - 1
					for(j=0;j<i; j++){
						a[i-1][j] = (a[i-1][j] - a[i-2][j])*factor;
					}
					il=il-1;
					i=i-1;
				}
				
		//  c
		//  c        for i=1,...,k, combine b-coeffs a(.,i) with b-spline values
		//  c        stored in dbiatx(.,m) to get value of  (m-1)st  derivative of
		//  c        i-th b-spline (of interest here) at  x , and store in
		//  c        dbiatx(i,m). storage of this value over the value of a b-spline
		//  c        of order m there is safe since the remaining b-spline derivat-
		//  c        ives of the same order do not use this value due to the fact
		//  c        that  a(j,i) = 0  for j .lt. i .
		//		   do 40 i=1,k
		//			  sum = 0.
		//			  jlow = max0(i,m)
		//			  do 35 j=jlow,k
		//	 35          sum = a(j,i)*dbiatx(j,m) + sum
		//	 40       dbiatx(i,m) = sum
				for(i=0; i<k; i++){
					sum=0.000;
					jlow=Math.max(i,m-1); //check the -1...do we need it for the java arrays
					for(j=jlow;j<=k;j++){
						sum=a[j-1][i]*dbiatx[j-1][m-1] + sum;
					}
					dbiatx[i][m-1] = sum;
				}
			}
			return;
		//	 99                                   return
//		end

	}
	
	
	
	public static void main(String[] args) {
		int jhigh=5;
		int left;
		int jmax=20;
		int maxderiv=3;
		double biatx[] = new double[jhigh];
		double dbiatx[][] = new double[jhigh][maxderiv];
		int ia;
//		double t[]={0,0,0,0,1,3,4,6,6,6,6};
//double t[]={0,0,0, 1, 1, 3, 4, 6,6,6};
//    	double t[]={0,0,0,0,.25,.5,.75,1,1,1,1};
		double t[]= new double[100];
		double deltal[] =new double[jmax];
		double deltar[] = new double[jmax];
		int j_it[]= new int[1];
		double a[][] =new double[jhigh][jhigh];
		int nderiv=0;		

	    double values[]=new double[5];
		double d1values[]=new double[5];
		double d2values[]=new double[5];
	    double x;
	    j_it[0]=1;
	    int n=7;
	    int k=3;
	    int npoint=31;
//	    double xl=t[k-1];
//	    double dx=(t[n]-t[k-1])/(double)(npoint-1);
	    
	    
//		for(x=0; x<180; x=x+1){
//		   for(left=3;left<=6; left++){
//			//System.out.println(x+"\t"+left+"\t"+ t[left] + "\t"+ t[left+1]);
//				if(t[left]<=x && x<t[left+1])
//				{
//					//bsplvb( t, 5, 1, x, left, biatx, deltal, deltar, j_it);
//					bsplvd( t,  4,   x,   left,   a,   dbiatx,  3);
//					System.out.println("&&\t"+x+"\t"+dbiatx[0][1]+"\t"+dbiatx[1][1]+"\t"+dbiatx[2][1]);
//		
//				}
//				//values[left-3] = biatx[7-left];
//
//				//d1values[left-4] = dbiatx[7-left][0];
//				//d2values[left-4] = dbiatx[7-left][1];
//		   }
		  // for(left=0;left<4; left++){
			//System.out.println(x+"\t"+left+"\t"+ t[left] + "\t"+ t[left+1]);
			//System.out.println("&&\t"+x+"\t"+dbiatx[0][1]+"\t"+dbiatx[1][1]+"\t"+dbiatx[2][1]);

		//		values[left-3] = dbiatx[6-left];
		//left=0;
		//System.out.println("&&\t"+x+"\t"+dbiatx[left][0]+"\t"+dbiatx[left][1]+"\t"+dbiatx[left][2]);
		//	System.out.println(x+"\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[left][2]+"\t"+biatx[left][3]);
			
				//d1values[left-4] = dbiatx[7-left][0];
				//d2values[left-4] = dbiatx[7-left][1];
		   //}
//x=0.6;
//		bsplvb( t, 3, 1, x, 6, biatx, deltal, deltar, j_it);
//		for(int   mmm = 1 ; mmm < 5; mmm++){
//		values[mmm ] = biatx[mmm];
//		}

	  //System.out.println(x+"\t"+biatx[0]+"\t"+biatx[1]+"\t"+biatx[2]+"\t"+biatx[3]);
	  
		 // System.out.println(x+"\t"+values[0]+"\t"+values[1]+"\t"+values[2]+"\t"+values[3]);
		//  System.out.println();

	//	double t[]= new double[100];
		
		
	//	}
		double data[][] = new double[100][100];
		int ntau2=100;
		double weight[] = new double[100];
		l2a_main ( data, ntau2,   weight,  t);
//		chapter x. example 2. plotting the pol,s which make up a b-spline
//		c  from  * a practical guide to splines *  by c. de boor    
//		calls  bsplvb
//		c
//			  integer ia,left
//			  real biatx(4),t(11),values(4),x
//		c                 knot sequence set here....
//			  data t / 4*0.,1.,3.,4.,4*6. /
//			  do 20 ia=1,40
//				 x = float(ia)/5. - 1.
//				 do 10 left=4,7
//					call bsplvb ( t, 4, 1, x, left, biatx )
//		c
//		c           according to  bsplvb  listing,  biatx(.) now contains value
//		c           at  x  of polynomial which agrees on the interval  (t(left),
//		c           t(left+1) )  with the b-spline  b(left-4 + . ,4,t) . hence,
//		c           biatx(8-left)  now contains value of that polynomial for
//		c           b(left-4 +(8-left) ,4,t) = b(4,4,t) . since this b-spline
//		c           has support  (t(4),t(8)), it makes sense to run  left = 4,
//		c           ...,7, storing the resulting values in  values(1),...,
//		c           values(4)  for later printing.
//		c
//		   10       values(left-3) = biatx(8-left)
//		   20    print 620, x, values
//		  620 format(f10.1,4f20.8)
//												stop
//			  end


		
	}
}
