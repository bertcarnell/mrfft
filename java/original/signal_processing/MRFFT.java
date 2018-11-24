/*
 * MRFFT.java	$Revision: 1.18 $ $Date: 2004/10/14 14:31:31EDT $ $Author: PIRASW $
 *
 * Copyright 2004 Battelle Memorial Institute. All rights reserved.
 * This computer program contains the confidential information of Battelle.
 * Use or duplication in any form is subject to conditions set forth in a
 * separate license agreement. Unauthorized use or duplication may violate
 * the license agreement, state, federal and/or international laws including
 * the Copyright Laws of the United States and of other international
 * jurisdictions.
 * 
 * Warning:  This document contains sensitive security information that is 
 * controlled under the provisions of 49 CFR Part 1520.  No part of this 
 * document may be released without the written permission of the Under 
 * Secretary of Transportation for Security, Washington, DC  20590.  
 * Unauthorized release may result in civil penalty or other action.  
 * For U.S. Government agencies, public availability to be determined under 
 * 5 U.S.C. 552.
 * 
 * Author: mooneyd
 * Date: 19 October 2004
 * 
 */
package signal_processing;

import ctt.api.mathlib.Complex;

/**
 * 
 * This class implements a java version of Singleton's mixed Raduis FFT
 * as presented in "An Algorithm for Computing the Mixed Radix FFT Fourier Transform,"
 * Richard C. Singleton, IEEE Transactions on Audio and Electroacoutics,
 * Vol AU-17, No. 2, 1969.
 * 
 * @version $Revision:   $ $Date:   $
 * @author: Doug Mooney and Rob Carnell
 */

public class MRFFT {
//	subroutine fft(a,b,ntot,n,nspan,isn)
//	c  multivariate complex fourier transform, computed in place
//	c    using mixed-radix fast fourier transform algorithm.
//	c  by r. c. singleton, stanford research institute, sept. 1968
//	c  arrays a and b originally hold the real and imaginary
//	c    components of the data, and return the real and
//	c    imaginary components of the resulting fourier coefficients.
//	c  multivariate data is indexed according to the fortran
//	c    array element successor function, without limit
//	c    on the number of implied multiple subscripts.
//	c    the subroutine is called once for each variate.
//	c    the calls for a multivariate transform may be in any order.
//	c  ntot is the total number of complex data values.
//	c  n is the dimension of the current variable.
//	c  nspan/n is the spacing of consecutive data values
//	c    while indexing the current variable.
//	c  the sign of isn determines the sign of the complex
//	c    exponential, and the magnitude of isn is normally one.
//	c  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
//	c    is computed by
//	c      call fft(a,b,n1*n2*n3,n1,n1,1)
//	c      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
//	c      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
//	c  for a single-variate transform,
//	c    ntot = n = nspan = (number of complex data values), e.g.
//	c      call fft(a,b,n,n,n,1)
//	c  the data can alternatively be stored in a single complex array c
//	c    in standard fortran fashion, i.e. alternating real and imaginary
//	c    parts. then with most fortran compilers, the complex array c can
//	c    be equivalenced to a real array a, the magnitude of isn changed
//	c    to two to give correct indexing increment, and a(1) and a(2) used
//	c    to pass the initial addresses for the sequences of real and
//	c    imaginary values, e.g.
//	c       complex c(ntot)
//	c       real    a(2*ntot)
//	c       equivalence (c(1),a(1))
//	c       call fft(a(1),a(2),ntot,n,nspan,2)
//	c  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
//	c    are used for temporary storage.  if the available storage
//	c    is insufficient, the program is terminated by a stop.
//	c    maxf must be .ge. the maximum prime factor of n.
//	c    maxp must be .gt. the number of prime factors of n.
//	c    in addition, if the square-free portion k of n has two or
//	c    more prime factors, then maxp must be .ge. k-1.
//	dimension a(1),b(1)
//	c  array storage in nfac for a maximum of 15 prime factors of n.
//	c  if n has more than one square-free factor, the product of the
//	c    square-free factors must be .le. 210
//	dimension nfac(11),np(209)
//	c  array storage for maximum prime factor of 23
//	dimension at(23),ck(23),bt(23),sk(23)
//	equivalence (i,ii)
//	c  the following two constants should agree with the array dimensions.
//	maxp=209
//	C
//	maxf=23
//	C
	protected final static int MAX_N_FACTORS = 11;
	protected final static int MAX_FACTOR = 23;
	protected final static int MAX_PRODUCT = 209;
	double [] at = new double[MAX_FACTOR];
	double [] ck = new double[MAX_FACTOR];
	double [] bt = new double[MAX_FACTOR];
	double [] sk = new double[MAX_FACTOR];
	
	Complex[] x; 
	int maxf = MAX_FACTOR;
	int nTot; 
	int n; 
	int nSpan; 
	int isn;
	int inc;
	int arrayLength;
	
	int nFac[] = new int[MAX_N_FACTORS];
	int nP[] = new int[MAX_PRODUCT];
	int nt;
	int ks;
	int kSpan;
	int nn;
	int jc;
	int jf;
	double radf;
	double [] a;
	double [] b;
	
	double rad = 2*Math.PI;
//  rad=6.2831853071796
	double s72;
	double c72;
	double s120;
	
	double sd;
	double cd;
	double s1;
	double c1;
	double ak;
	double bk;
	int i;
	int jj;
	int kk;
	int k1;
	int k2;
	int k3;
	boolean cont;
	int numFacs = 0;
	int kt;
	int m;
	int k;
	int j;
	int jsquared;
	int kSpnn;
	// fftOdd and fft3
	double aj;
	double bj;
	// fftOdd and fft5
	double aa;
	double bb;
	
	

	/**
	 *  
	 * Constructor for the Mixed Radix Fast Fourier Transform<p>
	 * For most univariate transforms, nTot = n = nSpan and isn = 1<p>
	 * 
	 */
	public MRFFT(Complex[] x, int nTot, int n,int nSpan, int isn ){
		System.out.println("Constructor");
		this.x = x;
		this.arrayLength = x.length;
		this.nTot = nTot;
		this.n = n;
		this.nSpan = nSpan;
		this.isn = isn;
		
		a = new double[arrayLength+1];
		b = new double[arrayLength+1];
		for(int i=1; i<arrayLength; i++){
			a[i] = x[i].getReal();
			b[i] = x[i].getImaginary();
		}
		inc = isn;
//      inc=isn
		c72 = Math.cos(rad / 5.0);
//      c72=0.30901699437494742
		// cosine of 72 degrees
		s72 = Math.sin(rad / 5.0);
//      s72=0.95105651629515357
		// sine of 72 degrees
		s120 = Math.sqrt(0.75);
//      s120=0.86602540378443865
		// s120 is the sine of 120 degrees, the sqrt of 3/4 is a more exact calculation
		if(isn < 0){
//      if(isn .ge. 0) go to 10
			s72 = -s72;
//	        s72=-s72
			s120 = -s120;
//          s120=-s120
			rad = -rad;
//          rad=-rad
			inc = -inc;
//          inc=-inc
		} // end if(isn < 0)
		nt = inc * nTot; // 10
//      10 nt=inc*ntot
		ks = inc * nSpan;
//      ks=inc*nspan
		kSpan = ks;
//      kspan=ks
		nn = nt - inc;
//      nn=nt-inc
		jc = ks / n;
//      jc=ks/n
		radf = rad * ((double)jc) * 0.5;
//      radf=rad*float(jc)*0.5
		jf = 0;
//      jf=0
		i = 0;
//      i=0
		m = factors(nFac, n); // added by Rob
		numFacs = m; // added by Rob
	} // end MRFFT constructor
	
	
	// compute the CTTFFT of x[], assuming its length is a power of 2
//	public void mrfft() {
//		if(arrayLength < 2) return;
//		fft2();
//		return;
//	}
	
	
	
	/**
	 * factors 
	 * @param nFac
	 * @param n
	 * @return
	 * Factors small integers into primes (and powers of 4) 
	 * as part of Singleton's FFT algorithm
	 * The factorization is such that square primes 
	 * bracket the odd primes
	 * 
	 */
	public int factors(int[] nFac, int n){
		/*
		 * Here are factorization examples:
		 * 
		 *6561:  3 3 3 3 3 3 3 3
		 *34574400:  4 3 5 7 7 4 7 7 5 3 4
		 *5625:  3 5 5 5 5 3
		 *99:  3 11 3
		 *270:  3 2 3 5 3
		 *30030:  2 3 5 7 11 13
		 *390390:  13 2 3 5 7 11 13
		 */
		// Determine the factors of n
		m = 0;
//      m=0
		k = n;
//      k=n
		// first check for factors of 4^2
//		go to 20
		while (k % 16 == 0){
			m++; // 15
//  		15 m=m+1
			nFac[m] = 4;
//	        nfac(m)=4
			k = k / 16;
//      	k=k/16
		} // end while(k % 16 == 0) 	
//		20 if(k-(k/16)*16 .eq. 0) go to 15
		// Next check for other squared primes
		j = 3;
//      j=3
		jsquared = 9;	
//      jj=9
//      go to 30
		while(jsquared <= k){
			while(k % jsquared == 0) {
				m++;
//			    25 m=m+1
				nFac[m] = j;
//		        nfac(m)=j
				k = k / jsquared;
//      		k=k/jj
			} // end while(k % jsquared == 0)
//   		30 if(mod(k,jj) .eq. 0) go to 25
			j = j + 2; 
			//steps up through integers...only primes will register
//          j=j+2
			jsquared=j*j;
//		    jj=j**2
		} // end while(jsquard <= k)
//      if(jj .le. k) go to 30
		
		// Now fill in the odd primes
		if(k <= 4){
//      if(k .gt. 4) go to 40
			kt = m;
//          kt=m
			nFac[m+1] = k;
//          nfac(m+1)=k
			if(k != 1){
				m++;
//              if(k .ne. 1) m=m+1
//              go to 80
			} // end if(k !=1)
		} // end if(k <= 4)
		else{//40
			if(k % 4 == 0){
//   			40 if(k-(k/4)*4 .ne. 0) go to 50
				m++;
//		        m=m+1
				nFac[m] = 2;
//		        nfac(m)=2
				k = k / 4;
//		        k=k/4
			} // end(k % 4 == 0)
			kt = m;//50
//		    50 kt=m
			j = 2;
//		    j=2
			while(j <= k){
				if(k % j == 0){//60
//			    60 if(mod(k,j) .ne. 0) go to 70
					m++;
//			        m=m+1
					nFac[m] = j;
//			        nfac(m)=j
					k = k / j;
//			        k=k/j
				} // end if(k % j == 0)
				j=((j+1)/2) * 2 + 1;
//			    70 j=((j+1)/2)*2+1
			} // end while(j <= k)
//	        if(j .le. k) go to 60
		} // end else block
		if(kt == 0 ){//80
//      80 if(kt .eq. 0) go to 100 // 100 is the fft()
			return m;
		} // end if(kt == 0)
		else{
			j = kt;
//	        j=kt
			while(j != 0){
				m++;
//			    90 m=m+1
				nFac[m] = nFac[j];
//		        nfac(m)=nfac(j)
				j = j - 1;
//			    j=j-1
			} // end while(j != 0)
//		    if(j .ne. 0) go to 90
			return m;
		}// end  else block
	} // end factors method
	
	
	
	/**
	 * fft 
	 * Calculates the Fast Fourier Transform  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft(){	
		System.out.println("Start method fft()");
		i = 0;
		cont = true;
//		c  compute fourier transform
		if(n >= 2) {
//      if(n .lt. 2) return
			while(cont){
				sd = radf/((double)kSpan);
//				100 sd=radf/float(kspan)
				cd = 2.0 * Math.pow(Math.sin(sd), 2.0);
//		        cd=2.0*sin(sd)**2
				sd = Math.sin(sd + sd);
//		        sd=sin(sd+sd)
				kk = 1;
//	            kk=1
				i++;
//		        i=i+1
				if(this.nFac[i] != 2) fft4();  //go to 400;
//		        if(nfac(i) .ne. 2) go to 400
				else fft2();
				if (i >= numFacs) cont = false;
				System.out.println("i: " + i  +" nFac[i]: " + nFac[i] + " cont:"+ cont);
			} // end while(cont)
		} // end if(n >= 2)
		for(int i=1; i<=n; i++){
			x[i]= new Complex(a[i],b[i]);
		} // end for block
	} // end fft method
	
	
	/**
	 * fft4() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 4  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft4(){
		System.out.println("Start Method fft4()");
//		c  transform for factor of 4
		if(nFac[i] != 4) fftOdd();
//      400 if(nfac(i) .ne. 4) go to 600 // 600 is fftOdd
		else{
			double akp;
			double akm;
			double ajp;
			double ajm;
			double bkp;
			double bkm;
			double bjp;
			double bjm;
			double c2 = 0;
			double c3 = 0;
			double s2 = 0;
			double s3 = 0;
			kSpnn = kSpan;
//	        kspnn=kspan
			kSpan = kSpan / 4;
//	        kspan=kspan/4
			do {//410 
				c1 = 1.0;
//			    410 c1=1.0
				s1 = 0;
//		        s1=0
				do{// 420
					do{//420
						k1 = kk + kSpan;
//					    420 k1=kk+kspan
						k2 = k1 + kSpan;
//				        k2=k1+kspan
						k3 = k2 + kSpan;
//				        k3=k2+kspan
						akp = a[kk] + a[k2];
//				        akp=a(kk)+a(k2)
						akm = a[kk] - a[k2];
//				        akm=a(kk)-a(k2)
						ajp = a[k1] + a[k3];
//				        ajp=a(k1)+a(k3)
						ajm = a[k1] - a[k3];
//				        ajm=a(k1)-a(k3)
						a[kk] = akp + ajp;
//				        a(kk)=akp+ajp
						ajp = akp - ajp;
//				        ajp=akp-ajp
						bkp = b[kk] + b[k2];
//				        bkp=b(kk)+b(k2)
						bkm = b[kk] - b[k2];
//				        bkm=b(kk)-b(k2)
						bjp = b[k1] + b[k3];
//				        bjp=b(k1)+b(k3)
						bjm = b[k1] - b[k3];
//				        bjm=b(k1)-b(k3)
						b[kk] = bkp + bjp;
//				        b(kk)=bkp+bjp
						bjp = bkp - bjp;
//				        bjp=bkp-bjp
						if(isn < 0){// go to 450
//					    if(isn .lt. 0) go to 450
							akp = akm + bjm;
//							450 akp=akm+bjm
							akm = akm - bjm;
//					        akm=akm-bjm
							bkp = bkm - ajm;
//					        bkp=bkm-ajm
							bkm = bkm + ajm;
//					        bkm=bkm+ajm
						} // end if(isn < 0 )
						else{
							akp = akm - bjm;
//					        akp=akm-bjm
							akm = akm + bjm;
//					        akm=akm+bjm
							bkp = bkm + ajm;
//					        bkp=bkm+ajm
							bkm = bkm - ajm;
//					        bkm=bkm-ajm
						} // end else block
						if(s1 == 0){ // 460
//						     if(s1 .eq. 0) go to 460
							a[k1] = akp;
//						    460 a(k1)=akp
							b[k1] = bkp;
//					        b(k1)=bkp
							a[k2] = ajp;
//					        a(k2)=ajp
							b[k2] = bjp;
//					        b(k2)=bjp
							a[k3] = akm;
//					        a(k3)=akm
							b[k3] = bkm;
//					        b(k3)=bkm
							kk = k3 + kSpan;
//					        kk=k3+kspan
//					        if(kk .le. nt) go to 420
//					        else go to 440
						} // end if(s1 == 0)
						else{//430 
//					    if(s1 .ne. 0) go to 430
							a[k1] = akp*c1 - bkp*s1;
//						    430 a(k1)=akp*c1-bkp*s1
							b[k1] = akp*s1 + bkp*c1;
//					        b(k1)=akp*s1+bkp*c1
							a[k2] = ajp*c2 - bjp*s2;
//					        a(k2)=ajp*c2-bjp*s2
							b[k2] = ajp*s2 + bjp*c2;
//					        b(k2)=ajp*s2+bjp*c2
							a[k3] = akm*c3 - bkm*s3;
//					        a(k3)=akm*c3-bkm*s3
							b[k3] = akm*s3 + bkm*c3;
//					        b(k3)=akm*s3+bkm*c3
							kk = k3 + kSpan;
//					        kk=k3+kspan
						} // end else block
					}while(kk <= nt); // go to 420
//			        if(kk .le. nt) go to 420
					//440 
					c2 = c1 - (cd*c1 + sd*s1);
//				    440 c2=c1-(cd*c1+sd*s1)
					s1 = (sd*c1 - cd*s1) + s1;
//			        s1=(sd*c1-cd*s1)+s1
					c1 = 2.0 - (c2*c2 + s1*s1);
//			        c1=2.0-(c2**2+s1**2)
					s1 = c1 * s1;
//			        s1=c1*s1
					c1 = c1 * c2;
//			        c1=c1*c2
					c2 = c1*c1 - s1*s1;
//			        c2=c1**2-s1**2
					s2 = 2.0 * c1 * s1;
//			        s2=2.0*c1*s1
					c3 = c2*c1 - s2*s1;
//			        c3=c2*c1-s2*s1
					s3 = c2*s1 + s2*c1;
//			        s3=c2*s1+s2*c1
					kk = kk - nt + jc;
//			        kk=kk-nt+jc
				}while(kk <= kSpan) ;//go to 420
//		        if(kk .le. kspan) go to 420
				kk = kk - kSpan + inc;
//		        kk=kk-kspan+inc
			}while(kk <= jc); //go to 410
//	        if(kk .le. jc) go to 410
			if(kSpan == jc) permute();
//	        if(kspan .eq. jc) go to 800
		} // end else block
		return;
//      go to 100
	} // end fft4 method
	
	
	
	/**
	 * fftOdd() 
	 * Calculates the Fast Fourier Transform
	 * for odd factors  
	 * as part of Singleton's FFT algorithm.
	 * fft3() and fft5() are called from fftOdd
	 * 
	 */
	public void fftOdd(){
//		c  transform for odd factors
		//600
		System.out.println("Start Method fftOdd()");
		double aa;
		double bb;
		k = nFac[i];
//      600 k=nfac(i)
		kSpnn = kSpan;
//      kspnn=kspan
		kSpan = kSpan / k;
//	    kspan=kspan/k
		if(k == 3) {fft3(); System.out.println("Leave fftOdd after fft3"); return;} // go to 320 // Rob added
//      if(k .eq. 3) go to 320
		if(k == 5) {fft5(); System.out.println("Leave fftOdd after fft5"); return;} // go to 510 // Rob added
//      if(k .eq. 5) go to 510
		if(k != jf) {
//		if(k .eq. jf) go to 640
			jf = k;
//		    jf=k
			s1 = rad / ((double)k);
//		    s1=rad/float(k)
			c1 = Math.cos(s1);
//		    c1=cos(s1)
			s1 = Math.sin(s1);
//		    s1=sin(s1)
			if(jf > maxf) {
				// error finish, insufficient array storage
				isn = 0;
				throw new Error("array bounds exceeded within subroutine fftOdd");
			} // end if(jf > maxf)
//		    if(jf .gt. maxf) go to 998
			ck[jf] = 1.0;
//		    ck(jf)=1.0
			sk[jf] = 0.0;
//		    sk(jf)=0.0
			j = 1;
//		    j=1
			do{// 630
				ck[j] = ck[k]*c1 + sk[k]*s1;
//				630 ck(j)=ck(k)*c1+sk(k)*s1
				sk[j] = ck[k]*s1 - sk[k]*c1;
//				sk(j)=ck(k)*s1-sk(k)*c1
				k = k - 1;
//				k=k-1
				ck[k] = ck[j];
//				ck(k)=ck(j)
				sk[k] = -sk[j];
//				sk(k)=-sk(j)
				j = j + 1;
//				j=j+1
			}while(j < k); //end do block
//			if(j .lt. k) go to 630
		} // end if(k!=jf)
		do{//640
			do{//640
				k1 = kk;
//				640 k1=kk
				k2 = kk + kSpnn;
//				k2=kk+kspnn
				aa = a[kk];
//				aa=a(kk)
				bb = b[kk];
//				bb=b(kk)
				ak = aa;
//				ak=aa
				bk = bb;
//				bk=bb
				j = 1;
//				j=1
				k1 = k1 + kSpan;
//				k1=k1+kSpan
				do{//650
					k2 = k2 - kSpan;
//					650 k2=k2-kSpan
					j = j + 1;
//					j=j+1
					at[j] = a[k1] + a[k2];
//					at(j)=a(k1)+a(k2)
					ak = at[j] + ak;
//					ak=at(j)+ak
					bt[j] = b[k1] + b[k2];
//					bt(j)=b(k1)+b(k2)
					bk = bt[j] + bk;
//					bk=bt(j)+bk
					j = j + 1;
//					j=j+1
					at[j] = a[k1] - a[k2];
//					at(j)=a(k1)-a(k2)
					bt[j] = b[k1] - b[k2];
//					bt(j)=b(k1)-b(k2)
					k1 = k1 + kSpan;
//					k1=k1+kSpan
				}while(k1 < k2); //end do block
//				if(k1 .lt. k2) go to 650
				a[kk] = ak;
//				a(kk)=ak
				b[kk] = bk;
//				b(kk)=bk
				k1 = kk;
//				k1=kk
				k2 = kk + kSpnn;
//				k2=kk+kspnn
				j = 1;
//				j=1
				do{//660
					k1 = k1 + kSpan;
//					660 k1=k1+kSpan
					k2 = k2 - kSpan;
//					k2=k2-kSpan
					jj = j;
//					jj=j
					ak = aa;
//					ak=aa
					bk = bb;
//					bk=bb
					aj = 0.0;
//					aj=0.0
					bj = 0.0;
//					bj=0.0
					k = 1;
//					k=1
					do{//670
						k = k + 1;
//						670 k=k+1
						ak = at[k]*ck[jj] + ak;
//						ak=at(k)*ck(jj)+ak
						bk = bt[k]*ck[jj] + bk;
//						bk=bt(k)*ck(jj)+bk
						k = k + 1;
//						k=k+1
						aj = at[k]*sk[jj] + aj;
//						aj=at(k)*sk(jj)+aj
						bj = bt[k]*sk[jj] + bj;
//						bj=bt(k)*sk(jj)+bj
						jj = jj + j;
//						jj=jj+j
						if(jj > jf) jj = jj - jf;
//						if(jj .gt. jf) jj=jj-jf
					}while(k < jf); // end do block
//					if(k .lt. jf) go to 670
					k = jf - j;
//					k=jf-j
					a[k1] = ak - bj;
//					a(k1)=ak-bj
					b[k1] = bk + aj;
//					b(k1)=bk+aj
					a[k2] = ak + bj;
//					a(k2)=ak+bj
					b[k2] = bk - aj;
//					b(k2)=bk-aj
					j = j + 1;
//					j=j+1
				}while(j < k); // end do block
//				if(j .lt. k) go to 660
				kk = kk + kSpnn;
//				kk=kk+kspnn
			}while(kk <= nn); // end do block
//			if(kk .le. nn) go to 640
			kk = kk - nn;
//			kk=kk-nn
		}while(kk <= kSpan); //end do block
//		if(kk .le. kSpan) go to 640
		rotate(); // rotate factor applied to all factors except 2 and 4
		System.out.println("fftOdd Complete");
		return;
	} // end fftOdd method
	
	
	/**
	 * fft3() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 3  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft3(){
//		c  transform for factor of 3 (optional code)
		System.out.println("Start Method fft3()");
		do{ // 320
			do{ // 320 
				k1 = kk + kSpan;
//  		    320 k1=kk+kspan
				k2 = k1 + kSpan;
//		        k2=k1+kspan
				ak = a[kk];
//			    ak=a(kk)
				bk = b[kk];
//			    bk=b(kk)
				aj = a[k1] + a[k2];
//		        aj=a(k1)+a(k2)
				bj = b[k1] + b[k2];
//		        bj=b(k1)+b(k2)
				a[kk] = ak + aj;
//		        a(kk)=ak+aj
				b[kk] = bk + bj;
//		        b(kk)=bk+bj
				ak = -0.5*aj + ak;
//		        ak=-0.5*aj+ak
				bk = -0.5*bj + bk;
//		        bk=-0.5*bj+bk
				aj = (a[k1] - a[k2])*s120;
//		        aj=(a(k1)-a(k2))*s120
				bj = (b[k1] - b[k2])*s120;
//		        bj=(b(k1)-b(k2))*s120
				a[k1] = ak - bj;
//		        a(k1)=ak-bj
				b[k1] = bk + aj;
//		        b(k1)=bk+aj
				a[k2] = ak + bj;
//		        a(k2)=ak+bj
				b[k2] = bk - aj;
//		        b(k2)=bk-aj
				kk = k2 + kSpan;
//		        kk=k2+kspan
			}while(kk < nn); // go to 320
//		    if(kk .lt. nn) go to 320
			kk=kk-nn;
//		    kk=kk-nn
		}while(kk <= kSpan);// go to 320
//      if(kk .le. kspan) go to 320
		rotate(); //go to 700
		System.out.println("fft3 Complete");
		return;
	}// end fft3 method
	
	
	/**
	 * fft5() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 5  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft5(){
//		c  transform for factor of 5 (optional code)
//		510 
		System.out.println("Start Method fft5()");
		double c2;
		double s2;
		double akp;
		double akm;
		double bkp;
		double bkm;
		double ajp;
		double ajm;
		double bjp;
		double bjm;
		int k4;
		c2 = c72*c72 - s72*s72;
//      510 c2=c72**2-s72**2
		s2 = 2.0*c72*s72;
//      s2=2.0*c72*s72
		do{ //520
			do{ //520
				k1 = kk + kSpan;
//			    520 k1=kk+kspan
				k2 = k1 + kSpan;
//		        k2=k1+kspan
				k3 = k2 + kSpan;
//		        k3=k2+kspan
				k4 = k3 + kSpan;
//		        k4=k3+kspan
				akp = a[k1] + a[k4];
//		        akp=a(k1)+a(k4)
				akm = a[k1] - a[k4];
//		        akm=a(k1)-a(k4)
				bkp = b[k1] + b[k4];
//		        bkp=b(k1)+b(k4)
				bkm = b[k1] - b[k4];
//		        bkm=b(k1)-b(k4)
				ajp = a[k2] + a[k3];
//		        ajp=a(k2)+a(k3)
				ajm = a[k2] - a[k3];
//		        ajm=a(k2)-a(k3)
				bjp = b[k2] + b[k3];
//		        bjp=b(k2)+b(k3)
				bjm = b[k2] - b[k3];
//		        bjm=b(k2)-b(k3)
				aa = a[kk];
//		        aa=a(kk)
				bb = b[kk];
//		        bb=b(kk)
				a[kk] = aa + akp + ajp;
//		        a(kk)=aa+akp+ajp
				b[kk] = bb + bkp + bjp;
//		        b(kk)=bb+bkp+bjp
				ak = akp*c72 + ajp*c2 + aa;
//		        ak=akp*c72+ajp*c2+aa
				bk = bkp*c72 + bjp*c2 + bb;
//		        bk=bkp*c72+bjp*c2+bb
				aj = akm*s72 + ajm*s2;
//		        aj=akm*s72+ajm*s2
				bj = bkm*s72 + bjm*s2;
//		        bj=bkm*s72+bjm*s2
				a[k1] = ak - bj;
//		        a(k1)=ak-bj
				a[k4] = ak + bj;
//		        a(k4)=ak+bj
				b[k1] = bk + aj;
//		        b(k1)=bk+aj
				b[k4] = bk - aj;
//		        b(k4)=bk-aj
				ak = akp*c2 + ajp*c72 + aa;
//		        ak=akp*c2+ajp*c72+aa
				bk = bkp*c2 + bjp*c72 + bb;
//		        bk=bkp*c2+bjp*c72+bb
				aj = akm*s2 - ajm*s72;
//		        aj=akm*s2-ajm*s72
				bj = bkm*s2 - bjm*s72;
//		        bj=bkm*s2-bjm*s72
				a[k2] = ak - bj;
//		        a(k2)=ak-bj
				a[k3] = ak + bj;
//		        a(k3)=ak+bj
				b[k2] = bk + aj;
//		        b(k2)=bk+aj
				b[k3] = bk - aj;
//		        b(k3)=bk-aj
				kk = k4 + kSpan;
//			    kk=k4+kspan
			}while(kk < nn);// go to 520
//	        if(kk .lt. nn) go to 520
			kk = kk - nn;
//	        kk=kk-nn
		}while(kk <= kSpan);// go to 520
//      if(kk .le. kspan) go to 520
//		go to 700 // 700 is rotate
		rotate();
		System.out.println("fft5 Complete");
		return;
	} // end fft5 method

	/**
	 * fft2() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 2 (not 4)  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft2(){
//		c  transform for factor of 2 (including rotation factor)
		System.out.println("Start Method fft2()");
		kSpan = kSpan / 2;
//      kspan=kspan/2
		k1 = kSpan + 2;
//      k1=kspan+2
		do { // 210
			do { // 210
				//System.out.println("fft2, kk = " + kk);
				k2 = kk + kSpan;
//			    210 k2=kk+kspan
				ak = a[k2];
//      		ak=a(k2)
				bk = b[k2];
//			    bk=b(k2)
				a[k2] = a[kk] - ak;
//			    a(k2)=a(kk)-ak
				b[k2] = b[kk] - bk;
//			    b(k2)=b(kk)-bk
				a[kk] = a[kk] + ak;
//			    a(kk)=a(kk)+ak
				b[kk] = b[kk] + bk;
//			    b(kk)=b(kk)+bk
				kk = k2 + kSpan;
//			    kk=k2+kspan
			} while (kk <= nn);
//	        if(kk .le. nn) go to 210
			kk = kk - nn;
//	        kk=kk-nn
		} while (kk <= jc);
//      if(kk .le. jc) go to 210
		if(kk > kSpan) {
			permute();  // go to 800 // 800 is permute
//      	if(kk .gt. kspan) go to 800
			return;
		} // end if(kk > kSpan)
		do{ // 220
			c1 = 1.0 - cd;
//			220 c1=1.0-cd
			s1 = sd;
//	        s1=sd
			do{ // 230
				do{ // 230
					do{ // 230
						k2 = kk + kSpan;
//  				    230 k2=kk+kspan
						//System.out.println("fft2, kk = " + kk + " k2 = " + k2);
						ak = a[kk] - a[k2];
//	   				    ak=a(kk)-a(k2)
						bk = b[kk] - b[k2];
//					    bk=b(kk)-b(k2)
						a[kk] = a[kk] + a[k2];
//					    a(kk)=a(kk)+a(k2)
						b[kk] = b[kk] + b[k2];
//					    b(kk)=b(kk)+b(k2)
						a[k2] = c1*ak - s1*bk;
//					    a(k2)=c1*ak-s1*bk
						b[k2] = s1*ak + c1*bk;
//					    b(k2)=s1*ak+c1*bk
						kk= k2 + kSpan;
//					    kk=k2+kspan
					}while(kk < nt) ;
//			        if(kk .lt. nt) go to 230
					k2 = kk - nt;
//				    k2=kk-nt
					c1 = -c1;
//			        c1=-c1
					kk = k1 - k2;
//			        kk=k1-k2
				}while(kk > k2) ;
//		        if(kk .gt. k2) go to 230
				ak = c1 - (cd*c1 + sd*s1);
//			    ak=c1-(cd*c1+sd*s1)
				s1 = (sd*c1 - cd*s1) + s1;
//			    s1=(sd*c1-cd*s1)+s1
				c1 = 2.0 - (Math.pow(ak, 2) + Math.pow(s1, 2));
//			    c1=2.0-(ak**2+s1**2)
				s1 = c1*s1;
//			    s1=c1*s1
				c1 = c1*ak;
//			    c1=c1*ak
				kk = kk + jc;
//			    kk=kk+jc
			}while(kk < k2) ;
//		    if(kk .lt. k2) go to 230
			k1 = k1 + inc + inc;
//	        k1=k1+inc+inc
			kk = (k1 - kSpan) / 2 + jc;
//    	    kk=(k1-kspan)/2+jc
		}while(kk <= (jc + jc)) ;
//      if(kk .le. jc+jc) go to 220
//      go to 100 // return to fft()
		System.out.println("fft2() Complete");
	} // end fft2 method


	/**
	 * Multiplies by a rotation factor as part of Singleton's FFT algorithm.<p>
	 * This method is not used with factors of 2 and 4 because the rotation is
	 * included in those methods.<p>
	 * The FFT of x is given by Tx = PF[1]F[2]..F[m] where P is the permutation step
	 * and each F[i] corresponds to the ith factor.  N=n[1]n[2]...n[m].  Each F[i]=R[i]W[i] where W[i] is the ith
	 * transform step and R[i] is the ith rotation step.
	 * 
	 */
	public void rotate() {
//		c  multiply by rotation factor (except for factors of 2 and 4)
		
		//700 
		System.out.println("Start Method rotate");
		double c2;
		double s2;
		if(i != m){
//		    700 if(i .eq. m) go to 800
			kk = jc + 1;
//	        kk=jc+1
			do{//710 
				c2 = 1.0 - cd;
//			    710 c2=1.0-cd
				s1 = sd;
//		        s1=sd
				do{// 720 
					c1 = c2;
//				    720 c1=c2
					s2 = s1;
//			        s2=s1
					kk = kk + kSpan;
//			        kk=kk+kspan
					do{ // 730
						do{ //730 
							ak = a[kk];
//						    730 ak=a(kk)
							a[kk] = c2*ak - s2*b[kk];
//					        a(kk)=c2*ak-s2*b(kk)
							b[kk] = s2*ak + c2*b[kk];
//					        b(kk)=s2*ak+c2*b(kk)
							kk = kk + kSpnn;
//					        kk=kk+kspnn
						}while(kk <= nt);// go to 730
//				        if(kk .le. nt) go to 730
						ak = s1*s2;
//				        ak=s1*s2
						s2 = s1*c2 + c1*s2;
//				        s2=s1*c2+c1*s2
						c2 = c1*c2 - ak;
//				        c2=c1*c2-ak
						kk = kk - nt + kSpan;
//				        kk=kk-nt+kspan
					}while(kk <= kSpnn);// go to 730
//			        if(kk .le. kspnn) go to 730
					c2= c1 - (cd*c1 + sd*s1);
//			        c2=c1-(cd*c1+sd*s1)
					s1= s1 + (sd*c1 - cd*s1);
//			        s1=s1+(sd*c1-cd*s1)
					c1 = 2.0 - (c2*c2 + s1*s1);
//			        c1=2.0-(c2**2+s1**2)
					s1 = c1*s1;
//			        s1=c1*s1
					c2 = c1*c2;
//			        c2=c1*c2
					kk = kk - kSpnn + jc;
//			        kk=kk-kspnn+jc
				}while(kk <= kSpan);// go to 720
//		        if(kk .le. kspan) go to 720
				kk = kk - kSpan + jc + inc;
//		        kk=kk-kspan+jc+inc
			}while(kk <= (jc + jc)); //go to 710
//	        if(kk .le. jc+jc) go to 710
			System.out.println("End method rotate, without permute");
			return; //go to 100
		} // end if (i != m)
		else{
			permute(); // 800
			System.out.println("End method rotate, with permute");

//			for(int zz=0; zz<arrayLength; zz++){
//				System.out.println(a[zz] + "\t" + b[zz]);
//			}
			
			return;
		} // end else block
	} // end method rotate
	
	
	
	/**
	 *  
	 * Re-arranges the results of the Fast Fourier Transform
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void permute() {
//		c  permute the results to normal order---done in two stages
//		c  permutation for square factors of n
		System.out.println("Start Method permute");
		
//		for(int zz=0; zz<arrayLength; zz++){
//			System.out.println(a[zz] + "\t" + b[zz]);
//		}
		
		int ii;
//		800 np(1)=ks
		nP[1] = ks;
//		if(kt .eq. 0) go to 890
		if(kt != 0){ 
			k = kt + kt + 1;
//          k=kt+kt+1
			if(m < k) k = k - 1;
//          if(m .lt. k) k=k-1
			j = 1;
//          j=1
			nP[k+1] = jc;
//          np(k+1)=jc
			do { // 810
				nP[j+1] = nP[j] / nFac[j];
//			    810 np(j+1)=np(j)/nfac(j)
				nP[k] = nP[k+1] * nFac[j];
//		        np(k)=np(k+1)*nfac(j)
				j = j + 1;
//		        j=j+1
				k = k - 1;
//		        k=k-1
			}while(j < k); //go to 810
//	        if(j .lt. k) go to 810
			k3 = nP[k+1];
//	        k3=np(k+1)
			kSpan = nP[2];
//	        kspan=np(2)
			kk = jc + 1;
//	        kk=jc+1
			k2 = kSpan + 1;
//	        k2=kspan+1
			j = 1;
//	        j=1
			if(n == nTot) {
//		    if(n .ne. ntot) go to 850 // multivariate transform
//				c  permutation for single-variate transform (optional code)
				boolean b820 = true;
				boolean loop = true;
				start:
				  while(loop){
				    if(b820) {
				      do{//820
				        ak = a[kk];
//				        820 ak=a(kk)
				        a[kk] = a[k2];
//				        a(kk)=a(k2)
				        a[k2] = ak;
//				        a(k2)=ak
				        bk = b[kk];
//				        bk=b(kk)
				        b[kk] = b[k2];
//				        b(kk)=b(k2)
				        b[k2] = bk;
//				        b(k2)=bk
				        kk = kk + inc;
//				        kk=kk+inc
				        k2 = kSpan + k2;
//				        k2=kspan+k2
				      }while(k2 < ks);
//				      if(k2 .lt. ks) go to 820
				    } // end if(b820)
				    do{//830
				      k2 = k2 - nP[j];
//				      k2=k2-np(j)
				      j = j + 1;
//				      j=j+1
				      k2 = nP[j+1] + k2;
//				      k2=np(j+1)+k2
				    }while(k2 > nP[j]);
//				    if(k2 .gt. np(j)) go to 830
				    j = 1;
				    do{//840
				      if(kk < k2) {b820 = true; continue start;}
//				      if(kk .lt. k2) go to 820
				      kk = kk + inc;
//				      kk=kk+inc
				      k2 = kSpan + k2;
//				      k2=kspan+k2
				    }while(k2 < ks);
//				    if(k2 .lt. ks) go to 840
				    if(kk < ks) {b820=false; continue start;}
//				    if(kk .lt. ks) go to 830
				    jc = k3;
//				    jc=k3
				    break start;
				    //go to 890
				  }// end while(loop)
			}// end if(n == nTot)	
				
////////////  doug's code 820 through 840 ///////////////				
//			    boolean b830 = false;
//				boolean b840 = false;
//				//c  permutation for single-variate transform (optional code)
//				do{
//					do{
//						if(b840 == false){
//							do{ //820  
//								ak = a[kk];
//								a[kk] = a[k2];
//								a[k2] = ak;
//								bk = b[kk];
//								b[kk] = b[k2];
//								b[k2] = bk;
//								kk = kk + inc;
//								k2 = kSpan + k2;
//							}while(k2 < ks); //go to 820
//						} // end if(b840 == false)
//						if(b830 == false){
//							do{ //830
//								k2 = k2 - nP[j];
//								j = j + 1;
//								k2 = nP[j+1] + k2;
//							}while(k2 > nP[j]); //go to 830
//							j = 1;
//						} // end if(b830 == false)
//						b840 = false;
//						b830 = false;						
//					}while(kk < k2);
//					kk = kk + inc;
//					k2 = kSpan + k2;
//					if(k2 < ks) {
//						b830 = true;
//						b840 = true;
//					} // end if(k2 <- ks)
//					else if(kk < ks) {
//						b840 = true;
//						b830 = false;
//					} // end else if(kk < ks)
//				}while( (k2 < ks ) || (kk < ks) ); // go to 840
//			} // end if(N == nTot)
//			jc = k3;
//			//go to 890
//////////////////////////////////////////////////				
			else throw new Error("n != nTot, a multivariate transform has been called");
			
//////////////////////////////////////////////////
//			c  permutation for multivariate transform
//			850 k=kk+jc
//			860 ak=a(kk)
//			a(kk)=a(k2)
//			a(k2)=ak
//			bk=b(kk)
//			b(kk)=b(k2)
//			b(k2)=bk
//			kk=kk+inc
//			k2=k2+inc
//			if(kk .lt. k) go to 860
//			kk=kk+ks-jc
//			k2=k2+ks-jc
//			if(kk .lt. nt) go to 850
//			k2=k2-nt+kSpan
//			kk=kk-nt+jc
//			if(k2 .lt. ks) go to 850
//			870 k2=k2-np(j)
//			j=j+1
//			k2=np(j+1)+k2
//			if(k2 .gt. np(j)) go to 870
//			j=1
//			880 if(kk .lt. k2) go to 850
//			kk=kk+jc
//			k2=kSpan+k2
//			if(k2 .lt. ks) go to 880
//			if(kk .lt. ks) go to 870
//			jc=k3
////////////////////////////////////////////////////
				
		} // end if(kt != 0)
		
		if(2*kt + 1 >= m) {System.out.println("End method permute, first position"); return;}
//		890 if(2*kt+1 .ge. m) return
		kSpnn = nP[kt+1];
//		kspnn=np(kt+1)

//		c  permutation for square-free factors of n
		j = m - kt;
//      j=m-kt
		nFac[j+1] = 1;
//      nfac(j+1)=1
		do { //900 
			nFac[j] = nFac[j]*nFac[j+1];
//          900 nfac(j)=nfac(j)*nfac(j+1)
			j = j - 1;
//          j=j-1
		}while(j != kt); 
//      if(j .ne. kt) go to 900
		kt = kt + 1;
//      kt=kt+1
		nn = nFac[kt] - 1;
//      nn=nfac(kt)-1
		if(nn > MAX_PRODUCT) {
//      if(nn .gt. maxp) go to 998
			throw new Error("array bounds exceeded within subroutine permute");
		}
		jj = 0;
//      jj=0
		j = 0;
//      j=0

		boolean b902 = false;
		boolean b904 = false;
		boolean loop = true;

//		go to 906
		start:
		  while(loop) {
		    if(b902) {
		      jj = jj - k2;
//		      jj=jj-k2
		      k2 = kk;
//		      k2=kk
		      k = k + 1;
//		      k=k+1
		      kk = nFac[k];
//		      kk=nfac(k)
		    } // end if(b902)
		    if(b904) {
		      jj = kk + jj;
//		      jj=kk+jj
		      if(jj >= k2) {b902 = true; b904 = true; continue start;}
//		      if(jj .ge. k2) go to 902
		      nP[j] = jj;
//		      np(j)=jj
		    } // end if(b904)
		    k2 = nFac[kt];
//		    906 k2=nfac(kt)
		    k = kt + 1;
//		    k=kt+1
		    kk = nFac[k];
//		    kk=nfac(k)
		    j = j + 1;
//		    j=j+1
		    if(j <= nn) {b902 = false; b904 = true; continue start;}
		    else break start;
//		    if(j .le. nn) go to 904
		  }// end while(loop)


////////////////////////////////////////////////////
//		906 
//		boolean first = true;
//		boolean b904 = false;
//		
//		//go to 906
//		do{
//			if(!first){
//				do{
//					
//					// 902 
//					if(b904){
//						jj = jj-k2;
//						k2 = kk;
//						k = k + 1;
//						kk = nFac[k];}
//					// 904 
//					jj = kk + jj;
//					b904 = true;
//				}while(jj >= k2);// go to 902
//				nP[j] = jj;
//			}
//			//906 
//			first = false;
//			b904 = false;
//			k2 = nFac[kt];
//			k = kt + 1;
//			kk = nFac[k];
//			j = j + 1;
//		}while(j <= nn);// go to 904
////////////////////////////////////////////////////		
		
		
//		c  determine the permutation cycles of length greater than
		j = 0;
//		j=0
		boolean b910 = false;
		loop = true;

		start:
//		go to 914
		  while(loop) {
		    if(b910) {
		      do{
		        k = kk;
//		        910 k=kk
		        kk = nP[k];
//		        kk=np(k)
		        nP[k] = -kk;
//		        np(k)=-kk
		      }while(kk != j);
//		      if(kk .ne. j) go to 910
		      k3 = kk;
//		      k3=kk
		    } // end if(b910)
		    do{// 914
		      do{ // 914
		        j = j + 1;
//		        914 j=j+1
		        kk = nP[j];
//		        kk=np(j)
		      }while(kk < 0);
//		      if(kk .lt. 0) go to 914
		      if(kk != j) {b910 = true; continue start;}
//		      if(kk .ne. j) go to 910
		      nP[j] = -j;
// 		      np(j)=-j
		    }while(j != nn);
//		    if(j .ne. nn) go to 914
		    maxf = inc*maxf;
//		    maxf=inc*maxf
		    break start;
		  }// end while(loop)
		  

//		c  reorder a and b, following the permutation cycles
//		go to 950
		j = k3 + 1;
//		950 j=k3+1
		nt = nt - kSpnn;
//		nt=nt-kspnn
		ii = nt - inc + 1;
//		ii=nt-inc+1
		if(nt >= 0){ //go to 924
//		if(nt .ge. 0) go to 924 else return
			do{ // 924
				do{ // 924	
					do{ // 924	
						j = j - 1;
//						924 j=j-1
					}while(nP[j] < 0);// go to 924
//					if(np(j) .lt. 0) go to 924
					jj = jc;
//					jj=jc
					do{ // 926
						kSpan = jj;
//						926 kSpan=jj
						if(jj > maxf) kSpan = maxf;
//						if(jj .gt. maxf) kSpan=maxf
						jj = jj - kSpan;
//						jj=jj-kSpan
						k = nP[j];
//						k=np(j)
						kk = jc*k + ii + jj;
//						kk=jc*k+ii+jj
						k1 = kk + kSpan;
//						k1=kk+kSpan
						k2 = 0;
//						k2=0
						do{ // 928
							k2 = k2 + 1;
//							928 k2=k2+1
							at[k2] = a[k1];
//							at(k2)=a(k1)
							bt[k2] = b[k1];
//							bt(k2)=b(k1)
							k1 = k1 - inc;
//							k1=k1-inc
						}while(k1 != kk);// go to 928
//						if(k1 .ne. kk) go to 928
						do{ // 932
							k1 = kk + kSpan;
//							932 k1=kk+kSpan
							k2 = k1 - jc*(k + nP[k]);
//							k2=k1-jc*(k+np(k))
							k = -nP[k];
//							k=-np(k)
							do{ // 936
								a[k1] = a[k2];
//								936 a(k1)=a(k2)
								b[k1] = b[k2];
//								b(k1)=b(k2)
								k1 = k1 - inc;
//								k1=k1-inc
								k2 = k2 - inc;
//								k2=k2-inc
							}while(k1 != kk);// go to 936
//							if(k1 .ne. kk) go to 936
							kk = k2;
//							kk=k2
						}while(k != j);// go to 932
//						if(k .ne. j) go to 932
						k1 = kk + kSpan;
//						k1=kk+kSpan
						k2 = 0;
//						k2=0
						do{ // 940
							k2 = k2 + 1;
//							940 k2=k2+1
							a[k1] = at[k2];
//							a(k1)=at(k2)
							b[k1] = bt[k2];
//							b(k1)=bt(k2)
							k1 = k1 - inc;
//							k1=k1-inc
						}while(k1 != kk);// go to 940
//						if(k1 .ne. kk) go to 940
					}while(jj != 0);// go to 926
//					if(jj .ne. 0) go to 926
				} while(j != 1); // go to 924
//				if(j .ne. 1) go to 924
				//950 
				j = k3 + 1;
				nt = nt - kSpnn;
				ii = nt - inc + 1;
			}while(nt >= 0);// go to 924
		} // end if(nt >=0)
		System.out.println("Method permute Complete, final lines");
		return;
	} // end method permute

	
} // end class MRFFT