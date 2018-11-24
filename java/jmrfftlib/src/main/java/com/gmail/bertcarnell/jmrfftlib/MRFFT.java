package com.gmail.bertcarnell.jmrfftlib;

import org.apache.commons.math3.complex.Complex;
import java.util.ArrayList;
import java.util.List;

/**
 * This class implements a java version of Singleton's mixed Radix FFT
 * as presented in "An Algorithm for Computing the Mixed Radix FFT Fourier Transform,"
 * Richard C. Singleton, IEEE Transactions on Audio and Electroacoutics,
 * Vol AU-17, No. 2, 1969.
 */
public class MRFFT 
{
	protected final static int MAX_N_FACTORS = 11;
	protected final static int MAX_FACTOR = 23;
	protected final static int MAX_PRODUCT = 209;
        /**
         * arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
         * are used for temporary storage.  if the available storage
         * is insufficient, the program is terminated by a stop.
         * maxf must be .ge. the maximum prime factor of n.
         * maxp must be .gt. the number of prime factors of n.
         * in addition, if the square-free portion k of n has two or
         * more prime factors, then maxp must be .ge. k-1.
         */
	private double[] at = new double[MAX_FACTOR];
	private double[] ck = new double[MAX_FACTOR];
	private double[] bt = new double[MAX_FACTOR];
	private double[] sk = new double[MAX_FACTOR];
	private int nP[] = new int[MAX_PRODUCT];
	
	private ArrayList<Complex> y; 
        /**
         * ntot is the total number of complex data values.
         */
	private int nTot; 
        /**
         * n is the dimension of the current variable.
         */
	private int n;
        /**
         * nspan/n is the spacing of consecutive data values while indexing the current variable.
         */
	//private int nSpan; 
        /**
         * the sign of isn determines the sign of the complex exponential, and the magnitude of isn is normally one.
         * used in the constructor and fft4()
         */
	private int isn;
	private int inc;
	
        /**
         * array storage in nfac for a maximum of 15 prime factors of n.
         * if n has more than one square-free factor, the product of the
         * square-free factors must be .le. 210
         */
	private int nFac[] = new int[MAX_N_FACTORS];
	private int nt;
	private int ks;
	private int kSpan;
	private int nn;
	private int jc;
	//private int jf;
	private double radf;
        // arrays a and b originally hold the real and imaginary
        // components of the data, and return the real and
        // imaginary components of the resulting fourier coefficients.
        
	private double [] a;
	private double [] b;
	
        // variables that are set by the constructor and used in the fft's
	private double rad;
	private double s72;
	private double c72;
	private double s120;
	
	private double sd;
	private double cd;
	private double s1;
	private double c1;
	private double ak;
	private double bk;
	//private int jj;
	private int kk;
	private int k1;
	private int k2;
	private int k3;
	//private boolean cont;
	private int numFacs;
	private int kt;
	private int m;
	private int k;
	private int j;
	//private int jsquared;
	private int kSpnn;

	/**
	 * Constructor for the Mixed Radix Fast Fourier Transform<p>
	 * For most univariate transforms, nTot = n = nSpan and isn = 1<p>
	 */
	public MRFFT(List<Complex> x, int nTot, int n, int nSpan, int isn)
        {
		//System.out.println("Constructor");
                int arrayLength = x.size();
                y = new ArrayList<Complex>(arrayLength);
                for (int i = 0; i < arrayLength; i++)
                {
                    y.add(Complex.ZERO);
                }
		this.nTot = nTot;
		this.n = n;
		//this.nSpan = nSpan;
		this.isn = isn;
		
		a = new double[arrayLength+1];
		b = new double[arrayLength+1];
		for (int i = 1; i < arrayLength; i++){
			a[i] = x.get(i).getReal();
			b[i] = x.get(i).getImaginary();
		}
		inc = isn;
                rad = 2*Math.PI;
		// cosine of 72 degrees
		c72 = Math.cos(rad / 5.0);
		// sine of 72 degrees
		s72 = Math.sin(rad / 5.0);
		// s120 is the sine of 120 degrees, the sqrt of 3/4 is a more exact calculation
		s120 = Math.sqrt(0.75);
		if (isn < 0)
                {
			s72 = -s72;
			s120 = -s120;
			rad = -rad;
			inc = -inc;
		} // end if (isn < 0)
		nt = inc * nTot; // 10
		ks = inc * nSpan;
		kSpan = ks;
		nn = nt - inc;
		jc = ks / n;
		radf = rad * ((double)jc) * 0.5;
		//jf = 0;
		//i = 0;
		m = factors(nFac, n);
		numFacs = m;
	} // end MRFFT constructor
	
	// compute the CTTFFT of x[], assuming its length is a power of 2
//	public void mrfft() {
//		if (arrayLength < 2) return;
//		fft2();
//		return;
//	}
        
        public List<Complex> getSignal()
        {
            return y;
        }
        
        public Complex getSignal(int i)
        {
            return y.get(i);
        }
	
	/**
	 * Factors small integers into primes (and powers of 4) as part of Singleton's FFT algorithm
	 * The factorization is such that square primes bracket the odd primes
	 * @param nFac
	 * @param n
	 * @return
	 */
	public int factors(int[] nFac, int n)
        {
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
            int jsquared;
		// Determine the factors of n
		m = 0;
		k = n;
		// first check for factors of 4^2
		while (k % 16 == 0){
			m++; // 15
			nFac[m] = 4;
			k = k / 16;
		} // end while(k % 16 == 0) 	
		// Next check for other squared primes
		j = 3;
		jsquared = 9;	
		while (jsquared <= k)
                {
			while (k % jsquared == 0) 
                        {
				m++;
				nFac[m] = j;
				k = k / jsquared;
			} // end while(k % jsquared == 0)
			j = j + 2; 
			//steps up through integers...only primes will register
			jsquared = j * j;
		} // end while(jsquared <= k)
		
		// Now fill in the odd primes
		if (k <= 4)
                {
			kt = m;
			nFac[m+1] = k;
			if (k != 1)
                        {
				m++;
			} // end if (k !=1)
		} // end if (k <= 4)
		else
                {//40
			if (k % 4 == 0)
                        {
				m++;
				nFac[m] = 2;
				k = k / 4;
			} // end(k % 4 == 0)
			kt = m;//50
			j = 2;
			while (j <= k)
                        {
				if (k % j == 0)
                                {//60
					m++;
					nFac[m] = j;
					k = k / j;
				} // end if (k % j == 0)
                                // this is not equivalent to j + 2 because of intger division
				j = ((j+1)/2) * 2 + 1;
			} // end while(j <= k)
		} // end else block
		if (kt == 0 )
                {//80
			return m;
		} // end if (kt == 0)
		else
                {
			j = kt;
			while (j != 0)
                        {
				m++;
				nFac[m] = nFac[j];
				j = j - 1;
			} // end while(j != 0)
			return m;
		}// end  else block
	} // end factors method
	
	/**
	 * Calculates the Fast Fourier Transform as part of Singleton's FFT algorithm
        * subroutine fft(a,b,ntot,n,nspan,isn)
        * multivariate complex fourier transform, computed in place using mixed-radix 
        * fast fourier transform algorithm.
        * by r. c. singleton, stanford research institute, sept. 1968
        * 
        * multivariate data is indexed according to the fortran array element successor function, without limit
        * on the number of implied multiple subscripts.  the subroutine is called once for each variate.
        * the calls for a multivariate transform may be in any order.
        * 
        * a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3) is computed by
        * call fft(a,b,n1*n2*n3,n1,n1,1)
        * call fft(a,b,n1*n2*n3,n2,n1*n2,1)
        * call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
        * 
        * for a single-variate transform, ntot = n = nspan = (number of complex data values), e.g.
        * call fft(a,b,n,n,n,1)
	 */
	public void fft()
        {	
		//System.out.println("Start method fft()");
		int factorNumber = 0;
		boolean cont = true;
		if (n >= 2) 
                {
			while (cont)
                        {
				sd = radf/((double)kSpan);
				cd = 2.0 * Math.pow(Math.sin(sd), 2.0);
				sd = Math.sin(sd + sd);
				kk = 1;
				factorNumber++;
				if (this.nFac[factorNumber] != 2)
                                {
                                    if (nFac[factorNumber] != 4)
                                    {
                                        fftOdd(nFac[factorNumber], factorNumber);
                                    }
                                    else
                                    {
                                        fft4();  //go to 400;
                                    }
                                }
                                else
                                {
                                    fft2();
                                }
				if (factorNumber >= numFacs)
                                {
                                    cont = false;
                                }
				//System.out.println("i: " + i  +" nFac[i]: " + nFac[i] + " cont:"+ cont);
			} // end while(cont)
		} // end if (n >= 2)
		for (int iOutput = 1; iOutput <= n; iOutput++)
                {
			y.set(iOutput, Complex.valueOf(a[iOutput], b[iOutput]));
		} // end for block
	} // end fft method
	
	
	/**
	 * fft4() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 4  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft4()
        {
		//System.out.println("Start Method fft4()");
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
                kSpan = kSpan / 4;
                do
                {//410 
                        c1 = 1.0;
                        s1 = 0;
                        do
                        {// 420
                                do
                                {//420
                                        k1 = kk + kSpan;
                                        k2 = k1 + kSpan;
                                        k3 = k2 + kSpan;
                                        akp = a[kk] + a[k2];
                                        akm = a[kk] - a[k2];
                                        ajp = a[k1] + a[k3];
                                        ajm = a[k1] - a[k3];
                                        a[kk] = akp + ajp;
                                        ajp = akp - ajp;
                                        bkp = b[kk] + b[k2];
                                        bkm = b[kk] - b[k2];
                                        bjp = b[k1] + b[k3];
                                        bjm = b[k1] - b[k3];
                                        b[kk] = bkp + bjp;
                                        bjp = bkp - bjp;
                                        if (isn < 0)
                                        {// go to 450
                                                akp = akm + bjm;
                                                akm = akm - bjm;
                                                bkp = bkm - ajm;
                                                bkm = bkm + ajm;
                                        } // end if (isn < 0 )
                                        else
                                        {
                                                akp = akm - bjm;
                                                akm = akm + bjm;
                                                bkp = bkm + ajm;
                                                bkm = bkm - ajm;
                                        } // end else block
                                        if (s1 == 0)
                                        { // 460
                                                a[k1] = akp;
                                                b[k1] = bkp;
                                                a[k2] = ajp;
                                                b[k2] = bjp;
                                                a[k3] = akm;
                                                b[k3] = bkm;
                                                kk = k3 + kSpan;
                                        } // end if (s1 == 0)
                                        else
                                        {//430 
                                                a[k1] = akp*c1 - bkp*s1;
                                                b[k1] = akp*s1 + bkp*c1;
                                                a[k2] = ajp*c2 - bjp*s2;
                                                b[k2] = ajp*s2 + bjp*c2;
                                                a[k3] = akm*c3 - bkm*s3;
                                                b[k3] = akm*s3 + bkm*c3;
                                                kk = k3 + kSpan;
                                        } // end else block
                                } while (kk <= nt); // go to 420
                                c2 = c1 - (cd*c1 + sd*s1);
                                s1 = (sd*c1 - cd*s1) + s1;
                                c1 = 2.0 - (c2*c2 + s1*s1);
                                s1 = c1 * s1;
                                c1 = c1 * c2;
                                c2 = c1*c1 - s1*s1;
                                s2 = 2.0 * c1 * s1;
                                c3 = c2*c1 - s2*s1;
                                s3 = c2*s1 + s2*c1;
                                kk = kk - nt + jc;
                        } while (kk <= kSpan) ;//go to 420
                        kk = kk - kSpan + inc;
                } while (kk <= jc); //go to 410
                if (kSpan == jc)
                {
                    permute();
                }
	} // end fft4 method
	
	/**
	 * fftOdd() 
	 * Calculates the Fast Fourier Transform
	 * for odd factors  
	 * as part of Singleton's FFT algorithm.
	 * fft3() and fft5() are called from fftOdd
	 * 
	 */
	public void fftOdd(int fac, int factorNumber)
        {
		//System.out.println("Start Method fftOdd()");
		double aa;
		double bb;
                double aj;
                double bj;
                int jf = 0;
                int jj;
		k = fac;
		kSpnn = kSpan;
		kSpan = kSpan / k;
		if (k == 3) 
                {
                    fft3(factorNumber); 
                    //System.out.println("Leave fftOdd after fft3"); 
                    return;
                } // go to 320 // Rob added
		if (k == 5) 
                {
                    fft5(factorNumber); 
                    //System.out.println("Leave fftOdd after fft5"); 
                    return;
                } // go to 510 // Rob added
		if (k != jf) 
                {
			jf = k;
			s1 = rad / ((double)k);
			c1 = Math.cos(s1);
			s1 = Math.sin(s1);
			if (jf > MAX_FACTOR) 
                        {
				// error finish, insufficient array storage
				//isn = 0; // TODO: Why was this set right before throwing an error?
				throw new Error("array bounds exceeded within subroutine fftOdd");
			} // end if (jf > MAX_FACTOR)
			ck[jf] = 1.0;
			sk[jf] = 0.0;
			j = 1;
			do
                        {// 630
				ck[j] = ck[k]*c1 + sk[k]*s1;
				sk[j] = ck[k]*s1 - sk[k]*c1;
				k = k - 1;
				ck[k] = ck[j];
				sk[k] = -sk[j];
				j = j + 1;
			} while (j < k); //end do block
		} // end if (k!=jf)
		do
                {//640
			do
                        {//640
				k1 = kk;
				k2 = kk + kSpnn;
				aa = a[kk];
				bb = b[kk];
				ak = aa;
				bk = bb;
				j = 1;
				k1 = k1 + kSpan;
				do
                                {//650
					k2 = k2 - kSpan;
					j = j + 1;
					at[j] = a[k1] + a[k2];
					ak = at[j] + ak;
					bt[j] = b[k1] + b[k2];
					bk = bt[j] + bk;
					j = j + 1;
					at[j] = a[k1] - a[k2];
					bt[j] = b[k1] - b[k2];
					k1 = k1 + kSpan;
				} while (k1 < k2); //end do block
				a[kk] = ak;
				b[kk] = bk;
				k1 = kk;
				k2 = kk + kSpnn;
				j = 1;
				do
                                {//660
					k1 = k1 + kSpan;
					k2 = k2 - kSpan;
					jj = j;
					ak = aa;
					bk = bb;
					aj = 0.0;
					bj = 0.0;
					k = 1;
					do
                                        {//670
						k = k + 1;
						ak = at[k]*ck[jj] + ak;
						bk = bt[k]*ck[jj] + bk;
						k = k + 1;
						aj = at[k]*sk[jj] + aj;
						bj = bt[k]*sk[jj] + bj;
						jj = jj + j;
						if (jj > jf)
                                                {
                                                    jj = jj - jf;
                                                }
					} while (k < jf); // end do block
					k = jf - j;
					a[k1] = ak - bj;
					b[k1] = bk + aj;
					a[k2] = ak + bj;
					b[k2] = bk - aj;
					j = j + 1;
				} while (j < k); // end do block
				kk = kk + kSpnn;
			} while (kk <= nn); // end do block
			kk = kk - nn;
		} while (kk <= kSpan); //end do block
		rotate(factorNumber); // rotate factor applied to all factors except 2 and 4
		//System.out.println("fftOdd Complete");
	} // end fftOdd method
	
	
	/**
	 * fft3() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 3  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft3(int factorNumber)
        {
            double aj;
            double bj;
		//System.out.println("Start Method fft3()");
		do
                { // 320
			do
                        { // 320 
				k1 = kk + kSpan;
				k2 = k1 + kSpan;
				ak = a[kk];
				bk = b[kk];
				aj = a[k1] + a[k2];
				bj = b[k1] + b[k2];
				a[kk] = ak + aj;
				b[kk] = bk + bj;
				ak = -0.5*aj + ak;
				bk = -0.5*bj + bk;
				aj = (a[k1] - a[k2])*s120;
				bj = (b[k1] - b[k2])*s120;
				a[k1] = ak - bj;
				b[k1] = bk + aj;
				a[k2] = ak + bj;
				b[k2] = bk - aj;
				kk = k2 + kSpan;
			} while (kk < nn); // go to 320
			kk = kk - nn;
		} while (kk <= kSpan);// go to 320
		rotate(factorNumber); //go to 700
		//System.out.println("fft3 Complete");
	}// end fft3 method
	
	
	/**
	 * fft5() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 5  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft5(int factorNumber)
        {
		//System.out.println("Start Method fft5()");
		double c2;
		double s2;
                double aa;
                double bb;
                double aj;
                double bj;
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
		s2 = 2.0*c72*s72;
		do
                { //520
			do
                        { //520
				k1 = kk + kSpan;
				k2 = k1 + kSpan;
				k3 = k2 + kSpan;
				k4 = k3 + kSpan;
				akp = a[k1] + a[k4];
				akm = a[k1] - a[k4];
				bkp = b[k1] + b[k4];
				bkm = b[k1] - b[k4];
				ajp = a[k2] + a[k3];
				ajm = a[k2] - a[k3];
				bjp = b[k2] + b[k3];
				bjm = b[k2] - b[k3];
				aa = a[kk];
				bb = b[kk];
				a[kk] = aa + akp + ajp;
				b[kk] = bb + bkp + bjp;
				ak = akp*c72 + ajp*c2 + aa;
				bk = bkp*c72 + bjp*c2 + bb;
				aj = akm*s72 + ajm*s2;
				bj = bkm*s72 + bjm*s2;
				a[k1] = ak - bj;
				a[k4] = ak + bj;
				b[k1] = bk + aj;
				b[k4] = bk - aj;
				ak = akp*c2 + ajp*c72 + aa;
				bk = bkp*c2 + bjp*c72 + bb;
				aj = akm*s2 - ajm*s72;
				bj = bkm*s2 - bjm*s72;
				a[k2] = ak - bj;
				a[k3] = ak + bj;
				b[k2] = bk + aj;
				b[k3] = bk - aj;
				kk = k4 + kSpan;
			} while (kk < nn);// go to 520
			kk = kk - nn;
		} while (kk <= kSpan);// go to 520
		rotate(factorNumber);
		//System.out.println("fft5 Complete");
	} // end fft5 method

	/**
	 * fft2() 
	 * Calculates the Fast Fourier Transform
	 * for factors of 2 (not 4)  
	 * as part of Singleton's FFT algorithm
	 * 
	 */
	public void fft2(){
		//System.out.println("Start Method fft2()");
		kSpan = kSpan / 2;
		k1 = kSpan + 2;
		do
                { // 210
			do
                        { // 210
				k2 = kk + kSpan;
				ak = a[k2];
				bk = b[k2];
				a[k2] = a[kk] - ak;
				b[k2] = b[kk] - bk;
				a[kk] = a[kk] + ak;
				b[kk] = b[kk] + bk;
				kk = k2 + kSpan;
			} while (kk <= nn);
			kk = kk - nn;
		} while (kk <= jc);
		if (kk > kSpan) 
                {
			permute();  // go to 800 // 800 is permute
			return;
		} // end if (kk > kSpan)
		do
                { // 220
			c1 = 1.0 - cd;
			s1 = sd;
			do
                        { // 230
				do
                                { // 230
					do
                                        { // 230
						k2 = kk + kSpan;
						ak = a[kk] - a[k2];
						bk = b[kk] - b[k2];
						a[kk] = a[kk] + a[k2];
						b[kk] = b[kk] + b[k2];
						a[k2] = c1*ak - s1*bk;
						b[k2] = s1*ak + c1*bk;
						kk= k2 + kSpan;
					} while (kk < nt) ;
					k2 = kk - nt;
					c1 = -c1;
					kk = k1 - k2;
				} while (kk > k2) ;
				ak = c1 - (cd*c1 + sd*s1);
				s1 = (sd*c1 - cd*s1) + s1;
				c1 = 2.0 - (Math.pow(ak, 2) + Math.pow(s1, 2));
				s1 = c1*s1;
				c1 = c1*ak;
				kk = kk + jc;
			} while (kk < k2) ;
			k1 = k1 + inc + inc;
			kk = (k1 - kSpan) / 2 + jc;
		} while (kk <= (jc + jc)) ;
		//System.out.println("fft2() Complete");
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
        // unused assignment is not correct in the do while loop
        @SuppressWarnings("UnusedAssignment")
	public void rotate(int factorNumber) {
		//System.out.println("Start Method rotate");
		double c2;
		double s2;
		if (factorNumber != m)
                {
			kk = jc + 1;
			do
                        {//710 
				c2 = 1.0 - cd;
				s1 = sd;
				do
                                {// 720 
					c1 = c2;
					s2 = s1;
					kk = kk + kSpan;
					do
                                        { // 730
						do
                                                { //730 
							ak = a[kk];
							a[kk] = c2*ak - s2*b[kk];
							b[kk] = s2*ak + c2*b[kk];
							kk = kk + kSpnn;
						} while (kk <= nt);// go to 730
						ak = s1*s2;
						s2 = s1*c2 + c1*s2;
						c2 = c1*c2 - ak;
						kk = kk - nt + kSpan;
					} while (kk <= kSpnn);// go to 730
					c2= c1 - (cd*c1 + sd*s1);
					s1= s1 + (sd*c1 - cd*s1);
					c1 = 2.0 - (c2*c2 + s1*s1);
					s1 = c1*s1;
					c2 = c1*c2;
					kk = kk - kSpnn + jc;
				} while (kk <= kSpan);// go to 720
				kk = kk - kSpan + jc + inc;
			} while (kk <= (jc + jc)); //go to 710
			//System.out.println("End method rotate, without permute");
		} // end if (i != m)
		else
                {
			permute(); // 800
			//System.out.println("End method rotate, with permute");
		} // end else block
	} // end method rotate
	
	/**
	 * Re-arranges the results of the Fast Fourier Transform
	 * as part of Singleton's FFT algorithm
	 */
	private void permute() 
        {
            int maxf = MAX_FACTOR;
            int jj;
		//System.out.println("Start Method permute");
		int ii;
		nP[1] = ks;
		if (kt != 0)
                { 
			k = kt + kt + 1;
			if (m < k)
                        {
                            k = k - 1;
                        }
			j = 1;
			nP[k+1] = jc;
			do
                        { // 810
				nP[j+1] = nP[j] / nFac[j];
				nP[k] = nP[k+1] * nFac[j];
				j = j + 1;
				k = k - 1;
			} while (j < k); //go to 810
			k3 = nP[k+1];
			kSpan = nP[2];
			kk = jc + 1;
			k2 = kSpan + 1;
			j = 1;
			if (n == nTot)
                        {
				boolean b820 = true;
				boolean loop = true;
				start:
				  while (loop)
                                  {
				    if (b820) 
                                    {
				      do
                                      {//820
				        ak = a[kk];
				        a[kk] = a[k2];
				        a[k2] = ak;
				        bk = b[kk];
				        b[kk] = b[k2];
				        b[k2] = bk;
				        kk = kk + inc;
				        k2 = kSpan + k2;
				      } while (k2 < ks);
				    } // end if (b820)
				    do
                                    {//830
				      k2 = k2 - nP[j];
				      j = j + 1;
				      k2 = nP[j+1] + k2;
				    } while (k2 > nP[j]);
				    j = 1;
				    do
                                    {//840
				      if (kk < k2)
                                      {
                                          b820 = true; 
                                          continue start;
                                      }
				      kk = kk + inc;
				      k2 = kSpan + k2;
				    } while (k2 < ks);
				    if (kk < ks) 
                                    {
                                        b820=false; 
                                        continue start;
                                    }
				    jc = k3;
				    break start;
				  }// end while(loop)
			}// end if (n == nTot)	
			else 
                        {
                            throw new Error("n != nTot, a multivariate transform has been called");
                        }
		} // end if (kt != 0)
		
		if (2*kt + 1 >= m) 
                {
                    //System.out.println("End method permute, first position"); 
                    return;
                }
		kSpnn = nP[kt+1];
		j = m - kt;
		nFac[j+1] = 1;
		do
                { //900 
			nFac[j] = nFac[j]*nFac[j+1];
			j = j - 1;
		} while (j != kt); 
		kt = kt + 1;
		nn = nFac[kt] - 1;
		if (nn > MAX_PRODUCT) 
                {
			throw new Error("array bounds exceeded within subroutine permute");
		}
		jj = 0;
		j = 0;
		boolean b902 = false;
		boolean b904 = false;
		boolean loop = true;

		start:
		  while (loop) 
                  {
		    if (b902) 
                    {
		      jj = jj - k2;
		      k2 = kk;
		      k = k + 1;
		      kk = nFac[k];
		    } // end if (b902)
		    if (b904) 
                    {
		      jj = kk + jj;
		      if (jj >= k2) 
                      {
                          b902 = true; 
                          b904 = true; 
                          continue start;
                      }
		      nP[j] = jj;
		    } // end if (b904)
		    k2 = nFac[kt];
		    k = kt + 1;
		    kk = nFac[k];
		    j = j + 1;
		    if (j <= nn) 
                    {
                        b902 = false; 
                        b904 = true; 
                        continue start;
                    }
		    else break start;
		  }// end while(loop)
		j = 0;
		boolean b910 = false;
		loop = true;

		start:
		  while (loop) 
                  {
		    if (b910) 
                    {
		      do
                      {
		        k = kk;
		        kk = nP[k];
		        nP[k] = -kk;
		      } while (kk != j);
		      k3 = kk;
		    } // end if (b910)
		    do
                    {// 914
		      do
                      { // 914
		        j = j + 1;
		        kk = nP[j];
		      } while (kk < 0);
		      if (kk != j) 
                      {
                          b910 = true; 
                          continue start;
                      }
		      nP[j] = -j;
		    } while (j != nn);
		    maxf = inc*maxf;
		    break start;
		  }// end while(loop)
                // reorder a and b, following the permutation cycles
		j = k3 + 1;
		nt = nt - kSpnn;
		ii = nt - inc + 1;
		if (nt >= 0)
                { //go to 924
			do
                        { // 924
				do
                                { // 924	
					do
                                        { // 924	
						j = j - 1;
					} while (nP[j] < 0);// go to 924
					jj = jc;
					do
                                        { // 926
						kSpan = jj;
						if (jj > maxf)
                                                {
                                                    kSpan = maxf;
                                                }
						jj = jj - kSpan;
						k = nP[j];
						kk = jc*k + ii + jj;
						k1 = kk + kSpan;
						k2 = 0;
						do
                                                { // 928
							k2 = k2 + 1;
							at[k2] = a[k1];
							bt[k2] = b[k1];
							k1 = k1 - inc;
						} while (k1 != kk);// go to 928
						do
                                                { // 932
							k1 = kk + kSpan;
							k2 = k1 - jc*(k + nP[k]);
							k = -nP[k];
							do
                                                        { // 936
								a[k1] = a[k2];
								b[k1] = b[k2];
								k1 = k1 - inc;
								k2 = k2 - inc;
							} while (k1 != kk);// go to 936
							kk = k2;
						} while (k != j);// go to 932
						k1 = kk + kSpan;
						k2 = 0;
						do
                                                { // 940
							k2 = k2 + 1;
							a[k1] = at[k2];
							b[k1] = bt[k2];
							k1 = k1 - inc;
						} while (k1 != kk);// go to 940
					} while (jj != 0);// go to 926
				} while (j != 1); // go to 924
				j = k3 + 1;
				nt = nt - kSpnn;
				ii = nt - inc + 1;
			} while (nt >= 0);// go to 924
		} // end if (nt >=0)
		//System.out.println("Method permute Complete, final lines");
	} // end method permute
} // end class MRFFT
