package signal_processing;

import java.util.Arrays;
import ctt.api.mathlib.Complex;
import junit.framework.TestCase;

public class MRFFTTest extends TestCase {

	/*
	 * Test method for 'signal_processing.MRFFT.MRFFT(Complex[], int, int, int, int)'
	 */
	public void testMRFFT() {

	}

	/*
	 * Test method for 'signal_processing.MRFFT.factors(int[], int)'
	 */
	public void testFactors() {
	   	// create a signal for factorization
		int num;
    	num=6561;
		Complex[] x = new Complex[num+1];
		for(int i=0; i<num; i++){
			// Start writing the array at the 1st position instead of the 0th
			x[i+1]=new Complex(1,1);
		}

		// Test that the method factors() can factor the signal properly
		MRFFT y = new MRFFT( x, num, num, num, 1);
		assertTrue(y.factors(y.nFac, y.n)==8);
     	int[] test = {0,3,3,3,3,3,3,3,3,0,0};
     	assertTrue(Arrays.equals(test, y.nFac));

     	// Test the method factors() independently of the MRFFT constructor
     	// by sending in n
     	int[] testArray = {6561, 34574400, 5625, 99, 270, 30030, 390390};
     	int[][] resultArray = {
     		{0,3,3,3,3,3,3,3,3,0,0,0},
     		{0,4,3,5,7,7,4,7,7,5,3,4},
     		{0,3,5,5,5,5,3,0,0,0,0,0},
     		{0,3,11,3,0,0,0,0,0,0,0,0},
     		{0,3,2,3,5,3,0,0,0,0,0,0},
     		{0,2,3,5,7,11,13,0,0,0,0,0},
     		{0,13,2,3,5,7,11,13,0,0,0,0}
     	};
     	int[] resultArray2 = {8,11,6,3,5,6,7};
     	// loop throught the test numbers to factor
     	for(int i = 0; i<testArray.length; i++) {
     		// max number of factors is 11 but need one more since the zero'th element not used by factors()
         	int[] nFac = new int[11+1];
         	// use the ith element of testArray
         	int n = testArray[i];
         	// Use the factors method in class MRFFT but pass in new values other than what exist in y
          	int numberOfFactors = y.factors(nFac, n);
          	assertTrue(numberOfFactors==resultArray2[i]);
          	assertTrue(Arrays.equals(nFac,resultArray[i]));
          	// Print the result
         	// System.out.println(n + " factors to " + Arrays.toString(nFac));
     	}
	}


	/*
	 * Test method for 'signal_processing.MRFFT.fft()'
	 */
	public void testFft() {
		for(int w=0; w<5; w++){
			// maximum prime factor is 23
			int[] sample = {6, 8, 10, 12, 18, 20, 24, 30, 32, 36, 40, 48, 50, 54, 60, 72, 80, 90, 96, 100};
			int num = sample[w];
			System.out.println(num);
		
			Complex[] x1 = new Complex[num + 1];
			Complex[] x2 = new Complex[num];
			for(int i=0; i<num; i++){
				x1[i+1] = new Complex(Math.exp(-i*0.5), Math.exp(-i*0.5));
				x2[i] = x1[i + 1];
			}
 
			MRFFT y = new MRFFT( x1, num, num, num, 1);
			y.fft();

			Complex[] a = new Complex[num];
			a = DFT.dft(x2);
		
			Complex [] diff = new Complex[num];
			double realDiff[] = new double[num];
			double imagDiff[] = new double[num];
			double lengthDiff[] = new double[num];
			for(int j=0; j<num; j++) {
				diff[j] = a[j].minus(x1[j+1]);
				realDiff[j] = Math.abs(diff[j].getReal());
				imagDiff[j] = Math.abs(diff[j].getImaginary());
				lengthDiff[j] = Math.abs(a[j].modulus()-x1[j+1].modulus());
			}
			Arrays.sort(realDiff);
			Arrays.sort(imagDiff);
			Arrays.sort(lengthDiff);
			
			assertTrue(realDiff[num-1] < 1.0E-14);
			assertTrue(imagDiff[num-1] < 1.0E-14);
			assertTrue(lengthDiff[num-1] < 1.0E-14);
		}

	}

	/*
	 * Test method for 'signal_processing.MRFFT.fft4()'
	 */
	public void testFft4() {
		for(int w=1; w<=4; w++){
			int num = (int) Math.pow(4, w);
			System.out.println(num);
		
			Complex[] x1 = new Complex[num + 1];
			Complex[] x2 = new Complex[num];
			for(int i=0; i<num; i++){
				x1[i+1] = new Complex(Math.exp(-i*0.5), Math.exp(-i*0.5));
				x2[i] = x1[i + 1];
			}
 
			MRFFT y = new MRFFT( x1, num, num, num, 1);
			y.fft();

			Complex[] a = new Complex[num];
			a = DFT.dft(x2);
		
			Complex [] diff = new Complex[num];
			double realDiff[] = new double[num];
			double imagDiff[] = new double[num];
			double lengthDiff[] = new double[num];
			for(int j=0; j<num; j++) {
				diff[j] = a[j].minus(x1[j+1]);
				realDiff[j] = Math.abs(diff[j].getReal());
				imagDiff[j] = Math.abs(diff[j].getImaginary());
				lengthDiff[j] = Math.abs(a[j].modulus()-x1[j+1].modulus());
			}
			Arrays.sort(realDiff);
			Arrays.sort(imagDiff);
			Arrays.sort(lengthDiff);
			
			assertTrue(realDiff[num-1] < 1.0E-14);
			assertTrue(imagDiff[num-1] < 1.0E-14);
			assertTrue(lengthDiff[num-1] < 1.0E-14);
		}
	}

	/*
	 * Test method for 'signal_processing.MRFFT.fftOdd()'
	 */
	public void testFftOdd() {
		for(int w=0; w<5; w++){
			// maximum prime factor is 23
			int[] sample = {7, 11, 13, 17, 19, 23};
			int num = sample[w];
			System.out.println(num);
		
			Complex[] x1 = new Complex[num + 1];
			Complex[] x2 = new Complex[num];
			for(int i=0; i<num; i++){
				x1[i+1] = new Complex(Math.exp(-i*0.5), Math.exp(-i*0.5));
				x2[i] = x1[i + 1];
			}
 
			MRFFT y = new MRFFT( x1, num, num, num, 1);
			y.fft();

			Complex[] a = new Complex[num];
			a = DFT.dft(x2);
		
			Complex [] diff = new Complex[num];
			double realDiff[] = new double[num];
			double imagDiff[] = new double[num];
			double lengthDiff[] = new double[num];
			for(int j=0; j<num; j++) {
				diff[j] = a[j].minus(x1[j+1]);
				realDiff[j] = Math.abs(diff[j].getReal());
				imagDiff[j] = Math.abs(diff[j].getImaginary());
				lengthDiff[j] = Math.abs(a[j].modulus()-x1[j+1].modulus());
			}
			Arrays.sort(realDiff);
			Arrays.sort(imagDiff);
			Arrays.sort(lengthDiff);
			
			assertTrue(realDiff[num-1] < 1.0E-14);
			assertTrue(imagDiff[num-1] < 1.0E-14);
			assertTrue(lengthDiff[num-1] < 1.0E-14);
		}

	}

	/*
	 * Test method for 'signal_processing.MRFFT.fft3()'
	 */
	public void testFft3() {
		for(int w=1; w<=5; w++){
			int num = (int) Math.pow(3, w);
			System.out.println(num);
		
			Complex[] x1 = new Complex[num + 1];
			Complex[] x2 = new Complex[num];
			for(int i=0; i<num; i++){
				x1[i+1] = new Complex(Math.exp(-i*0.5), Math.exp(-i*0.5));
				x2[i] = x1[i + 1];
			}
 
			MRFFT y = new MRFFT( x1, num, num, num, 1);
			y.fft();

			Complex[] a = new Complex[num];
			a = DFT.dft(x2);
		
			Complex [] diff = new Complex[num];
			double realDiff[] = new double[num];
			double imagDiff[] = new double[num];
			double lengthDiff[] = new double[num];
			for(int j=0; j<num; j++) {
				diff[j] = a[j].minus(x1[j+1]);
				realDiff[j] = Math.abs(diff[j].getReal());
				imagDiff[j] = Math.abs(diff[j].getImaginary());
				lengthDiff[j] = Math.abs(a[j].modulus()-x1[j+1].modulus());
			}
			Arrays.sort(realDiff);
			Arrays.sort(imagDiff);
			Arrays.sort(lengthDiff);
			
			assertTrue(realDiff[num-1] < 1.0E-14);
			assertTrue(imagDiff[num-1] < 1.0E-14);
			assertTrue(lengthDiff[num-1] < 1.0E-14);
		}

	}

	/*
	 * Test method for 'signal_processing.MRFFT.fft5()'
	 */
	public void testFft5() {
		for(int w=1; w<=3; w++){
			int num = (int) Math.pow(5, w);
			System.out.println(num);
		
			Complex[] x1 = new Complex[num + 1];
			Complex[] x2 = new Complex[num];
			for(int i=0; i<num; i++){
				x1[i+1] = new Complex(Math.exp(-i*0.5), Math.exp(-i*0.5));
				x2[i] = x1[i + 1];
			}
 
			MRFFT y = new MRFFT( x1, num, num, num, 1);
			y.fft();

			Complex[] a = new Complex[num];
			a = DFT.dft(x2);
		
			Complex [] diff = new Complex[num];
			double realDiff[] = new double[num];
			double imagDiff[] = new double[num];
			double lengthDiff[] = new double[num];
			for(int j=0; j<num; j++) {
				diff[j] = a[j].minus(x1[j+1]);
				realDiff[j] = Math.abs(diff[j].getReal());
				imagDiff[j] = Math.abs(diff[j].getImaginary());
				lengthDiff[j] = Math.abs(a[j].modulus()-x1[j+1].modulus());
			}
			Arrays.sort(realDiff);
			Arrays.sort(imagDiff);
			Arrays.sort(lengthDiff);
			
			assertTrue(realDiff[num-1] < 1.0E-14);
			assertTrue(imagDiff[num-1] < 1.0E-14);
			assertTrue(lengthDiff[num-1] < 1.0E-14);
		}

	}

	/*
	 * Test method for 'signal_processing.MRFFT.fft2()'
	 */
	public void testFft2() {
		for(int w=1; w<=7; w++){
			int num = (int) Math.pow(2, w);
			System.out.println(num);
		
			Complex[] x1 = new Complex[num + 1];
			Complex[] x2 = new Complex[num];
			for(int i=0; i<num; i++){
				x1[i+1] = new Complex(Math.exp(-i*0.5), Math.exp(-i*0.5));
				x2[i] = x1[i + 1];
			}
 
			MRFFT y = new MRFFT( x1, num, num, num, 1);
			y.fft();

			Complex[] a = new Complex[num];
			a = DFT.dft(x2);
		
			Complex [] diff = new Complex[num];
			double realDiff[] = new double[num];
			double imagDiff[] = new double[num];
			double lengthDiff[] = new double[num];
			for(int j=0; j<num; j++) {
				diff[j] = a[j].minus(x1[j+1]);
				realDiff[j] = Math.abs(diff[j].getReal());
				imagDiff[j] = Math.abs(diff[j].getImaginary());
				lengthDiff[j] = Math.abs(a[j].modulus()-x1[j+1].modulus());
			}
			Arrays.sort(realDiff);
			Arrays.sort(imagDiff);
			Arrays.sort(lengthDiff);
			
			assertTrue(realDiff[num-1] < 1.0E-14);
			assertTrue(imagDiff[num-1] < 1.0E-14);
			assertTrue(lengthDiff[num-1] < 1.0E-14);
		}

	}

	/*
	 * Test method for 'signal_processing.MRFFT.rotate()'
	 */
	public void testRotate() {
		int N = 15;
		Complex[] x1 = new Complex[N+1];
		for(int k=0; k<N; k++) {
			x1[k+1]= new Complex(1, 0);
			//System.out.println(x1[k+1].toString());
		}
		MRFFT y1 = new MRFFT( x1, N, N, N, 1);
		y1.i = 1;
		y1.m=2;
		y1.jc=1;
		y1.cd=0.5;
		y1.s1=0.0;
		y1.sd=0.5;
		y1.c1=0.0;
		y1.kSpan=5;
		y1.ak=.5;
		Arrays.fill(y1.a, 1.0);
		y1.a[0] = 0.0;
		y1.a[N+1] = 0.0;
		Arrays.fill(y1.b, 0.0);
		y1.b[0] = 0.0;
		y1.b[N+1] = 0.0;
		//System.out.println(Arrays.toString(y1.a));
		y1.kSpnn=15;
		y1.rotate();
		//System.out.println(Arrays.toString(y1.a));
		//System.out.println(Arrays.toString(y1.b));
		for(int zz=1; zz<=15; zz++){
			Complex yy = new Complex(y1.a[zz], y1.b[zz]);
			double xx = Math.toDegrees(yy.toTheta());
			//System.out.println(xx);
			if(zz<=5){
				assertTrue(xx==0.0);
			}
			else if(zz>5 && zz<=10){
				assertTrue(xx==0.0+45.0*(zz-6));
			}
			else if(zz>10 && zz<=13) {
				assertTrue(xx==0.0 + 90.0*(zz-11));
			}
			else if(zz==14){
				assertTrue(xx==-90.0);
			}
			else if(zz==15){
				assertTrue(xx==0.0);
			}
		}

	}

	/*
	 * Test method for 'signal_processing.MRFFT.permute()'
	 */
	public void testPermute() {
		int N = 15;
		Complex[] x1 = new Complex[N+1];
		for(int k=0; k<N; k++) {
			x1[k+1]= new Complex(1, 0);
			//System.out.println(x1[k+1].toString());
		}
		MRFFT y1 = new MRFFT( x1, N, N, N, 1);
		y1.i = 1;
		y1.m=2;
		y1.jc=1;
		y1.cd=0.5;
		y1.s1=0.0;
		y1.sd=0.5;
		y1.c1=0.0;
		y1.kSpan=5;
		y1.ak=.5;
		for(int f=0; f<N; f++){
			y1.a[f]=f;
		}
		//Arrays.fill(y1.a, 1.0);
		y1.a[N+1] = 0.0;
		Arrays.fill(y1.b, 0.0);
		y1.b[0] = 0.0;
		y1.b[N+1] = 0.0;
		//System.out.println(Arrays.toString(y1.a));
		y1.kSpnn=15;
		y1.permute();
		//System.out.println(Arrays.toString(y1.a));
		//System.out.println(Arrays.toString(y1.b));
		double [] result = {0.0, 1.0, 6.0, 11.0, 2.0, 7.0, 12.0, 3.0, 8.0, 13.0, 4.0, 9.0, 14.0, 5.0, 10.0, 1.0, 0};
		assertTrue(Arrays.equals(y1.a, result));
	}

}
