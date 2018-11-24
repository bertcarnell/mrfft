/*
 * Created on Oct 19, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package signal_processing;

import ctt.api.mathlib.Complex;
import ctt.api.mathlib.CTTFFT;
import java.util.Arrays;



/**
 * @author MOONEYD
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class Driver {

	public static void main(String[] args) {
    	int num;
    	num = 4;
    	
		Complex[] x = new Complex[num+1];
		for(int i=0; i<num; i++){
			x[i+1]=new Complex(Math.exp(-i*0.5),Math.exp(-i*0.5));
			//System.out.println(Math.exp(-i*0.5));
		}
 
		long time1 = System.nanoTime(); 
		
		MRFFT y = new MRFFT( x, num, num, num, 1);
		y.fft();
		
		double deltaTime1 = System.nanoTime() - time1;
		
		//for(int i=1; i<=num; i++){
		//	System.out.println(x[i].getReal() + "\t" + x[i].getImaginary());
		//}

		// Testing the CTT fft
		//Complex[] x2 = new Complex[num];
		//for(int i=0; i<num; i++){
		//	x2[i]=new Complex(Math.exp(-i*0.5),Math.exp(-i*0.5));
		//}

		//System.out.println(x2.length);
		//Complex [] z;
		//z = CTTFFT.fft(x2);
		//for(int i=0; i<num; i++){
		//	System.out.println(z[i].getReal() + "\t" + z[i].getImaginary());
		//}

		
		// a = T x2
		int N1 = num;

		Complex[] x2 = new Complex[N1];
		for(int j=0; j<N1; j++) {
			x2[j] = new Complex(Math.exp(-j*0.5), Math.exp(-j*0.5));
		}

		time1 = System.nanoTime();
		
		Complex[] a = new Complex[N1];
		a = DFT.dft(x2);
		
		double deltaTime2 = System.nanoTime() - time1;

		//for(int i=0; i<num; i++){
		//	System.out.println(a[i].getReal() + "\t" + a[i].getImaginary());
		//}
		
		System.out.println("\n\nThe MRFFT method took " + deltaTime1/1E9 + " s");
		System.out.println("The DFT code took " + deltaTime2/1E9 + " s");
		
		System.out.println("\nDifference between FFT and DFT abs(DFT-FFT)");
		Complex [] diff = new Complex[N1];
		double realDiff[] = new double[N1];
		double imagDiff[] = new double[N1];
		double lengthDiff[] = new double[N1];
		for(int j=0; j<N1; j++) {
			diff[j] = a[j].minus(x[j+1]);
			realDiff[j] = Math.abs(diff[j].getReal());
			imagDiff[j] = Math.abs(diff[j].getImaginary());
			lengthDiff[j] = a[j].modulus()-x[j+1].modulus();
		}
		Arrays.sort(realDiff);
		Arrays.sort(imagDiff);
		Arrays.sort(lengthDiff);
		
		System.out.println("Maximum absolute difference between real parts " + realDiff[N1-1]);
		System.out.println("Maximum absolute difference between imaginary parts " + imagDiff[N1-1]);
		System.out.println("Maximum absolute difference between moduli " + lengthDiff[N1-1]);

	
	
	}
	
}
