package signal_processing;

import ctt.api.mathlib.Complex;

public class DFT {

	/**Method dft of class DFT
	 * <p>
	 * The Discrete Fourier Transform (DFT) of x is found by premultiplying by a matrix T <p>
	 * where t(j,k) = exp(i*2*pi*j*k/N), length(x) = N, j:[0,N], k:[0,N]<p>
	 * a = T x<p>
	 * <p>
	 * Singleton, R.C. "An algorithm for computing the mixed radix fast Fourier transform."<p>
	 * IEEE Trans. Audio and Electroacoustics, Vol. AU-17, Issue 2, June 1969, pg. 93-103.
	 */
	public static Complex[] dft(Complex[] x){
		int N = x.length;
		
		Complex[] a = new Complex[N];
		double temp1;
		Complex temp2;
		for(int j=0; j<N; j++) {
			a[j] = new Complex(0,0);
			for(int k=0; k<N; k++) {
				temp1 = 2*Math.PI*j*k/N;
				temp2 = new Complex(0, temp1).exp().times(x[k]);
				a[j] = a[j].plus(temp2);
			} // end for k loop
			//System.out.println(a[j].toString());
		} // end for j loop

		return a;
		
	} // end dft method

} // end DFT class
