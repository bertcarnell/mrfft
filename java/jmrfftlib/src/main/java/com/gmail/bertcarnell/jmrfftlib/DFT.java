package com.gmail.bertcarnell.jmrfftlib;

import org.apache.commons.math3.complex.Complex;
import java.util.List;
import java.util.ArrayList;

/**
 * Discrete Fourier Transform Class
 * @author carnellr
 */
public class DFT 
{
        /**
         * <pre>
	 * The Discrete Fourier Transform (DFT) of x is found by premultiplying by a matrix T
	 * where t(j,k) = exp(i*2*pi*j*k/N), length(x) = N, j:[0,N], k:[0,N]
	 * a = T x
	 * 
	 * Singleton, R.C. "An algorithm for computing the mixed radix fast Fourier transform."<p>
	 * IEEE Trans. Audio and Electroacoustics, Vol. AU-17, Issue 2, June 1969, pg. 93-103.
         * </pre>
         * @param x Complex signal
         * @return the discrete fourier transform of the signal
         */ 
	public static List<Complex> dft(List<Complex> x)
        {
            int N = x.size();
            double Nd = (double) N;
            List<Complex> a = new ArrayList<Complex>(N);
            double temp1;
            Complex temp2;
            for (int j = 0; j < N; j++) 
            {
                a.add(new Complex(0.0, 0.0));
                for (int k = 0; k < N; k++)
                {
                    temp1 = 2.0 * Math.PI * (double) j * (double) k / Nd;
                    temp2 = a.get(j).add(new Complex(0, temp1).exp().multiply(x.get(k)));
                    a.set(j, new Complex(temp2.getReal(), temp2.getImaginary()));
                } // end for k loop
            } // end for j loop
            return a;
	} // end dft method
} // end DFT class
