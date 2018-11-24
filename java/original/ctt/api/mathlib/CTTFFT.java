/*--------------------------------------------------------------------------*
 * CTTFFT
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
 * Date: 5/11/2004
 *
 * Class to compute CTTFFT (radix 2 only currently)
 * References:  Complex 
 * 
 * Revision: 
 *
 * History: 
 * 
 *
\*--------------------------------------------------------------------------*/
package  ctt.api.mathlib;
import ctt.api.mathlib.Complex;
public class CTTFFT {

    // compute the CTTFFT of x[], assuming its length is a power of 2
    public static Complex[] fft(Complex[] x) {
        int N = x.length;

        Complex[] y = new Complex[N];
        if (N == 1) {
            y[0] = x[0];
            return y;
        }
        if (N % 2 != 0) throw new RuntimeException("N is not a power of 2");

        Complex[] even = new Complex[N/2];
        Complex[] odd  = new Complex[N/2];
        for (int k = 0; k < N/2; k++) even[k] = x[2*k];
        for (int k = 0; k < N/2; k++) odd[k]  = x[2*k + 1];

        Complex[] q = fft(even);
        Complex[] r = fft(odd);

        for (int k = 0; k < N/2; k++) {
            Complex wk = new Complex(Math.cos(-2 * k * Math.PI / N), Math.sin(-2 * k * Math.PI / N));
            y[k]       = q[k].plus(wk.times(r[k]));
            y[k + N/2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }


}
