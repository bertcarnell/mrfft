package com.gmail.bertcarnell.jmrfft;

import com.gmail.bertcarnell.jmrfftlib.MRFFT;
import org.apache.commons.math3.complex.Complex;
import java.util.List;
import java.util.ArrayList;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
    	int num = 2*3*5;
    	
        List<Complex> x = new ArrayList<Complex>(num+1);
        x.add(new Complex(0.0, 0.0));
        for (int i = 0; i < num; i++)
        {
            x.add(new Complex(Math.exp(-i*0.5),Math.exp(-i*0.5)));
        }

        long time1 = System.nanoTime(); 

        MRFFT y = new MRFFT( x, num, num, num, 1);
        y.fft();

        double deltaTime1 = System.nanoTime() - time1;
		
	System.out.println("\n\nThe MRFFT method took " + deltaTime1/1E9 + " s");
    }
}
