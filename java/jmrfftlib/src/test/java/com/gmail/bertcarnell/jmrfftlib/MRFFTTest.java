/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.gmail.bertcarnell.jmrfftlib;

import junit.framework.TestCase;
import org.apache.commons.math3.complex.Complex;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author carnellr
 */
public class MRFFTTest extends TestCase {
    
    public MRFFTTest(String testName) {
        super(testName);
    }
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void basicFFTTest(int num)
    {
        List<Complex> x1 = new ArrayList<Complex>(num + 1);
        List<Complex> x2 = new ArrayList<Complex>(num);
        x1.add(new Complex(0.0, 0.0));
        Complex temp;
        for (int i = 0; i < num; i++)
        {
            temp = new Complex(Math.exp(-1.0 * (double) i * 0.5), Math.exp(-1.0 * (double) i * 0.5));
            x1.add(temp);
            x2.add(temp);
        }

        MRFFT y = new MRFFT(x1, num, num, num, 1);
        y.fft();

        List<Complex> a = DFT.dft(x2);

        Complex [] diff = new Complex[num];
        double realDiff[] = new double[num];
        double imagDiff[] = new double[num];
        double lengthDiff[] = new double[num];
        for (int j = 0; j < num; j++) 
        {
            diff[j] = a.get(j).subtract(y.x.get(j+1));
            realDiff[j] = Math.abs(diff[j].getReal());
            imagDiff[j] = Math.abs(diff[j].getImaginary());
            lengthDiff[j] = Math.abs(a.get(j).abs()-y.x.get(j+1).abs());
        }
        Arrays.sort(realDiff);
        Arrays.sort(imagDiff);
        Arrays.sort(lengthDiff);

        assertTrue(realDiff[num-1] < 1.0E-14);
        assertTrue(imagDiff[num-1] < 1.0E-14);
        assertTrue(lengthDiff[num-1] < 1.0E-14);
    }
    /**
     * Test of fft method, of class MRFFT.
     */
    public void testFft() {
        for (int w = 0; w < 5; w++)
        {
            // maximum prime factor is 23
            int[] sample = {6, 8, 10, 12, 18, 20, 24, 30, 32, 36, 40, 48, 50, 54, 60, 72, 80, 90, 96, 100};
            int num = sample[w];
            System.out.println(num);
            basicFFTTest(num);
        }

    }

    /**
     * Test of fft4 method, of class MRFFT.
     */
    public void testFft4() {
        for (int w = 1; w <= 4; w++)
        {
            int num = (int) Math.pow(4, w);
            System.out.println(num);
            basicFFTTest(num);
        }
    }

    /**
     * Test of fftOdd method, of class MRFFT.
     */
    public void testFftOdd() {
        for (int w = 0; w < 5; w++)
        {
            // maximum prime factor is 23
            int[] sample = {7, 11, 13, 17, 19, 23};
            int num = sample[w];
            System.out.println(num);
            basicFFTTest(num);
        }
    }

    /**
     * Test of fft3 method, of class MRFFT.
     */
    public void testFft3() {
        for (int w = 1; w <= 5; w++)
        {
            int num = (int) Math.pow(3, w);
            System.out.println(num);
            basicFFTTest(num);
        }
    }

    /**
     * Test of fft5 method, of class MRFFT.
     */
    public void testFft5() {
        for (int w = 1; w <= 3; w++)
        {
            int num = (int) Math.pow(5, w);
            System.out.println(num);
            basicFFTTest(num);
        }
    }

    /**
     * Test of fft2 method, of class MRFFT.
     */
    public void testFft2() {
        for (int w = 1; w <= 7; w++)
        {
            int num = (int) Math.pow(2, w);
            System.out.println(num);
            basicFFTTest(num);
        }
    }

    /**
     * Test of rotate method, of class MRFFT.
     */
    /*public void testRotate() {
        int N = 15;
        Complex[] x1 = new Complex[N+1];
        for (int k = 0; k < N; k++) 
        {
            x1[k+1]= new Complex(1, 0);
        }
        MRFFT y1 = new MRFFT( Arrays.asList(x1), N, N, N, 1);
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

    }*/
    
    public void testFactors() {
        // create a signal for factorization
        int num = 6561;
	List<Complex> x = new ArrayList<Complex>(num+1);
        x.add(new Complex(0.0, 0.0));
        for (int i = 0; i < num; i++){
            // Start writing the array at the 1st position instead of the 0th
            x.add(new Complex(1.0, 1.0));
        }
        // Test that the method factors() can factor the signal properly
        MRFFT y = new MRFFT( x, num, num, num, 1);
        int[] nFac = new int[11];
        assertEquals(8, y.factors(nFac, num));
        
     	int[] test = {0,3,3,3,3,3,3,3,3,0,0};
     	assertTrue(Arrays.equals(test, nFac));

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
     	for (int i = 0; i < testArray.length; i++) 
        {
     		// max number of factors is 11 but need one more since the zero'th element not used by factors()
         	nFac = new int[11+1];
         	// use the ith element of testArray
         	int n = testArray[i];
         	// Use the factors method in class MRFFT but pass in new values other than what exist in y
          	int numberOfFactors = y.factors(nFac, n);
          	assertEquals(numberOfFactors, resultArray2[i]);
                assertTrue(Arrays.equals(nFac, resultArray[i]));
     	}
    }
}
