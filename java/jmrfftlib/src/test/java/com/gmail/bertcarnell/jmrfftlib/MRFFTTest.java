/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.gmail.bertcarnell.jmrfftlib;

import static com.gmail.bertcarnell.assertextensions.AssertExtensions.*;
import junit.framework.TestCase;
import org.apache.commons.math3.complex.Complex;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import static junit.framework.TestCase.assertEquals;

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
            diff[j] = a.get(j).subtract(y.getSignal(j+1));
            realDiff[j] = Math.abs(diff[j].getReal());
            imagDiff[j] = Math.abs(diff[j].getImaginary());
            lengthDiff[j] = Math.abs(a.get(j).abs()-y.getSignal(j+1).abs());
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
        // maximum prime factor is 23
        int[] sample = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 20, 24, 30, 32, 36, 40, 48, 50, 54, 60, 72, 80, 90, 96, 100};
        for (int w = 0; w < sample.length; w++)
        {
            System.out.println("Testing: " + sample[w]);
            basicFFTTest(sample[w]);
        }
        
        // Test an impulse function at 1
        int num = 8;
        List<Complex> testImpulse = new ArrayList<Complex>(num + 1);
        List<Complex> testImpulse2 = new ArrayList<Complex>(num);
        testImpulse.add(Complex.ZERO);
        testImpulse.add(Complex.ONE);
        testImpulse2.add(Complex.ONE);
        for (int i = 2; i < 9; i++)
        {
            testImpulse.add(Complex.ZERO);
            testImpulse2.add(Complex.ZERO);
        }
        MRFFT y = new MRFFT(testImpulse, 8, 8, 8, 1);
        y.fft();
        List<Complex> a = DFT.dft(testImpulse2);
        for (int i = 1; i < 9; i++)
        {
            assertEqualsLRE(1.0, y.getSignal(i).getReal(), 12);
            assertEqualsLRE(1.0, a.get(i-1).getReal(), 12);
            assertEqualsLRE(0.0, y.getSignal(i).getImaginary(), 12);
            assertEqualsLRE(0.0, a.get(i-1).getImaginary(), 12);
        }

        num = 8;
        List<Complex> testImpulse3 = new ArrayList<Complex>(num + 1);
        List<Complex> testImpulse4 = new ArrayList<Complex>(num);
        testImpulse3.add(Complex.ZERO);
        testImpulse3.add(Complex.ZERO);
        testImpulse4.add(Complex.ZERO);
        testImpulse3.add(Complex.ONE);
        testImpulse4.add(Complex.ONE);
        for (int i = 3; i < 9; i++)
        {
            testImpulse3.add(Complex.ZERO);
            testImpulse4.add(Complex.ZERO);
        }
        y = new MRFFT(testImpulse3, 8, 8, 8, 1);
        y.fft();
        a = DFT.dft(testImpulse4);
        for (int i = 1; i < 9; i++)
        {
            assertEqualsLRE(1.0, y.getSignal(i).abs(), 12);
            assertEqualsLRE(1.0, a.get(i-1).abs(), 12);
        }

        // test for linearity 
        num = 8;
        List<Complex> testImpulse13 = new ArrayList<Complex>(num + 1);
        List<Complex> testImpulse24 = new ArrayList<Complex>(num);
        double c1 = 3.0;
        double c2 = 7.0;
        testImpulse13.add(Complex.ZERO);
        for (int i = 0; i < num; i++)
        {
            testImpulse13.add(testImpulse.get(i+1).multiply(c1).add(testImpulse3.get(i+1).multiply(c2)));
            testImpulse24.add(testImpulse2.get(i).multiply(c1).add(testImpulse4.get(i).multiply(c2)));
        }
        MRFFT y1 = new MRFFT(testImpulse, 8, 8, 8, 1);
        MRFFT y3 = new MRFFT(testImpulse3, 8, 8, 8, 1);
        MRFFT y13 = new MRFFT(testImpulse13, 8, 8, 8, 1);
        y1.fft();
        y3.fft();
        y13.fft();
        Complex temp;
        for (int i = 1; i < 9; i++)
        {
            temp = y1.getSignal(i).multiply(c1).add(y3.getSignal(i).multiply(c2));
            assertEqualsLRE(temp.getReal(), y13.getSignal(i).getReal(), 12);
            assertEqualsLRE(temp.getImaginary(), y13.getSignal(i).getImaginary(), 12);
        }
        List<Complex> y2 = DFT.dft(testImpulse2);
        List<Complex> y4 = DFT.dft(testImpulse4);
        List<Complex> y24 = DFT.dft(testImpulse24);
        for (int i = 0; i < 8; i++)
        {
            temp = y2.get(i).multiply(c1).add(y4.get(i).multiply(c2));
            assertEqualsLRE(temp.getReal(), y24.get(i).getReal(), 12);
            assertEqualsLRE(temp.getImaginary(), y24.get(i).getImaginary(), 12);
        }
        for (int i = 0; i < 8; i++)
        {
            assertEqualsLRE(y13.getSignal(i+1).abs(), y24.get(i).abs(), 12);
            assertEqualsLRE(y1.getSignal(i+1).abs(), y2.get(i).abs(), 12);
            assertEqualsLRE(y3.getSignal(i+1).abs(), y4.get(i).abs(), 12);
        }
        
        // Input random data
        num = 64;
        List<Complex> randomSignal = new ArrayList<Complex>(num + 1);
        randomSignal.add(Complex.ZERO);
        for (int i = 0; i < num; i++)
        {
            randomSignal.add(new Complex(Math.random()*10.0, Math.random()*10.0));
        }
        // test that it does not rail
        MRFFT yrand = new MRFFT(randomSignal, num, num, num, 1);
        yrand.fft();
        
        // Input all zeros
        num = 64;
        List<Complex> zeroSignal = new ArrayList<Complex>(num + 1);
        for (int i = 0; i < num + 1; i++)
        {
            zeroSignal.add(Complex.ZERO);
        }
        MRFFT yzero = new MRFFT(zeroSignal, num, num, num, 1);
        yzero.fft();
        for (int i = 0; i < num; i++)
        {
            assertEquals(0.0, yzero.getSignal(i+1).getReal());
            assertEquals(0.0, yzero.getSignal(i+1).getImaginary());
        }
        
        // inputs are all ones
        num = 64;
        List<Complex> oneSignal = new ArrayList<Complex>(num + 1);
        for (int i = 0; i < num + 1; i++)
        {
            oneSignal.add(Complex.ONE);
        }
        MRFFT yone = new MRFFT(oneSignal, num, num, num, 1);
        yone.fft();
        assertEquals((double) num, yone.getSignal(1).getReal(), 1E-22);
        assertEquals(0.0, yone.getSignal(1).getImaginary(), 1E-22);
        for (int i = 1; i < num; i++)
        {
            assertEquals(0.0, yone.getSignal(i+1).getReal(), 1E-22);
            assertEquals(0.0, yone.getSignal(i+1).getImaginary(), 1E-22);
        }
        
        // Inputs alternate between +1 and -1
        num = 64;
        List<Complex> altSignal = new ArrayList<Complex>(num + 1);
        altSignal.add(Complex.ZERO);
        for (int i = 1; i < num / 2 + 1; i++)
        {
            altSignal.add(Complex.ONE);
            altSignal.add(new Complex(-1.0, 0.0));
        }
        MRFFT yalt = new MRFFT(altSignal, num, num, num, 1);
        yalt.fft();
        for (int i = 0; i < num / 2; i++)
        {
            assertEquals(0.0, yalt.getSignal(i+1).getReal(), 1E-22);
            assertEquals(0.0, yalt.getSignal(i+1).getImaginary(), 1E-22);
        }
        assertEquals((double) num, yalt.getSignal(num / 2 + 1).getReal(), 1E-22);
        assertEquals(0.0, yalt.getSignal(num / 2 + 1).getImaginary(), 1E-22);
        for (int i = num / 2 + 1; i < num; i++)
        {
            assertEquals(0.0, yalt.getSignal(i+1).getReal(), 1E-22);
            assertEquals(0.0, yalt.getSignal(i+1).getImaginary(), 1E-22);
        }
        
        // e^(8*j*2*pi*i/N) for i = 0,1,2, ...,N-1. (j = sqrt(-1))
        num = 64;
        List<Complex> testSignal = new ArrayList<Complex>(num + 1);
        testSignal.add(Complex.ZERO);
        for (int i = 0; i < num; i++)
        {
            testSignal.add(new Complex(0.0, (double) (8*2*i) * Math.PI / (double) num).exp());
        }
        MRFFT ytest1 = new MRFFT(testSignal, num, num, num, 1);
        ytest1.fft();
        for (int i = 0; i < num; i++)
        {
            assertEquals(1.0, ytest1.getSignal(i+1).abs(), 1E-22);
        }
        
        /*
A.Single FFT tests - N inputs and N outputs
    5.Input is e^(8*j*2*pi*i/N) for i = 0,1,2, ...,N-1. (j = sqrt(-1))
    6.Input is cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1.
    7.Input is e^((43/7)*j*2*pi*i/N) for i = 0,1,2, ...,N-1. (j sqrt(-1))
    8.Input is cos((43/7)*2*pi*i/N) for i = 0,1,2, ...,N-1.        
    * */
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
    public void testFftOdd() 
    {
        int[] sample = {7, 11, 13, 17, 19};
        for (int w = 0; w < sample.length; w++)
        {
            // maximum prime factor is 23
            System.out.println(sample[w]);
            basicFFTTest(sample[w]);
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
