package ctt.api.mathlib;

import junit.framework.TestCase;

public class ComplexTest extends TestCase {

	/*
	 * Test method for 'ctt.api.mathlib.Complex.Complex(double, double)'
	 */
	public void testComplex() {

	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.toString()'
	 */
	public void testToString() {
		Complex p1 = new Complex(2,3);
		String p2 = p1.toString();
		//System.out.println(p2);
		assertTrue(p2.equalsIgnoreCase("2.0+3.0i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.modulus()'
	 */
	public void testModulus() {
		Complex p1 = new Complex(2,3);
		double p2 = p1.modulus();
		//System.out.println(p2);
		assertTrue(p2 == Math.sqrt(13));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.plus(Complex)'
	 */
	public void testPlus() {
		Complex p1 = new Complex(2,3);
		Complex p2 = new Complex(5,7);
		Complex p3 = p1.plus(p2);
		//System.out.println(p3);
		assertTrue(p3.toString().equalsIgnoreCase("7.0+10.0i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.minus()'
	 */
	public void testMinus() {
		Complex p1 = new Complex(2,3);
		Complex p2 = p1.minus();
		//System.out.println(p2);
		assertTrue(p2.toString().equalsIgnoreCase("-2.0+-3.0i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.minus(Complex)'
	 */
	public void testMinusComplex() {
		Complex p1 = new Complex(2,3);
		Complex p2 = new Complex(5,7);
		Complex p3 = p1.minus(p2);
		//System.out.println(p3);
		assertTrue(p3.toString().equalsIgnoreCase("-3.0+-4.0i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.times(Complex)'
	 */
	public void testTimesComplex() {
		Complex p1 = new Complex(2,3);
		Complex p2 = new Complex(4,5);
		Complex p3 = p1.times(p2);
		//System.out.println(p3);
		assertTrue(p3.toString().equalsIgnoreCase("-7.0+22.0i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.times(double)'
	 */
	public void testTimesDouble() {
		Complex p1 = new Complex(2,3);
		Complex p2 = p1.times(2.5);
		//System.out.println(p2);
		assertTrue(p2.toString().equalsIgnoreCase("5.0+7.5i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.conjugate()'
	 */
	public void testConjugate() {
		Complex p1 = new Complex(2,3);
		Complex p2 = p1.conjugate();
		//System.out.println(p2);
		assertTrue(p2.toString().equalsIgnoreCase("2.0+-3.0i"));
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.getReal()'
	 */
	public void testGetReal() {
		Complex p1 = new Complex(2,3);
		double p2 = p1.getReal();
		//System.out.println(p2);
		assertTrue(p2 == 2.0);
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.getImaginary()'
	 */
	public void testGetImaginary() {
		Complex p1 = new Complex(2,3);
		double p2 = p1.getImaginary();
		//System.out.println(p2);
		assertTrue(p2 == 3.0);
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.exp()'
	 */
	public void testExp() {
		Complex p1 = new Complex(0.0,Math.PI/2);
		Complex p2 = p1.exp();
		//System.out.println(p2);
		assertTrue(Math.abs(p2.getReal() - 0.0)<1.0E-15 && Math.abs(p2.getImaginary() - 1.0)<1E-15);
		p1 = new Complex(2.0, 0.0);
		p2 = p1.exp();
		//System.out.println(p2);
		assertTrue(Math.abs(p2.getReal() - Math.exp(2))<1.0E-15 && Math.abs(p2.getImaginary() - 0.0)<1E-15);
		p1 = new Complex(1.0, 1.0);
		p2 = p1.exp();
		//System.out.println(p2);
		assertTrue(Math.abs(p2.getReal() - Math.exp(1)*Math.cos(1))<1.0E-15 && Math.abs(p2.getImaginary() - Math.exp(1)*Math.sin(1))<1.0E-15);
	}

	/*
	 * Test method for 'ctt.api.mathlib.Complex.toTheta()'
	 */
	public void testToTheta() {
		Complex p1 = new Complex(1.0, 0.0);
		double p2 = p1.toTheta();
		//System.out.println(p2);
		assertTrue(Math.abs(p2-0.0)<1.0E-15);
		p1 = new Complex(0.0, 1.0);
		p2 = p1.toTheta();
		//System.out.println(p2);
		assertTrue(Math.abs(p2-Math.PI/2)<1.0E-15);
		p1 = new Complex(0.0, -1.0);
		p2 = p1.toTheta();
		//System.out.println(p2);
		assertTrue(Math.abs(p2+Math.PI/2)<1.0E-15);
		p1 = new Complex(-1.0, 0.0);
		p2 = p1.toTheta();
		//System.out.println(p2);
		assertTrue(Math.abs(p2-Math.PI)<1.0E-15);
		p1 = new Complex(1.0, 1.0);
		p2 = p1.toTheta();
		//System.out.println(p2);
		assertTrue(Math.abs(p2-Math.PI/4)<1.0E-15);
		p1 = new Complex(0.0, 0.0);
		p2 = p1.toTheta();
		//System.out.println(p2);
		assertTrue(Math.abs(p2-0.0)<1.0E-15);
	}
	
}
