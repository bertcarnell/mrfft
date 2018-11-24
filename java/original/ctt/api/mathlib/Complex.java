/*--------------------------------------------------------------------------*
 * Complex
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
 * A complex numbers class 
 * 
 * Revision: 
 *
 * History: 
 * 
 *
\*--------------------------------------------------------------------------*/
package ctt.api.mathlib;
public class Complex {
    private final double re;   // the real part
    private final double im;   // the imaginary part

    // create a new object with the given real and imaginary parts
    public Complex(double real, double imag) {
        this.re = real;
        this.im = imag;
    }
    
     // return a string representation of the invoking object
    public String toString()  { return re + "+" + im + "i"; }

    // return the modulus of this
    public double modulus() { return Math.sqrt(re*re + im*im);  }
	
    // Addition
    // return a new object whose value is (this + b)
    public Complex plus(Complex b) { 
        Complex a = this;             // invoking object
        double real = a.re + b.re;
        double imag = a.im + b.im;
        Complex sum = new Complex(real, imag);
        return sum;
    }
	
    // Additive Inverse
	// return a new object whose value is (-this)
	public Complex minus() { 
		Complex a = this;             // invoking object
		double real = -a.re ;
		double imag = -a.im ;
		Complex negative = new Complex(real, imag);
		return negative;
	}
    
    //Subtraction
    // returns an object whose value is (this - b)
	public Complex minus(Complex b) {
		Complex a = this;             // invoking object
		return a.plus(b.minus());
	}
    
    // multiplication by complex number (this * b)
    public Complex times(Complex b) {
        Complex a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        Complex prod = new Complex(real, imag);
        return prod;
    }
    
    // Multiplication by a real
    // returns an object which the this multiplied by a real (double) (this * scalar)
    public Complex times(double scalar) {
        return new Complex(scalar * re, scalar * im);
    }
    
    // Complex conjugate
    //  returns the complex conjugate of this
    public Complex conjugate() {  return new Complex(re, -im); }
    
    // Real Part
    public double getReal(){
    	return this.re;
    }
    
    // Imaginary Part
    public double getImaginary(){
    	return this.im;
    }
    
    // Exponential of a complex number
    // e^(x+i*y) = e^x * (cos(y)+i*sin(y))
    public Complex exp() {
    	Complex a = this;
    	return new Complex(Math.exp(a.re)*Math.cos(a.im), Math.exp(a.re)*Math.sin(a.im));
    } // end method exp()
    
    // converts a complex number into polar coordinates
    // theta = arctan(im / re)
    public double toTheta(){
    	Complex a = this;
    	return Math.atan2(a.im, a.re);
    } // end method toTheta()
    
} // end class Complex