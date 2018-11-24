/*
 * Created on Apr 1, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/**
 * @author MOONEYD
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package ctt.api.mathlib;

public class Util {

	public static int hsprmx(double[] a, int n, int[] indx){
		int ier=0;
		double b[] = new double[n];
		if (n<1){
			ier=-1;
		}
		for(int i=0; i<n; i++){
			b[i]=a[indx[i]];
			
		}
		for(int i=0; i<n; i++){
			a[i]=b[i];
		}
		for(int i=0; i<n; i++){
			//System.out.println( "a " + a[i] + " indx " + indx[i]+"\n");
		}	
		return ier;
	}
	public static int ismax(int n, double[][] x, int incx, int slice){
		/*
		 * Returns index of max value in indicated slice of 2-d arrary
		 */
		int currentMaxIndex=0;
		double currentMaxValue;
		
		currentMaxValue=x[0][slice];
		for(int i=0;i<n;i=i+incx){
			if (x[i][slice] > currentMaxValue){
				currentMaxValue=x[i][slice];
				currentMaxIndex=i;
			}
			 
		}
		return currentMaxIndex;
	}
	public static int isamax(int n, double[][] x, int iBegin, int incx, int slice){
		/*
		 * Returns index of max absolute value in indicated slice of 2-d arrary
		 */
		int currentMaxIndex=0;
		double currentMaxValue;
		
		currentMaxValue=Math.abs(x[iBegin][slice]);
		
		for(int i=iBegin;i<n;i=i+incx){
			if (Math.abs(x[i][slice]) > currentMaxValue){
				currentMaxValue=Math.abs(x[i][slice]);
				currentMaxIndex=i;
			}
		}
		return currentMaxIndex;
	}	
	
	public static double smax(int n, double[] x, int incx){
		/*
		 * Returns maximum value of a single dimensional array
		 */
		int currentMaxIndex=0;
		double currentMaxValue;
		
		currentMaxValue=x[0];
		for(int i=0;i<n;i=i+incx){
			if (x[i] > currentMaxValue){
				currentMaxValue=x[i];
				currentMaxIndex=i;
			} 
		}
		return currentMaxValue;
	}
	public static double[] smaxData(int n, double[] x, int incx){
		/*
		 * Returns maximum value and index of value for a single dimensional array 
		 */
		int currentMaxIndex=0;
		double currentMaxValue;
		double maxdata[]= new double[2];
		
		currentMaxValue=x[0];
		for(int i=0;i<n;i=i+incx){
			if (x[i] > currentMaxValue){
				currentMaxValue=x[i];
				currentMaxIndex=i;
			} 
		}
		maxdata[0]=(double)currentMaxIndex;
		maxdata[1]=currentMaxValue;
		return maxdata;
	}	
	public static double[] smaxData(int n, double[] x, boolean[] y, int incx){
		/*
		 * Returns maximum value and index of value for a single dimensional array 
		 * Only considers indices for which for which y is true;
		 * 
		 */
		int currentMaxIndex=0;
		double currentMaxValue=0;
		for(int i=0;i<n;i=i+incx){
			if(y[i]==true){
				currentMaxIndex=i;
				currentMaxValue=x[i];
				break;
			}
			
		}
		double maxdata[]= new double[2];
		
		currentMaxValue=x[0];
		for(int i=0;i<n;i=i+incx){
			if(y[i]==true){
				if (x[i] > currentMaxValue){
					currentMaxValue=x[i];
					currentMaxIndex=i;
				} 
			}
		}
		maxdata[0]=(double)currentMaxIndex;
		maxdata[1]=currentMaxValue;
		return maxdata;
	}	
	
	public static int sscale_slice(int n, double scale, double[][] x, int iBegin, int slice){
		/*
		 * Upon return the indicated slice of x is scaled by "scale" between iMin 
		 * and iMin + lenBin
		 */
		int ierr = 0;
		
		for(int i=iBegin; i< n; i++){
			x[i][slice]=scale*x[i][slice];
		}
		
		return ierr;
		
	}	
	public static int sscale(int n, double scale, double[] x, int iBegin){
		/*
		 * Upon return x is scaled by "scale" between iMin 
		 * and iMin + lenBin
		 */
		int ierr = 0;
		
		for(int i=iBegin; i< n; i++){
			x[i]=scale*x[i];
		}
		
		return ierr;
		
	}

	public static double smin(int n, double[] x, int incx){
		/*
		 * Returns minimum value of a single dimensional array
		 */
		int currentMinIndex=0;
		double currentMinValue;
		
		currentMinValue=x[0];
		for(int i=0;i<n;i=i+incx){
			if (x[i] < currentMinValue){
				currentMinValue=x[i];
				currentMinIndex=i;
			} 
		}
		return currentMinValue;
	}
	
	public static int scopy(int n, double[][] x, int incx, int xStart, double[][] y, int incy, int yStart, int slice){
		/*
		 * copies n values of x into y respecting start positions
		 * this application uses 2-d arrays and this method copies only the indicated slice
		 */
		int currentMaxIndex=0;
		int ierr=0;
		double currentMaxValue;
		
		int xPosition=xStart;
		int yPosition=yStart;
	
		for(int idx=0;idx<n;idx++){
				y[yPosition][slice]=x[xPosition][slice];
				yPosition=yPosition+incy;
				xPosition=xPosition+incx;		
		}
		return ierr;
	}
	
	public static int sfill(int n, double  x, double[][] y, int incy, int yStart, int slice){
		/*
		 * copys n copies of x into y starting at yStart in the indicated slice of y
		 */
		int currentMaxIndex=0;
		int ierr=0;
		double currentMaxValue;
		
		 
		int yPosition=yStart;
	
		for(int idx=0;idx<n;idx++){
				y[yPosition][slice]=x ;
				yPosition=yPosition+incy;
				 
		}
		return ierr;
	}
}
