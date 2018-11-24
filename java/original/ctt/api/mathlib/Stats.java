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

import java.util.Arrays;

public class Stats {
	public static double median(double[] list, int min_idx, int max_idx) {
		/* This method returns the median of the values of an array between min_idx and max_idx
		 * 	array x with nx values
		 * 
		 */	
		Arrays.sort(list,min_idx,max_idx);
		double median=0;
		int idx;
		idx = max_idx/2 + min_idx;
		if(max_idx % 2 ==0){
			
			median = (list[idx] + list[idx+1])*0.5;
		}
		else if (max_idx % 2 ==1){
			median=list[idx+1];
		}
		else{
			// Error action
		}
		return median;
	}

	public static double[] mode(double[] list, int min_idx, int max_idx, int numBins) {
		/* This method returns the mode of the values of an array between min_idx and max_idx
		 * 	mode[0] = location of mode
		 *  mode[1] = count
		 * 
		 */	
		Arrays.sort(list,min_idx,max_idx);
		double minVal = list[min_idx];
		double maxVal = list[max_idx-1];
		double binWidth = (maxVal - minVal)/numBins;
		double binArray[] = new double[numBins];
		double mode[]= new double[2];
		for(int idx=0; idx< numBins; idx++){
			binArray[idx]=0.0;
		}
		

		int idx=0;
		for(int binIndex=0; binIndex<numBins; binIndex++){
			while(((minVal+binIndex*binWidth) <= list[idx]) && (list[idx]< (minVal +(binIndex+1)*binWidth))&& idx<max_idx){
				binArray[binIndex]=binArray[binIndex]+1;
				idx++;
			}
		}
		mode = Util.smaxData(numBins, binArray, 1);
		
		mode[0]=minVal+mode[0]*binWidth;
		return mode;
	}
}
