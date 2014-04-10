package model;
public class SumLogValueTable{
	/**
	 * sumLogValue[i] = \sum_{j=0}^i log(j + bias); Consequently, \sum_{j=a}^b log(j + bias) = sumLogValue[b] - sumLogValue[a-1];
	 */
	double[] sumLogValues;
	double bias;
	int maxN;//length of logValues is maxN + 1
	public SumLogValueTable(double bias, int maxN){
		this.sumLogValues = new double[maxN + 1];
		this.bias = bias;
		this.maxN = maxN;
		this.sumLogValues[0] = Math.log(bias);
		for(int i = 1; i < maxN + 1; i++){
			this.sumLogValues[i] = this.sumLogValues[i-1] + Math.log(i + bias);
		}
	}
	/**
	 * \sum_{j=a}^b log(bias + j)
	 * @param a
	 * @param b
	 * @return
	 */
	public double getSumLog(int a, int b){
		if(a < 0)
			System.err.println("SumLogValueTable.java: lower a is negative: a < 0 !");
		if(b > this.maxN)
			System.err.println("SumLogValueTable.java: upper b is larger than maximum value allowed: b > maxN !");
		
		if(a == 0)
			return this.sumLogValues[b];
		else
			return this.sumLogValues[b] - this.sumLogValues[a - 1];
			
	}
}