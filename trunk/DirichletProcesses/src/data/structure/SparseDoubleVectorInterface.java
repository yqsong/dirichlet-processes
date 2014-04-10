package data.structure;
public interface SparseDoubleVectorInterface{
	
	
	/**
	 * @author Jianwen Zhang
	 */
	public SparseDoubleVectorInterface Clone();
	
	/**
	 * 
	 * @param nzIdx, the idx of the indices value, rather than the idx of the value
	 * @return
	 */
	public double getQuickViaNzIdx(int nzIdx);//note! it's different with that method in other type of SparseVectors
	
	public void setQuickViaNzIdx(int nzIdx, double s);
	
	//v = v + this; v will be overwritten
	public void PlusToDense(double[] v);
	
	//v = v - this
	public void MinusFromDense(double[] v);
	
	
	//vi *= s
	public void TimesEqualScalar(double s);
	//vo = vi * s
	public SparseDoubleVectorInterface TimesScalar(double s);
	
	public double InnerWithDense(final double[] v);
	//s = norm(vi)
	public double Norm2();
	
	//
	public double Sum();
	
	
	//res = Min(v), res[0] = minv, res[1] = min idx
	public double[] Min();
	//res = Max(v), res[0] = maxv, res[1] = max idx
	public double[] Max();
	
	public String ToString();
	
	public double[] ToDouble();
	
	public SparseDoubleVectorInterface FromDouble(double[] v);
	
	public int[] getIndices();
	public double[] getValues();
	
	public int getDim();
}