package data.structure;
public interface SparseIntVectorInterface{
	
	/**
	 * @author Jianwen Zhang
	 */	
	public SparseIntVectorInterface Clone();
	
	/**
	 * 
	 * @param nzIdx, the idx of the indices value, rather than the idx of the value
	 * @return
	 */
	public int getQuickViaNzIdx(int nzIdx);//note! it's different with that method in other type of SparseVectors
	
	public void setQuickViaNzIdx(int nzIdx, int s);
	
	//v = v + this; v will be overwritten
	public void PlusToDense(int[] v);
	
	//v = v - this
	public void MinusFromDense(int[] v);
	
	
	//vi *= s
	public void TimesEqualScalar(int s);
	//vo = vi * s
	public SparseIntVectorInterface TimesScalar(int s);
	
	public double InnerWithDense(final int[] v);
//	//s = norm(vi)
//	public long Norm2();
//	
	//
	public int Sum();
	
	
	//res = Min(v), res[0] = minv, res[1] = min idx
	public int[] Min();
	//res = Max(v), res[0] = maxv, res[1] = max idx
	public int[] Max();
	
	public String ToString();
	
	public int[] ToDense();
	
	public SparseIntVectorInterface FromDense(int[] v);
	
	public int[] getIndices();
	public int[] getValues();
	
	public int getDim();
}