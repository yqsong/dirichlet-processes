package data.structure;

import util.VectorOper;

/**
 * @author Jianwen Zhang
 */
public class SparseIntVector implements SparseIntVectorInterface{

	int[] m_Indices;
	int [] m_Vals;
	int m_dim;
	public SparseIntVector(int[] indices, int[] vals, int dim, boolean copy){
		if(copy == false){
			m_Indices = indices;
			m_Vals = vals;
		}
		else{
			m_Indices = new int[indices.length];
			System.arraycopy(indices, 0, m_Indices, 0, indices.length);
			m_Vals = new int[vals.length];
			System.arraycopy(vals, 0, m_Vals, 0, vals.length);
		}
		m_dim = dim;
	}
	public SparseIntVector(int[] indices, int[] vals, int dim){
		this(indices, vals, dim, true);
	}
	
	public SparseIntVector(int[] dense, boolean copy){
		m_Indices = new int[dense.length];
		if(copy == true)
			m_Vals = dense;
		else{
			m_Vals = VectorOper.VectorClone(dense);
		}
		for(int i = 0; i < dense.length; i++){
			m_Indices[i] = i;
		}
		m_dim = dense.length;
	}
	public SparseIntVector(int[] dense){
		this(dense, true);
	}
	
	@Override
	public SparseIntVectorInterface Clone() {
		// TODO Auto-generated method stub
		return new SparseIntVector(m_Indices, m_Vals, m_dim, true);
	}


	@Override
	public SparseIntVectorInterface FromDense(int[] v) {
		// TODO Auto-generated method stub
		
		int nz = 0;
		for(int i = 0; i < v.length; i++){
			if(Math.abs(v[i]) > 0){
				nz++;
			}
		}
		int[] ind = new int[nz];
		int[] val = new int[nz];
		int j = 0;
		for(int i = 0; i < v.length; i++){
			if(Math.abs(v[i]) > 0){
				ind[j] = i;
				val[j] = v[i];
				j++;
			}
		}
		return new SparseIntVector(ind, val, v.length, false);
	}



	@Override
	/**
	 * return res[0], the index of the indices, rather than the index of the vals!
	 * res[1], the maximum values
	 */
	public int[] Max() {
		// TODO Auto-generated method stub
		int mi = 0;
		for(int i = 0; i < m_Indices.length; i++){
			if(m_Vals[i] > m_Vals[mi])
				mi = i;
		}
		int[] res = new int[2];
		res[0] = mi;
		res[1] = m_Vals[mi];
		return res;
	}



	@Override
	public int[] Min() {
		// TODO Auto-generated method stub
		int mi = 0;
		for(int i = 0; i < m_Indices.length; i++){
			if(m_Vals[i] < m_Vals[mi])
				mi = i;
		}
		int[] res = new int[2];
		res[0] = mi;
		res[1] = m_Vals[mi];
		return res;
	}



	@Override
	public void MinusFromDense(int[] v) {
		// TODO Auto-generated method stub
		int idx;
		for(int i = 0; i < m_Indices.length; i++){
			idx = m_Indices[i];
			v[idx] -= m_Vals[i];
		}
	}



//	@Override
//	public long Norm2() {
//		// TODO Auto-generated method stub
//		long res = 0;
//		for(int i = 0; i < m_Indices.length; i++){
//			res += m_Vals[i] * m_Vals[i];
//		}
//		return res;
//	}



	@Override
	public void PlusToDense(int[] v) {
		// TODO Auto-generated method stub
		int idx;
		for(int i = 0; i < m_Indices.length; i++){
			idx = m_Indices[i];
			v[idx] += m_Vals[i];
		}
	}


	@Override
	public int Sum() {
		// TODO Auto-generated method stub
		int res = 0;
		for(int i = 0; i < m_Indices.length; i++){
			res += m_Vals[i];
		}
		return res;
	}



	@Override
	public void TimesEqualScalar(int s) {
		// TODO Auto-generated method stub
		for(int i = 0; i < m_Indices.length; i++){
			m_Vals[i] *= s;
		}
	}



	@Override
	public SparseIntVectorInterface TimesScalar(int s) {
		// TODO Auto-generated method stub
		SparseIntVectorInterface res = this.Clone();
		res.TimesEqualScalar(s);
		return res;
	}



	@Override
	public int[] ToDense() {
		// TODO Auto-generated method stub
//		int midx = 0;
//		for(int i = 0; i < m_Indices.length; i++){
//			if(midx < m_Indices[i]){
//				midx = i;
//			}
//		}
//		int n = m_Indices[midx] + 1;
		int[] res = new int[m_dim];
		for(int i = 0; i < m_dim; i++){
			res[i] = 0;
		}
		for(int i = 0; i < m_Indices.length; i++){
			res[m_Indices[i]] = m_Vals[i];
		}
		return res;
	}



	@Override
	public String ToString() {
		// TODO Auto-generated method stub
		return VectorOper.ToString(ToDense());
	}

//
//	private int getMaxDim(){
//		int midx = 0;
//		for(int i = 0; i < m_Indices.length; i++){
//			if(midx < m_Indices[i]){
//				midx = i;
//			}
//		}
//		return m_Indices[midx] + 1;
//	}

	@Override
	public final int[] getIndices() {
		// TODO Auto-generated method stub
		return m_Indices;
	}



	@Override
	public int getQuickViaNzIdx(int nzIdx) {
		// TODO Auto-generated method stub
		return m_Vals[nzIdx];
	}



	@Override
	public final int[] getValues() {
		// TODO Auto-generated method stub
		return m_Vals;
	}



	@Override
	public void setQuickViaNzIdx(int nzIdx, int s) {
		// TODO Auto-generated method stub
		m_Vals[nzIdx] = s;
	}
	

	@Override
	public double InnerWithDense(final int[] v) {
		// TODO Auto-generated method stub
		double res = 0;
		for(int i = 0; i < m_Indices.length; i++){
			res += m_Vals[i] * v[m_Indices[i]];
		}
		return res;
	}
	@Override
	public int getDim() {
		// TODO Auto-generated method stub
		return m_dim;
	}

	public static void main(String[] arg){
		SparseIntVectorInterface vi0 = new SparseIntVector(new int[]{1, 3, 4}, new int[]{1, -5, 2}, 6);
		System.out.println("vi0 = " + vi0.ToString());
		
		SparseIntVectorInterface vt = (SparseIntVector) vi0.FromDense(new int[]{1, 2, 0, 4});
		System.out.println("vt = " + vt.ToString());
		
		int[] maxres = vi0.Max();
		System.out.println("vi max inx = " + maxres[0] + "max val = " + maxres[1]);
		
		int[] minres = vi0.Min();
		System.out.println("vi min inx = " + minres[0] + "min val = " + minres[1]);
		
		SparseIntVectorInterface vi = vi0.Clone();
		vi.MinusFromDense(new int[]{1, 2, 3, 4, 5});
		System.out.println("vi = " + vi.ToString());
		
		
		vi = vi0.Clone();
		vi.PlusToDense(new int[]{1, 2, 3, 4, 5});
		System.out.println("vi = " + vi.ToString());

		System.out.println("Sum(vi0) = " + vi0.Sum());
		
		vi = vi0.Clone();
		vi.TimesEqualScalar(5);
		System.out.println("vi * 5 = " + vi.ToString());

		System.out.println("indices = " + vi.getIndices().toString());
		System.out.println("vals = " + vi.getValues().toString());
			
	}
	
}