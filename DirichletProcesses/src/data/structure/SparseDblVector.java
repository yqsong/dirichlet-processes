package data.structure;

import util.VectorOper;

/**
 * @author Jianwen Zhang
 */
public class SparseDblVector implements SparseDoubleVectorInterface{

	int[] m_Indices;
	double [] m_Vals;
	int m_dim;
	public SparseDblVector(int[] indices, double[] vals, int dim, boolean copy){
		if(copy == false){
			m_Indices = indices;
			m_Vals = vals;
		}
		else{
			m_Indices = new int[indices.length];
			System.arraycopy(indices, 0, m_Indices, 0, indices.length);
			m_Vals = new double[vals.length];
			System.arraycopy(vals, 0, m_Vals, 0, vals.length);
		}
		m_dim = dim;
	}
	public SparseDblVector(int[] indices, double[] vals, int dim){
		this(indices, vals, dim, true);
	}
	
	public SparseDblVector(double[] dense, boolean copy){
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
	public SparseDblVector(double[] dense){
		this(dense, true);
	}
	public SparseDblVector(final cern.colt.matrix.impl.SparseDoubleMatrix1D x){
		m_dim = x.size();
		cern.colt.list.IntArrayList indexList = new cern.colt.list.IntArrayList();
		cern.colt.list.DoubleArrayList valueList = new cern.colt.list.DoubleArrayList();
		x.getNonZeros(indexList, valueList);
		int nzsize = indexList.size();
		m_Indices = new int[nzsize];
		m_Vals = new double[nzsize];
		System.arraycopy(indexList.elements(), 0, m_Indices, 0, nzsize);
		System.arraycopy(valueList.elements(), 0, m_Vals, 0, nzsize);
	}
	@Override
	public SparseDoubleVectorInterface Clone() {
		// TODO Auto-generated method stub
		return new SparseDblVector(m_Indices, m_Vals, m_dim, true);
	}


	@Override
	public SparseDoubleVectorInterface FromDouble(double[] v) {
		// TODO Auto-generated method stub
		final double EPS = 10E-50;
		int nz = 0;
		for(int i = 0; i < v.length; i++){
			if(Math.abs(v[i]) > EPS){
				nz++;
			}
		}
		int[] ind = new int[nz];
		double[] val = new double[nz];
		int j = 0;
		for(int i = 0; i < v.length; i++){
			if(Math.abs(v[i]) > EPS){
				ind[j] = i;
				val[j] = v[i];
				j++;
			}
		}
		return new SparseDblVector(ind, val, v.length, false);
	}



	@Override
	/**
	 * return res[0], the index of the indices, rather than the index of the vals!
	 * res[1], the maximum values
	 */
	public double[] Max() {
		// TODO Auto-generated method stub
		int mi = 0;
		for(int i = 0; i < m_Indices.length; i++){
			if(m_Vals[i] > m_Vals[mi])
				mi = i;
		}
		double[] res = new double[2];
		res[0] = mi;
		res[1] = m_Vals[mi];
		return res;
	}



	@Override
	public double[] Min() {
		// TODO Auto-generated method stub
		int mi = 0;
		for(int i = 0; i < m_Indices.length; i++){
			if(m_Vals[i] < m_Vals[mi])
				mi = i;
		}
		double[] res = new double[2];
		res[0] = mi;
		res[1] = m_Vals[mi];
		return res;
	}



	@Override
	public void MinusFromDense(double[] v) {
		// TODO Auto-generated method stub
		int idx;
		for(int i = 0; i < m_Indices.length; i++){
			idx = m_Indices[i];
			v[idx] -= m_Vals[i];
		}
	}



	@Override
	public double Norm2() {
		// TODO Auto-generated method stub
		double res = 0;
		for(int i = 0; i < m_Indices.length; i++){
			res += m_Vals[i] * m_Vals[i];
		}
		return res;
	}



	@Override
	public void PlusToDense(double[] v) {
		// TODO Auto-generated method stub
		int idx;
		for(int i = 0; i < m_Indices.length; i++){
			idx = m_Indices[i];
			v[idx] += m_Vals[i];
		}
	}


	@Override
	public double Sum() {
		// TODO Auto-generated method stub
		double res = 0;
		for(int i = 0; i < m_Indices.length; i++){
			res += m_Vals[i];
		}
		return res;
	}



	@Override
	public void TimesEqualScalar(double s) {
		// TODO Auto-generated method stub
		for(int i = 0; i < m_Indices.length; i++){
			m_Vals[i] *= s;
		}
	}



	@Override
	public SparseDoubleVectorInterface TimesScalar(double s) {
		// TODO Auto-generated method stub
		SparseDoubleVectorInterface res = Clone();
		res.TimesEqualScalar(s);
		return res;
	}



	@Override
	public double[] ToDouble() {
		// TODO Auto-generated method stub
//		int midx = 0;
//		for(int i = 0; i < m_Indices.length; i++){
//			if(midx < m_Indices[i]){
//				midx = i;
//			}
//		}
//		int n = m_Indices[midx] + 1;
		double[] res = new double[m_dim];
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
		return VectorOper.ToString(ToDouble());
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
	public double getQuickViaNzIdx(int nzIdx) {
		// TODO Auto-generated method stub
		return m_Vals[nzIdx];
	}



	@Override
	public final double[] getValues() {
		// TODO Auto-generated method stub
		return m_Vals;
	}



	@Override
	public void setQuickViaNzIdx(int nzIdx, double s) {
		// TODO Auto-generated method stub
		m_Vals[nzIdx] = s;
	}
	

	@Override
	public double InnerWithDense(final double[] v) {
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
		SparseDoubleVectorInterface vi0 = new SparseDblVector(new int[]{1, 3, 4}, new double[]{0.1, -0.5, 0.2}, 6);
		System.out.println("vi0 = " + vi0.ToString());
		
		SparseDoubleVectorInterface vt = (SparseDblVector) vi0.FromDouble(new double[]{1, 2, 0, 4});
		System.out.println("vt = " + vt.ToString());
		
		double[] maxres = vi0.Max();
		System.out.println("vi max inx = " + maxres[0] + "max val = " + maxres[1]);
		
		double[] minres = vi0.Min();
		System.out.println("vi min inx = " + minres[0] + "min val = " + minres[1]);
		
		SparseDoubleVectorInterface vi = vi0.Clone();
		vi.MinusFromDense(new double[]{1, 2, 3, 4, 5});
		System.out.println("vi = " + vi.ToString());
		
		System.out.println("Norm(vi0) = " + vi0.Norm2());
		
		vi = vi0.Clone();
		vi.PlusToDense(new double[]{1, 2, 3, 4, 5});
		System.out.println("vi = " + vi.ToString());

		System.out.println("Sum(vi0) = " + vi0.Sum());
		
		vi = vi0.Clone();
		vi.TimesEqualScalar(0.5);
		System.out.println("vi * 0.5 = " + vi.ToString());

		System.out.println("indices = " + vi.getIndices().toString());
		System.out.println("vals = " + vi.getValues().toString());
		

		cern.colt.matrix.impl.SparseDoubleMatrix1D coltv = new cern.colt.matrix.impl.SparseDoubleMatrix1D(5);
		coltv.setQuick(1, 2);
		coltv.setQuick(3, 1);
		System.out.println("coltv = " + coltv.toString());
		
		SparseDblVector sv = new SparseDblVector(coltv);
		System.out.println("sv = " + sv.ToString());
		
		
		
		

		
		
		
	}
	
}