package util;

import java.util.List;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntDoubleHashMap;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import data.structure.SparseDoubleVectorInterface;

public class VectorOper{
	public static void VectorMinusEqual(OpenIntIntHashMap v1, final OpenIntIntHashMap v2){
		IntArrayList id2 = v2.keys();
		IntArrayList va2 = v2.values();
		int id = 0, v = 0;
		for(int i = 0; i < id2.size(); i++){
			id = id2.get(i);
			v = v1.get(id);
			v1.put(id, v - va2.get(i));
		}
	}
	public static double[] VectorAdd(double[] v1, double[] v2){
		double[] res = new double[v1.length];
		for(int i = 0; i < v1.length; i++){
			res[i] = v1[i] + v2[i];
		}
		return res;
	}
	public static OpenIntDoubleHashMap VectorAdd(OpenIntDoubleHashMap v1, OpenIntIntHashMap v2){
		OpenIntDoubleHashMap ret = (OpenIntDoubleHashMap) v1.copy();
		
		IntArrayList id2 = v2.keys();
		IntArrayList va2 = v2.values();
		int id = 0;
		double v = 0;
		for(int i = 0; i < id2.size(); i++){
			id = id2.get(i);
			v = ret.get(id);
			ret.put(id, v + va2.get(i));
		}
		return ret;
	}
	public static OpenIntIntHashMap VectorAdd(OpenIntIntHashMap v1, OpenIntIntHashMap v2){
		OpenIntIntHashMap ret = (OpenIntIntHashMap) v1.copy();
		
		IntArrayList id2 = v2.keys();
		IntArrayList va2 = v2.values();
		int id = 0, v = 0;
		for(int i = 0; i < id2.size(); i++){
			id = id2.get(i);
			v = ret.get(id);
			ret.put(id, v + va2.get(i));
		}
		return ret;
	}
	public static void VectorAddEqual(OpenIntIntHashMap v1, final OpenIntIntHashMap v2){
		IntArrayList id2 = v2.keys();
		IntArrayList va2 = v2.values();
		int id = 0, v = 0;
		for(int i = 0; i < id2.size(); i++){
			id = id2.get(i);
			v = v1.get(id);
			v1.put(id, v + va2.get(i));
		}
	}
	public static OpenIntDoubleHashMap VectorAdd(OpenIntDoubleHashMap v1
			, OpenIntDoubleHashMap v2){
		OpenIntDoubleHashMap ret = (OpenIntDoubleHashMap) v1.copy();
		
		IntArrayList id2 = v2.keys();
		DoubleArrayList va2 = v2.values();
		int id = 0;
		double v = 0;
		for(int i = 0; i < id2.size(); i++){
			id = id2.get(i);
			v = ret.get(i);
			ret.put(id, v + va2.get(i));
		}
		return ret;
	}
	public static void VectorAddEqual(double[] v, double k){
		for(int i = 0; i < v.length; i++){
			v[i] += k;
		}
	}
	public static void VectorAddEqual(double[] v1, final double[] v2){
		for(int i = 0; i < v1.length; i++){
			v1[i] += v2[i];
		}
	}
	public static void VectorAddEqual(OpenIntDoubleHashMap v1, OpenIntDoubleHashMap v2){
		IntArrayList indices2 = v2.keys();
		DoubleArrayList vals2 = v2.values();
		int len = indices2.size();
		double v = 0;
		int k = 0;
		for(int ki = 0; ki < len; ki++){
			k = indices2.get(ki);
			v = vals2.get(ki) + v1.get(k);
			v1.put(k, v);
		}
	}
	public static double[] VectorAdd(double[] v, double k){
		double[] res = new double[v.length];
		for(int i = 0; i < v.length; i++){
			res[i] = v[i] + k;
		}		
		return res;
	}
	public static void VectorExpEqual(double[] v){
		for(int i = 0; i < v.length; i++){
			v[i] = Math.exp(v[i]);
		}
	}
	public static void VectorTimesEqual(double[] v1, double[] v2){
		for(int i = 0; i < v1.length; i++){
			v1[i] *= v2[i];
		}
	}
	public static double VectorSum(final double[] v){
		double sum = 0;
		for(int i = 0; i < v.length; i++){
			sum += v[i];
		}
		return sum;
	}
	public static double VectorSum(List<Integer> v){
		double sum = 0;
		for(int i = 0; i < v.size(); i++){
			sum += v.get(i);
		}
		return sum;
	}
	public static double VectorSum(List<Double> v, boolean flag){
		double sum = 0;
		for(int i = 0; i < v.size(); i++){
			sum += v.get(i);
		}
		return sum;
	}
	public static int VectorSum(OpenIntIntHashMap v){
		int ret = 0;
		IntArrayList val = v.values();
		for(int i = 0; i < val.size(); i++){
			ret += val.get(i);
		}
		return ret;
	}
	public static double VectorSum(OpenIntDoubleHashMap v){
		double ret = 0;
		DoubleArrayList val = v.values();
		for(int i = 0; i < val.size(); i++){
			ret += val.get(i);
		}
		return ret;
	}
	public static long VectorSum(final int[] v){
		long sum = 0;
		for(int i = 0; i < v.length; i++){
			sum += v[i];
		}
		return sum;
	}
	public static double[] VectorTimes(double[] v, double k){
		double[] res = new double[v.length];
		for(int i = 0; i < v.length; i++){
			res[i] = v[i] * k;
		}
		return res;
	}
	public static OpenIntDoubleHashMap VectorTimes(OpenIntDoubleHashMap v, double k){
		OpenIntDoubleHashMap ret = (OpenIntDoubleHashMap) v.copy();
		IntArrayList ks = ret.keys();
		DoubleArrayList vs = ret.values();
		for(int ki = 0; ki < ks.size(); ki++){
			ret.put(ks.get(ki), vs.get(ki) * k);
		}
		return ret;
	}
	public static void VectorTimesEqual(double[] v, double k){
		for(int i = 0; i < v.length; i++){
			v[i] *= k;
		}
	}
	public static void VectorSqrtEqual(double[] v){
		for(int i = 0; i < v.length; i++){
			v[i] = Math.sqrt(v[i]);
		}
	}
	public static void VectorTimesEqual(OpenIntDoubleHashMap v, double k){
//		OpenIntDoubleHashMap ret = (OpenIntDoubleHashMap) v.copy();
		IntArrayList ks = v.keys();
		DoubleArrayList vs = v.values();
		for(int ki = 0; ki < ks.size(); ki++){
			v.put(ks.get(ki), vs.get(ki) * k);
		}
//		return ret;
	}
	public static double VectorInnerProduct(double[] v1, double[] v2){
		double res = 0;
		for(int i = 0; i < v1.length; i++){
			res += v1[i] * v2[i];
		}
		return res;
	}
	
	public static double VectorInnerProduct(OpenIntDoubleHashMap v1, OpenIntDoubleHashMap v2){
		double ret = 0;
		IntArrayList indices1 = v1.keys();
		DoubleArrayList vals1 = v1.values();
		int len = indices1.size();
		for(int ki = 0; ki < len; ki++){
			ret += vals1.get(ki) * v2.get(indices1.get(ki));
		}
		return ret;
	}
//	public static double[] int2double(int[] vi){
//		double[] vo = new double[vi.length];
//		for(int i = 0; i < vi.length; i++){
//			vo[i] = (double)(vi[i]);
//		}
//		return vo;
//	}
	/**
	 * 
	 * @param v, the vector
	 * @return res[0] is the minmal value, res[1] is the position (index) of the minimum
	 */
	public static double[] min(double[] v){
		double[] res = new double[2];
		int mi = 0;
		for(int i = 1; i < v.length; i++){
			if(v[i] < v[mi])
				mi = i;
		}
		res[0] = v[mi];
		res[1] = mi;
		return res;
	}
	/**
	 * 
	 * @param v, the vector
	 * @return res[0] is the minmal value, res[1] is the position (index) of the minimum
	 */
	public static double[] min(DoubleArrayList v){
		double[] res = new double[2];
		int mi = 0;
		for(int i = 1; i < v.size(); i++){
			if(v.get(i) < v.get(mi))
				mi = i;
		}
		res[0] = v.get(mi);
		res[1] = mi;
		return res;
	}
	
	/**
	 * 
	 * @param v, the vector
	 * @return res[0] is the maximal value, res[1] is the position (index) of the maximum
	 */
	public static double[] max(double[] v){
		double[] res = new double[2];
		int mi = 0;
		for(int i = 1; i < v.length; i++){
			if(v[i] > v[mi])
				mi = i;
		}
		res[0] = v[mi];
		res[1] = mi;
		return res;
	}
	/**
	 * Found the minimum along row or column in a matrix expressed by double[][]
	 * @param mvals, output, the minimum values
	 * @param midxs, output, the index of the minimum in original matrix
	 * @param mat, input, matrix, R x C
	 * @param direction, input, 1: find a minimum along a column, then mvals and midxs are
	 * with length C; 2: find a minimum along a row, then mvals and midxs are with
	 * length R.
	 */
	public static void min(double[] mvals, int[] midxs, final double[][] mat, int direction){
		int mi;
		if(direction == 1){//"min" is taken along a column
			for(int j = 0; j < mat[0].length; j++){//traverse all the columns
				mi = 0;
				for(int i = 0; i < mat.length; i++){
					if(mat[i][j] < mat[mi][j])
					{
						mi = i;
					}
				}
				if(midxs != null){
					midxs[j] = mi;
				}
				if(mvals != null){
					mvals[j] = mat[mi][j];
				}
			}
		}
		else{
			for(int i = 0; i < mat.length; i++){
				mi = 0;
				for(int j = 0; j < mat[0].length; j++){
					if(mat[i][j] < mat[i][mi])
					{
						mi = j;
					}
				}
				if(midxs != null){
					midxs[i] = mi;
				}
				if(mvals != null){
					mvals[i] = mat[i][mi];
				}
			}	
		}
	}
	
	/**
	 * Found the maximum along row or column in a matrix expressed by double[][]
	 * @param mvals, output, the maximum values
	 * @param midxs, output, the index of the maximum in original matrix
	 * @param mat, input, matrix, R x C
	 * @param direction, input, 1: find a maximum along a column, then mvals and midxs are
	 * with length C; 2: find a maximum along a row, then mvals and midxs are with
	 * length R.
	 */
	public static void max(double[] mvals, int[] midxs, final double[][] mat, int direction){
		int mi;
		if(direction == 1){//"min" is taken along a column
			for(int j = 0; j < mat[0].length; j++){//traverse all the columns
				mi = 0;
				for(int i = 0; i < mat.length; i++){
					if(mat[i][j] > mat[mi][j])
					{
						mi = i;
					}
				}
				if(midxs != null){
					midxs[j] = mi;
				}
				if(mvals != null){
					mvals[j] = mat[mi][j];
				}
			}
		}
		else{
			for(int i = 0; i < mat.length; i++){
				mi = 0;
				for(int j = 0; j < mat[0].length; j++){
					if(mat[i][j] > mat[i][mi])
					{
						mi = j;
					}
				}
				if(midxs != null){
					midxs[i] = mi;
				}
				if(mvals != null){
					mvals[i] = mat[i][mi];
				}
			}	
		}
	}
	

	/**
	 * Sum along a dimension of a matrix expressed as an array
	 * @param mat, the matrix
	 * @param direction, 1: along a column, 2: along a row
	 * @return direction=1: res[j] is the summation of the j-th column's elements,
	 * 			direction = 2: res[i] is the summation of the i-th row's elements
	 */
	public static double[] sum(final double[][] mat, int direction){
		double[] res = null;
		int nr = mat.length;
		int nc = mat[0].length;
		if(direction == 1){
			res = new double[nc];
			for(int j = 0; j < nc; j++){
				res[j] = 0;
				for(int i = 0; i < nr; i++){
					res[j] += mat[i][j];
				}
			}
			return res;
		}
		else{
			res = new double[nr];
			for(int i = 0; i < nr; i++){
				res[i] = 0;
				for(int j = 0; j < nc; j++){
					res[i] += mat[i][j];
				}
			}
			return res;
		}
	}
	
	public static double[] sum(final DoubleMatrix2D mat, int direction){
		int nr = mat.rows();
		int nc = mat.columns();
		double[] ret = null;
		if(direction == 1){
			ret = new double[nc];
			DoubleMatrix1D col = null;
			for(int i = 0; i < nc; i++){
				col = mat.viewColumn(i);
				ret[i] = col.zSum();
				//ret.setQuick(i, col.zSum());
			}
			return ret;
		}
		else{
			ret = new double[nr];
			DoubleMatrix1D row = null;
			for(int i = 0; i < nr; i++){
				row = mat.viewRow(i);
				ret[i] = row.zSum();
				//ret.setQuick(i, row.zSum());
			}
			return ret;
		}
	}

	public static double[] sum(final List<SparseDoubleVectorInterface>mat, int direction){
		int nr = mat.size();
		SparseDoubleVectorInterface[] smat = new SparseDoubleVectorInterface[nr];
		for(int i = 0; i < nr; i++){
			smat[i] = mat.get(i);
		}
		return sum(smat, direction);
	}
	public static double[] sum(final SparseDoubleVectorInterface[] mat, int direction){
		int nr = mat.length;
		int nc = mat[0].getDim();
		double[] ret = null;
		if(direction == 1){
			ret = new double[nc];
			for(int i = 0; i < nc; i++){
				ret[i] = 0;
			}
			for(int j = 0; j < nr; j++){
				mat[j].PlusToDense(ret);
			}
			return ret;
		}
		else{
			ret = new double[nr];
			
			for(int i = 0; i < nr; i++){
				ret[i] = mat[i].Sum();
				//ret.setQuick(i, row.zSum());
			}
			return ret;
		}
	}
	public static OpenIntIntHashMap sum(OpenIntIntHashMap[] mat){
		OpenIntIntHashMap ret = (OpenIntIntHashMap) mat[0].copy();
		for(int i = 1; i < mat.length; i++){
			VectorAddEqual(ret, mat[i]);
		}
		return ret;
	}
	public static OpenIntDoubleHashMap Normalize(OpenIntIntHashMap v){
		OpenIntDoubleHashMap ret = new OpenIntDoubleHashMap();
		IntArrayList ks = v.keys();
		IntArrayList vals = v.values();
		double ss = VectorSum(v);
		double val;
		int k = 0;
		for(int ki = 0; ki < ks.size(); ki++){
			k = ks.get(ki);
			val = vals.get(ki);
			ret.put(k, val / ss);
		}
		return ret;
	}
	
	public static void NormalizeEqual(OpenIntDoubleHashMap v, double sum){
		IntArrayList ks = v.keys();
		DoubleArrayList vals = v.values();
		double ss = sum / VectorSum(v);
		double val;
		int k = 0;
		for(int ki = 0; ki < ks.size(); ki++){
			k = ks.get(ki);
			val = vals.get(ki);
			v.put(k, val * ss);
		}
	}
	public static void NormalizeEqual(OpenIntDoubleHashMap[] vs, double sum){
		for(int i = 0; i < vs.length; i++){
			NormalizeEqual(vs[i], sum);
		}
	}
	public static void Normalize(OpenIntDoubleHashMap v){
//		OpenIntDoubleHashMap ret = new OpenIntDoubleHashMap();
		IntArrayList ks = v.keys();
		DoubleArrayList vals = v.values();
		double ss = VectorSum(v);
		if(ss == 0)
			return;
		double val;
		int k = 0;
		for(int ki = 0; ki < ks.size(); ki++){
			k = ks.get(ki);
			val = vals.get(ki);
			v.put(k, val / ss);
		}
	}
	public static double[][] Transpose(double[][] X){
		double[][] XT = new double[X[0].length][X.length];
		for(int i = 0; i < X.length; i++){
			for(int j = 0; j < X[0].length; j++){
				XT[j][i] = X[i][j];
			}
		}
		return XT;
	}
	public static String ToString(double [][] mat){
		String str = "";
		for(int i = 0; i < mat.length; i++){
			for(int j = 0; j < mat[i].length; j++){
				str += mat[i][j] + ",";
			}
			str += "\n";
		}
		return str;
	}
	public static String ToString(double[] v){
		String str = "[";
		for(int i = 0; i < v.length; i++){
			str += v[i] + ",";
		}
		return str + "]";
	}
	public static String ToString(double[] v, String format){
		String str = "[";
		for(int i = 0; i < v.length; i++){
			str += String.format(format, v[i]) + ",";
		}
		return str + "]";
	}
	public static String ToString(double[] v, String pre, String post, String format){
		String str = "[";
		for(int i = 0; i < v.length; i++){
			str += pre + String.format(format, v[i]) + post;
		}
		return str + "]";
	}
	
	public static String ToString(int[] v, String pre, String post){
		String str = "[";
		for(int i = 0; i < v.length; i++){
			str += pre + v[i] + post;
		}
		return str + "]";
	}
	public static String ToString(int[] v, String pre, String post, String format){
		String str = "[";
		for(int i = 0; i < v.length; i++){
			str += pre + String.format(format, v[i]) + post;
		}
		return str + "]";
	}
	public static String ToString(List<Double> v, String format){
		String str = "[";
		for(int i = 0; i < v.size(); i++){
			str += String.format(format, v.get(i)) + ",";
		}
		return str + "]";
	}
	public static String ToString(List<Integer> v, String format, boolean flag){
		//flag no use...
		String str = "[";
		for(int i = 0; i < v.size(); i++){
			str += String.format(format, v.get(i)) + ",";
		}
		return str + "]";
	}
	public static String ToString(List<Double> v, String pre, String post, String format){
		String str = "[";
		for(int i = 0; i < v.size(); i++){
			str += pre + String.format(format, v.get(i)) + post;
		}
		return str + "]";
	}
	public static String ToString(List<Integer> v, String pre, String post, String format, boolean flag){
		//flag no use...
		String str = "[";
		for(int i = 0; i < v.size(); i++){
			str += pre + String.format(format, v.get(i)) + post;
		}
		return str + "]";
	}
	public static String ToString(int [][] mat){
		String str = "";
		for(int i = 0; i < mat.length; i++){
			for(int j = 0; j < mat[i].length; j++){
				str += mat[i][j] + ",";
			}
			str += "\n";
		}
		return str;
	}
	public static String ToString(int[] v){
		String str = "[";
		for(int i = 0; i < v.length; i++){
			str += v[i] + ",";
		}
		return str + "]";
	}
	public static double[] VectorClone(double[] v){
		double[] res = new double[v.length];
		System.arraycopy(v, 0, res, 0, v.length);
		return res;
	}
	public static int[] VectorClone(final int[] v){
		int[] res = new int[v.length];
		System.arraycopy(v, 0, res, 0, v.length);
		return res;
	}
	public static double VectorNorm2(final double[] v){
		double res = 0;
		for(int i = 0; i < v.length; i++){
			res += v[i] * v[i];
		}
		return res;
	}
	public static void VectorSet(double[] v, double s){
		for(int i = 0; i < v.length; i++){
			v[i] = s;
		}
	}
	public static void VectorSet(double[][]v, double s){
		for(int i = 0; i < v.length; i++){
			VectorSet(v[i], s);
		}
	}
	public static void VectorSet(int[] v, int s){
		for(int i = 0; i < v.length; i++){
			v[i] = s;
		}
	}
	public static double[] newVector(int len, double s){
		double[] ret = new double[len];
		VectorSet(ret, s);
		return ret;
	}
	public static double[][] newArray(int row, int col, double s){
		double[][] ret = new double[row][col];
		for(int i = 0; i < row; i++){
			for(int j = 0; j < col; j++){
				ret[i][j] = s;
			}
		}
		return ret;
	}
	public static int[] newVector(int len, int s){
		int[] ret = new int[len];
		VectorSet(ret, s);
		return ret;
	}
	
	
	/**
	 * return int[0][i] the unique value, [1][i] the count of unique value
	 */
	public static int[][] Unique(final int[] v){
		OpenIntIntHashMap imap = new OpenIntIntHashMap();
		int cnt = 0;
		for(int i = 0; i < v.length; i++){
			cnt = imap.get(v[i]);
			imap.put(v[i], cnt+1);
		}
		IntArrayList allv = imap.keys();
		int[][] ret = new int[2][allv.size()];
		int key = -1;
		for(int i = 0; i < allv.size(); i++){
			key = allv.get(i);
			ret[0][i] = key;
			ret[1][i] = imap.get(key);
		}
		return ret;
	}
	
	public static IntArrayList[] Find(final int[] v, final int[] alphabet){
		IntArrayList[] idxs = new IntArrayList[alphabet.length];
		for(int i = 0; i < idxs.length; i++){
			idxs[i] = new IntArrayList();
		}
		for(int i = 0; i < v.length; i++){
			for(int c = 0; c < alphabet.length; c++){
				if(v[i] == alphabet[c]){
					idxs[c].add(i);
				}
			}
		}
		return idxs;
	}
	
	public static double[] Posterior(double[] prior, double[] loglike){
		double minloglike = VectorOper.min(loglike)[0];
		double[] post = VectorOper.VectorAdd(loglike, -minloglike);
		VectorOper.VectorExpEqual(post);
		VectorOper.VectorTimesEqual(post, prior);
		double s = VectorOper.VectorSum(post);
		
		VectorOper.VectorTimesEqual(post, 1/s);
		return post;
	}
	public static double[] DirectPosterior(double[] prior, double[] like){
		int K = prior.length;
		double[] join = new double[K];
		double ss = 0;
		for(int k = 0; k < K; k++){
			join[k] = prior[k] * like[k];
			ss += join[k];
		}
		for(int k = 0; k < K; k++){
			join[k] /= ss;
		}
		return join;
	}
	
	public static double[] MeanStd(double[] v){
		int n = v.length;
		double ret[] = new double[2];
		ret[0] = VectorSum(v)/n;
		double var = 0;
		double dif = 0;
		for(int i = 0; i < n; i++){
			dif = v[i] - ret[0];
			var += dif * dif;
		}
		ret[1] = Math.sqrt(var / (n - 1));
		return ret;
	}
	/**
	 * Find the maximal maxK elements in a vector, suitable for small maxK.
	 * @param v
	 * @param maxK
	 * @return
	 */
	public static int[] MaxK(double[] v, int maxK){
		int L = v.length;
		maxK = (maxK < L)? maxK : L;
		int[] ret = new int[maxK];
		double[] vcopy = VectorOper.VectorClone(v);
		int n = v.length;
		for(int k = 0; k < maxK; k++){
			int mi = 0;
			double mx = vcopy[mi];
			for(int i = 1; i < n; i++){
				double vi = vcopy[i];
				if(vi > mx){
					mi = i;
					mx = vi;
				}
			}
			ret[k] = mi;
			vcopy[mi] = Double.MIN_VALUE;
		}
		return ret;
	}
	public static int[] MaxK(int[] v, int maxK){
		int L = v.length;
		maxK = (maxK < L)? maxK : L;
		int[] ret = new int[maxK];
		int[] vcopy = VectorOper.VectorClone(v);
		int n = v.length;
		for(int k = 0; k < maxK; k++){
			int mi = 0;
			double mx = vcopy[mi];
			for(int i = 1; i < n; i++){
				double vi = vcopy[i];
				if(vi > mx){
					mi = i;
					mx = vi;
				}
			}
			ret[k] = mi;
			vcopy[mi] = Integer.MIN_VALUE;
		}
		return ret;
	}
	static double KL_EPS = 10E-5;
	public static double KL(double[] p, double[] q){
		int n = p.length;
		double ret = 0;
		double pi, qi;
		for(int i = 0; i < n; i++){
			pi = p[i] + KL_EPS;
			qi = q[i] + KL_EPS;
			ret += pi * Math.log(pi / qi);
		}
		return ret;
	}
	public static double SymKL(double[] p, double[] q){
		return 0.5 * (KL(p,q) + KL(q, p));
	}
	public static double KL(OpenIntDoubleHashMap p, OpenIntDoubleHashMap q){
		OpenIntIntHashMap atoms = new OpenIntIntHashMap();
		IntArrayList patoms = p.keys();
		IntArrayList qatoms = q.keys();
		for(int i = 0; i < patoms.size(); i++){
			atoms.put(patoms.get(i), 0);
		}
		for(int i = 0; i < qatoms.size(); i++){
			atoms.put(qatoms.get(i), 0);
		}
		IntArrayList iatoms = atoms.keys();
		int n = iatoms.size();
		double[] dp = new double[n];
		double[] dq = new double[n];
		for(int i = 0; i < n; i++){
			dp[i] = p.get(iatoms.get(i));
			dq[i] = q.get(iatoms.get(i));
		}
		return KL(dp, dq);
	}
	public static double SymKL(OpenIntDoubleHashMap p, OpenIntDoubleHashMap q){
		return 0.5 *  (KL(p,q) + KL(q, p));
	}
		
	public static void main(String[] arg){
		double[][] mat = new double[2][];
		mat[0] = new double[]{1, 2, 3, 4};
		mat[1] = new double[]{4, 3, 2, 1};
		System.out.println(ToString(mat));
		double[] rowmin = new double[2];
		double[] colmin = new double[4];
		int[] rowminIdx = new int[2];
		int[] colminIdx = new int[4];
		min(rowmin, rowminIdx, mat, 2);
		
		System.out.println(ToString(rowmin));
		System.out.println(ToString(rowminIdx));
		
		min(colmin, colminIdx, mat, 1);
		
		System.out.println(ToString(colmin));
		System.out.println(ToString(colminIdx));
	
		int[] maxK = VectorOper.MaxK(mat[0], 3);
		System.out.println(VectorOper.ToString(maxK));
	}

}