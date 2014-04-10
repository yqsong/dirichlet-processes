package util;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntDoubleHashMap;

public class RandomSampler extends cc.mallet.util.Randoms{

	/**
	 * 
	 */
	public static double m_EPS_PrK = 10E-10;
	private static final long serialVersionUID = -2708396183850779371L;

	public RandomSampler(){
		super();
	}
	final static int MAX_SAMPLE_SIZE = 1000000;
	static int[] RAND_PERM_INTS = new int[MAX_SAMPLE_SIZE];
//	private static Random m_Rander;
	static{
		for(int i = 0; i < MAX_SAMPLE_SIZE; i++){
			RAND_PERM_INTS[i] = i;
		}
	}
	public synchronized int[] RandPerm(int n){
		int[] vas = new int[n];
		System.arraycopy(RAND_PERM_INTS, 0, vas, 0, n);
		for(int i = n - 1; i >= 0; i--){
			int idx = (int)((i+1) * Math.random());
			if(idx != i){
				int temp = vas[idx];
				vas[idx] = vas[i];
				vas[i] = temp;
			}
		}
		return vas;
	}
	public static IntArrayList[] K_Folder(int n, int K){
		IntArrayList[] ret = new IntArrayList[K];
		for(int k = 0; k < K; k++){
			ret[k] = new IntArrayList();
		}
		for(int i = 0; i < n; i++){
			int k = (int)(Math.random() * K);
			if(k == K)
				k -= 1;
			ret[k].add(i);
		}
		return ret;
	}
	public double[] nextGamma(double[] alpha){
		//sampling from a gamma distribution with shape alpha and scale 1
		double[] res = new double[alpha.length];
		for(int i = 0; i < alpha.length; i++){
			res[i] = nextGamma(alpha[i], 1);
		}
		return res;
	}
	
	public OpenIntDoubleHashMap nextDirichlet(OpenIntDoubleHashMap alpha){
		OpenIntDoubleHashMap ret = new OpenIntDoubleHashMap();
		IntArrayList ks = alpha.keys();
		DoubleArrayList vals = alpha.values();
		double[] da = new double[vals.size()];
		int K = ks.size();
		for(int ki = 0; ki < K; ki++){
			da[ki] = vals.getQuick(ki);
		}
		Smooth(da);
		double[] dv = nextDirichlet(da);
		Smooth(dv);
		for(int ki = 0; ki < K; ki++){
			ret.put(ks.getQuick(ki), dv[ki]);
		}
		return ret;
		
	}
	private void Smooth(double[] v){
		boolean bsmooth = false;
		for(int ki = 0; ki < v.length; ki++){
			if(v[ki] == 0){
				bsmooth = true;
				break;
			}
		}
		if(bsmooth){
			VectorOper.VectorAddEqual(v, m_EPS_PrK);
			VectorOper.VectorTimesEqual(v, 1/VectorOper.VectorSum(v));
		}
	}
	public double[] nextDirichlet(double[] alpha) {
		/**
		 * draw a multinomial from a dirichlet with parameter aa
		 * according to the property that the mean of a set of i.i.d. gamma samples 
		 * conforms to a dirichlet distribution
		 */
//		Smooth(alpha);
		double[] gmspl = nextGamma(alpha);
        double sum = 0;
        for (int i = 0; i < gmspl.length; i++) {
        	sum += gmspl[i];
        }
        for (int i = 0; i < gmspl.length; i++) {
        	gmspl[i] /= sum;
        }
        return gmspl;
    }
	public static int nextMultnomial(double[] pi){
		double sum = 0;
		double rd = Math.random();
		for(int i = 0; i < pi.length; i++){
			sum += pi[i];
			if(rd < sum)
				return i;
		}
		return pi.length - 1;
	}
	public static int[] nextMultnomial(double[] pi, int n){
		/* sampling a multinomial, repeat n times, return a vector with counts in each atom
		 * n: repeat time
		 */
		int[] res = new int[pi.length];
		double cumpp[] = new double[pi.length];
		cumpp[0] = pi[0];
		res[0] = 0;
		for(int i = 1; i < pi.length; i++)
		{
			res[i] = 0;
			cumpp[i] = cumpp[i-1] + pi[i];
		}
		for(int r = 0; r < n; r++){
			double rd = Math.random();
			for(int i = 0; i < cumpp.length; i++){
				if(rd < cumpp[i]){
					res[i] += 1;
					break;
				}
			}
		}
		return res;
		
	}
	public static int[] nextMultinomialSeq(double[] pi, int n){
		int[] res = new int[n];
		
		if (pi == null || pi.length == 0 || n == 0) {
			return res;
		}
		
		double cumpp[] = new double[pi.length];
		cumpp[0] = pi[0];
		
		for(int i = 1; i < pi.length; i++)
		{
			cumpp[i] = cumpp[i-1] + pi[i];
		}
		for(int r = 0; r < n; r++){
			double rd = Math.random();
			for(int i = 0; i < cumpp.length; i++){
				if(rd < cumpp[i]){
					res[r] =  i;
					break;
				}
			}
		}
		return res;
	}
	
	public double[] nextSphericalGaussian(double[] mu, double sigma2){
		double[] res = new double[mu.length];
		for(int i = 0; i < mu.length; i++){
			res[i] = nextGaussian(mu[i], sigma2);
		}
		return res;
	}
	
	public static int CRPTableNumber(int nd, double alpha){
		if(nd < 1){
			return 0;
		}
		int t = 1;
		for(int i = 1; i < nd; i++){
			if(Math.random() < (alpha / (nd + alpha))){
				t++;
			}
		}
		return t;
	}
	
	@Override
	public double nextBeta(double alpha, double beta){
		if((alpha <= 0) || (beta <= 0))
		{
			alpha += 10E-10;
			beta += 10E-10;
		}
		double ret = super.nextBeta(alpha, beta);
		if(ret <= 0)
			ret = 10E-10;
		if(ret >= 1)
			ret = 1 - 10E-10;
		
		if(Double.isNaN(ret)){
			ret = alpha / (alpha + beta);
		}
		return ret;
	}
	public static void main(String[] arg){
		RandomSampler rdm = new RandomSampler();
		int[] rpm = rdm.RandPerm(10);
		System.out.println(VectorOper.ToString(rpm));
	}
}