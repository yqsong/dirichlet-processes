package model;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import util.RandomSampler;
import util.VectorOper;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntDoubleHashMap;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.DocCorpora;
import data.structure.DocTimeCorpora;
import data.structure.Keywords;
import description.CorporaClusteringPostInfo;

public class MultiCorporaHDPClusteringSplitMerge {
	DocTimeCorpora docTimeCorpora;
	/**
	 * The top measure's parameter (Dirichlet distribution with symmetric parameters of \b_0)
	 */
	double b0;

	/** [t] is time
	 * Z[t][j][i] is z_{ji}, i.e., the cluster assignment of document i of corpus j
	 */
	List<List<List<Integer>>> Z;
	/**
	 * clusterSize[j].get(k) is number of documents in corpus j assigned to cluster k, i.e., n_{jk} 
	 */
	List<List<OpenIntIntHashMap>> localClusterSize;
	OpenIntIntHashMap globalClusterSize;
	/**
	 * clusterTermFreq.get(k) is an OpenIntIntHashMap mp, and mp.get(w) is the frequency
	 * of term w appearing in cluster k. 
	 */
	OpenIntObjectHashMap clusterTermFreq;
	/**
	 * clusterTotalTermFreq.get(k) = \sum_{w=1}^W clusterTermFreq.get(k).get(w),
	 * preserve this summation to speed up sampling.
	 */
	OpenIntIntHashMap clusterTotalTermFreq;
	/**
	 * the probabilities for measure G_0, not containing \beta_u
	 */
	OpenIntDoubleHashMap beta;
	
	/**
	 * tableNumbers.get(k) is the number of tables (sub-clusters of documents) with 
	 * dish (cluster) k in all restaurants (corpus)
	 */
	OpenIntIntHashMap tableNumbers;
	/**
	 * probability of G_0 on a new component
	 */
	double beta_u;
	
	//data structures for accelerating the computing the marginal of a Dirichlet distribution
	/**
	 * b0SumLogTable storing the summation of log values, \sum_{j=0}^a log(b0 + j)
	 */
	SumLogValueTable b0SumLogTable;
	/**
	 * Storing the summation of log values, \sum_{j=0}^a log(bN + j), bN = b0 * V
	 */
	SumLogValueTable bNSumLogTable;
	/**
	 * index for a cluster
	 */
	int clusterID;
	/**
	 * the index indicating the new cluster
	 */
	static int NEW_CLUSTER_ID = -1;//actually not used but for result report
	/**
	 * concentration parameter for $G_0 \sim DP(\gamma, H)$
	 */
	double gamma;
	
	RandomSampler randomSampler;
	/**
	 * concentration parameter for $G_j \sim DP(\alpha, G_0)$
	 */
	double alpha;
	
	/**
	 * Hyper parameters for vague gamma prior if hyper 
	 * parameters \alpha and \gamma are also required to be sampled
	 * from a vague gamma prior.
	 */
	double alpha_a, alpha_b, gamma_a, gamma_b;
	
	/**
	 * iteration number
	 */
	int maxIter;
	/**
	 * burn-in time
	 */
	int burnIn;
	/**
	 * Whether sampling hyper parameters $\alpha$ and $\gamma$, i.e., the concentration
	 * parameters of two levels' DPs.
	 */
	boolean bSamplingHyper;
	/**
	 * 0: release; 1: iter + time; 2: iter + K + time; 3: iter + K + time + like; 
	 * 4: iter + K + Ks + like + time; 5: iter + K + Ks + like + time + beta + topicsizes
	 */
	int logLevel;
	/**
	 * number of documents selected for convergence diagnosis in each corpus
	 */
	List<List<Integer>> CVG_SEL_NUM;
	/**
	 * documents selected for convergence diagnosis
	 * cvgDocs[j][] is the indices of documents selected in corpus j.  
	 */
	List<List<List<Integer>>> cvgDocs;//the set of documents selected for convergence diagnosis
	
	
	/**
	 * Constructor for the case not need to sample hyper parameters \alpha or \gamma
	 * @param docCorpora
	 * @param b0
	 * @param gamma
	 * @param alpha
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCorporaHDPClusteringSplitMerge(
			DocTimeCorpora docTimeCorpora
			, double b0
			, double gamma, double alpha
			, int maxIter
			, int burnIn
			){
		this(docTimeCorpora
			, b0
			, 0.0, 0.0
			, 0.0, 0.0
			, gamma, alpha
			, false//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
	}
	/**
	 * Constructor for the case when sampling hyper parameters are required.
	 * @param docCorpora
	 * @param b0
	 * @param gamma_a
	 * @param gamma_b
	 * @param alpha_a
	 * @param alpha_b
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCorporaHDPClusteringSplitMerge(
			DocTimeCorpora docTimeCorpora
			, double b0
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, int maxIter
			, int burnIn
			){
		this(docTimeCorpora
			, b0
			, gamma_a, gamma_b
			, alpha_a, alpha_b
			, 0.0, 0.0
			, true//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
	}
	/**
	 * Private constructor.
	 * @param docCorpora, a DocCorpora structure containing the document corpora
	 * @param b0, the parameter of top measure, i.e., the symmetric Dirichlet distribution 
	 * @param gamma_a, shape parameter of the vague gamma distribution for \gamma
	 * @param gamma_b, inverse scale parameter of the vague gamma distribution for \gamma
	 * @param alpha_a, shape parameter of the vague gamma distribution for \alpha
	 * @param alpha_b, inverse scale parameter of the vague gamma distribution for \alpha
	 * @param gamma, if not need to sample hyper-parameters, \gamma is specified.
	 * @param alpha, if not need to sample hyper-parameters, \alpha is specified.
	 * @param bSamplingHyper, whether need to sample hyper parameters. 
	 * @param maxIter, maximal iteration steps
	 * @param burnIn, the burn-in time.
	 */
	protected MultiCorporaHDPClusteringSplitMerge(
			DocTimeCorpora docTimeCorpora
			, double b0
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, double gamma, double alpha
			, boolean bSamplingHyper//whether sampling hyper parameters \alpha and \gamma
			, int maxIter
			, int burnIn
			){
		this.docTimeCorpora = docTimeCorpora;
		this.b0 = b0;
		this.maxIter = maxIter;
		this.burnIn = burnIn;
		this.gamma_a = gamma_a;
		this.gamma_b = gamma_b;
		this.alpha_a = alpha_a;
		this.alpha_b = alpha_b;
		this.bSamplingHyper = bSamplingHyper;
		
		this.randomSampler = new RandomSampler();
		
		//initialize hyperparameter
		this.alpha = (bSamplingHyper)? (this.alpha_a / this.alpha_b) : alpha;
		this.gamma = (bSamplingHyper)? (this.gamma_a / this.gamma_b) : gamma;
		
		//allocate assignments
		int J = this.docTimeCorpora.corporaNumber;
		int T = this.docTimeCorpora.timestampNumber;
		this.Z = new ArrayList<List<List<Integer>>>();
		for (int t = 0; t < T; ++t) {
			this.Z.add(new ArrayList<List<Integer>>());
			for(int j = 0; j < J; j++){
				this.Z.get(t).add(new ArrayList<Integer>());
				for(int i = 0; i < this.docTimeCorpora.docNumbers.get(t).get(j); ++i) {
					this.Z.get(t).get(j).add(0);
				}
			}
		}
		
		//allocate clusterSize
		this.globalClusterSize = new OpenIntIntHashMap();
		this.localClusterSize = new ArrayList<List<OpenIntIntHashMap>>();
		for (int t = 0; t < T; ++t) {
			this.localClusterSize.add(new ArrayList<OpenIntIntHashMap>());
			for(int j = 0; j < J; j++){
				this.localClusterSize.get(t).add(new OpenIntIntHashMap());
			}
			
		}
		//allocate clusterTermFreq
		this.clusterTermFreq = new OpenIntObjectHashMap();
		this.clusterTotalTermFreq = new OpenIntIntHashMap();
		
		
		//allocate beta
		this.beta = new OpenIntDoubleHashMap();
		this.beta_u = 1.0;
		
		//allocate the two log tables
		//the extreme case is that all documents are assigned to a same cluster
		this.b0SumLogTable = new SumLogValueTable(this.b0, this.docTimeCorpora.totalTokenNumber);
		this.bNSumLogTable = new SumLogValueTable(this.b0 * this.docTimeCorpora.dictSize, this.docTimeCorpora.totalTokenNumber);
		
		//allocate tableNumbers
		this.tableNumbers = new OpenIntIntHashMap();
		
		this.logLevel = 0;
		
		int max_sel = Integer.MAX_VALUE;
		this.CVG_SEL_NUM = new ArrayList<List<Integer>>();
		this.cvgDocs = new ArrayList<List<List<Integer>>>();
		for (int t = 0; t < T; ++t) {
			this.CVG_SEL_NUM.add(new ArrayList<Integer>());
			this.cvgDocs.add(new ArrayList<List<Integer>>());
			for(int j = 0; j < J; j++){
				this.CVG_SEL_NUM.get(t).add(
						(this.docTimeCorpora.docNumbers.get(t).get(j) < max_sel)? this.docTimeCorpora.docNumbers.get(t).get(j) : max_sel);
				this.cvgDocs.get(t).add(new ArrayList<Integer>());
				int[] rpm = this.randomSampler.RandPerm(this.docTimeCorpora.docNumbers.get(t).get(j));
				for(int i = 0; i < CVG_SEL_NUM.get(t).get(j); i++){
					this.cvgDocs.get(t).get(j).add(rpm[i]);
				}
			}
		}
	}

	
	/**
	 * @param corpus_j, indicating where this operation comes from
	 */
	private int AddNewCluster(int corpus_j){
		int id = this.clusterID++;
		//sampling beta_u, beta_K
		double nu = this.randomSampler.nextBeta(this.gamma, 1.0);
		double beta_k = this.beta_u * (1 - nu);
		this.beta_u = this.beta_u * nu;
		this.beta.put(id, beta_k);
		this.clusterTermFreq.put(id, new OpenIntIntHashMap());
		return id;
	}
	
	public int getLocalClusterSize(int t, int j, int kid){
		return this.localClusterSize.get(t).get(j).get(kid);
	}
	
	private void RemoveCluster(int clusterID){
		this.beta.removeKey(clusterID);
		this.globalClusterSize.removeKey(clusterID);
		this.clusterTermFreq.removeKey(clusterID);
		this.clusterTotalTermFreq.removeKey(clusterID);
		this.tableNumbers.removeKey(clusterID);
		for (int t = 0; t < this.docTimeCorpora.timestampNumber; ++t) { 
			for(int j = 0; j < this.docTimeCorpora.corporaNumber; j++){
				this.localClusterSize.get(t).get(j).removeKey(clusterID);
			}
		}
	}
	
	public void setLogLevel(int ll){
		this.logLevel = ll;
	}
	private int[] randInitLabel(int len, int[] labelDict){
		int[] ret = new int[len];
		double[] pr = new double[len];
		VectorOper.VectorSet(pr, 1.0 / labelDict.length);
		int[] indices =  RandomSampler.nextMultinomialSeq(pr, len);
		for(int i = 0; i < ret.length; i++){
			ret[i] = labelDict[indices[i]];
		}
		return ret;
	}
	
	private void ClusterAddDoc(int kid, int time_t, int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(time_t).get(corpus_j).get(doc_i);
		int docLen = this.docTimeCorpora.docLengths.get(time_t).get(corpus_j).get(doc_i);
		OpenIntIntHashMap cluster = (OpenIntIntHashMap) this.clusterTermFreq.get(kid);
		VectorOper.VectorAddEqual(cluster, docTf);
		int totalTf = this.clusterTotalTermFreq.get(kid);
		this.clusterTotalTermFreq.put(kid, totalTf + docLen);
		this.Z.get(time_t).get(corpus_j).set(doc_i, kid);
		int njk = this.localClusterSize.get(time_t).get(corpus_j).get(kid);
		this.localClusterSize.get(time_t).get(corpus_j).put(kid, njk + 1);
		int nk = this.globalClusterSize.get(kid);
		this.globalClusterSize.put(kid, nk + 1);
	}
	private void ClusterRemoveDoc(int kid, int time_t, int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(time_t).get(corpus_j).get(doc_i);
		int docLen = this.docTimeCorpora.docLengths.get(time_t).get(corpus_j).get(doc_i);
		OpenIntIntHashMap cluster = (OpenIntIntHashMap) this.clusterTermFreq.get(kid);
		VectorOper.VectorMinusEqual(cluster, docTf);
		int totalTf = this.clusterTotalTermFreq.get(kid);
		this.clusterTotalTermFreq.put(kid, totalTf - docLen);
		this.Z.get(time_t).get(corpus_j).set(doc_i, MultiCorporaHDPClustering.NEW_CLUSTER_ID);//not assigned yet
		int njk = this.localClusterSize.get(time_t).get(corpus_j).get(kid);
		this.localClusterSize.get(time_t).get(corpus_j).put(kid, njk - 1);
		if(njk == 1){
			this.localClusterSize.get(time_t).get(corpus_j).removeKey(kid);
		}
		int nk = this.globalClusterSize.get(kid);
		this.globalClusterSize.put(kid, nk - 1);
		if(nk == 1)
			this.RemoveCluster(kid);
	}
	public int getLocalClusterNumber(int t, int j){
		return this.localClusterSize.get(t).get(j).size();
	}
	public int[] getLocalClusterNumbers(int t){
		int[] ret = new int[this.docTimeCorpora.corporaNumber];
		for(int j = 0; j < ret.length; j++){
			ret[j] = this.localClusterSize.get(t).get(j).size();
		}
		return ret;
	}
	public int getClusterNumber(){
		return this.globalClusterSize.size();
	}
	/**
	 * Calculating the log p(x) of the marginal of Dirichlet distribution, NOTICE that the
	 * constant items having nothing to do with b0 are ignored.
	 * p(x | b) = \int Mult(x | m) Dirichlet(m | b)d m
	 * @param tfVector, the term frequency vector of a sample
	 * @param sumtf is the total count of tfVector
	 * @param tf0, tf0.get(w) + b0 = b_w
	 * @param sumtf0, the total count of tf0
	 * @param b0
	 * @return
	 */
	private double DirichletLogMarginal(OpenIntIntHashMap tfVector, int sumtf, OpenIntIntHashMap tf0, int sumtf0){
		double part1 = 0;
		IntArrayList wids = tfVector.keys();
		IntArrayList vals = tfVector.values();
		if (vals.size() == 0) { // TODO: Potential problems..
			return 0;
		}
		for(int w = 0; w < wids.size(); w++){
			int wid = wids.get(w);
			int tfw0 = tf0.get(wid);
			int tfw = vals.get(w);
			part1 += this.b0SumLogTable.getSumLog(tfw0, tfw + tfw0 - 1);
		}
		double part2 = this.bNSumLogTable.getSumLog(sumtf0, sumtf + sumtf0 - 1);
		return part1- part2;
	}
	private double ClusterLogMarginal(int kid, OpenIntIntHashMap tfVector, int sumtf){
		OpenIntIntHashMap tf0 = (OpenIntIntHashMap) this.clusterTermFreq.get(kid);
		int sumtf0 = this.clusterTotalTermFreq.get(kid);
		return this.DirichletLogMarginal(tfVector, sumtf, tf0, sumtf0);
	}
	private double EmptyDirichletLogMarginal(OpenIntIntHashMap tfVector, int sumtf){
		double part1 = 0;
		IntArrayList vals = tfVector.values();
		if (vals.size() == 0) { // TODO: Potential problems..
			return 0;
		}
		for(int w = 0; w < vals.size(); w++){
			int tfw = vals.get(w);
			part1 += this.b0SumLogTable.getSumLog(0, tfw - 1);
		}
		double part2 = this.bNSumLogTable.getSumLog(0, sumtf  - 1);
		return part1 - part2;
	}
	
	protected int SamplingDocCluster(int current_t, int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(current_t).get(corpus_j).get(doc_i);
		int sumTf = this.docTimeCorpora.docLengths.get(current_t).get(corpus_j).get(doc_i);
		return this.SamplingDocCluster(current_t, corpus_j, docTf, sumTf, true);
	}
	
	protected int PredictDocCluster(int current_t, int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(current_t).get(corpus_j).get(doc_i);
		int sumTf = this.docTimeCorpora.docLengths.get(current_t).get(corpus_j).get(doc_i);
		return this.SamplingDocCluster(current_t, corpus_j, docTf, sumTf, false);
	}
	
	protected int SamplingDocCluster(int current_t, int corpus_j, OpenIntIntHashMap docTf, int sumtf, boolean includeNew){
		int K = this.globalClusterSize.size();
		if(Math.abs(this.alpha) < 10E-8)
			includeNew = false;
		double[] logpost = includeNew? new double[K + 1] : new double[K];
		IntArrayList kids = this.globalClusterSize.keys();
		for(int ki = 0; ki < K; ki++){
			int k = kids.get(ki);
			double beta_k = this.beta.get(k);
			int njk = 0;
			for (int t = 0; t <= current_t; ++t) {
				njk += this.getLocalClusterSize(t, corpus_j, k);
			}
			double logmarg = this.ClusterLogMarginal(k, docTf, sumtf);
			logpost[ki] = Math.log(njk + this.alpha * beta_k) + logmarg;
		}
		if(includeNew){
			logpost[K] = Math.log(this.alpha * beta_u) + this.EmptyDirichletLogMarginal(docTf, sumtf);
		}
		//find a max logpost[i] and return i
		double[] mmax = VectorOper.max(logpost);
		int midx = (int) Math.round(mmax[1]);
		if(midx < K){
			return kids.get(midx);
		}
		else//sampled a new cluster
			return this.AddNewCluster(corpus_j);
	}
	
	public void InitializeT0(int K){
		int[] alphabet = new int[K];
		for(int i = 0; i < K; i++){
			alphabet[i] = i;
		}
		int J = this.docTimeCorpora.corporaNumber;
		int[][] initlabel = new int[J][];
		for (int j = 0; j < J; ++j) {
			initlabel[j] = this.randInitLabel(this.docTimeCorpora.docNumbers.get(0).get(j), alphabet);
		}
		OpenIntIntHashMap labelmap = new OpenIntIntHashMap();
		int ilabel;
		int clusterid = -1;
		for(int j = 0; j < J; j++){
			for(int i = 0; i < initlabel[j].length; i++){
				ilabel = initlabel[j][i];
				if(labelmap.containsKey(ilabel) == false){
					clusterid = AddNewCluster(j);
					labelmap.put(ilabel, clusterid);
				}
				else{
					clusterid = labelmap.get(ilabel);///getComponent(labelmap.get(ilabel));
				}
				
				this.ClusterAddDoc(clusterid, 0, j, i);
			}
		}
	}
	
	public void InitializeT (int time_t) {
		int J = this.docTimeCorpora.corporaNumber;
		int D, z, znew;
		
		for(int j = 0; j < J; j++){
			D = this.docTimeCorpora.docNumbers.get(time_t).get(j);
			for(int i = 0; i < D; i++){
//				z = this.Z.get(time_t).get(j).get(i);

				//sampling a new cluster
				znew = this.PredictDocCluster(time_t, j, i);
				
				//assign the new cluster
				this.ClusterAddDoc(znew, time_t, j, i);
				
			}//for i
		}//for j
	}
	
	public void ResamplingZ (int time_t) {
		int J = this.docTimeCorpora.corporaNumber;
		int D, z, znew;
		
		for(int j = 0; j < J; j++){
			D = this.docTimeCorpora.docNumbers.get(time_t).get(j);
			for(int i = 0; i < D; i++){
				z = this.Z.get(time_t).get(j).get(i);

				//remove the doc from current cluster
				this.ClusterRemoveDoc(z, time_t, j, i);
				
				//sampling a new cluster
				znew = this.PredictDocCluster(time_t, j, i);
				
				//assign the new cluster
				this.ClusterAddDoc(znew, time_t, j, i);
				
			}//for i
		}//for j
	}
	
	public void SamplingZ(int time_t, boolean shuffle){
		int J = this.docTimeCorpora.corporaNumber;
		int[] rpmj = shuffle? this.randomSampler.RandPerm(J) : null;
		int D, z, znew;
		
		for(int rj = 0; rj < J; rj++){
			int j = shuffle? rpmj[rj] : rj;
			D = this.docTimeCorpora.docNumbers.get(time_t).get(j);
			int[] rpmd = shuffle? this.randomSampler.RandPerm(D) : null;
			for(int ri = 0; ri < D; ri++){
				int i = shuffle? rpmd[ri] : ri;
				z = this.Z.get(time_t).get(j).get(i);

				//remove the doc from current cluster
				this.ClusterRemoveDoc(z, time_t, j, i);
				
				//sampling a new cluster
				znew = this.SamplingDocCluster(time_t, j, i);
				
				//assign the new cluster
				this.ClusterAddDoc(znew, time_t, j, i);
				
			}//for i
		}//for j
	}
	
	/**
	 * Sampling the number of sub-clusters for each cluster k in each corpus j
	 */
	protected void SamplingMandBeta(int current_T){
		int J = this.docTimeCorpora.corporaNumber;
		int njk = 0, kid = -1, mk = 0;
		this.tableNumbers.clear();
		
		IntArrayList clusterIDs = this.globalClusterSize.keys();
		int K = clusterIDs.size();
		double[] ms = new double[K + 1];
		ms[K] = this.gamma;
		double bb = 0;
		for(int k = 0; k < K; k++){
			kid = clusterIDs.get(k);
			mk = 0;
			bb = this.alpha * this.beta.get(kid);
			for(int j = 0; j < J; j++){
				njk = 0;
				for (int t = 0; t <= current_T; ++t) {
					njk += this.localClusterSize.get(t).get(j).get(kid);
				}
				//simulate a CRP to sample the table number assigned with 
				//a same component in each restaurant
				if(njk > 0){
					mk += RandomSampler.CRPTableNumber(njk, bb);
				}
			}
			ms[k] = mk;
			if(mk == 0)
				System.err.println("SamplingMandBeta: mk = 0 !");
			this.tableNumbers.put(kid, mk);
		}
//		Smooth(ms);
		double[] beta = this.randomSampler.nextDirichlet(ms);
//		Smooth(beta);
		for(int k = 0; k < K; k++){
			kid = clusterIDs.get(k);
			this.beta.put(kid, beta[k]);
		}
		this.beta_u = beta[K];
	}
	
	/**
	 * Sampling the concentration parameter \gamma, using the auxiliary variable method of 
	 * 		Escobar, M. & West, M. 
	 * 		Bayesian Density Estimation and Inference Using Mixtures. 
	 * 		JASA, 1995, 90, 577-588
	 * @param iter, iteration number
	 */
	protected void SamplingGamma(int iter){
		double ga, eta = 0;
		double m = VectorOper.VectorSum(this.tableNumbers);
		
		int K = this.getClusterNumber();
		double kai = 0, beta = 0;//parameter for gamma distribution
		ga = this.gamma;//_a / m_gamma_b;
		for(int i = 0; i < iter; i++){
			eta = this.randomSampler.nextBeta(ga + 1, m);
			kai = this.gamma_a + K - 1;
			beta = this.gamma_b - Math.log(eta);
			if(Math.random() < (kai / (kai + m * beta)))
				kai += 1;
			ga = this.randomSampler.nextGamma(kai, 1/beta);
			//notice that here 'beta' is the inverse scale parameter,
			//while in Mallet's Randoms class, the nextGamma(alpha, beta) , 
			//'beta' there is the scale parameter
		}
		this.gamma = ga;
	}
	
	/**
	 * Sampling the concentration parameter alpha, using the auxiliary variable method of
	 * 		Teh, Y.; Jordan, M.; Beal, M. & Blei, D. 
	 * 		Hierarchical Dirichlet Processes. 
	 * 		JASA, 2006, 101, 1566-1581
	 * @param iter
	 */
	protected void SamplingAlpha(int iter, int current_T){
		int J = this.docTimeCorpora.corporaNumber;
		
		List<Integer> nj = new ArrayList<Integer>();
		for (int j = 0; j < this.docTimeCorpora.corporaNumber; ++j) {
			int value = 0;
			for (int t = 0; t <= current_T; ++t) {
				value +=  this.docTimeCorpora.docNumbers.get(t).get(j);
			}
			nj.add(value);
		}
		
		double m = VectorOper.VectorSum(this.tableNumbers);
		double w, s;
		
		double t = 0;
		double a = this.alpha_a + m, b = this.alpha_b;//parameters for gamma distribution
		double al = this.alpha;
		for(int i = 0; i < iter; i++){
			a = this.alpha_a + m;
			b = this.alpha_b;
			//sampling wj and sj
			for(int j = 0; j < J; j++){
				w = this.randomSampler.nextBeta(al + 1, nj.get(j));
				t = nj.get(j) / al;
				s = (Math.random() < t / (t + 1))? 1:0;
				a -= s;
				b -= Math.log(w);
			}
			//sampling alpha
			al = this.randomSampler.nextGamma(a, 1 / b);
		}
		this.alpha = al;
	}
	
	public void Sampling(int initK, String outputFolder){
		
		int conpiter = 20;//default iteration for sampling hyperparameters \alpha and \gamma
		long starttime, time0, usedtime;
		int topN = 30;
		
		InitializeT0 (initK);
		Log(0, -1, 0, false);
		
		for(int iter = 0; iter < this.maxIter; iter++){
			starttime = System.currentTimeMillis();
			time0 = starttime;
			SamplingZ(0, true);
			SamplingMandBeta(0);
		
			if(this.bSamplingHyper){
				SamplingAlpha(conpiter, 0);
				SamplingGamma(conpiter);
			}
			usedtime = System.currentTimeMillis() - time0;
			if (iter % 10 == 0) {
				Log(0, iter, usedtime, false);	
			}
		}
		SummarySplitMergeResults splitmergeSum1 = new SummarySplitMergeResults(this.docTimeCorpora, this.globalClusterSize);
		String strTimeT = splitmergeSum1.SummarizeTimestampT (0, this.Z.get(0), topN);
		System.out.println(strTimeT);
		
		File outputFile = null;
		FileWriter writer;
		if (outputFolder != null) {
			try {
				outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(0) + "_cluster.txt" );
				writer = new FileWriter(outputFile, true);
				writer.write(strTimeT);
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		for (int t = 1; t < this.docTimeCorpora.timestampNumber; ++t) {

			List<List<List<Integer>>> passedLabels = PredictPassedTimestamps(t);
//			List<List<Integer>> predictLabels = PredictTimestamp(t);
			
			System.out.println("LogInfo on Single Timestamp. \n");
			
			InitializeT (t);
			for(int iter = 0; iter < this.maxIter; iter++){
				starttime = System.currentTimeMillis();
				time0 = starttime;
				SamplingZ(t, true);
				SamplingMandBeta(t);
				
				if(this.bSamplingHyper){
					SamplingAlpha(conpiter, t);
					SamplingGamma(conpiter);
				}
				usedtime = System.currentTimeMillis() - time0;
				if (iter % 10 == 0) {
					Log(t, iter, usedtime, false);	
				}
			}
			
			System.out.println("LogInfo on Multiple Timestamp. \n");
			for(int iter = 0; iter < this.maxIter; iter++){
				starttime = System.currentTimeMillis();
				time0 = starttime;
				for (int tao = 0; tao <= t; ++tao) {
					ResamplingZ(tao);
				}
				SamplingMandBeta(t);
				
				if(this.bSamplingHyper){
					SamplingAlpha(conpiter, t);
					SamplingGamma(conpiter);
				}
				usedtime = System.currentTimeMillis() - time0;
				if (iter % 10 == 0) {
					Log(t, iter, usedtime, false);	
				}
			}
			
			SummarySplitMergeResults splitmergeSum = new SummarySplitMergeResults(this.docTimeCorpora, this.globalClusterSize);

			strTimeT = splitmergeSum.SummarizeTimestampT (t, this.Z.get(t), topN);
			System.out.println(strTimeT);
			if (outputFolder != null) {
				try {
					outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(t) + "_cluster.txt" );
					writer = new FileWriter(outputFile, true);
					writer.write(strTimeT);
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			String strTimePassedOutput = splitmergeSum.SummarizeTimestampOutput (t, passedLabels, this.Z.subList(0, t + 1), topN);
//			System.out.println(strTimePassedOutput);
			if (outputFolder != null) {
				try {
					outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(t - 1) + "_output.txt" );
					writer = new FileWriter(outputFile, true);
					writer.write(strTimePassedOutput);
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			String strTimeNewInput = splitmergeSum.SummarizeTimestampInput (t, passedLabels.get(t), this.Z.get(t), topN);
//			System.out.println(strTimeNewInput);
			if (outputFolder != null) {
				try {
					outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(t) + "_input.txt" );
					writer = new FileWriter(outputFile, true);
					writer.write(strTimeNewInput);
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			
			System.out.println("\n" + "Training end of Timestamp " + t + ". \n");
		}
		this.docTimeCorpora.docLabels = this.Z;
	}
	
//	List<List<Integer>> PredictTimestamp(int time_t) {
//		List<List<Integer>> predictLabels = new ArrayList<List<Integer>>();
//		
//		List<List<OpenIntIntHashMap>> data_time_t = this.docTimeCorpora.tfVectors.get(time_t);
//		int znew;
//		
//		for(int j = 0; j < data_time_t.size(); j++){
//			predictLabels.add(new ArrayList<Integer>());
//			for(int i = 0; i < data_time_t.get(j).size(); i++){
//				//sampling a new cluster
//				znew = this.PredictDocCluster(time_t, j, i);
//				predictLabels.get(j).add(znew);
//				
//			}//for i
//		}//for j
//		
//		return predictLabels;
//	}
	
	List<List<List<Integer>>> PredictPassedTimestamps(int time_t) {
		List<List<List<Integer>>> predictPassedLabels = new ArrayList<List<List<Integer>>>();
		
		List<List<List<OpenIntIntHashMap>>> data_passed_time_t = 
			this.docTimeCorpora.tfVectors.subList(0, time_t + 1);
		int znew;
		for (int t = 0; t < data_passed_time_t.size(); ++t) {
			predictPassedLabels.add(new ArrayList<List<Integer>>());
			for(int j = 0; j < data_passed_time_t.get(t).size(); j++){
				predictPassedLabels.get(t).add(new ArrayList<Integer>());
				for(int i = 0; i < data_passed_time_t.get(t).get(j).size(); i++){
					//sampling a new cluster
					znew = this.PredictDocCluster(t, j, i);
					predictPassedLabels.get(t).get(j).add(znew);
				}//for i
			}//for j
		}
		
		return predictPassedLabels;
	}
	

	
	
	protected void Log(int timestamp, int iter, long usedTime, boolean isPast){
		int K = -1;
		if(this.logLevel > 1)
			K = this.getClusterNumber();
		
		double loglike = 0;
		if(this.logLevel >= 3){
			if (isPast == true) {
				loglike = CalLogLikePastT(this.cvgDocs, timestamp);
			} else {
				loglike = CalLogLikeT(this.cvgDocs, timestamp);
			}
		}
		
		String log = "[LOG]>>";
		switch(this.logLevel){
		case 0:
			return;//release
		case 1://iter + time
			log += "TimeStamp: " + this.docTimeCorpora.corporaIndex2Time.get(timestamp) + ", Iteration=" + iter + "; Time="+ usedTime + "ms;";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 2:// iter + K + time
			log += "TimeStamp: " + this.docTimeCorpora.corporaIndex2Time.get(timestamp) + ", Iteration=" + iter + "; Time="+ usedTime + "ms;" 
			+ " K=" + K +";";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 3://iter + K + like + time
			log += "TimeStamp: " + this.docTimeCorpora.corporaIndex2Time.get(timestamp) + ", Iteration=" + iter + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		default:
			break;
		}
		if(this.logLevel > 0)
			System.out.println(log);
	}
	
	/**
	 * not strictly the likelihood, but just cumulation of the likelihood of each word
	 * from its topic,  for convergence diagnosis
	 * @return
	 */
	private double CalLogLikeT(List<List<List<Integer>>> docs, int current_t){
		double ret = 0;
		int J = docs.get(current_t).size();
		int kid = 0;
		OpenIntIntHashMap tfDoc = null;
		for(int j = 0; j < J; j++){
			for(int i = 0; i < docs.get(current_t).get(j).size(); i++){
				int d = docs.get(current_t).get(j).get(i);
				tfDoc = this.docTimeCorpora.tfVectors.get(current_t).get(j).get(d);
				int tfsum = this.docTimeCorpora.docLengths.get(current_t).get(j).get(d);
				kid = this.Z.get(current_t).get(j).get(d);
				ret += this.ClusterLogMarginal(kid, tfDoc, tfsum);
			}
		}
		return ret;
	}
	
	private double CalLogLikePastT(List<List<List<Integer>>> docs, int current_t){
		double ret = 0;
		int kid = 0;
		OpenIntIntHashMap tfDoc = null;
		for (int t = 0; t <= current_t; ++t) {
			for(int j = 0; j < docs.get(t).size(); j++){
				for(int i = 0; i < docs.get(t).get(j).size(); i++){
					int d = docs.get(t).get(j).get(i);
					tfDoc = this.docTimeCorpora.tfVectors.get(t).get(j).get(d);
					int tfsum = this.docTimeCorpora.docLengths.get(t).get(j).get(d);
					kid = this.Z.get(t).get(j).get(d);
					ret += this.ClusterLogMarginal(kid, tfDoc, tfsum);
				}
			}
		}
		return ret;
	}
	
	
	public DocTimeCorpora getDocTimeCorpora() {
		return docTimeCorpora;
	}
	public void setDocTimeCorpora(DocTimeCorpora docTimeCorpora) {
		this.docTimeCorpora = docTimeCorpora;
	}
	public List<List<OpenIntIntHashMap>> getLocalClusterSize() {
		return localClusterSize;
	}
	public void setLocalClusterSize(List<List<OpenIntIntHashMap>> localClusterSize) {
		this.localClusterSize = localClusterSize;
	}
	public OpenIntIntHashMap getGlobalClusterSize() {
		return globalClusterSize;
	}
	public void setGlobalClusterSize(OpenIntIntHashMap globalClusterSize) {
		this.globalClusterSize = globalClusterSize;
	}
	public OpenIntObjectHashMap getClusterTermFreq() {
		return clusterTermFreq;
	}
	public void setClusterTermFreq(OpenIntObjectHashMap clusterTermFreq) {
		this.clusterTermFreq = clusterTermFreq;
	}
	public OpenIntIntHashMap getClusterTotalTermFreq() {
		return clusterTotalTermFreq;
	}
	public void setClusterTotalTermFreq(OpenIntIntHashMap clusterTotalTermFreq) {
		this.clusterTotalTermFreq = clusterTotalTermFreq;
	}
	
	public static void main(String[] args){
		
//		double b0 = 0.5;
//		double gamma = 5, alpha = 1;
////		double gamma_a = 2, gamma_b = 1;
////		double alpha_a = 2, alpha_b = 1;
//		
//		int maxIter = 200, burnin = 10;
//		int initK = 5;
//		int loglevel = 3;
//		MultiCorporaHDPClusteringSplitMerge hdpCluster = new MultiCorporaHDPClusteringSplitMerge(
//				docTimeCorpora
//				, b0
////				, gamma_a, gamma_b
////				, alpha_a, alpha_b
//				, gamma, alpha //if sampling hyper, comment this line and uncomment above 2 lines, vice versa
//				, maxIter
//				, burnin);
//		hdpCluster.setLogLevel(loglevel);
//		hdpCluster.Sampling(initK, outputFolder);
		
	}
}
