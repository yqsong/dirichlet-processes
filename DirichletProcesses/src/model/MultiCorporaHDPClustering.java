/**
 * HDPMultinomialClustering.java
 * @author Jianwen Zhang (jw-zhang06@mails.tsinghua.edu.cn, jianwenzh@gmail.com)
 * Feb. 21, 2010.
 * 
 * HDP is implemented as an clustering model for multiple document corpora,
 * with multinomial-dirichlet conjugate pair, 
 * using direct posterior sampling scheme.
 * 
 */
package model;


import java.io.File;
import java.util.ArrayList;
import java.util.List;

import util.RandomSampler;
import util.VectorOper;

import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntDoubleHashMap;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.DocCorpora;
import data.structure.DocCorpus;
import data.structure.DocTimeCorpora;
import description.CorporaClusteringPostInfo;

public class MultiCorporaHDPClustering{
	
	DocCorpora docCorpora;
	/**
	 * The top measure's parameter (Dirichlet distribution with symmetric parameters of \b_0)
	 */
	double b0;

	/**
	 * Z[j][i] is z_{ji}, i.e., the cluster assignment of document i of corpus j
	 */
	List<List<Integer>> Z;
	/**
	 * clusterSize[j].get(k) is number of documents in corpus j assigned to cluster k, i.e., n_{jk} 
	 */
	List<OpenIntIntHashMap> localClusterSize;
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
	List<Integer> CVG_SEL_NUM;
	/**
	 * documents selected for convergence diagnosis
	 * cvgDocs[j][] is the indices of documents selected in corpus j.  
	 */
	List<List<Integer>> cvgDocs;//the set of documents selected for convergence diagnosis
	
	
	/**
	 * Constructor for the case not need to sample hyper parameters \alpha or \gamma
	 * @param docCorpora
	 * @param b0
	 * @param gamma
	 * @param alpha
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCorporaHDPClustering(
			DocCorpora docCorpora
			, double b0
			, double gamma, double alpha
			, int maxIter
			, int burnIn
			){
		this(docCorpora
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
	public MultiCorporaHDPClustering(
			DocCorpora docCorpora
			, double b0
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, int maxIter
			, int burnIn
			){
		this(docCorpora
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
	protected MultiCorporaHDPClustering(
			DocCorpora docCorpora
			, double b0
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, double gamma, double alpha
			, boolean bSamplingHyper//whether sampling hyper parameters \alpha and \gamma
			, int maxIter
			, int burnIn
			){
		this.docCorpora = docCorpora;
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
		int J = this.docCorpora.corporaNumber;
		this.Z = new ArrayList<List<Integer>>();
		for(int j = 0; j < J; j++){
			this.Z.add(new ArrayList<Integer>());
			for(int i = 0; i < this.docCorpora.docNumbers.get(j); ++i) {
				this.Z.get(j).add(0);
			}
		}
		
		//allocate clusterSize
		this.globalClusterSize = new OpenIntIntHashMap();
		this.localClusterSize = new ArrayList<OpenIntIntHashMap>();
		for(int j = 0; j < J; j++){
			this.localClusterSize.add(new OpenIntIntHashMap());
		}
		//allocate clusterTermFreq
		this.clusterTermFreq = new OpenIntObjectHashMap();
		this.clusterTotalTermFreq = new OpenIntIntHashMap();
		
		
		//allocate beta
		this.beta = new OpenIntDoubleHashMap();
		this.beta_u = 1.0;
		
		//allocate the two log tables
		//the extreme case is that all documents are assigned to a same cluster
		this.b0SumLogTable = new SumLogValueTable(this.b0, this.docCorpora.totalTokenNumber);
		this.bNSumLogTable = new SumLogValueTable(this.b0 * this.docCorpora.dictSize, this.docCorpora.totalTokenNumber);
		
		//allocate tableNumbers
		this.tableNumbers = new OpenIntIntHashMap();
		
		this.logLevel = 0;
		
		int max_sel = Integer.MAX_VALUE;
		this.CVG_SEL_NUM = new ArrayList<Integer>();
		this.cvgDocs = new ArrayList<List<Integer>>();
		for(int j = 0; j < J; j++){
			this.CVG_SEL_NUM.add(
					(this.docCorpora.docNumbers.get(j) < max_sel)? this.docCorpora.docNumbers.get(j) : max_sel);
			this.cvgDocs.add(new ArrayList<Integer>());
			int[] rpm = this.randomSampler.RandPerm(this.docCorpora.docNumbers.get(j));
			for(int i = 0; i < CVG_SEL_NUM.get(j); i++){
				this.cvgDocs.get(j).add(rpm[i]);
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
	
	public int getLocalClusterSize(int j, int kid){
		return this.localClusterSize.get(j).get(kid);
	}
	
	private void RemoveCluster(int clusterID){
		this.beta.removeKey(clusterID);
		this.globalClusterSize.removeKey(clusterID);
		this.clusterTermFreq.removeKey(clusterID);
		this.clusterTotalTermFreq.removeKey(clusterID);
		this.tableNumbers.removeKey(clusterID);
		for(int j = 0; j < this.docCorpora.corporaNumber; j++){
			this.localClusterSize.get(j).removeKey(clusterID);
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
	public void Initialize(int K){
		int[] alphabet = new int[K];
		for(int i = 0; i < K; i++){
			alphabet[i] = i;
		}
		int J = this.docCorpora.corporaNumber;
		int[][] initlabel = new int[J][];
		for(int j = 0; j < J; j++){
			initlabel[j] = this.randInitLabel(this.docCorpora.docNumbers.get(j), alphabet);
		}
		Initialize(initlabel);
	}
	
	private void ClusterAddDoc(int kid, int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docCorpora.tfVectors.get(corpus_j).get(doc_i);
		int docLen = this.docCorpora.docLengths.get(corpus_j).get(doc_i);
		OpenIntIntHashMap cluster = (OpenIntIntHashMap) this.clusterTermFreq.get(kid);
		VectorOper.VectorAddEqual(cluster, docTf);
		int totalTf = this.clusterTotalTermFreq.get(kid);
		this.clusterTotalTermFreq.put(kid, totalTf + docLen);
		this.Z.get(corpus_j).set(doc_i, kid);
		int njk = this.localClusterSize.get(corpus_j).get(kid);
		this.localClusterSize.get(corpus_j).put(kid, njk + 1);
		int nk = this.globalClusterSize.get(kid);
		this.globalClusterSize.put(kid, nk + 1);
	}
	private void ClusterRemoveDoc(int kid, int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docCorpora.tfVectors.get(corpus_j).get(doc_i);
		int docLen = this.docCorpora.docLengths.get(corpus_j).get(doc_i);
		OpenIntIntHashMap cluster = (OpenIntIntHashMap) this.clusterTermFreq.get(kid);
		VectorOper.VectorMinusEqual(cluster, docTf);
		int totalTf = this.clusterTotalTermFreq.get(kid);
		this.clusterTotalTermFreq.put(kid, totalTf - docLen);
		this.Z.get(corpus_j).set(doc_i, MultiCorporaHDPClustering.NEW_CLUSTER_ID);//not assigned yet
		int njk = this.localClusterSize.get(corpus_j).get(kid);
		this.localClusterSize.get(corpus_j).put(kid, njk - 1);
		if(njk == 1){
			this.localClusterSize.get(corpus_j).removeKey(kid);
		}
		int nk = this.globalClusterSize.get(kid);
		this.globalClusterSize.put(kid, nk - 1);
		if(nk == 1)
			this.RemoveCluster(kid);
	}
	public int getLocalClusterNumber(int j){
		return this.localClusterSize.get(j).size();
	}
	public int[] getLocalClusterNumbers(){
		int[] ret = new int[this.docCorpora.corporaNumber];
		for(int j = 0; j < ret.length; j++){
			ret[j] = this.localClusterSize.get(j).size();
		}
		return ret;
	}
	public void Initialize(int[][] initZ){
		OpenIntIntHashMap labelmap = new OpenIntIntHashMap();
		int J = initZ.length;
		int ilabel;
		int clusterid = -1;
		for(int j = 0; j < J; j++){
			for(int i = 0; i < initZ[j].length; i++){
				ilabel = initZ[j][i];
				if(labelmap.containsKey(ilabel) == false){
					clusterid = AddNewCluster(j);
					labelmap.put(ilabel, clusterid);
				}
				else{
					clusterid = labelmap.get(ilabel);///getComponent(labelmap.get(ilabel));
				}
				
				this.ClusterAddDoc(clusterid, j, i);
			}
		}
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
	
	protected int SamplingDocCluster(int corpus_j, int doc_i){
		OpenIntIntHashMap docTf = this.docCorpora.tfVectors.get(corpus_j).get(doc_i);
		int sumTf = this.docCorpora.docLengths.get(corpus_j).get(doc_i);
		return this.SamplingDocCluster(corpus_j, docTf, sumTf, true);
	}
	protected int SamplingDocCluster(int corpus_j, OpenIntIntHashMap docTf, int sumtf, boolean includeNew){
		int K = this.globalClusterSize.size();
		if(Math.abs(this.alpha) < 10E-8)
			includeNew = false;
		double[] logpost = includeNew? new double[K + 1] : new double[K];
		IntArrayList kids = this.globalClusterSize.keys();
		for(int ki = 0; ki < K; ki++){
			int k = kids.get(ki);
			double beta_k = this.beta.get(k);
			int njk = this.getLocalClusterSize(corpus_j, k);
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
	
	public void SamplingZ(boolean shuffle){
		int J = this.docCorpora.corporaNumber;
		int[] rpmj = shuffle? this.randomSampler.RandPerm(J) : null;
		int D, z, znew;
		
		for(int rj = 0; rj < J; rj++){
			int j = shuffle? rpmj[rj] : rj;
			D = this.docCorpora.docNumbers.get(j);
			int[] rpmd = shuffle? this.randomSampler.RandPerm(D) : null;
			for(int ri = 0; ri < D; ri++){
				int i = shuffle? rpmd[ri] : ri;
				z = this.Z.get(j).get(i);

				//remove the doc from current cluster
				this.ClusterRemoveDoc(z, j, i);
				
				//sampling a new cluster
				znew = this.SamplingDocCluster(j, i);
				
				//assign the new cluster
				this.ClusterAddDoc(znew, j, i);
				
			}//for i
		}//for j
	}
	
	/**
	 * Sampling the number of sub-clusters for each cluster k in each corpus j
	 */
	protected void SamplingMandBeta(){
		int J = this.docCorpora.corporaNumber;
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
				njk = this.localClusterSize.get(j).get(kid);
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
	protected void SamplingAlpha(int iter){
		int J = this.docCorpora.corporaNumber;
		
		List<Integer> nj = this.docCorpora.docNumbers;
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
	
	public void Sampling(){
		
		int conpiter = 20;//default iteration for sampling hyperparameters \alpha and \gamma
		long starttime, time0, usedtime;
		Log(-1, 0);

		for(int iter = 0; iter < this.maxIter; iter++){
			starttime = System.currentTimeMillis();
			time0 = starttime;
			SamplingZ(true);
//			System.out.println("	Sampling Z: " + (System.currentTimeMillis() - starttime) + "ms");
//			starttime = System.currentTimeMillis();
			SamplingMandBeta();
//			System.out.println("	Sampling M and Beta: " + (System.currentTimeMillis() - starttime) + "ms");
		
			if(this.bSamplingHyper){
//				starttime = System.currentTimeMillis();
				SamplingAlpha(conpiter);
//				System.out.println("	Sampling alpha: " + (System.currentTimeMillis() - starttime) + "ms");
//				
//				starttime = System.currentTimeMillis();
				SamplingGamma(conpiter);
//				System.out.println("	Sampling gamma: " + (System.currentTimeMillis() - starttime) + "ms");
			}
			
			usedtime = System.currentTimeMillis() - time0;
			Log(iter, usedtime);	
			
			if (iter % 10 == 0 || iter == (this.maxIter - 1) ) {
				CorporaClusteringPostInfo topicInfo1 = new CorporaClusteringPostInfo(
						this.docCorpora.dictionary
						, this.docCorpora.forwardHash
						, this.docCorpora.backwardHash
						, this.docCorpora.tfVectors
						, null
						, this.Z
						, null
						, null
						, null
						, null);
				int topN = 40;
				String msg = topicInfo1.PrintTopicsKeywords(topN);
				for(int j = 0; j < this.docCorpora.corporaNumber; j++){
					msg += "\n corpus " + j + "---------------------------\n";
					msg += topicInfo1.getLocalInfo(j).PrintTopicsKeywords(topN);
				}
				System.out.println(msg);
			}
		}
		this.docCorpora.docLabels = this.Z;
	}
	
	protected void Log(int iter, long usedTime){
		int K = -1;
		int[] Ks = null;
		if(this.logLevel > 1)
			K = this.getClusterNumber();
		
		if(this.logLevel > 2)
			Ks = this.getLocalClusterNumbers();
		
		double loglike = 0;
		if(this.logLevel >= 3){
			loglike = CalLogLike(this.cvgDocs);
		}
		
		String log = "[LOG]>>";
		switch(this.logLevel){
		case 0:
			return;//release
		case 1://iter + time
			log += "Iteration=" + iter + "; Time="+ usedTime + "ms;";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 2:// iter + K + time
			log += "Iteration=" + iter + "; Time="+ usedTime + "ms;" 
			+ " K=" + K +";";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 3://iter + K + like + time
			log += "Iteration=" + iter + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 4://iter + K + Ks + like + time
			log += "Iteration=" + iter  + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; Ks=" + VectorOper.ToString(Ks) + 
			"; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 5://iter + K + Ks + like + time + nu + beta + topicsize
			log += "Iteration=" + iter  + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; Ks=" + VectorOper.ToString(Ks) + 
			"; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += ";gamma=" + this.gamma + "; alpha = " + this.alpha;
			
			double[] beta = new double[K + 1];
			int[] kids, sizes;
			//global info
			kids = new int[K + 1];
			sizes = new int[K + 1];
			int[] tableNum = new int[K + 1];
			tableNum[0] = 0;
			beta[0] = this.beta_u;
			kids[0] = MultiCorporaHDPClustering.NEW_CLUSTER_ID;
			sizes[0] = 0;
			
			IntArrayList clusterIDs = this.globalClusterSize.keys();
			for(int k = 1; k < K + 1; k++){
				kids[k] = clusterIDs.get(k - 1);
				sizes[k] = this.globalClusterSize.get(kids[k]);
				beta[k] = this.beta.get(kids[k]);
				tableNum[k] = this.tableNumbers.get(kids[k]);
			}
			log += "\nID		|" + VectorOper.ToString(kids, "", "	|", "%6d") + "\n";
			log += "beta(1E-2)	|" + VectorOper.ToString(VectorOper.VectorTimes(beta, 100)
					, "", "	|", "%6.3f") + "\n";
			log += "sizes		|" + VectorOper.ToString(sizes, "", "	|", "%6d") + "\n";
			log += "table		|" + VectorOper.ToString(tableNum, "", "	|", "%6d") + "\n";
			log += "totalsize = " + VectorOper.VectorSum(sizes);
			
			
			
			for(int j = 0; j < this.docCorpora.corporaNumber; j++){
				
				IntArrayList jTopicIds = this.localClusterSize.get(j).keys();
				int Kj = jTopicIds.size();
				int[] jkids = new int[Kj + 1];
				jkids[0] = MultiCorporaHDPClustering.NEW_CLUSTER_ID;
				sizes = new int[Kj + 1];
				sizes[0] = 0;
				for(int ki = 1; ki < Kj + 1; ki++){
					jkids[ki] = jTopicIds.get(ki - 1);
					sizes[ki] = this.getLocalClusterSize(j, jkids[ki]);
				}
				log += "\n	corpus " + j + "-------------------------------------";
				log += "\n	ID		|" + VectorOper.ToString(jkids, "", "	|", "%6d") + "\n";
				log += "	sizes		|" + VectorOper.ToString(sizes, "", "	|", "%6d") + "\n";
				log += "	totalsize = " + VectorOper.VectorSum(sizes);
			}
			
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
	private double CalLogLike(List<List<Integer>> docs){
		double ret = 0;
		int J = docs.size();
		int kid = 0;
		OpenIntIntHashMap tfDoc = null;
		for(int j = 0; j < J; j++){
			for(int i = 0; i < docs.get(j).size(); i++){
				int d = docs.get(j).get(i);
				tfDoc = this.docCorpora.tfVectors.get(j).get(d);
				int tfsum = this.docCorpora.docLengths.get(j).get(d);
				kid = this.Z.get(j).get(d);
				ret += this.ClusterLogMarginal(kid, tfDoc, tfsum);
			}
		}
		return ret;
	}
	

	public DocCorpora getDocTimeCorpora() {
		return docCorpora;
	}
	public void setDocTimeCorpora(DocTimeCorpora docTimeCorpora) {
		this.docCorpora = docCorpora;
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
//		double gamma = 10, alpha = 1;
//		double gamma_a = 2, gamma_b = 1;
//		double alpha_a = 0.1, alpha_b = 0.1;
//		
//		int maxIter = 50, burnin = 10;
//		int initK = 15;
//		int loglevel = 3;
//		MultiCorporaHDPClustering hdpCluster = new MultiCorporaHDPClustering(
//				docCorpora
//				, b0
//				, gamma_a, gamma_b
//				, alpha_a, alpha_b
////				, gamma, alpha //if sampling hyper, comment this line and uncomment above 2 lines, vice versa
//				, maxIter
//				, burnin);
//		hdpCluster.Initialize(initK);
//		hdpCluster.setLogLevel(loglevel);
//		hdpCluster.Sampling();

	}
	
	
}