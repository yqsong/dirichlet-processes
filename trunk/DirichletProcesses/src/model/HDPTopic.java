/**
 * HDPTopic.java
 * @author Jianwen Zhang (jw-zhang06@mails.tsinghua.edu.cn, jianwenzh@gmail.com)
 * Feb. 21, 2010.
 * 
 * HDP is implemented as an infinite LDA topic model, using direct posterior sampling
 * scheme.
 * 
 */
package model;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import util.RandomSampler;
import util.VectorOper;

import cc.mallet.types.Alphabet;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntDoubleHashMap;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.DocCorpora;
import data.structure.DocCorpus;
import description.CorpusTopicPostInfo;

public class HDPTopic{
	
	DocCorpus docCorpus;
	/**
	 * The top measure's parameter (Dirichlet distribution with symmetric parameters of \b_0)
	 */
	double b0;
	/**
	 * Z[j][i] is z_{ji}, i.e., the topic assignment of token i of document j
	 */
	List<int[]> Z;
	/**
	 * docTopicSize[j].get(k) is number of tokens assigned to topic k, i.e., n_{jk} 
	 */
	List<OpenIntIntHashMap> docTopicSize;
	/**
	 * topicTermFreq.get(k) is a OpenIntIntHashMap mp, and mp.get(w) is the frequence
	 * of term w appearing in topic k. 
	 * In fact, topicTermFreq describes the matrix \Phi, i.e., topic's word distribution
	 */
	OpenIntObjectHashMap topicTermFreq;
	/**
	 * topicTotalTermFreq.get(k) = \sum_{w=1}^W topicTermFreq.get(k).get(w),
	 * preserve this summation to speed up sampling.
	 */
	OpenIntIntHashMap topicTotalTermFreq;
	/**
	 * the probabilities for measure G_0, not containing \beta_u
	 */
	OpenIntDoubleHashMap beta;
	
	/**
	 * tableNumbers.get(k) is the number of tables (sub-clusters) with 
	 * dish (topic) k in all restaurants (documents)
	 */
	OpenIntIntHashMap tableNumbers;
	/**
	 * probability of G_0 on a new component
	 */
	double beta_u;
	double newTopicLike;//1 / W;
	double b0TimesW;//b0 * W
	/**
	 * index for a topic
	 */
	int topicID;
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
	 * 0: release; 1: iter + time; 2: iter + K + time; 3: iter + K + Ks + time; 
	 * 4: iter + K + Ks + like + time; 5: iter + K + Ks + like + time + beta + topicsizes
	 */
	int logLevel;
	int CVG_SEL_NUM;
	int[] cvgDocs;//the set of documents selected for convergence diagnosis
	
	/**
	 * Constructor for the case not need to sample hyper parameters \alpha or \gamma
	 * @param docCorpus
	 * @param b0
	 * @param alpha
	 * @param gamma
	 * @param maxIter
	 * @param burnIn
	 */
	public HDPTopic(
			DocCorpus docCorpus
			, double b0
			, double alpha, double gamma
			, int maxIter
			, int burnIn
			){
		this(docCorpus
			, b0
			, 0.0, 0.0
			, 0.0, 0.0
			, alpha, gamma
			, false//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
	}
	/**
	 * Constructor for the case when sampling hyper parameters are required.
	 * @param docCorpus
	 * @param b0
	 * @param alpha_a
	 * @param alpha_b
	 * @param gamma_a
	 * @param gamma_b
	 * @param maxIter
	 * @param burnIn
	 */
	public HDPTopic(
			DocCorpus docCorpus
			, double b0
			, double alpha_a, double alpha_b
			, double gamma_a, double gamma_b
			, int maxIter
			, int burnIn
			){
		this(docCorpus
			, b0
			, alpha_a, alpha_b
			, gamma_a, gamma_b
			, 0.0, 0.0
			, true//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
	}
	/**
	 * Private constructor.
	 * @param docCorpus, a DocCorpus structure containing the feature sequences
	 * @param b0, the parameter of top measure, i.e., the symmetric Dirichlet distribution 
	 * @param alpha_a, shape parameter of the vague gamma distribution for \alpha
	 * @param alpha_b, inverse scale parameter of the vague gamma distribution for \alpha
	 * @param gamma_a, shape parameter of the vague gamma distribution for \gamma
	 * @param gamma_b, inverse scale parameter of the vague gamma distribution for \gamma
	 * @param alpha, if not need to sample hyper-parameters, \alpha is specified.
	 * @param gamma, if not need to sample hyper-parameters, \gamma is specified.
	 * @param bSamplingHyper, whether need to sample hyper parameters. 
	 * @param maxIter, maximal iteration steps
	 * @param burnIn, the burn-in time.
	 */
	protected HDPTopic(
			DocCorpus docCorpus
			, double b0
			, double alpha_a, double alpha_b
			, double gamma_a, double gamma_b
			, double alpha, double gamma
			, boolean bSamplingHyper//whether sampling hyper parameters \alpha and \gamma
			, int maxIter
			, int burnIn
			){
		this.docCorpus = docCorpus;
		this.b0 = b0;
		this.maxIter = maxIter;
		this.burnIn = burnIn;
		this.alpha_a = alpha_a;
		this.alpha_b = alpha_b;
		this.gamma_a = gamma_a;
		this.gamma_b = gamma_b;
		this.bSamplingHyper = bSamplingHyper;
		
		this.randomSampler = new RandomSampler();
		
		//initialize hyperparameter
		this.alpha = (bSamplingHyper)? (this.alpha_a / this.alpha_b) : alpha;
		this.gamma = (bSamplingHyper)? (this.gamma_a / this.gamma_b) : gamma;
		
		//allocate assignments
		int J = this.docCorpus.docNumber;
		this.Z = new ArrayList<int[]>();
		for(int j = 0; j < J; j++){
			this.Z.add( new int[this.docCorpus.docLengths.get(j)] );
		}
		
		//allocate docTopicSize
		this.docTopicSize = new ArrayList<OpenIntIntHashMap>();
		for(int j = 0; j < J; j++){
			this.docTopicSize.add( new OpenIntIntHashMap() );
		}
		//allocate topicTermFreq
		this.topicTermFreq = new OpenIntObjectHashMap();
		this.topicTotalTermFreq = new OpenIntIntHashMap();
		
		//allocate beta
		this.beta = new OpenIntDoubleHashMap();
		this.beta_u = 1.0;
		this.newTopicLike = 1.0 / this.docCorpus.dictSize;
		this.b0TimesW = this.b0 * this.docCorpus.dictSize;
		this.topicID = 1;//leave 0 as the ID of a new topic
		
		//allocate tableNumbers
		this.tableNumbers = new OpenIntIntHashMap();
		
		this.logLevel = 0;
		int max_sel = Integer.MAX_VALUE;
		this.CVG_SEL_NUM = (this.docCorpus.docNumber < max_sel)? this.docCorpus.docNumber : max_sel;
		this.cvgDocs = new int[this.CVG_SEL_NUM];
		int[] rpm = this.randomSampler.RandPerm(this.docCorpus.docNumber);
		for(int i = 0; i < CVG_SEL_NUM; i++){
			this.cvgDocs[i] = rpm[i];
		}
	}
//	public int getDocNumber(){
//		return this.fetSequences.length;
//	}
//	public int getDocLen(int j){
//		return this.fetSequences[j].length;
//	}
	
	/**
	 * predict the topic assignment for a held out token
	 */
	protected int PredictTokenTopic(int wIdx){
		IntArrayList topicIDs = this.topicTermFreq.keys();
		int K = topicIDs.size();
		double[] post = new double[K];
		for(int k = 0; k < K; k++){
			int kid = topicIDs.get(k);
			OpenIntIntHashMap topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(kid);
			double beta_k = this.beta.get(kid);
			int	nkw = topicDist.get(wIdx);
			int	nk = this.topicTotalTermFreq.get(kid);
			
			post[k] = (this.alpha * beta_k) * (nkw + b0) / (nk + this.b0TimesW);
		}
		double ss = VectorOper.VectorSum(post);
		VectorOper.VectorTimesEqual(post, 1/ss);
		int zkid = this.randomSampler.nextDiscrete(post);
		int z = topicIDs.get(zkid);
		return z;
		
	}
	/**
	 * Predict the topic assignments for a held out document
	 * @param docFtSeq, the feature sequence of the document
	 * @return the predicted topic sequence of the document
	 */
	public int[] PredictDocTopics(int[] docFtSeq){
		int[] zs = new int[docFtSeq.length];
		for(int i = 0; i < docFtSeq.length; i++){
			zs[i] = this.PredictTokenTopic(docFtSeq[i]);
		}
		return zs;
	}
	public int AddNewTopic(){
		int id = this.topicID++;
		//sampling beta_u, beta_K
		double nu = this.randomSampler.nextBeta(this.gamma, 1.0);
		double beta_k = this.beta_u * (1 - nu);
		this.beta_u = this.beta_u * nu;
		this.beta.put(id, beta_k);
		this.topicTermFreq.put(id, new OpenIntIntHashMap());
		return id;
	}
	public void RemoveTopic(int topicID){
		this.beta.removeKey(topicID);
		for(int j = 0; j < this.docCorpus.docNumber; j++){
			this.docTopicSize.get(j).removeKey(topicID);
		}
		this.topicTermFreq.removeKey(topicID);
		this.topicTotalTermFreq.removeKey(topicID);
		this.tableNumbers.removeKey(topicID);
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
		int J = this.docCorpus.docNumber;
		int[][] initlabel = new int[J][];
		for(int j = 0; j < J; j++){
			initlabel[j] = this.randInitLabel(this.docCorpus.docLengths.get(j), alphabet);
		}
		Initialize(initlabel);
	}
	
//	private void assignTopic(int doc_j, int word_i, int topicid){
//		int wji = this.fetSequences[doc_j][word_i];
//		this.Z[doc_j][word_i] = topicid;
//		
//		int njk = this.docTopicSize[doc_j].get(topicid);
//		this.docTopicSize[doc_j].put(topicid, njk + 1);
//		
//		OpenIntIntHashMap topicDist = (OpenIntIntHashMap)this.topicTermFreq.get(topicid);
//		int nkw = topicDist.get(wji);
//		topicDist.put(wji, nkw + 1);
//	}
//	private void removeAssign(int doc_j, int word_i){
//		this.Z[doc_j]
//	}
	public void Initialize(int[][] initZ){
		OpenIntIntHashMap labelmap = new OpenIntIntHashMap();
		int J = initZ.length;
		int ilabel;
		int topicid = -1;
		int wji = -1;
		int njk, nkw, nk;
		OpenIntIntHashMap topicDist = null;
		for(int j = 0; j < J; j++){
			for(int i = 0; i < initZ[j].length; i++){
				ilabel = initZ[j][i];
				wji = this.docCorpus.ftSequences.get(j)[i];
				if(labelmap.containsKey(ilabel) == false){
					topicid = AddNewTopic();
					labelmap.put(ilabel, topicid);
				}
				else{
					topicid = labelmap.get(ilabel);///getComponent(labelmap.get(ilabel));
				}
				
				//assign the word to the topic
				this.Z.get(j)[i] = topicid;
				//update count variables
				//update the count n_{jk} of words assigned to a topic k in document j
				njk = this.docTopicSize.get(j).get(topicid);
				this.docTopicSize.get(j).put(topicid, njk + 1);
				//update the count n_{kw} of word w assigned to topic k in all documents
				topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(topicid);
				nkw = topicDist.get(wji);
				topicDist.put(wji, nkw + 1);
				//update the summation variable of n_k = \sum_w n_{kw}
				nk = this.topicTotalTermFreq.get(topicid);
				this.topicTotalTermFreq.put(topicid, nk + 1);
			}
		}
	}
	
	public int getTopicNumber(){
		return this.topicTermFreq.size();
	}
	
	public void SamplingZ(boolean shuffle){
		int J = this.docCorpus.docNumber;
		int[] rpm = shuffle? this.randomSampler.RandPerm(J) : null;
		int d, w, L, z, znew, K;
		int njk, nkw, nk;
		int kid;
		double beta_k, ss;
		OpenIntIntHashMap topicDist = null;
		OpenIntIntHashMap docTopicSize = null;
		for(int j = 0; j < J; j++){
			d = shuffle? rpm[j] : j;
			L = this.docCorpus.docLengths.get(d);
			for(int i = 0; i < L; i++){
				w = this.docCorpus.ftSequences.get(d)[i];
				z = this.Z.get(d)[i];
				
				//sampling a new topic
				docTopicSize = this.docTopicSize.get(d);
				
				IntArrayList topicIDs = this.topicTermFreq.keys();
				K = topicIDs.size();
				double[] post = new double[K + 1];
				for(int k = 0; k < K; k++){
					kid = topicIDs.get(k);
					topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(kid);
					beta_k = this.beta.get(kid);
					if(kid != z){
						njk = docTopicSize.get(kid);
						nkw = topicDist.get(w);
						nk = this.topicTotalTermFreq.get(kid);
					}
					else{
						njk = docTopicSize.get(kid) - 1;
						nkw = topicDist.get(w) - 1;
						nk = this.topicTotalTermFreq.get(kid) - 1;
					}
					//nk = 0 means the topic k should be removed
					post[k] = (nk == 0)? 0 : ((njk + this.alpha * beta_k) * (nkw + this.b0) / (nk + this.b0TimesW));
				}
				post[K] = this.alpha * this.beta_u * this.newTopicLike;
				ss = VectorOper.VectorSum(post);
				VectorOper.VectorTimesEqual(post, 1/ss);
				znew = this.randomSampler.nextDiscrete(post);
				znew = (znew < K)? topicIDs.get(znew) : this.AddNewTopic();
				//if znew == z, nothing need to update
				if(znew != z){
					this.Z.get(d)[i] = znew;
					//remove i from original topic
					njk = docTopicSize.get(z);
					if(njk == 1){
						docTopicSize.removeKey(z);
					}
					else{
						docTopicSize.put(z, njk - 1);
					}
					topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(z);
					nkw = topicDist.get(w);
					if(nkw == 1){
						topicDist.removeKey(w);
					}
					else{
						topicDist.put(w, nkw - 1);
					}
					nk = this.topicTotalTermFreq.get(z);
					if(nk == 1){
						this.RemoveTopic(z);
					}
					else{
						this.topicTotalTermFreq.put(z, nk - 1);
					}
					
					//assign the new topic
					njk = docTopicSize.get(znew);
					docTopicSize.put(znew, njk + 1);
					
					topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(znew);
					nkw = topicDist.get(w);
					topicDist.put(w, nkw + 1);
					
					nk = this.topicTotalTermFreq.get(znew);
					this.topicTotalTermFreq.put(znew, nk + 1);
					
				}//if
				
			}//for i
		}//for j
	}
	
	/**
	 * Sampling the number of sub-clusters for each topic k in each document j
	 */
	protected void SamplingMandBeta(){
		int J = this.docCorpus.docNumber;
		int njk = 0, kid = -1, mk = 0;
		this.tableNumbers.clear();
		
		IntArrayList topicIDs = this.topicTermFreq.keys();
		int K = topicIDs.size();
		double[] ms = new double[K + 1];
		ms[K] = this.gamma;
		double bb = 0;
		for(int k = 0; k < K; k++){
			kid = topicIDs.get(k);
			mk = 0;
			bb = this.alpha * this.beta.get(kid);
			for(int j = 0; j < J; j++){
				njk = this.docTopicSize.get(j).get(kid);
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
		
//		// test smooth beta
//		double meanBeta = 0.0;
//		for (int k = 0; k < K; ++k) {
//			meanBeta += beta[k];
//		}
//		for (int k = 0; k < K; ++k) {
//			beta[k] = (beta[k] + meanBeta/K) / (meanBeta * 2);
//		}
//		// end test smooth bets
		
		for(int k = 0; k < K; k++){
			kid = topicIDs.get(k);
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
		
		int K = this.getTopicNumber();
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
		int J = this.docCorpus.docNumber;
		
		List<Integer> nj = this.docCorpus.docLengths;
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
				CorpusTopicPostInfo topicInfo = new CorpusTopicPostInfo(
						this.docCorpus.dictionary
						, this.docCorpus.forwardHash
						, this.docCorpus.backwardHash
						, this.docCorpus.ftSequences, this.Z, null, null, null);
				String msg = topicInfo.PrintTopicsKeywords(40);
				System.out.println(msg);
			}
		}
	}
	
	protected void Log(int iter, long usedTime){
		int K = -1;
		int[] Ks = null;
		if(this.logLevel > 1)
			K = this.getTopicNumber();
		
		if(this.logLevel > 2)
			Ks = getDocTopicNums(this.cvgDocs);
		
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
				log += "gamma = " + this.gamma + "; alpha = " + this.alpha; 
			break;
			
		case 2:// iter + K + time
			log += "Iteration=" + iter + "; Time="+ usedTime + "ms;" 
			+ " K=" + K +";";
			if(this.bSamplingHyper)
				log +=  "  gamma = " + this.gamma + "; alpha = " + this.alpha  ;
			break;
			
		case 3://iter + K + like + time
			log += "Iteration=" + iter + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log +="  gamma = " + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 4://iter + K + Ks + like + time
			log += "Iteration=" + iter  + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; Ks=" + VectorOper.ToString(Ks) + 
			"; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += "  gamma = " + this.gamma + "; alpha = " + this.alpha;
			break;
			
		case 5://iter + K + Ks + like + time + beta + topicsize
			log += "Iteration=" + iter  + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; Ks=" + VectorOper.ToString(Ks) + 
			"; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += "  gamma = " + this.gamma + "; alpha = " + this.alpha + "\n";
			
			double[] betas = new double[K + 1];
			int[] kids = new int[K + 1];
			int[] sizes = new int[K + 1];
			betas[0] = this.beta_u;
			kids[0] = 0;
			sizes[0] = 0;
			
			IntArrayList topicIDs = this.topicTermFreq.keys();
			for(int k = 1; k < K + 1; k++){
				kids[k] = topicIDs.get(k - 1);
				sizes[k] = this.topicTotalTermFreq.get(kids[k]);
				betas[k] = this.beta.get(kids[k]);
			}
			log += "\nID	|" + VectorOper.ToString(kids, "", "	|", "%6d") + "\n";
			log += "beta	|" + VectorOper.ToString(betas, "", "	|", "%.4f") + "\n";
			log += "sizes	|" + VectorOper.ToString(sizes, "", "	|", "%6d") + "\n";
			log += "totalsize = " + VectorOper.VectorSum(sizes);
			break;
			
		default:
			break;
		}
		if(this.logLevel > 0)
			System.out.println(log);
	}
	
	int[] getDocTopicNums(){
		int J = this.docCorpus.docNumber;
		int[] ret = new int[J];
		for(int j = 0; j < J; j++){
			ret[j] = this.docTopicSize.get(j).size();
		}
		return ret;
	}
	
	int[] getDocTopicNums(int[] docs){
		int[] ret = new int[docs.length];
		for(int j = 0; j < docs.length; j++){
			ret[j] = this.docTopicSize.get(docs[j]).size();
		}
		return ret;
	}
	
	/**
	 * not strictly the likelihood, but just cumulation of the likelihood of each word
	 * from its topic,  for convergence diagnosis
	 * @return
	 */
	private double CalLogLike(int[] docs){
		double ret = 0;
		int J = docs.length;
		int kid = 0, nkw, nk;
		int d,w;
		for(int j = 0; j < J; j++){
			d = docs[j];
			for(int i = 0; i < this.docCorpus.ftSequences.get(d).length; i++){
				w = this.docCorpus.ftSequences.get(d)[i];
				kid = this.Z.get(d)[i];
				nkw = ((OpenIntIntHashMap) this.topicTermFreq.get(kid)).get(w);
				nk = this.topicTotalTermFreq.get(w);
				ret += Math.log((nkw + this.b0) / (nk + this.b0TimesW));
			}
		}
		return ret;
	}
	
	
	
	/**
	 * doc-topic count matrix. 
	 * @return ret.
	 * ret[d][k] is number of words assigned to topic k in document d, i.e., the empirical unormalized
	 * estimation of \theta in LDA
	 */
	public List<OpenIntIntHashMap> getDocTopicMatrix(){
		return this.docTopicSize;
	}
	/**
	 * topic-word empirical count matrix
	 * @return ret
	 * ret.get(k) is an OpenIntIntHashMap data, e.g., dist, then dist.get(w) is the term frequency
	 * of w appearing in topic k.
	 */
	public OpenIntObjectHashMap getTopicDists(){
		return this.topicTermFreq;
	}
	/**
	 * topic size: the number of words assigned to a topic
	 * @return ret
	 * ret.get(k) is the size of topic k.
	 */
	public OpenIntIntHashMap getTopicSizes(){
		return this.topicTotalTermFreq;
	}
	public static void main(String[] args){
		
//		double b0 = 0.5;
//		double gamma = 5, alpha = 1;
//		double gamma_a = 10, gamma_b = 1;
//		double alpha_a = 10, alpha_b = 1;
//		
//		int maxIter = 1000, burnin = 10;
//		int initK = 20;
//		int loglevel = 3;
//		HDPTopic hdptopic = new HDPTopic(
//				docCorpus
//				, b0
////				, alpha, gamma   //if sampling hyper, comment this line and uncomment following two lines, vice versa
//				, alpha_a, alpha_b
//				, gamma_a, gamma_b
//				, maxIter
//				, burnin);
//		hdptopic.Initialize(initK);
//		hdptopic.setLogLevel(loglevel);
//		hdptopic.Sampling();
		

	}

}