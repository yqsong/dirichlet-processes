/**
 * MultiHDPTopic.java
 * @author Jianwen Zhang (jw-zhang06@mails.tsinghua.edu.cn, jianwenzh@gmail.com)
 * HDP-LDA for multiple text corpora.
 */
package model;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import util.RandomSampler;
import util.VectorOper;

import bsh.This;

import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntDoubleHashMap;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.DocCorpora;
import description.CorporaTopicPostInfo;

public class MultiCoporaHDPTopic{
	DocCorpora docCorpora;
	
	/**
	 * The top measure's parameter (Dirichlet distribution with symmetric parameters of \b_0)
	 */
	double b0;
	
	/**
	 * Z[j][d][i] is z^j_{di}, i.e., the topic assignment of token i of document d of corpus j
	 */
	List<List<int[]>> Z;
	/**
	 * docTopicSize[j][d].get(k) is number of tokens assigned to topic k, i.e., n^j_{dk} 
	 */
	List<List<OpenIntIntHashMap>> docTopicSize;
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
	 * the probabilities for measure G_0, not containing \nu_u
	 */
	OpenIntDoubleHashMap nu;
	double nu_u;
	/**
	 * betas[j] are the probabilities for measure G^j, not containing \beta^j_u
	 */
	List<OpenIntDoubleHashMap> betas;
	List<Double> beta_us;
	
	/**
	 * tableNumbers[j].get(k) is the number of tables (sub-clusters) with 
	 * dish (topic) k in all (documents) of corpus j. Tables in corpus j are drawn from G^j.
	 * These count variables are used to sample G^j
	 */
	List<OpenIntIntHashMap> tableNumbers;
	
	/**
	 * metaTableNumbers.get(k) is the number of meta-tables (sub-clusters of tables) with 
	 * dish (topic) k in all corpora. All metatables are drawn from G_0.
	 * These count variables are used to sample G_0.
	 */
	OpenIntIntHashMap metaTableNumbers;
	
	double newTopicLike;//1 / W;
	double b0TimesW;//b0 * W
	/**
	 * index for a topic
	 */
	int topicID;
	static int NEW_TOPIC_ID = -1;//actually not used but for result report
	/**
	 * concentration parameter for $G_0 \sim DP(\xi, H)$
	 */
	double xi;
	/**
	 * concentration parameter for $G^j \ sim DP(\gamma, G_0)$
	 */
	double gamma;
	/**
	 * concentration parameter for $G^j_d \ sim DP(\alpha^j, G^j)$
	 */
	List<Double> alphas;
	/**
	 * parameters of the vague gamma prior for \xi \sim Gamma(\xi_a, \xi_b), where
	 * \xi_a is the shape parameter and \xi_b is the inverse scale parameter.
	 */
	double xi_a, xi_b;
	/**
	 * parameters of the vague gamma prior for \gamma \sim Gamma(\gamma_a, \gamma_b), where
	 * \gamma_a is the shape parameter and \gamma_b is the inverse scale parameter.
	 */
	double gamma_a, gamma_b;
	/**
	 * parameters of the vague gamma prior for \alpha^d \sim Gamma(\alpha^d_a, \alpha^d_b), where
	 * \alpha^d_a is the shape parameter and \alpha^d_b is the inverse scale parameter.
	 */
	List<Double> alphas_a, alphas_b;
	
	RandomSampler randomSampler;
	
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
	 * @param xi
	 * @param gamma
	 * @param alpha
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCoporaHDPTopic(
			DocCorpora docCorpora
			, double b0
			, double xi, double gamma, double alpha
			, int maxIter
			, int burnIn
			){
		this(
			docCorpora
			, b0
			, 0.0, 0.0
			, 0.0, 0.0
			, 0.0, 0.0
			, xi, gamma, alpha
			, false//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
	}
	/**
	 * Constructor for the case when sampling hyper parameters are required.
	 * @param docCorpora
	 * @param b0
	 * @param xi_a, xi_b
	 * @param gamma_a, gamma_b
	 * @param alpha_a, alpha_b
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCoporaHDPTopic(
			DocCorpora docCorpora
			, double b0
			, double xi_a, double xi_b
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, int maxIter
			, int burnIn
			){
		this(docCorpora
			, b0
			, xi_a, xi_b
			, gamma_a, gamma_b
			, alpha_a, alpha_b
			, 0.0, 0.0, 0.0
			, true//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
	}
	/**
	 * Private constructor.
	 * @param docCorpora, the DocCorpora structure containing the feature sequences of all
	 * documents
	 * @param b0, the parameter of top measure, i.e., the symmetric Dirichlet distribution
	 * @param xi_a, xi_b, the vague gamma hyper parameters for xi~Gamma(xi_a, xi_b),
	 * where xi_a is the shape parameter and xi_b is the inverse scale parameter
	 * @param gamma_a, gamma_b, hyper parameters of the vague gamma prior for \gamma,
	 * gamma_a is the shape parameter and gamma_b is the inverse scale parameter.
	 * @param alpha_a, alpha_b, hyper parameters of the vague gamma prior for \alpha^j,
	 * gamma_a is the shape parameter and gamma_b is the inverse scale parameter.
	 * @param xi, if not need to sample hyper parameters, xi is specified.
	 * @param gamma, if not need to sample hyper-parameters, \gamma is specified.
	 * @param alpha, if not need to sample hyper-parameters, \alpha^j is specified.
	 * @param bSamplingHyper, whether need to sample hyper parameters. 
	 * @param maxIter, maximal iteration steps
	 * @param burnIn, the burn-in time.
	 */
	protected MultiCoporaHDPTopic(
			DocCorpora docCorpora
			, double b0
			, double xi_a, double xi_b
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, double xi, double gamma, double alpha
			, boolean bSamplingHyper//whether sampling hyper parameters \alpha and \gamma
			, int maxIter
			, int burnIn
			){
		this.docCorpora = docCorpora;
		this.b0 = b0;
		this.maxIter = maxIter;
		this.burnIn = burnIn;
		this.bSamplingHyper = bSamplingHyper;
		this.xi_a = xi_a;
		this.xi_b = xi_b;
		this.xi = (bSamplingHyper)? this.xi_a / this.xi_b : xi;
		this.gamma_a = gamma_a;
		this.gamma_b = gamma_b;
		this.gamma = (bSamplingHyper)? this.gamma_a / this.gamma_b : gamma;
		
		int J = this.docCorpora.corporaNumber;
		this.alphas_a = new ArrayList<Double>();
		this.alphas_b = new ArrayList<Double>();
		this.alphas = new ArrayList<Double>();
		for(int j = 0; j < J; j++){
			this.alphas_a.add(alpha_a);
			this.alphas_b.add(alpha_b);
			this.alphas.add((bSamplingHyper)? alpha_a / alpha_b : alpha);;
		}
		
		this.randomSampler = new RandomSampler();
		
		//allocate assignments
		this.Z = new ArrayList<List<int[]>>();
		for(int j = 0; j < J; j++){
			this.Z.add(new ArrayList<int[]>());
			for(int d = 0; d < this.docCorpora.docNumbers.get(j); d++){
				this.Z.get(j).add(new int[this.docCorpora.docLengths.get(j).get(d)]);
			}
		}
		
		//allocate docTopicSize
		this.docTopicSize = new ArrayList<List<OpenIntIntHashMap>>();
		for(int j = 0; j < J; j++){
			this.docTopicSize.add(new ArrayList<OpenIntIntHashMap>());
			for(int d = 0; d < this.docCorpora.docNumbers.get(j); d++){
				this.docTopicSize.get(j).add(new OpenIntIntHashMap());
			}
		}
		//allocate topicTermFreq
		this.topicTermFreq = new OpenIntObjectHashMap();
		this.topicTotalTermFreq = new OpenIntIntHashMap();
		
		//allocate nu
		this.nu = new OpenIntDoubleHashMap();
		this.nu_u = 1.0;
		//allocate beta
		this.betas = new ArrayList<OpenIntDoubleHashMap>();
		this.beta_us = new ArrayList<Double>();
		for(int j = 0; j < J; j++){
			this.betas.add(new OpenIntDoubleHashMap());
			this.beta_us.add(1.0);
		}
		this.newTopicLike = 1.0 / this.docCorpora.dictSize;
		this.b0TimesW = this.b0 * this.docCorpora.dictSize;
		this.topicID = 1;//leave 0 as the ID of a new topic
		
		//allocate tableNumbers
		this.metaTableNumbers = new OpenIntIntHashMap();
		this.tableNumbers = new ArrayList<OpenIntIntHashMap>();
		for(int j = 0; j < J; j++){
			this.tableNumbers.add(new OpenIntIntHashMap());
		}
		
		this.logLevel = 0;
		
		//select a subset of documents to convergence diagnosis
		int max_sel = 50;
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
//	public int getDocNumber(){
//		return this.fetSequences.length;
//	}
//	public int getDocLen(int j){
//		return this.fetSequences[j].length;
//	}
	
	/**
	 * predict the topic assignment for a held out token
	 */
	protected int PredictTokenTopic(int corpusIdx, int wIdx){
		double alpha = this.alphas.get(corpusIdx);
		
		IntArrayList topicIDs = this.topicTermFreq.keys();
		int K = topicIDs.size();
		double[] post = new double[K];
		for(int k = 0; k < K; k++){
			int kid = topicIDs.get(k);
			OpenIntIntHashMap topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(kid);
			double beta_k = this.betas.get(corpusIdx).get(kid);
			int	nkw = topicDist.get(wIdx);
			int	nk = this.topicTotalTermFreq.get(kid);
			
			post[k] = (alpha * beta_k) * (nkw + b0) / (nk + this.b0TimesW);
		}
		double ss = VectorOper.VectorSum(post);
		VectorOper.VectorTimesEqual(post, 1/ss);
		int zkid = this.randomSampler.nextDiscrete(post);
		int z = topicIDs.get(zkid);
		return z;
		
	}
	/**
	 * Predict the topic assignments for a held out document
	 * @param corpusIdx, the index of the corpus
	 * @param docFtSeq, the feature sequence of the document
	 * @return the predicted topic sequence of the document
	 */
	public int[] PredictDocTopics(int corpusIdx, int[] docFtSeq){
		int[] zs = new int[docFtSeq.length];
		for(int i = 0; i < docFtSeq.length; i++){
			zs[i] = this.PredictTokenTopic(corpusIdx, docFtSeq[i]);
		}
		return zs;
	}
	public int AddNewTopic(){
		int id = this.topicID++;
		//sampling nu_u, nu_K
		double nu_split = this.randomSampler.nextBeta(this.xi, 1.0);
		double nu_K = this.nu_u * (1 - nu_split);
		this.nu.put(id, nu_K);
		this.nu_u = this.nu_u * nu_split;
		//sampling beta_u^j, beta_K^j
		for(int j = 0; j < this.docCorpora.corporaNumber; j++){
			double beta_split = this.randomSampler.nextBeta(this.gamma * this.nu_u, this.gamma * nu_K);
			double beta_K = this.beta_us.get(j) * (1 - beta_split);
			double beta_u = this.beta_us.get(j) * beta_split;
			this.beta_us.set(j, beta_u);
			this.betas.get(j).put(id, beta_K);
		}
		this.topicTermFreq.put(id, new OpenIntIntHashMap());
		return id;
	}
	public void RemoveTopic(int topicID){
		this.nu.removeKey(topicID);
		for(int j = 0; j < this.docCorpora.corporaNumber; j++){
			this.betas.get(j).removeKey(topicID);
		}
		for(int j = 0; j < this.docCorpora.corporaNumber; j++){
			int docnum = this.docCorpora.docNumbers.get(j);
			for(int d = 0; d < docnum; d++){
				this.docTopicSize.get(j).get(d).removeKey(topicID);
			}
		}
		
		this.topicTermFreq.removeKey(topicID);
		this.topicTotalTermFreq.removeKey(topicID);
		
		this.metaTableNumbers.removeKey(topicID);
		for(int j = 0; j < this.docCorpora.corporaNumber; j++){
			this.tableNumbers.get(j).removeKey(topicID);
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
		List<List<int[]>> initlabel = new ArrayList<List<int[]>>();
		for(int j = 0; j < J; j++){
			int docnum = this.docCorpora.docNumbers.get(j);
			initlabel.add(new ArrayList<int[]>());
			for(int d = 0; d < docnum; d++){
				initlabel.get(j).add(this.randInitLabel(this.docCorpora.docLengths.get(j).get(d), alphabet));
			}
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
	public void Initialize(List<List<int[]>> initZ){
		OpenIntIntHashMap labelmap = new OpenIntIntHashMap();
		int J = initZ.size();
		int ilabel;
		int topicid = -1;
		int wjdi = -1;
		int njdk, nkw, nk;
		OpenIntIntHashMap topicDist = null;
		for(int j = 0; j < J; j++){
			for(int d = 0; d < initZ.get(j).size(); d++){
				for(int i = 0; i < initZ.get(j).get(d).length; i++){
					ilabel = initZ.get(j).get(d)[i];
					wjdi = this.docCorpora.ftSequences.get(j).get(d)[i];
					if(labelmap.containsKey(ilabel) == false){
						topicid = AddNewTopic();
						labelmap.put(ilabel, topicid);
					}
					else{
						topicid = labelmap.get(ilabel);///getComponent(labelmap.get(ilabel));
					}
					
					//assign the word to the topic
					this.Z.get(j).get(d)[i] = topicid;
					//update count variables
					//update the count n_{jk} of words assigned to a topic k in document j
					njdk = this.docTopicSize.get(j).get(d).get(topicid);
					this.docTopicSize.get(j).get(d).put(topicid, njdk + 1);
					//update the count n_{kw} of word w assigned to topic k in all documents
					topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(topicid);
					nkw = topicDist.get(wjdi);
					topicDist.put(wjdi, nkw + 1);
					//update the summation variable of n_k = \sum_w n_{kw}
					nk = this.topicTotalTermFreq.get(topicid);
					this.topicTotalTermFreq.put(topicid, nk + 1);
				}
			}
		}
	}
	
	public int getTopicNumber(){
		return this.topicTermFreq.size();
	}
	
	public void SamplingZ(boolean shuffle){
		int J = this.docCorpora.corporaNumber;
		int[] rpmj = shuffle? this.randomSampler.RandPerm(J) : null;
		
		int d, w, L, z, znew, K;
		int njdk, nkw, nk;
		int kid;
		double beta_k, ss;
		OpenIntIntHashMap topicDist = null;
		OpenIntIntHashMap docTopicSize = null;
		for(int rj = 0; rj < J; rj++){
			int j = shuffle? rpmj[rj] : rj;
			int docnum = this.docCorpora.docNumbers.get(j);
			int[] rpmd = shuffle? this.randomSampler.RandPerm(docnum) : null;
			for(int rd = 0; rd < docnum; rd++){
				d = shuffle? rpmd[rd] : rd;
				L = this.docCorpora.docLengths.get(j).get(d);
				for(int i = 0; i < L; i++){
					w = this.docCorpora.ftSequences.get(j).get(d)[i];
					z = this.Z.get(j).get(d)[i];
					
					//sampling a new topic
					docTopicSize = this.docTopicSize.get(j).get(d);
					
					IntArrayList topicIDs = this.topicTermFreq.keys();
					K = topicIDs.size();
					double[] post = new double[K + 1];
					for(int k = 0; k < K; k++){
						kid = topicIDs.get(k);
						topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(kid);
						beta_k = this.betas.get(j).get(kid);
						if(kid != z){
							njdk = docTopicSize.get(kid);
							nkw = topicDist.get(w);
							nk = this.topicTotalTermFreq.get(kid);
						}
						else{
							njdk = docTopicSize.get(kid) - 1;
							nkw = topicDist.get(w) - 1;
							nk = this.topicTotalTermFreq.get(kid) - 1;
						}
						//nk = 0 means the topic k should be removed
						post[k] = (nk == 0)? 
								0 : ((njdk + this.alphas.get(j) * beta_k) * (nkw + b0) / (nk + this.b0TimesW));
					}
					post[K] = this.alphas.get(j) * this.beta_us.get(j) * this.newTopicLike;
					ss = VectorOper.VectorSum(post);
					VectorOper.VectorTimesEqual(post, 1/ss);
					znew = this.randomSampler.nextDiscrete(post);
					znew = (znew < K)? topicIDs.get(znew) : this.AddNewTopic();
					//if znew == z, nothing need to update
					if(znew != z){
						this.Z.get(j).get(d)[i] = znew;
						//remove i from original topic
						njdk = docTopicSize.get(z);
						if(njdk == 1){
							docTopicSize.removeKey(z);
						}
						else{
							docTopicSize.put(z, njdk - 1);
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
						njdk = docTopicSize.get(znew);
						docTopicSize.put(znew, njdk + 1);
						
						topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(znew);
						nkw = topicDist.get(w);
						topicDist.put(w, nkw + 1);
						
						nk = this.topicTotalTermFreq.get(znew);
						this.topicTotalTermFreq.put(znew, nk + 1);
						
					}//if
					
				}//for i
			}//for d
		}//for j
	}
	
	/**
	 * Sampling the number of sub-clusters for each topic k in each document j
	 */
	protected void SamplingTableNumbersAndBeta(){
		int J = this.docCorpora.corporaNumber;
		int njdk = 0, kid = -1, mk = 0;
		IntArrayList topicIDs = this.topicTermFreq.keys();
		int K = topicIDs.size();
		double[] ms = new double[K + 1];
		
		for(int j = 0; j < J; j++){
			int docnum = this.docCorpora.docNumbers.get(j);
			this.tableNumbers.get(j).clear();
			ms[K] = this.gamma * this.nu_u;
			double bb = 0;
			for(int k = 0; k < K; k++){
				kid = topicIDs.get(k);
				mk = 0;
				bb = this.alphas.get(j) * this.betas.get(j).get(kid);
				for(int d = 0; d < docnum; d++){
					njdk = this.docTopicSize.get(j).get(d).get(kid);
					//simulate a CRP to sample the table number assigned with 
					//a same component in each restaurant
					if(njdk > 0){
						mk += RandomSampler.CRPTableNumber(njdk, bb);
					}
				}
				ms[k] = this.gamma * this.nu.get(kid) + mk;
				
				this.tableNumbers.get(j).put(kid, mk);
			}
	//		Smooth(ms);
			double[] beta = this.randomSampler.nextDirichlet(ms);
	//		Smooth(beta);
			for(int k = 0; k < K; k++){
				kid = topicIDs.get(k);
				this.betas.get(j).put(kid, beta[k]);
			}
			this.beta_us.set(j, beta[K]);
		}
	}
	private void SamplingMetaTableNumbersAndNu(){
		int J = this.docCorpora.corporaNumber;
		int tjk = 0, kid = -1, mk = 0;
		this.metaTableNumbers.clear();
		
		IntArrayList topicIDs = this.topicTermFreq.keys();
		int K = topicIDs.size();
		double[] ms = new double[K + 1];
		ms[K] = this.xi;
		double bb = 0;
		for(int k = 0; k < K; k++){
			kid = topicIDs.get(k);
			mk = 0;
			bb = this.gamma * this.nu.get(kid);
			for(int j = 0; j < J; j++){
				tjk = this.tableNumbers.get(j).get(kid);
				//simulate a CRP to sample the table number assigned with 
				//a same component in each restaurant
				if(tjk > 0){
					mk += RandomSampler.CRPTableNumber(tjk, bb);
				}
			}
			ms[k] = mk;
			if(mk == 0)
				System.err.println("SamplingMandBeta: mk = 0 !");
			this.metaTableNumbers.put(kid, mk);
		}
//		Smooth(ms);
		double[] nu = this.randomSampler.nextDirichlet(ms);
//		Smooth(beta);
		for(int k = 0; k < K; k++){
			kid = topicIDs.get(k);
			this.nu.put(kid, nu[k]);
		}
		this.nu_u = nu[K];
	}
	
	/**
	 * Sampling the concentration parameter \xi, using the auxiliary variable method of 
	 * 		Escobar, M. & West, M. 
	 * 		Bayesian Density Estimation and Inference Using Mixtures. 
	 * 		JASA, 1995, 90, 577-588
	 * @param iter, iteration number
	 */
	protected void SamplingXi(int iter){
		double ga, eta = 0;
		double m = VectorOper.VectorSum(this.metaTableNumbers);
		
		int K = this.getTopicNumber();
		double kai = 0, beta = 0;//parameter for gamma distribution
		ga = this.xi;//_a / m_gamma_b;
		kai = this.xi_a + K - 1;

		for(int i = 0; i < iter; i++){
			eta = this.randomSampler.nextBeta(ga + 1, m);
			beta = this.xi_b - Math.log(eta);
			if(Math.random() < (kai / (kai + m * beta)))
				kai += 1;
			ga = this.randomSampler.nextGamma(kai, 1/beta);
			//notice that here 'beta' is the inverse scale parameter,
			//while in Mallet's Randoms class, the nextGamma(alpha, beta) , 
			//'beta' there is the scale parameter
		}
		this.xi = ga;
	}
	/**
	 * Sampling the concentration parameter \gamma, using the auxiliary variable method (for \alpha there) of
	 * 		Teh, Y.; Jordan, M.; Beal, M. & Blei, D. 
	 * 		Hierarchical Dirichlet Processes. 
	 * 		JASA, 2006, 101, 1566-1581
	 * @param iter
	 */
	protected void SamplingGamma(int iter){
		int J = this.docCorpora.corporaNumber;
		
		//tj[j] is the table numbers in all documents in corpus j, i.e., number of
		//samples drawn from G^j to compose the tables of all documents in corpus j
		int[] tj = new int[J];
		for(int j = 0; j < J; j++){
			tj[j] = VectorOper.VectorSum(this.tableNumbers.get(j));
		}
		
		double m = VectorOper.VectorSum(this.metaTableNumbers);
		double w, s;
		
		double t = 0;
		double a = this.gamma_a + m, b = this.gamma_b;//parameters for gamma distribution
		double al = this.gamma;
		for(int i = 0; i < iter; i++){
			a = this.gamma_a + m;
			b = this.gamma_b;
			//sampling wj and sj
			for(int j = 0; j < J; j++){
				w = this.randomSampler.nextBeta(al + 1, tj[j]);
				t = tj[j] / al;
				s = (Math.random() < t / (t + 1))? 1:0;
				a -= s;
				b -= Math.log(w);
			}
			//sampling gamma
			al = this.randomSampler.nextGamma(a, 1 / b);
		}
		this.gamma = al;
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
		for(int j = 0; j < J; j++){
		
			int D = this.docCorpora.docNumbers.get(j);
			List<Integer> nd = this.docCorpora.docLengths.get(j);
			double m = VectorOper.VectorSum(this.tableNumbers.get(j));
			double w, s;
			
			double t = 0;
			double a = this.alphas_a.get(j) + m, b = this.alphas_b.get(j);//parameters for gamma distribution
			double al = this.alphas.get(j);
			for(int i = 0; i < iter; i++){
				a = this.alphas_a.get(j) + m;
				b = this.alphas_b.get(j);
				//sampling wj and sj
				for(int d = 0; d < D; d++){
					w = this.randomSampler.nextBeta(al + 1, nd.get(d));
					t = nd.get(d) / al;
					s = (Math.random() < t / (t + 1))? 1:0;
					a -= s;
					b -= Math.log(w);
				}
				//sampling alpha
				al = this.randomSampler.nextGamma(a, 1 / b);
			}
			this.alphas.set(j, al);
		
		}
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
			SamplingTableNumbersAndBeta();
//			System.out.println("	Sampling M and Beta: " + (System.currentTimeMillis() - starttime) + "ms");
		
			SamplingMetaTableNumbersAndNu();

			if(this.bSamplingHyper){
//				starttime = System.currentTimeMillis();
				SamplingAlpha(conpiter);
//				System.out.println("	Sampling alpha: " + (System.currentTimeMillis() - starttime) + "ms");
//				
//				starttime = System.currentTimeMillis();
				SamplingGamma(conpiter);
//				System.out.println("	Sampling gamma: " + (System.currentTimeMillis() - starttime) + "ms");
				
				SamplingXi(conpiter);
			}
			
			usedtime = System.currentTimeMillis() - time0;
			Log(iter, usedtime);			
		}
	}
	
	protected void Log(int iter, long usedTime){
		int K = -1;
		int[][] Ks = null;
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
				log += "xi=" + this.xi +";gamma=" + this.gamma + "; alpha = " 
				 + VectorOper.ToString(this.alphas, "%.2f");
			break;
			
		case 2:// iter + K + time
			log += "Iteration=" + iter + "; Time="+ usedTime + "ms;" 
			+ " K=" + K +";";
			if(this.bSamplingHyper)
				log += "xi=" + this.xi +";gamma=" + this.gamma + "; alpha = " 
				 + VectorOper.ToString(this.alphas, "%.2f");
			break;
			
		case 3://iter + K + like + time
			log += "Iteration=" + iter + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; loglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += "\nxi=" + this.xi +";gamma=" + this.gamma + "; alpha = " 
				 + VectorOper.ToString(this.alphas, "%.2f");
			break;
			
		case 4://iter + K + Ks + like + time
			log += "Iteration=" + iter  + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; \nKs=" + VectorOper.ToString(Ks) + 
			"; \nloglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += "xi=" + this.xi +";gamma=" + this.gamma + "; alpha = " 
				 + VectorOper.ToString(this.alphas, "%.2f");
			break;
			
		case 5://iter + K + Ks + like + time + nu + beta + topicsize
			log += "Iteration=" + iter  + 
			"; Time="+ usedTime + "ms;"
			+ " K=" + K + "; \nKs=" + VectorOper.ToString(Ks) + 
			"; \nloglike=" + loglike + ";";
			if(this.bSamplingHyper)
				log += "xi=" + this.xi +";gamma=" + this.gamma + "; alpha = " 
				 + VectorOper.ToString(this.alphas, "%.2f");
			
			double[] nus;
			double[][] betas;
			int[] kids, sizes;
			//global info
			nus = new double[K + 1];
			kids = new int[K + 1];
			sizes = new int[K + 1];
			int[] metaTableNum = new int[K + 1];
			metaTableNum[0] = 0;
			betas = new double[this.docCorpora.corporaNumber][K + 1];
			for(int j = 0; j < this.docCorpora.corporaNumber; j++){
				betas[j][0] = this.beta_us.get(j);
			}
			nus[0] = this.nu_u;
			kids[0] = MultiCoporaHDPTopic.NEW_TOPIC_ID;
			sizes[0] = 0;
			
			IntArrayList topicIDs = this.topicTermFreq.keys();
			for(int k = 1; k < K + 1; k++){
				kids[k] = topicIDs.get(k - 1);
				sizes[k] = this.topicTotalTermFreq.get(kids[k]);
				nus[k] = this.nu.get(kids[k]);
				metaTableNum[k] = this.metaTableNumbers.get(kids[k]);
				for(int j = 0; j < this.docCorpora.corporaNumber; j++){
					betas[j][k] = this.betas.get(j).get(kids[k]);
				}
				
			}
			log += "\nID		|" + VectorOper.ToString(kids, "", "	|", "%6d") + "\n";
			log += "nu(1E-2)	|" + VectorOper.ToString(VectorOper.VectorTimes(nus, 100)
					, "", "	|", "%6.3f") + "\n";
			for(int j = 0; j < this.docCorpora.corporaNumber; j++){
				log += "beta(1E-2)" + j + "	|" + VectorOper.ToString(VectorOper.VectorTimes(betas[j], 100)
						, "", "	|", "%6.3f") + "\n";
			}
			log += "sizes		|" + VectorOper.ToString(sizes, "", "	|", "%6d") + "\n";
			log += "metaT		|" + VectorOper.ToString(metaTableNum, "", "	|", "%6d") + "\n";
			log += "totalsize = " + VectorOper.VectorSum(sizes);
			
			
			CorporaTopicPostInfo topicInfo = new CorporaTopicPostInfo(
					this.docCorpora.dictionary
					, this.docCorpora.forwardHash
					, this.docCorpora.backwardHash
					, this.docCorpora.ftSequences
					, this.Z
					, this.docTopicSize
					, this.topicTermFreq
					, this.topicTotalTermFreq);
			for(int j = 0; j < this.docCorpora.corporaNumber; j++){
				
				IntArrayList jTopicIds = topicInfo.getCorpusTopicTermFreqs(j).keys();
				int Kj = jTopicIds.size();
				int[] jkids = new int[Kj + 1];
				jkids[0] = MultiCoporaHDPTopic.NEW_TOPIC_ID;
				int[] tableNum = new int[Kj + 1];
				tableNum[0] = 0;
				sizes = new int[Kj + 1];
				sizes[0] = 0;
				for(int ki = 1; ki < Kj + 1; ki++){
					jkids[ki] = jTopicIds.get(ki - 1);
					sizes[ki] = topicInfo.getCorpusTopicSize(j, jkids[ki]);
					tableNum[ki] = this.tableNumbers.get(j).get(jkids[ki]);
				}
				log += "\n	corpus " + j + "-------------------------------------";
				log += "\n	ID		|" + VectorOper.ToString(jkids, "", "	|", "%6d") + "\n";
				log += "	sizes		|" + VectorOper.ToString(sizes, "", "	|", "%6d") + "\n";
				log += "	TNum		|" + VectorOper.ToString(tableNum, "", "	|", "%6d") + "\n";
				log += "	totalsize = " + VectorOper.VectorSum(sizes);
			}
			
			break;
			
			
		default:
			break;
		}
		if(this.logLevel > 0)
			System.out.println(log);
	}
	
	int[][] getDocTopicNums(){
		int J = this.docCorpora.corporaNumber;
		int[][] ret = new int[J][];
		for(int j = 0; j < J; j++){
			ret[j] = new int[this.docCorpora.docNumbers.get(j)];
			for(int d = 0; d < this.docCorpora.docNumbers.get(j); d++){
				ret[j][d] = this.docTopicSize.get(j).get(d).size();
			}
		}
		return ret;
	}
	
	int[][] getDocTopicNums(List<List<Integer>> docs){
		int[][] ret = new int[docs.size()][];
		for(int j = 0; j < docs.size(); j++){
			ret[j] = new int[docs.get(j).size()];
				for(int d = 0; d < docs.get(j).size(); d++){
					ret[j][d] = this.docTopicSize.get(j).get(docs.get(j).get(d)).size();
			}
		}
		return ret;
	}
	
	/**
	 * not strictly the likelihood, but just cumulation of the likelihood of each word
	 * from its topic,  for convergence diagnosis
	 * @return
	 */
	private double CalLogLike(List<List<Integer>> docs){
		double ret = 0;
		int J = docs.size();
		int kid = 0, nkw, nk;
		int d,w;
		for(int j = 0; j < J; j++){
			for(int rd = 0; rd < docs.get(j).size(); rd++){
				d = docs.get(j).get(rd);
				for(int i = 0; i < this.docCorpora.ftSequences.get(j).get(d).length; i++){
					w = this.docCorpora.ftSequences.get(j).get(d)[i];
					kid = this.Z.get(j).get(d)[i];
					nkw = ((OpenIntIntHashMap) this.topicTermFreq.get(kid)).get(w);
					nk = this.topicTotalTermFreq.get(w);
					ret += Math.log((nkw + this.b0) / (nk + this.b0TimesW));
				}
			}
		}
		return ret;
	}
	
//	public OpenIntObjectHashMap getTopicKeywords(int topN) {
//		// TODO Auto-generated method stub
//		if(this.docCorpora.dictionary == null){
//			System.err.println("getTopicKeywords failed, as no dictionary provided!");
//			return null;
//		}
//		OpenIntObjectHashMap ret = new OpenIntObjectHashMap();
//		IntArrayList ks = this.topicTermFreq.keys();
//		OpenIntIntHashMap topicDist = null;
//		for(int ki = 0; ki < ks.size(); ki++){
//			int kid = ks.get(ki);
//			topicDist = (OpenIntIntHashMap) this.topicTermFreq.get(kid);
//			IntArrayList wids = topicDist.keys();
//			IntArrayList wfrqs = topicDist.values();
//			int[] iwfrqs = new int[wfrqs.size()];
//			for(int i = 0; i < iwfrqs.length; i++){
//				iwfrqs[i] = wfrqs.get(i);
//			}
//			int[] maxN = VectorOper.MaxK(iwfrqs, topN);
//			int realTopN = maxN.length;
//			String[] words = new String[realTopN];
//			int[] freqs = new int[realTopN];
//			for(int i = 0; i < realTopN; i++){
//				words[i] = (String) this.docCorpora.dictionary.lookupObject(wids.get(maxN[i]));
//				freqs[i] = iwfrqs[maxN[i]];
//			}
//			ret.put(kid, new Keywords(words, freqs));
//		}
//		return ret;
//	} 
//	
//	public String PrintTopicKeywords(){
//		String msg = "ID	|	keyword(frequency)\n";
//		OpenIntObjectHashMap topicKeywords = this.getTopicKeywords(40);
//		if(topicKeywords == null)
//			return msg;
//		IntArrayList topicIDs = topicKeywords.keys();
//		for(int k = 0; k < topicIDs.size(); k++){
//			int kid = topicIDs.get(k);
//			msg += kid + "	|";
//			Keywords keywds = (Keywords) topicKeywords.get(kid);
//			String[] words = keywds.getWords();
//			int[] freqs = keywds.getFreqs();
//			for(int i = 0; i < words.length; i++){
//				msg += words[i] + "(" + freqs[i] + "),";
//			}
//			msg += "\n";
//		}
//		return msg;
//	}
	
	/**
	 * doc-topic count matrix. 
	 * @return ret.
	 * ret[j][d][k] is number of words assigned to topic k in document d of 
	 * corpus j
	 */
	public List<List<OpenIntIntHashMap>> getDocTopicMatrix(){
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
	
	// test
	public static void main(String[] args){

//		double b0 = 0.1;
//		double xi = 5.0, gamma = 5.0, alpha = 1.0;
//		double xi_a = 5.0, xi_b = 1.0;
//		double gamma_a = 5.0, gamma_b = 1.0;
//		double alpha_a = 5.0, alpha_b = 1.0;
//		
//		int maxIter = 20, burnin = 10;
//		int initK = 5;
//		int loglevel = 3;
//		MultiCoporaHDPTopic multHdpTopic = new MultiCoporaHDPTopic(
//				docCorpora
//				, b0
////				, xi, gamma, alpha //if sampling hyper, comment this line and uncomment following 3 lines, vice versa
//				, xi_a, xi_b
//				, gamma_a, gamma_b
//				, alpha_a, alpha_b
//				, maxIter
//				, burnin);
//		multHdpTopic.Initialize(initK);
//		multHdpTopic.setLogLevel(loglevel);
//		multHdpTopic.Sampling();
//		CorporaTopicPostInfo topicInfo = new CorporaTopicPostInfo(
//				docCorpora.dictionary
//				, docCorpora.forwardHash
//				, docCorpora.backwardHash
//				, docCorpora.ftSequences
//				, multHdpTopic.Z
//				, null
//				, null
//				, null);
//		int topN = 40;
//		String msg = topicInfo.PrintTopicsKeywords(topN);
//		for(int j = 0; j < docCorpora.corporaNumber; j++){
//			msg += "\n corpus " + j + ": " + docCorpora.getTypeStr(j) + "---------------------------\n";
//			msg += topicInfo.getLocalInfo(j).PrintTopicsKeywords(topN);
//		}
//		System.out.println(msg);
	}
	
}