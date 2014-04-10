package description;

import java.util.ArrayList;
import java.util.List;

import util.VectorOper;
import cc.mallet.types.Alphabet;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.Keywords;

/**
 * @author Jianwen Zhang
 */
public class CorporaClusteringPostInfo{
	Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;

	List<CorpusClusteringPostInfo> corpusInfo;
	public List<List<OpenIntIntHashMap>> tfVectors;
	List<List<Integer>> clusterAssignments;
	List<List<Integer>> labels;
	OpenIntObjectHashMap clusterTermFreqs;
	OpenIntIntHashMap clusterTotalTermFreqs;

	OpenIntIntHashMap globalClusterSize;
	
	public CorpusClusteringPostInfo getLocalInfo(int j){
		return this.corpusInfo.get(j);
	}
	/**
	 * 
	 * @param dict
	 * @param tfs
	 * @param labels
	 * @param clusterAssignments
	 * @param globalClusterSize
	 * @param localClusterSizes
	 * @param clusterTermFreqs
	 * @param clusterTotalTermFreqs
	 */
	public CorporaClusteringPostInfo(
			Alphabet dict
			, OpenIntIntHashMap forwardHash
			, OpenIntIntHashMap backwardHash
			, List<List<OpenIntIntHashMap>> tfs
			, List<List<Integer>> labels
			, List<List<Integer>> clusterAssignments
			, OpenIntIntHashMap globalClusterSize
			, List<OpenIntIntHashMap> localClusterSizes
			, OpenIntObjectHashMap clusterTermFreqs
			, OpenIntIntHashMap clusterTotalTermFreqs
			){
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;

		this.tfVectors = tfs;
		this.labels = labels;
		this.clusterAssignments = clusterAssignments;
		this.globalClusterSize = globalClusterSize;
		
		int J = tfs.size();
		this.corpusInfo = new ArrayList<CorpusClusteringPostInfo>();
		for(int j = 0; j < J; j++){
			List<Integer> labelj = (this.labels == null)? null : this.labels.get(j);
			OpenIntIntHashMap localClusterSize = (localClusterSizes == null)? null : localClusterSizes.get(j);
			this.corpusInfo.add( new CorpusClusteringPostInfo(
					this.dictionary, this.forwardHash, this.backwardHash,
					this.tfVectors.get(j), this.clusterAssignments.get(j), labelj, 					
					localClusterSize, null, null) );
		}
		if((this.globalClusterSize == null) || (this.clusterTermFreqs == null) || (this.clusterTotalTermFreqs == null)){
			this.UpdateGlobal();
		}
	}

	public void UpdateGlobal(){
		this.clusterTermFreqs = new OpenIntObjectHashMap();
		this.clusterTotalTermFreqs = new OpenIntIntHashMap();
		this.globalClusterSize = new OpenIntIntHashMap();
		int J = this.tfVectors.size();
		for(int j = 0; j < J; j++){
			OpenIntObjectHashMap clusterTfs = this.corpusInfo.get(j).clusterTermFreqs;
			IntArrayList kids = clusterTfs.keys();
			for(int k = 0; k < kids.size(); k++){
				int ki = kids.get(k);
				OpenIntIntHashMap kTf = (OpenIntIntHashMap) clusterTfs.get(ki);
				if(this.clusterTermFreqs.containsKey(ki)){
					VectorOper.VectorAddEqual((OpenIntIntHashMap)this.clusterTermFreqs.get(ki), kTf);
				}
				else
					this.clusterTermFreqs.put(ki, kTf.clone());
			}
			VectorOper.VectorAddEqual(this.globalClusterSize, this.corpusInfo.get(j).clusterSize);
			VectorOper.VectorAddEqual(this.clusterTotalTermFreqs, this.corpusInfo.get(j).clusterTotalTermFreqs);
		}
	}
	
	public OpenIntObjectHashMap getTopicsKeywords(int topN) {
		// TODO Auto-generated method stub
		if(this.dictionary == null){
			System.err.println("getTopicsKeywords failed, as no dictionary provided!");
			return null;
		}
		OpenIntObjectHashMap ret = new OpenIntObjectHashMap();
		IntArrayList ks = this.globalClusterSize.keys();
		OpenIntIntHashMap topicDist = null;
		for(int ki = 0; ki < ks.size(); ki++){
			int kid = ks.get(ki);
			topicDist = (OpenIntIntHashMap) this.clusterTermFreqs.get(kid);
			IntArrayList wids = topicDist.keys();
			IntArrayList wfrqs = topicDist.values();
			int[] iwfrqs = new int[wfrqs.size()];
			for(int i = 0; i < iwfrqs.length; i++){
				iwfrqs[i] = wfrqs.get(i);
			}
			int[] maxN = VectorOper.MaxK(iwfrqs, topN);
			int realTopN = maxN.length;
			String[] words = new String[realTopN];
			int[] freqs = new int[realTopN];
			for(int i = 0; i < realTopN; i++){
				int wordID = wids.get(maxN[i]);
				if (this.backwardHash != null) {
					wordID = this.backwardHash.get(wordID);
				}
				words[i] = (String) this.dictionary.lookupObject(wordID);
				freqs[i] = iwfrqs[maxN[i]];
			}
			ret.put(kid, new Keywords(words, freqs));
		}
		return ret;
	} 
	
	public String PrintTopicsKeywords(int topN){
		String msg = "ID	|Size	|	keyword(frequency)\n";
		OpenIntObjectHashMap topicKeywords = this.getTopicsKeywords(topN);
		if(topicKeywords == null)
			return msg;
		IntArrayList topicIDs = topicKeywords.keys();
		for(int k = 0; k < topicIDs.size(); k++){
			int kid = topicIDs.get(k);
			msg += kid + "	|" + this.globalClusterSize.get(kid) + "	|";
			Keywords keywds = (Keywords) topicKeywords.get(kid);
			String[] words = keywds.getWords();
			int[] freqs = keywds.getFreqs();
			for(int i = 0; i < words.length; i++){
				msg += words[i] + "(" + freqs[i] + "),";
			}
			msg += "\n";
		}
		return msg;
	}
	
	public String PrintLocalTopicsKeywords(int j, int topN){
		String msg = this.corpusInfo.get(j).PrintTopicsKeywords(topN);
		return msg;
	}
}