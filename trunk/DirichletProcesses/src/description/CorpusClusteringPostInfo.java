package description;


import java.util.List;

import util.VectorOper;
import cc.mallet.types.Alphabet;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.Keywords;

public class CorpusClusteringPostInfo{
	Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;

	List<OpenIntIntHashMap> tfVectors;
	List<Integer> clusterAssignments;
	List<Integer> labels;
	OpenIntIntHashMap clusterSize;//number of documents in a cluster
	OpenIntObjectHashMap clusterTermFreqs;
	OpenIntIntHashMap clusterTotalTermFreqs;
	
	OpenIntObjectHashMap clusterMembers;//
	
	public CorpusClusteringPostInfo(
			Alphabet dict
			, OpenIntIntHashMap forwardHash
			, OpenIntIntHashMap backwardHash
			, List<OpenIntIntHashMap> tfs
			, List<Integer>  clusterAssignments
			, List<Integer>  labels
			, OpenIntIntHashMap clusterSize
			, OpenIntObjectHashMap clusterTermFreqs
			, OpenIntIntHashMap clusterTotalTermFreqs){
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;

		this.labels = labels;
		this.tfVectors = tfs;
		this.clusterAssignments = clusterAssignments;
		this.clusterTermFreqs = clusterTermFreqs;
		this.clusterTotalTermFreqs = clusterTotalTermFreqs;
		this.clusterSize = clusterSize;
		
		this.UpdateMember();
		
		if((clusterTermFreqs == null) || (clusterTotalTermFreqs == null) || (clusterSize == null))
			this.Update();
	}
	
	private void UpdateMember(){
		this.clusterMembers = new OpenIntObjectHashMap();
		IntArrayList member = null;
		for(int d = 0; d < this.clusterAssignments.size(); d++){
			int z = this.clusterAssignments.get(d);
			if(this.clusterMembers.containsKey(z)){
				member = (IntArrayList) this.clusterMembers.get(z);
				member.add(d);
			}
			else{
				member = new IntArrayList();
				member.add(d);
				this.clusterMembers.put(z, member);
			}
				
		}
	}
	
	public void Update(){
		int D = this.clusterAssignments.size();
	
		this.clusterTermFreqs = new OpenIntObjectHashMap();
		this.clusterTotalTermFreqs = new OpenIntIntHashMap();
		this.clusterSize = new OpenIntIntHashMap();
		OpenIntIntHashMap clusterTf = null;
		for(int d = 0; d < D; d++){
			int z = this.clusterAssignments.get(d);
			if(this.clusterTermFreqs.containsKey(z))
				clusterTf = (OpenIntIntHashMap) this.clusterTermFreqs.get(z);
			else{
				clusterTf = new OpenIntIntHashMap();
				this.clusterTermFreqs.put(z, clusterTf);
			}
			VectorOper.VectorAddEqual(clusterTf, this.tfVectors.get(d));

			int nz = this.clusterSize.get(z);
			this.clusterSize.put(z, nz + 1);
			
		}
		IntArrayList kids = this.clusterTermFreqs.keys();
		for(int k = 0; k < kids.size(); k++){
			int kid = kids.get(k);
			clusterTf = (OpenIntIntHashMap) this.clusterTermFreqs.get(kid);
			int tfsum = VectorOper.VectorSum(clusterTf);
			this.clusterTotalTermFreqs.put(kid, tfsum);
		}		
	}
	
	
	public OpenIntObjectHashMap getTopicsKeywords(int topN) {
		// TODO Auto-generated method stub
		if(this.dictionary == null){
			System.err.println("getTopicsKeywords failed, as no dictionary provided!");
			return null;
		}
		OpenIntObjectHashMap ret = new OpenIntObjectHashMap();
		IntArrayList ks = this.clusterSize.keys();
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
			msg += kid + "	|" + this.clusterSize.get(kid) + "	|";
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
}