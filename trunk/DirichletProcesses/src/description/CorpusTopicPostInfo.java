package description;

import java.util.ArrayList;
import java.util.List;

import util.VectorOper;
import cc.mallet.types.Alphabet;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.Keywords;

public class CorpusTopicPostInfo{
	Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;

	List<int[]> ftSequences;
	List<int[]> topicAssignments;
	List<OpenIntIntHashMap> docTopicSizes;
	OpenIntObjectHashMap topicTermFreqs;
	OpenIntIntHashMap topicTotalTermFreqs;
	
	public CorpusTopicPostInfo(
			Alphabet dict
			, OpenIntIntHashMap forwardHash
			, OpenIntIntHashMap backwardHash
			, List<int[]> ftSequences
			, List<int[]> topicAssignments
			, List<OpenIntIntHashMap> docTopicSizes
			, OpenIntObjectHashMap topicTermFreqs
			, OpenIntIntHashMap topicTotalTermFreqs){
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;
		
		this.ftSequences = ftSequences;
		this.topicAssignments = topicAssignments;
		this.docTopicSizes = docTopicSizes;
		this.topicTermFreqs = topicTermFreqs;
		this.topicTotalTermFreqs = topicTotalTermFreqs;
		
		if((docTopicSizes == null) || (topicTermFreqs == null) || (topicTotalTermFreqs == null))
			this.Update();
	}
	
	public void Update(){
		int D = this.topicAssignments.size();
		this.docTopicSizes = new ArrayList<OpenIntIntHashMap>();
		this.topicTermFreqs = new OpenIntObjectHashMap();
		this.topicTotalTermFreqs = new OpenIntIntHashMap();
		for(int d = 0; d < D; d++){
			this.docTopicSizes.add(new OpenIntIntHashMap());
			int nd = this.topicAssignments.get(d).length;
			for(int w = 0; w < nd; w++){
				int z = this.topicAssignments.get(d)[w];
				int ndz = this.docTopicSizes.get(d).get(z);
				this.docTopicSizes.get(d).put(z, ndz+1);
				OpenIntIntHashMap topicDist = null;
				if(this.topicTermFreqs.containsKey(z))
					topicDist = (OpenIntIntHashMap) this.topicTermFreqs.get(z);
				else{
					topicDist = new OpenIntIntHashMap();
					this.topicTermFreqs.put(z, topicDist);
				}
				int wi = this.ftSequences.get(d)[w];
				int nzw = topicDist.get(wi);
				topicDist.put(wi, nzw + 1);
				int nz = this.topicTotalTermFreqs.get(z);
				this.topicTotalTermFreqs.put(z, nz + 1);
			}
		}
		
	}
	
	public OpenIntObjectHashMap getTopicsKeywords(int topN) {
		// TODO Auto-generated method stub
		if(this.dictionary == null){
			System.err.println("getTopicsKeywords failed, as no dictionary provided!");
			return null;
		}
		OpenIntObjectHashMap ret = new OpenIntObjectHashMap();
		IntArrayList ks = this.topicTermFreqs.keys();
		OpenIntIntHashMap topicDist = null;
		for(int ki = 0; ki < ks.size(); ki++){
			int kid = ks.get(ki);
			topicDist = (OpenIntIntHashMap) this.topicTermFreqs.get(kid);
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
		String msg = "ID	|	keyword(frequency)\n";
		OpenIntObjectHashMap topicKeywords = this.getTopicsKeywords(topN);
		if(topicKeywords == null)
			return msg;
		IntArrayList topicIDs = topicKeywords.keys();
		for(int k = 0; k < topicIDs.size(); k++){
			int kid = topicIDs.get(k);
			msg += kid + "	|";
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