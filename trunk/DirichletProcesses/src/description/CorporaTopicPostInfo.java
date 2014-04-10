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
public class CorporaTopicPostInfo{
	Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;

	List<List<int[]>> ftSequences;
	List<List<int[]>> topicAssignments;
	
	//global info
	List<List<OpenIntIntHashMap>> docTopicSizes;
	OpenIntObjectHashMap topicTermFreqs;
	OpenIntIntHashMap topicTotalTermFreqs;
	
	//local info
	CorpusTopicPostInfo[] corpusTopicInfo;
	
	public CorporaTopicPostInfo(
			Alphabet dict
			, OpenIntIntHashMap forwardHash
			, OpenIntIntHashMap backwardHash
			, List<List<int[]>> ftSequences
			, List<List<int[]>> topicAssignments
			, List<List<OpenIntIntHashMap>> docTopicSizes
			, OpenIntObjectHashMap topicTermFreqs
			, OpenIntIntHashMap topicTotalTermFreqs
			){
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;
		
		this.ftSequences = ftSequences;
		this.topicAssignments = topicAssignments;
		this.docTopicSizes = docTopicSizes;
		this.topicTermFreqs = topicTermFreqs;
		this.topicTotalTermFreqs = topicTotalTermFreqs;
		
		int corporaNumber = this.topicAssignments.size();
		this.corpusTopicInfo = new CorpusTopicPostInfo[corporaNumber];
		for(int j = 0; j < corporaNumber; j++){
			this.corpusTopicInfo[j] = new CorpusTopicPostInfo(
					this.dictionary, this.forwardHash, this.backwardHash,
					this.ftSequences.get(j), this.topicAssignments.get(j), null, null, null);
		}
		
		if((this.docTopicSizes == null) || (this.topicTermFreqs == null) || (this.topicTotalTermFreqs == null)){
			UpdateGlobal();
		}
		
	}
	public CorpusTopicPostInfo getLocalInfo(int j){
		return this.corpusTopicInfo[j];
	}
	public int getCorpusTopicSize(int j, int kid){
		return this.corpusTopicInfo[j].topicTotalTermFreqs.get(kid);
	}
	public OpenIntObjectHashMap getCorpusTopicTermFreqs(int j){
		return this.corpusTopicInfo[j].topicTermFreqs;
	}
	public void UpdateGlobal(){
		int J = this.topicAssignments.size();
		this.docTopicSizes = new ArrayList<List<OpenIntIntHashMap>>();
		this.topicTermFreqs = new OpenIntObjectHashMap();
		this.topicTotalTermFreqs = new OpenIntIntHashMap();
		
		for(int j = 0; j < J; j++){
			CorpusTopicPostInfo corpusInfo = this.corpusTopicInfo[j];
			this.docTopicSizes.add(corpusInfo.docTopicSizes);
			IntArrayList topics = corpusInfo.topicTermFreqs.keys();
			for(int ki = 0; ki < topics.size(); ki++){
				int k = topics.get(ki);
				OpenIntIntHashMap topicDist = (OpenIntIntHashMap) corpusInfo.topicTermFreqs.get(k);
				if(this.topicTermFreqs.containsKey(k))
					VectorOper.VectorAddEqual((OpenIntIntHashMap)(this.topicTermFreqs.get(k)), topicDist);
				else
					this.topicTermFreqs.put(k, topicDist.clone());
				
				VectorOper.VectorAddEqual(this.topicTotalTermFreqs, corpusInfo.topicTotalTermFreqs);
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
	
	public String PrintLocalTopicsKeywords(int j, int topN){
		String msg = this.corpusTopicInfo[j].PrintTopicsKeywords(topN);
		return msg;
	}

}