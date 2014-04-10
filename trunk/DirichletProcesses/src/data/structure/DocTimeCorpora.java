package data.structure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import util.VectorOper;
import cc.mallet.types.Alphabet;
import cern.colt.map.OpenIntIntHashMap;

/**
 * @author Jianwen Zhang
 */
public class DocTimeCorpora {

	public Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;
	public HashMap<Integer, String> corporaIndex2Name;
	public HashMap<String, Integer> corporaName2Index;
	public HashMap<String, Integer> corporaTime2Index;
    public HashMap<Integer, String> corporaIndex2Time;
    
	/**
	 * ftSequences[t][j][d][] is the feature sequence of time t document d of corpus j
	 */
	public List<List<List<int[]>>> ftSequences;
	/**
	 * tfVectors[t][j][d] is the term frequences of time t document d of corpus j
	 */
	public List<List<List<OpenIntIntHashMap>>> tfVectors;
	
	/**
	 * docCategories[t][j][d] is the category of document 
	 */
	public List<List<List<Integer>>> docLabels;
	/**
	 * labelDescription is a bi-mapping between a label integer value and its description 
	 */
	public Alphabet labelDescription;
	
	/**
	 * size of dictionary
	 */
	public int dictSize;
	/**
	 * number of corpora
	 */
	public int corporaNumber;
	public int timestampNumber;
	/**
	 * docNumbers[t][j] is the number of documents in corpus j
	 */
	public List<List<Integer>> docNumbers;
	/**
	 * total number of documents in all corpora.
	 */
	public int totalDocNumber;
	/**
	 * docLenghts[t][j][d] is the length of document d of corpus j
	 */
	public List<List<List<Integer>>> docLengths;
	/**
	 * corpusTokenNumbers[t][j] is the number of tokens in corpus j
	 */
	public List<List<Double>> averageDocLengths;
	public double averageDocLength;
	public List<List<Integer>> corpusTokenNumbers;
	/**
	 * The total number of tokens in all corpora.
	 */
	public int totalTokenNumber;
	
	public String corporaDescription;
	
	public DocTimeCorpora(){
		this.dictionary = null;
		this.ftSequences = null;
		this.tfVectors = null;
		this.docLabels = null;
		this.labelDescription = null;
		this.dictSize = 0;
		this.corporaNumber = 0;
		this.timestampNumber = 0;
		this.docLengths = null;
		this.averageDocLength = 0;
		this.averageDocLengths = null;
		this.corpusTokenNumbers = null;
		this.totalDocNumber = 0;
		this.totalTokenNumber = 0;
		this.corporaDescription = null;
		this.forwardHash = null;
		this.backwardHash = null;
		this.corporaIndex2Name = null;
		this.corporaName2Index = null;
		this.corporaTime2Index = null;
	    this.corporaIndex2Time = null;
	}
	
	
//	public String MsgCorporaInfo(){
//		String msg = "";
//		if(this.corporaDescription != null){
//			msg += "Corpora description: " + this.corporaDescription + "\n";
//		}
//		msg += "Corpora number			:" + this.corporaNumber + "\n";
//		msg += "Dictionary size			:" + this.dictSize + "\n";
//		msg += "Total document number		:" + this.totalDocNumber + "\n";
//		msg += "Total token number		:" + this.totalTokenNumber + "\n";
//		msg += "Average document length		:" + String.format("%.2f", this.averageDocLength) + "\n";
//		msg += "Info of each corpus:";
//		int[] corpusid = new int[this.corporaNumber];
//		for(int j = 0; j < this.corporaNumber; j++){
//			corpusid[j] = j;
//		}
//		msg += "\nCorpus			|" + VectorOper.ToString(corpusid, "", "|", "%10d");
//		msg += "\nDoc number		|" + VectorOper.ToString(this.docNumbers, "", "|", "%10d", false);
//		msg += "\nToken number		|" + VectorOper.ToString(this.corpusTokenNumbers, "", "|", "%10d", false);
//		msg += "\nAvg doc length		|" + VectorOper.ToString(this.averageDocLengths, "", "|", "%10.2f");
//		
//		return msg;
//		
//	}
	
	public void Create(Alphabet dict,
			OpenIntIntHashMap forwardHash,
			OpenIntIntHashMap backwardHash,
			HashMap<Integer, String> corporaIndex2Name,
			HashMap<String, Integer> corporaName2Index,
			HashMap<String, Integer> corporaTime2Index,
		    HashMap<Integer, String> corporaIndex2Time,
			List<List<List<int[]>>> ftSequences,
			List<List<List<OpenIntIntHashMap>>> tfVectors,
			List<List<List<Integer>>> docLabels,
			Alphabet labelDescription,
			String corporaDescription
			){
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;
		this.corporaIndex2Name = corporaIndex2Name;
		this.corporaName2Index = corporaName2Index;
		this.corporaTime2Index = corporaTime2Index;
		this.corporaIndex2Time = corporaIndex2Time;
		this.ftSequences = ftSequences;
		this.tfVectors = tfVectors;
		this.docLabels = docLabels;
		this.labelDescription = labelDescription;
		this.corporaDescription = corporaDescription;
		
		if(dict == null)
			System.err.println("DocCorpora: there is no dictionary!");
		if((ftSequences == null) && (tfVectors == null))
			System.err.println("DocCorpora: there is not any feature provided!");
		
		this.dictSize = (this.forwardHash == null) ? this.dictionary.size() : this.forwardHash.size();
		this.timestampNumber = (this.ftSequences == null)? this.tfVectors.size() : this.ftSequences.size();
		this.corporaNumber = (this.ftSequences == null)? this.tfVectors.get(0).size() : this.ftSequences.get(0).size();
		this.docNumbers = new ArrayList<List<Integer>>();
		this.corpusTokenNumbers = new ArrayList<List<Integer>>();
		this.averageDocLengths = new ArrayList<List<Double>>();
		this.docLengths = new ArrayList<List<List<Integer>>>();
		
		for (int t = 0; t < this.timestampNumber; ++t) {
			this.docNumbers.add(new ArrayList<Integer>());
			this.corpusTokenNumbers.add(new ArrayList<Integer>());
			this.averageDocLengths.add(new ArrayList<Double>());
			this.docLengths.add(new ArrayList<List<Integer>>());
			for(int j = 0; j < this.corporaNumber; j++){
				this.docNumbers.get(t).add(
						(this.ftSequences == null)? this.tfVectors.get(t).get(j).size() : this.ftSequences.get(t).get(j).size());
				this.docLengths.get(t).add(new ArrayList<Integer>());
				this.corpusTokenNumbers.get(t).add(0);
				for(int d = 0; d < this.docNumbers.get(t).get(j); d++){
					this.docLengths.get(t).get(j).add((this.ftSequences == null)? 
							VectorOper.VectorSum(this.tfVectors.get(t).get(j).get(d)) : this.ftSequences.get(t).get(j).get(d).length);
					this.corpusTokenNumbers.get(t).set(j, this.corpusTokenNumbers.get(t).get(j) + this.docLengths.get(t).get(j).get(d));
				}
				this.averageDocLengths.get(t).add(((double)this.corpusTokenNumbers.get(t).get(j)) / this.docNumbers.get(t).get(j));
			}
			this.totalDocNumber += (int) VectorOper.VectorSum(this.docNumbers.get(t));
			this.totalTokenNumber += (int) VectorOper.VectorSum(this.corpusTokenNumbers.get(t));
		}
		this.averageDocLength = ((double)this.totalTokenNumber) / this.totalDocNumber;
		

	}
	public DocTimeCorpora(
			Alphabet dict,
			OpenIntIntHashMap forwardHash,
			OpenIntIntHashMap backwardHash,
			HashMap<Integer, String> corporaIndex2Name,
			HashMap<String, Integer> corporaName2Index,
			HashMap<String, Integer> corporaTime2Index,
		    HashMap<Integer, String> corporaIndex2Time,
		    List<List<List<int[]>>> ftSequences,
			List<List<List<OpenIntIntHashMap>>> tfVectors,
			List<List<List<Integer>>> docLabels,
			Alphabet labelDescription,
			String corporaDescription
			){
		this.Create(dict, forwardHash, backwardHash, 
				corporaIndex2Name, corporaName2Index, corporaTime2Index, corporaIndex2Time,
				ftSequences, tfVectors, docLabels, labelDescription, corporaDescription);
	}
	
	public String getTypeStr(int type){
		if (corporaIndex2Name != null) {
			return (String) this.corporaIndex2Name.get(type);
		} else {
			return "Unknown name.";
		}
	}
	public String getTimeStr(int time){
		if (corporaIndex2Time != null) {
			return (String) this.corporaIndex2Time.get(time);
		} else {
			return "Unknown name.";
		}
	}
}
