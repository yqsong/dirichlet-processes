package data.structure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import util.VectorOper;

import cc.mallet.types.Alphabet;
import cern.colt.map.OpenIntIntHashMap;

public class DocCorpora {

	public Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;
	public HashMap<Integer, String> corporaIndex2Name;
	public HashMap<String, Integer> corporaName2Index;
	/**
	 * ftSequences[j][d][] is the feature sequence of document d of corpus j
	 */
	public List<List<int[]>> ftSequences;
	/**
	 * tfVectors[j][d] is the term frequences of document d of corpus j
	 */
	public List<List<OpenIntIntHashMap>> tfVectors;
	
	/**
	 * docCategories[j][d] is the category of document 
	 */
	public List<List<Integer>> docLabels;
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
	/**
	 * docNumbers[j] is the number of documents in corpus j
	 */
	public List<Integer> docNumbers;
	/**
	 * total number of documents in all corpora.
	 */
	public int totalDocNumber;
	/**
	 * docLenghts[j][d] is the length of document d of corpus j
	 */
	public List<List<Integer>> docLengths;
	/**
	 * corpusTokenNumbers[j] is the number of tokens in corpus j
	 */
	public List<Double> averageDocLengths;
	public double averageDocLength;
	public List<Integer> corpusTokenNumbers;
	/**
	 * The total number of tokens in all corpora.
	 */
	public int totalTokenNumber;
	
	public String corporaDescription;
	
	public DocCorpora(){
		this.dictionary = null;
		this.ftSequences = null;
		this.tfVectors = null;
		this.docLabels = null;
		this.labelDescription = null;
		this.dictSize = 0;
		this.corporaNumber = 0;
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
	}
	
	
	public String MsgCorporaInfo(){
		String msg = "";
		if(this.corporaDescription != null){
			msg += "Corpora description: " + this.corporaDescription + "\n";
		}
		msg += "Corpora number			:" + this.corporaNumber + "\n";
		msg += "Dictionary size			:" + this.dictSize + "\n";
		msg += "Total document number		:" + this.totalDocNumber + "\n";
		msg += "Total token number		:" + this.totalTokenNumber + "\n";
		msg += "Average document length		:" + String.format("%.2f", this.averageDocLength) + "\n";
		msg += "Info of each corpus:";
		int[] corpusid = new int[this.corporaNumber];
		for(int j = 0; j < this.corporaNumber; j++){
			corpusid[j] = j;
		}
		msg += "\nCorpus			|" + VectorOper.ToString(corpusid, "", "|", "%10d");
		msg += "\nDoc number		|" + VectorOper.ToString(this.docNumbers, "", "|", "%10d", false);
		msg += "\nToken number		|" + VectorOper.ToString(this.corpusTokenNumbers, "", "|", "%10d", false);
		msg += "\nAvg doc length		|" + VectorOper.ToString(this.averageDocLengths, "", "|", "%10.2f");
		
		return msg;
		
	}
	
	public void Create(Alphabet dict,
			OpenIntIntHashMap forwardHash,
			OpenIntIntHashMap backwardHash,
			HashMap<Integer, String> corporaIndex2Name,
			HashMap<String, Integer> corporaName2Index,
			List<List<int[]>> ftSequences,
			List<List<OpenIntIntHashMap>> tfVectors,
			List<List<Integer>> docLabels,
			Alphabet labelDescription,
			String corporaDescription
			){
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;
		this.corporaIndex2Name = corporaIndex2Name;
		this.corporaName2Index = corporaName2Index;
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
		this.corporaNumber = (this.ftSequences == null)? this.tfVectors.size() : this.ftSequences.size();
		this.docNumbers = new ArrayList<Integer>();
		this.corpusTokenNumbers = new ArrayList<Integer>();
		this.averageDocLengths = new ArrayList<Double>();
		this.docLengths = new ArrayList<List<Integer>>();
		
		for(int j = 0; j < this.corporaNumber; j++){
			this.docNumbers.add(
					(this.ftSequences == null)? this.tfVectors.get(j).size() : this.ftSequences.get(j).size());
			this.docLengths.add(new ArrayList<Integer>());
			this.corpusTokenNumbers.add(0);
			for(int d = 0; d < this.docNumbers.get(j); d++){
				this.docLengths.get(j).add((this.ftSequences == null)? 
						VectorOper.VectorSum(this.tfVectors.get(j).get(d)) : this.ftSequences.get(j).get(d).length);
				this.corpusTokenNumbers.set(j, this.corpusTokenNumbers.get(j) + this.docLengths.get(j).get(d));
			}
			this.averageDocLengths.add(((double)this.corpusTokenNumbers.get(j)) / this.docNumbers.get(j));
		}
		this.totalDocNumber = (int) VectorOper.VectorSum(this.docNumbers);
		this.totalTokenNumber = (int) VectorOper.VectorSum(this.corpusTokenNumbers);
		this.averageDocLength = ((double)this.totalTokenNumber) / this.totalDocNumber;

	}
	public DocCorpora(
			Alphabet dict,
			OpenIntIntHashMap forwardHash,
			OpenIntIntHashMap backwardHash,
			HashMap<Integer, String> corporaIndex2Name,
			HashMap<String, Integer> corporaName2Index,
			List<List<int[]>> ftSequences,
			List<List<OpenIntIntHashMap>> tfVectors,
			List<List<Integer>> docLabels,
			Alphabet labelDescription,
			String corporaDescription
			){
		this.Create(dict, forwardHash, backwardHash, 
				corporaIndex2Name, corporaName2Index,
				ftSequences, tfVectors, docLabels, labelDescription, corporaDescription);
	}
	
	public String getTypeStr(int type){
		if (corporaIndex2Name != null) {
			return (String) this.corporaIndex2Name.get(type);
		} else {
			return "Unknown name.";
		}
	}
	
}
