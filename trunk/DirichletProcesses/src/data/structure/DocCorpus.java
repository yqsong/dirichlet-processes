package data.structure;

import java.util.ArrayList;
import java.util.List;

import util.VectorOper;

import cc.mallet.types.Alphabet;
import cern.colt.map.OpenIntIntHashMap;

/**
 * @author Jianwen Zhang
 */
public class DocCorpus{
	public Alphabet dictionary;
	public OpenIntIntHashMap forwardHash;
	public OpenIntIntHashMap backwardHash;

	public List<int[]> ftSequences;
	public List<OpenIntIntHashMap> tfVectors;//term frequencies
	
	public List<Integer> docLabels;
	public Alphabet labelDescription;
	
	public int dictSize;
	public List<Integer> docLengths;
	public double averageDocLength;
	public int docNumber;
	public int tokenNumber;
	
	String corpusDescription;
	public DocCorpus(Alphabet dict
			, OpenIntIntHashMap forwardHash
			, OpenIntIntHashMap backwardHash
			, List<int[]> ftSequences
			, List<OpenIntIntHashMap> tfVectors
			, List<Integer> docLabels
			, Alphabet labelDescription
			, String corpusDescription){
		if(dict == null)
			System.err.println("DocCorpus: there is no dictionary!");
		if((ftSequences == null) && (tfVectors == null))
			System.err.println("DocCorpus: there is not any feature provided!");
		
		this.dictionary = dict;
		this.forwardHash = forwardHash;
		this.backwardHash = backwardHash;
		
		this.ftSequences = ftSequences;
		this.tfVectors = tfVectors;
		this.docLabels = docLabels;
		this.labelDescription = labelDescription; 
		this.corpusDescription = corpusDescription;
		
//		this.dictSize = this.dictionary.size();
		this.dictSize = (this.forwardHash == null) ? this.dictionary.size() : this.forwardHash.size();
		this.docNumber = this.ftSequences.size();
		this.docLengths = new ArrayList<Integer>();
		this.tokenNumber = 0;
		for(int i = 0; i < this.docNumber; i++){
			this.docLengths.add( (this.ftSequences == null)? 
					VectorOper.VectorSum(this.tfVectors.get(i)) : this.ftSequences.get(i).length );
			this.tokenNumber += this.docLengths.get(i);
		}
		this.averageDocLength = ((double)this.tokenNumber) / this.docNumber;
	}
	
	public String MsgCorpusInfo(){
		String msg = "";
		if(this.corpusDescription != null){
			msg += "Corpus description: " + this.corpusDescription + "\n";
		}
		msg += "Dictionary size: 			" + this.dictSize + "\n";
		msg += "Document number: 			" + this.docNumber + "\n";
		msg += "Token number: 				" + this.tokenNumber + "\n";
		msg += "Average document length: 		" + String.format("%.2f", this.averageDocLength) + "\n";
		
		return msg;
	}
	
}