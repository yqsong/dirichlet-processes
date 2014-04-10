package model;

import java.util.List;

import util.VectorOper;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.DocTimeCorpora;
import data.structure.Keywords;

public class SummarySplitMergeResults {
	
	DocTimeCorpora docTimeCorpora;
	OpenIntIntHashMap globalClusterSize;
	
	SummarySplitMergeResults (DocTimeCorpora docTimeCorpora,
								OpenIntIntHashMap globalClusterSize) {
		this.docTimeCorpora = docTimeCorpora;
		this.globalClusterSize = globalClusterSize;
	}
	
	String SummarizeTimestampT (int time_t, List<List<Integer>> labels, int topN) {
		String msg = "";
		msg += this.docTimeCorpora.corporaIndex2Time.get(time_t) + " \n\n" ;
		OpenIntObjectHashMap clusterTermFreqNew = new OpenIntObjectHashMap();
		OpenIntIntHashMap globalClusterSizeAtT = new OpenIntIntHashMap();

		IntArrayList clusterIDs = this.globalClusterSize.keys();
		for (int i = 0; i < this.globalClusterSize.size(); ++i) {
			int kid = clusterIDs.get(i);
			clusterTermFreqNew.put(kid, new OpenIntIntHashMap());
			globalClusterSizeAtT.put(kid, 0);
		}
		
		for (int j = 0; j < labels.size(); ++j) {
			for (int i = 0; i < labels.get(j).size(); ++i) {
				int label = labels.get(j).get(i);
				int value = globalClusterSizeAtT.get(label); 
				globalClusterSizeAtT.put(label, value + 1);
				
//				if (globalClusterSizeAtT.get(label) != 0) {
//					int value = globalClusterSizeAtT.get(label); 
//					globalClusterSizeAtT.put(label, value + 1);
//				} else {
//					globalClusterSizeAtT.put(label, 1);
//				}
				
				OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(time_t).get(j).get(i);
				OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(label);
				VectorOper.VectorAddEqual(cluster, docTf);
				
//				if (clusterTermFreqNew.get(label) != null) {
//					OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(label);
//					VectorOper.VectorAddEqual(cluster, docTf);
//				} else {
//					clusterTermFreqNew.put(label, new OpenIntIntHashMap());
//					OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(label);
//					VectorOper.VectorAddEqual(cluster, docTf);
//				}
			}
		}
		
		msg += SummarizeClusterKeywords(topN, globalClusterSizeAtT, clusterTermFreqNew);
		return msg;
	}
	
	String SummarizeTimestampOutput (int time_t, List<List<List<Integer>>> orgLabels,  List<List<List<Integer>>> updateLabels, int topN) {
		OpenIntObjectHashMap matrix = new OpenIntObjectHashMap();
		String msg = "Output from ";
		msg += this.docTimeCorpora.corporaIndex2Time.get(time_t - 1) + " \n\n" ;
		
		OpenIntObjectHashMap clusterClusterTermFreqNew = new OpenIntObjectHashMap();
		for (int t = 0; t < orgLabels.size(); ++t) {
			for (int j = 0; j < orgLabels.get(t).size(); ++j) {
				for (int i = 0; i < orgLabels.get(t).get(j).size(); ++i) {
					int orgLabel = orgLabels.get(t).get(j).get(i);
					int updateLabel = updateLabels.get(t).get(j).get(i);
					OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(t).get(j).get(i);
					
					if (clusterClusterTermFreqNew.get(orgLabel) != null) {
						OpenIntObjectHashMap clusterTermFreqNew = 
							(OpenIntObjectHashMap) clusterClusterTermFreqNew.get(orgLabel);
						
						OpenIntIntHashMap row = (OpenIntIntHashMap) matrix.get(orgLabel);
						
						if (clusterTermFreqNew.get(updateLabel) != null) {
							OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(updateLabel);
							VectorOper.VectorAddEqual(cluster, docTf);
							
							int num = row.get(updateLabel);
							row.put(updateLabel, num + 1);
						} else {
							clusterTermFreqNew.put(updateLabel, new OpenIntIntHashMap());
							OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(updateLabel);
							VectorOper.VectorAddEqual(cluster, docTf);
							
							row.put(updateLabel, 1);
						}
					} else {
						clusterClusterTermFreqNew.put(orgLabel, new OpenIntObjectHashMap());
						OpenIntObjectHashMap clusterTermFreqNew = 
							(OpenIntObjectHashMap) clusterClusterTermFreqNew.get(orgLabel);
						clusterTermFreqNew.put(updateLabel, new OpenIntIntHashMap());
						OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(updateLabel);
						VectorOper.VectorAddEqual(cluster, docTf);
						
						matrix.put(orgLabel, new OpenIntIntHashMap());
						OpenIntIntHashMap row = (OpenIntIntHashMap) matrix.get(orgLabel);
						row.put(updateLabel, 1);
					}
				}
			}
		}
		
		IntArrayList rowKeys = this.globalClusterSize.keys();
//		IntArrayList rowKeys = matrix.keys();
		for (int i = 0; i < rowKeys.size(); ++i) {
			int rowID = rowKeys.get(i);
			OpenIntIntHashMap row = (OpenIntIntHashMap) matrix.get(rowID);
			
			if (row != null) {
				IntArrayList columnKeys = row.keys();
				
				int sum = 0;
				for (int j = 0; j < columnKeys.size(); ++j) {
					int columnID = columnKeys.get(j);
					int value = row.get(columnID);
					sum += value;
				}
				
				if (sum == 0) {
					OpenIntObjectHashMap temp = new OpenIntObjectHashMap();
					temp.put(rowID, new OpenIntIntHashMap());
					row = new OpenIntIntHashMap();
					row.put(rowID, 0);
					msg += "From Cluster: " + rowID + ", " + 1.0 + "\n";
					msg += "To Cluster: " + rowID + "\n";
					msg += SummarizeClusterKeywords(topN, row, temp);
					msg += "\n";
					msg += "\n";
				} else {
					for (int j = 0; j < columnKeys.size(); ++j) {
						int columnID = columnKeys.get(j);
						int value = row.get(columnID);
						if (value > 0.0 ) {
							OpenIntObjectHashMap clusterTermFreqNew = 
								(OpenIntObjectHashMap) clusterClusterTermFreqNew.get(rowID);
							OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(columnID);
							
							OpenIntObjectHashMap temp = new OpenIntObjectHashMap();
							temp.put(columnID, cluster);
							
							msg += "From Cluster: " + rowID + ", " + (double)value/(double)sum + "\n";
							msg += "To Cluster: " + columnID + "\n";
							msg += SummarizeClusterKeywords(topN, row, temp);
							msg += "\n";
							msg += "\n";
						}
					}
				}
			} else {
				OpenIntObjectHashMap temp = new OpenIntObjectHashMap();
				temp.put(rowID, new OpenIntIntHashMap());
				row = new OpenIntIntHashMap();
				row.put(rowID, 0);
				msg += "From Cluster: " + rowID + ", " + 1.0 + "\n";
				msg += "To Cluster: " + rowID + "\n";
				msg += SummarizeClusterKeywords(topN, row, temp);
				msg += "\n";
				msg += "\n";
			}
		}
		
		return msg;
	}
	
	String SummarizeTimestampInput (int time_t, List<List<Integer>> predictLabels, List<List<Integer>> updateLabels, int topN) {
		OpenIntObjectHashMap matrix = new OpenIntObjectHashMap();
		String msg = "Input to ";
		msg += this.docTimeCorpora.corporaIndex2Time.get(time_t) + " \n\n" ;
		
		OpenIntObjectHashMap clusterClusterTermFreqNew = new OpenIntObjectHashMap();
		for (int j = 0; j < predictLabels.size(); ++j) {
			for (int i = 0; i < predictLabels.get(j).size(); ++i) {
				int predictLabel = predictLabels.get(j).get(i);
				int updateLabel = updateLabels.get(j).get(i);
				OpenIntIntHashMap docTf = this.docTimeCorpora.tfVectors.get(time_t).get(j).get(i);
				
				if (clusterClusterTermFreqNew.get(updateLabel) != null) {
					OpenIntObjectHashMap clusterTermFreqNew = 
						(OpenIntObjectHashMap) clusterClusterTermFreqNew.get(updateLabel);
					
					OpenIntIntHashMap row = (OpenIntIntHashMap) matrix.get(updateLabel);
					
					if (clusterTermFreqNew.get(predictLabel) != null) {
						OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(predictLabel);
						VectorOper.VectorAddEqual(cluster, docTf);
						
						int num = row.get(predictLabel);
						row.put(predictLabel, num + 1);
					} else {
						clusterTermFreqNew.put(predictLabel, new OpenIntIntHashMap());
						OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(predictLabel);
						VectorOper.VectorAddEqual(cluster, docTf);
						
						row.put(predictLabel, 1);
					}
				} else {
					clusterClusterTermFreqNew.put(updateLabel, new OpenIntObjectHashMap());
					OpenIntObjectHashMap clusterTermFreqNew = 
						(OpenIntObjectHashMap) clusterClusterTermFreqNew.get(updateLabel);
					clusterTermFreqNew.put(predictLabel, new OpenIntIntHashMap());
					OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(predictLabel);
					VectorOper.VectorAddEqual(cluster, docTf);
					
					matrix.put(updateLabel, new OpenIntIntHashMap());
					OpenIntIntHashMap row = (OpenIntIntHashMap) matrix.get(updateLabel);
					row.put(predictLabel, 1);
				}	
					
			}
		}
		
		IntArrayList rowKeys = this.globalClusterSize.keys();
//		IntArrayList rowKeys = matrix.keys();
		for (int i = 0; i < rowKeys.size(); ++i) {
			int rowID = rowKeys.get(i);
			OpenIntIntHashMap row = (OpenIntIntHashMap) matrix.get(rowID);
			
			if (row != null) {
				IntArrayList columnKeys = row.keys();
				
				int sum = 0;
				for (int j = 0; j < columnKeys.size(); ++j) {
					int columnID = columnKeys.get(j);
					int value = row.get(columnID);
					sum += value;
				}
				
				if (sum == 0) {
					OpenIntObjectHashMap temp = new OpenIntObjectHashMap();
					temp.put(rowID, new OpenIntIntHashMap());
					row = new OpenIntIntHashMap();
					row.put(rowID, 0);
					msg += "To Cluster: " + rowID + ", " + 1.0 + "\n";
					msg += "From Cluster: " + rowID + "\n";
					msg += SummarizeClusterKeywords(topN, row, temp);
					msg += "\n";
					msg += "\n";
				} else {
					for (int j = 0; j < columnKeys.size(); ++j) {
						int columnID = columnKeys.get(j);
						int value = row.get(columnID);
						if (value > 0.0 ) {
							OpenIntObjectHashMap clusterTermFreqNew = 
								(OpenIntObjectHashMap) clusterClusterTermFreqNew.get(rowID);
							OpenIntIntHashMap cluster = (OpenIntIntHashMap) clusterTermFreqNew.get(columnID);
							
							OpenIntObjectHashMap temp = new OpenIntObjectHashMap();
							temp.put(columnID, cluster);
							
							msg += "To Cluster: " + rowID + ", " + (double)value/(double)sum + "\n";
							msg += "From Cluster: " + columnID + "\n";
							msg += SummarizeClusterKeywords(topN, row, temp);
							msg += "\n";
							msg += "\n";
						}
					}
				}
				
			} else {
				OpenIntObjectHashMap temp = new OpenIntObjectHashMap();
				temp.put(rowID, new OpenIntIntHashMap());
				row = new OpenIntIntHashMap();
				row.put(rowID, 0);
				msg += "To Cluster: " + rowID + ", " + 1.0 + "\n";
				msg += "From Cluster: " + rowID + "\n";
				msg += SummarizeClusterKeywords(topN, row, temp);
				msg += "\n";
				msg += "\n";
			}
			
		}
		
		return msg;
	}
	
	public String SummarizeClusterKeywords(int topN, OpenIntIntHashMap sizes, OpenIntObjectHashMap clusterTermFreqNew){
		String msg = "ID	|Size	|	keyword:frequency,\n";
		OpenIntObjectHashMap topicKeywords = GetClusterKeywords(topN, clusterTermFreqNew);
		if(topicKeywords == null)
			return msg;
		IntArrayList topicIDs = topicKeywords.keys();
		for(int k = 0; k < topicIDs.size(); k++){
			int kid = topicIDs.get(k);
			msg += kid + "	|" + sizes.get(kid) + "	|";
			Keywords keywds = (Keywords) topicKeywords.get(kid);
			String[] words = keywds.getWords();
			int[] freqs = keywds.getFreqs();
			for(int i = 0; i < words.length; i++){
				msg += words[i] + ":" + freqs[i] + ", ";
			}
			msg += "\n";
		}
		return msg;
	}
	
	public OpenIntObjectHashMap GetClusterKeywords(int topN, OpenIntObjectHashMap clusterTermFreqNew) {
		// TODO Auto-generated method stub
		if(this.docTimeCorpora.dictionary == null){
			System.err.println("getTopicsKeywords failed, as no dictionary provided!");
			return null;
		}
		OpenIntObjectHashMap ret = new OpenIntObjectHashMap();
		IntArrayList ks = clusterTermFreqNew.keys();
		OpenIntIntHashMap topicDist = null;
		for(int ki = 0; ki < ks.size(); ki++){
			int kid = ks.get(ki);
			topicDist = (OpenIntIntHashMap) clusterTermFreqNew.get(kid);
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
				if (this.docTimeCorpora.backwardHash != null) {
					wordID = this.docTimeCorpora.backwardHash.get(wordID);
				}
				words[i] = (String) this.docTimeCorpora.dictionary.lookupObject(wordID);
				freqs[i] = iwfrqs[maxN[i]];
			}
			ret.put(kid, new Keywords(words, freqs));
		}
		return ret;
	} 
	
}
