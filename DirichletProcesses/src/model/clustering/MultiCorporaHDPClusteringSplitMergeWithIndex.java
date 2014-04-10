package model.clustering;

import index.PrepareCorporaTimeDataForIndex;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import model.SummarySplitMergeResults;
import util.RandomSampler;
import cern.colt.map.OpenIntDoubleHashMap;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenIntObjectHashMap;
import data.structure.DocTimeCorpora;

public class MultiCorporaHDPClusteringSplitMergeWithIndex extends
		MultiCorporaHDPClusteringSplitMerge {
	
	public List<List<List<String>>> structuredURIs = null;

	/**
	 * Constructor for the case not need to sample hyper parameters \alpha or \gamma
	 * @param docCorpora
	 * @param b0
	 * @param gamma
	 * @param alpha
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCorporaHDPClusteringSplitMergeWithIndex(
			DocTimeCorpora docTimeCorpora
			, List<List<List<String>>> structuredURIs
			, double b0
			, double gamma, double alpha
			, int maxIter
			, int burnIn
			){
		super(docTimeCorpora
			, b0
			, 0.0, 0.0
			, 0.0, 0.0
			, gamma, alpha
			, false//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
		this.structuredURIs = structuredURIs;
	}
	/**
	 * Constructor for the case when sampling hyper parameters are required.
	 * @param docCorpora
	 * @param b0
	 * @param gamma_a
	 * @param gamma_b
	 * @param alpha_a
	 * @param alpha_b
	 * @param maxIter
	 * @param burnIn
	 */
	public MultiCorporaHDPClusteringSplitMergeWithIndex(
			DocTimeCorpora docTimeCorpora
			, List<List<List<String>>> structuredURIs
			, double b0
			, double gamma_a, double gamma_b
			, double alpha_a, double alpha_b
			, int maxIter
			, int burnIn
			){
		super(docTimeCorpora
			, b0
			, gamma_a, gamma_b
			, alpha_a, alpha_b
			, 0.0, 0.0
			, true//whether sampling hyper parameters \alpha and \gamma
			, maxIter
			, burnIn);
		this.structuredURIs = structuredURIs;
	}

	
	
	public void Sampling(int initK, String outputFolder){
		
		int conpiter = 20;//default iteration for sampling hyperparameters \alpha and \gamma
		long starttime, time0, usedtime;
		int topN = 30;
		
		InitializeT0 (initK);
		Log(0, -1, 0, false);
		
		for(int iter = 0; iter < this.maxIter; iter++){
			starttime = System.currentTimeMillis();
			time0 = starttime;
			SamplingZ(0, true);
			SamplingMandBeta(0);
		
			if(this.bSamplingHyper){
				SamplingAlpha(conpiter, 0);
				SamplingGamma(conpiter);
			}
			usedtime = System.currentTimeMillis() - time0;
			if (iter % 10 == 0) {
				Log(0, iter, usedtime, false);	
			}
		}
		SummarySplitMergeResults splitmergeSum1 = new SummarySplitMergeResults(this.docTimeCorpora, this.globalClusterSize);
		String strTimeT = splitmergeSum1.SummarizeTimestampT (0, this.Z.get(0), topN);
		System.out.println(strTimeT);
		
		File outputFile = null;
		FileWriter writer;
		if (outputFolder != null) {
			try {
				outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(0) + "_cluster.txt" );
				writer = new FileWriter(outputFile, true);
				writer.write(strTimeT);
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		for (int t = 1; t < this.docTimeCorpora.timestampNumber; ++t) {

			List<List<List<Integer>>> passedLabels = PredictPassedTimestamps(t);
//			List<List<Integer>> predictLabels = PredictTimestamp(t);
			
			System.out.println("LogInfo on Single Timestamp. \n");
			
			InitializeT (t);
			for(int iter = 0; iter < this.maxIter; iter++){
				starttime = System.currentTimeMillis();
				time0 = starttime;
				SamplingZ(t, true);
				SamplingMandBeta(t);
				
				if(this.bSamplingHyper){
					SamplingAlpha(conpiter, t);
					SamplingGamma(conpiter);
				}
				usedtime = System.currentTimeMillis() - time0;
				if (iter % 10 == 0) {
					Log(t, iter, usedtime, false);	
				}
			}
			
			System.out.println("LogInfo on Multiple Timestamp. \n");
			for(int iter = 0; iter < this.maxIter; iter++){
				starttime = System.currentTimeMillis();
				time0 = starttime;
				for (int tao = 0; tao <= t; ++tao) {
					ResamplingZ(tao);
				}
				SamplingMandBeta(t);
				
				if(this.bSamplingHyper){
					SamplingAlpha(conpiter, t);
					SamplingGamma(conpiter);
				}
				usedtime = System.currentTimeMillis() - time0;
				if (iter % 10 == 0) {
					Log(t, iter, usedtime, false);	
				}
			}
			
			SummarySplitMergeResults splitmergeSum = new SummarySplitMergeResults(this.docTimeCorpora, this.globalClusterSize);

			strTimeT = splitmergeSum.SummarizeTimestampT (t, this.Z.get(t), topN);
			System.out.println(strTimeT);
			if (outputFolder != null) {
				try {
					outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(t) + "_cluster.txt" );
					writer = new FileWriter(outputFile, true);
					writer.write(strTimeT);
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			String strTimePassedOutput = splitmergeSum.SummarizeTimestampOutput (t, passedLabels, this.Z.subList(0, t + 1), topN);
//			System.out.println(strTimePassedOutput);
			if (outputFolder != null) {
				try {
					outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(t - 1) + "_output.txt" );
					writer = new FileWriter(outputFile, true);
					writer.write(strTimePassedOutput);
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			String strTimeNewInput = splitmergeSum.SummarizeTimestampInput (t, passedLabels.get(t), this.Z.get(t), topN);
//			System.out.println(strTimeNewInput);
			if (outputFolder != null) {
				try {
					outputFile = new File(outputFolder + "\\" + this.docTimeCorpora.corporaIndex2Time.get(t) + "_input.txt" );
					writer = new FileWriter(outputFile, true);
					writer.write(strTimeNewInput);
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			
			System.out.println("\n" + "Training end of Timestamp " + t + ". \n");
		}
		this.docTimeCorpora.docLabels = this.Z;
	}
	
	public static void main(String[] args){

		String luceneSrc = "D:\\data_test\\vis\\text_index\\";
		File stopwordFile = new File("D:\\data_test\\vis\\stopwords.txt");
		String outputFolder = "D:\\data_test\\vis\\splitmerge_output";
//		String queryStr = "*:*";
		String queryStr = "confName:infovis";
		PrepareCorporaTimeDataForIndex	prepareTimeData = new PrepareCorporaTimeDataForIndex();
		DocTimeCorpora docTimeCorpora = prepareTimeData.initializeData (
				luceneSrc, 
				queryStr, 
				true, 
				1, 
				1, 
				stopwordFile,
				"link",
				"title",
				"abstract",
				"confName",
				"year",
				"year");

		
		
		double b0 = 0.5;
		double gamma = 5, alpha = 1;
//		double gamma_a = 2, gamma_b = 1;
//		double alpha_a = 2, alpha_b = 1;
		
		int maxIter = 200, burnin = 10;
		int initK = 5;
		int loglevel = 3;
		MultiCorporaHDPClusteringSplitMerge hdpCluster = new MultiCorporaHDPClusteringSplitMerge(
				docTimeCorpora
				, b0
//				, gamma_a, gamma_b
//				, alpha_a, alpha_b
				, gamma, alpha //if sampling hyper, comment this line and uncomment above 2 lines, vice versa
				, maxIter
				, burnin);
		hdpCluster.setLogLevel(loglevel);
		hdpCluster.Sampling(initK, outputFolder);
		
	}

}
