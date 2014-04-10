package index;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import cc.mallet.pipe.CharSequence2TokenSequence;
import cc.mallet.pipe.FeatureSequence2FeatureVector;
import cc.mallet.pipe.Input2CharSequence;
import cc.mallet.pipe.Noop;
import cc.mallet.pipe.Pipe;
import cc.mallet.pipe.SerialPipes;
import cc.mallet.pipe.Target2Label;
import cc.mallet.pipe.TokenSequence2FeatureSequence;
import cc.mallet.pipe.TokenSequenceLowercase;
import cc.mallet.pipe.TokenSequenceRemoveStopwords;
import cc.mallet.types.Alphabet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.types.SparseVector;
import cern.colt.map.OpenIntIntHashMap;

import data.structure.DocCorpora;
import data.structure.DocTimeCorpora;

public class PrepareCorporaDataForIndex {
	protected static double minDFCount = 2;
	protected static double maxDFCount = 10000000;
	protected static int seed = 0;
	
	public InstanceList tflist = null;
	public List<List<OpenIntIntHashMap>> tfVectors = null;
	public List<List<int[]>> featureSequences = null;
	public List<List<String>> structuredURIs = null;
    protected int wordNumber = 0;
    
    protected List<String> docStrings = null;
    protected List<String> docURIs = null;
    protected List<String> docCorporaName = null;
    protected List<Long> docTime = null;
    
    protected Alphabet alphabet = null;
    protected OpenIntIntHashMap forwardHash = null;
    protected OpenIntIntHashMap backwardHash = null;
    
    protected int corporaNumber = 0;
    protected int maxTimestampNumber = 0;
    protected int minTimestampNumber = 0;
    protected HashMap<String, Integer> corporaName2Index = null;
    protected HashMap<Integer, String> corporaIndex2Name = null;
    
    protected HashMap<String, Integer> corporaTime2Index = null;
    protected HashMap<Integer, String> corporaIndex2Time = null;
    
	private InstanceList loadInstanceSubsetList (File stopwordFile, boolean isUseTFVector) {
		Pipe instancePipe;
		instancePipe = new SerialPipes (new Pipe[] {
                new Target2Label(),
                new Input2CharSequence(),
                ((Pipe) new CharSequence2TokenSequence()),
                ((Pipe) new TokenSequenceLowercase()),
                ((stopwordFile == null) ? ((Pipe) new TokenSequenceRemoveStopwords(false, true)) : 
                	((Pipe) new TokenSequenceRemoveStopwords(false, true).addStopWords(stopwordFile))),
                ((Pipe) new TokenSequence2FeatureSequence()),
                (isUseTFVector ? ((Pipe) new FeatureSequence2FeatureVector()) : (Pipe) new Noop()),
            });
             
		tflist = new InstanceList(instancePipe); 
		
		corporaName2Index = new HashMap<String, Integer>();
        corporaIndex2Name = new HashMap<Integer, String>();
        corporaNumber = 0;
        
        maxTimestampNumber = 0;
        minTimestampNumber = Integer.MAX_VALUE;
        corporaTime2Index = new HashMap<String, Integer>();
        corporaIndex2Time = new HashMap<Integer, String>();
	        
        int docIndex = 0;
		 for (int i = 0; i < docStrings.size(); ++i) {
			 if (docIndex % 1000 == 0) {
	             	System.out.println(">>>[LOG]: Mallet parse " + docIndex + " documents.");
	             }           
			 docIndex++;
	             
			 String text = docStrings.get(i);
			 Instance carrier;
			 carrier = instancePipe.instanceFrom(new Instance (text, 0, null, null));
			 tflist.add (carrier);
			 if (!corporaName2Index.containsKey(docCorporaName.get(i))) {
				 corporaName2Index.put(docCorporaName.get(i), corporaNumber);
				 corporaIndex2Name.put(corporaNumber, docCorporaName.get(i));
				 corporaNumber++;
			 }
			 
			 long time = docTime.get(i);
			 Calendar calendar = Calendar.getInstance();
			 calendar.setTimeInMillis(time);
			 int year = calendar.get(Calendar.YEAR);
			 int month = calendar.get(Calendar.MONTH);
			 int timeIndex = (year - 2000) * 12 + month;
			 if (maxTimestampNumber < timeIndex) {
				 maxTimestampNumber = timeIndex;
			 }
			 if (minTimestampNumber > timeIndex) {
				 minTimestampNumber = timeIndex;
			 }
		 }
		 for (int i = minTimestampNumber; i < (maxTimestampNumber + 1); ++i){
			 int year = i / 12;
			 int month = (i - year * 12);
			 year += 2000;
			 SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMM");
			 Calendar cal = Calendar.getInstance();
			 cal.set(Calendar.YEAR, year);
			 cal.set(Calendar.MONTH, month);
			 cal.set(Calendar.DAY_OF_MONTH, 1);
			 Date date = cal.getTime();
			 String timeStr = dateFormat.format(date);
			 corporaTime2Index.put(timeStr, i - minTimestampNumber);
			 corporaIndex2Time.put(i - minTimestampNumber, timeStr);
		 }
		 
		 return tflist;
	}
	
	public void loadFeatureVector(File stopwordFile) {
		
		InstanceList tflist = loadInstanceSubsetList (stopwordFile, true);
        alphabet = tflist.getAlphabet();
        tfVectors = new ArrayList<List<OpenIntIntHashMap>>(); 
        structuredURIs = new ArrayList<List<String>>(); 
        for (int i = 0; i < corporaNumber; ++i) {
        	tfVectors.add(new ArrayList<OpenIntIntHashMap>());
        	structuredURIs.add(new ArrayList<String>());
        }
        
        int[] alphabetCount = new int[alphabet.size()];
        double[] dfCount = new double[alphabet.size()];
        double[] tfCount = new double[alphabet.size()];
        int docCount = 0;
        for (Instance carrier : tflist) {
            if (docCount % 1000 == 0) {
            	System.out.println(">>>[LOG]: Sum feature DF " + docCount + " documents.");
            }       
            SparseVector sv = (SparseVector) carrier.getData();
            int[] index = sv.getIndices();
            double[] value = sv.getValues();
            
            for (int j = 0; j < index.length; ++j) {
            	if (value[j] > 0.0) {
            		alphabetCount[index[j]] += value[j];
            		dfCount[index[j]] += 1;
            		tfCount[index[j]] += value[j];
            	}
            }
            docCount++;

        }
        for (int i = 0; i < dfCount.length; ++i) {
        	if (dfCount[i] > 0) {
        		dfCount[i] = Math.log(tflist.size()/dfCount[i]);
        	}
        }
        
        double[] tempTFCount = Arrays.copyOf(tfCount, tfCount.length);
        Arrays.sort(tempTFCount);
//        minDFCount = tempTFCount[tempTFCount.length - 5000];//temp[mutualInformation.length - 5000];
        
        // remove unusual words
        forwardHash = new OpenIntIntHashMap();
        backwardHash = new OpenIntIntHashMap();
        int wordIndex = 0;
        for (int i = 0; i < alphabetCount.length; ++i) {
        	if (alphabetCount[i] > minDFCount 
        			&& alphabetCount[i] < maxDFCount) {
        		forwardHash.put(i, wordIndex);
        		backwardHash.put(wordIndex, i);
        		wordIndex++;
        	}
        }
        wordNumber = wordIndex;

        int totalWord = 0;
		docCount = 0;
        for (int i = 0; i < tflist.size(); ++i) {
        	Instance carrier = tflist.get(i);
        	if (docCount % 1000 == 0) {
            	System.out.println(">>>[LOG]: Convert to TF vector " + docCount + " documents.");
            }       
            docCount++;
            SparseVector sv = (SparseVector) carrier.getData();
            int[] index = sv.getIndices();
            double[] value = sv.getValues();
            OpenIntIntHashMap vector = new OpenIntIntHashMap();
            for (int j = 0; j < index.length; ++j) {
            	if (value[j] > 0.0 
            			&& alphabetCount[index[j]] > minDFCount 
            			&& alphabetCount[index[j]] < maxDFCount) {
            		vector.put(forwardHash.get(index[j]), (int) value[j]);
            		totalWord += value[j];
            	}
            }
            String corporaName = docCorporaName.get(i);
            int corporaIndex = corporaName2Index.get(corporaName);
            tfVectors.get(corporaIndex).add(vector);
            structuredURIs.get(corporaIndex).add(this.docURIs.get(i));
//            carrier.unLock();
//            carrier.clearSource();
//            carrier.setData(null);
        }
        System.out.println(">>>[LOG]: Doc number: " + tflist.size());
        System.out.println(">>>[LOG]: Word number: " + totalWord);
        System.out.println(">>>[LOG]: Vocabulary size: " + forwardHash.size());
	}
	
	public void loadFeatureSequence(File stopwordFile) {
		
		InstanceList tflist = loadInstanceSubsetList (stopwordFile, false);
        alphabet = tflist.getAlphabet();
        featureSequences = new ArrayList<List<int[]>>(); 
        structuredURIs = new ArrayList<List<String>>(); 
        for (int i = 0; i < corporaNumber; ++i) {
        	featureSequences.add(new ArrayList<int[]>());
        	structuredURIs.add(new ArrayList<String>());
        }
        
        int totalWord = 0;
		int docCount = 0;
        for (int i = 0; i < tflist.size(); ++i) {
        	Instance carrier = tflist.get(i);
        	if (docCount % 1000 == 0) {
            	System.out.println(">>>[LOG]: Convert to sequence " + docCount + " documents.");
            }       
            docCount++;
            FeatureSequence fs = (FeatureSequence) carrier.getData();
            int[] indices = fs.toFeatureIndexSequence();
            String corporaName = docCorporaName.get(i);
            int corporaIndex = corporaName2Index.get(corporaName);
            featureSequences.get(corporaIndex).add(indices);
            structuredURIs.get(corporaIndex).add(this.docURIs.get(i));
        }
        System.out.println(">>>[LOG]: Doc number: " + tflist.size());
        System.out.println(">>>[LOG]: Word number: " + totalWord);
        System.out.println(">>>[LOG]: Vocabulary size: " + alphabet.size());
	}
	
	public DocCorpora initializeData (List<String> uri, 
			List<String> contents, 
			List<Long> timeList, 
			List<String> corporaList, 
			boolean isUseTFVector,
			File stopwordFile) {
    	
		docStrings = contents;
        docURIs = uri;
        docCorporaName = corporaList;
        docTime = timeList;
        
        DocCorpora corpora = null;
    	if (isUseTFVector) {
    		loadFeatureVector(stopwordFile);
    		corpora = new DocCorpora(alphabet, forwardHash, backwardHash,
        			corporaIndex2Name, corporaName2Index,
        			null, tfVectors, null, null, "TIARA Data Analysis");
    	} else {
    		loadFeatureSequence(stopwordFile);
    		corpora = new DocCorpora(alphabet, forwardHash, backwardHash,
        			corporaIndex2Name, corporaName2Index,
        			featureSequences, null, null, null, "TIARA Data Analysis");
    	}
    	
    	return corpora;
    }
}
