package index;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.queryParser.QueryParser;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;

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

public class PrepareVisTimeCorpora {
	protected static double minDFCount = 10;
	protected static double maxDFCount = 1000;
	protected static int seed = 0;
	
	public List<List<List<OpenIntIntHashMap>>> tfVectors = null;
	public List<List<List<int[]>>> featureSequences = null;
    protected int wordNumber = 0;
    
    protected List<String> docStringsAll = null;
    protected List<String> docURIsAll = null;
    protected List<String> docCorporaNameAll = null;
    protected List<Integer> docTimeAll = null;
    
    protected List<String> docStrings = null;
    protected List<String> docURIs = null;
    protected List<String> docCorporaName = null;
    protected List<Integer> docTime = null;
    
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
    
	public List<String> loadDataString (String inputDirectory, String queryStr, double dataPercentage) {
        IndexSearcher searcher = null;
        
        try{
            searcher = new IndexSearcher(inputDirectory);
        } catch (Exception e) {
            e.printStackTrace();
        }        
        
		QueryParser queryParser = new QueryParser("queryName", new StandardAnalyzer());
		
        docStringsAll = new ArrayList<String>();
        docURIsAll = new ArrayList<String>();
        docCorporaNameAll = new ArrayList<String>();
        docTimeAll = new ArrayList<Integer>();
        
        int docIndex = 0;
        
        try {
			Query query = queryParser.parse(queryStr);
			TopDocs hits = searcher.search(query, searcher.maxDoc());
			
			ScoreDoc[] docs = hits.scoreDocs;
			List<Integer> permIndex = new ArrayList<Integer>();
			 for (int i = 0; i < docs.length; ++i) {
				 permIndex.add(i);
			 }
			 Collections.shuffle(permIndex);
			 int subsetNum = (int)(docs.length * dataPercentage);
			 List<Integer> subIndices = permIndex.subList(0, subsetNum);
			
			
			for (int i = 0; i < subIndices.size(); ++i) {
				ScoreDoc doc = docs[subIndices.get(i)];
				Document document = searcher.doc(doc.doc);

				String uri = document.get("link");
                docURIsAll.add(uri);
                
                String content = document.get("title") + document.get("abstract");
                docStringsAll.add(content);
                
                String corporaName = document.get("confName");
                docCorporaNameAll.add(corporaName);
                
                String time = document.get("year");
                docTimeAll.add(Integer.parseInt(time));

                if (docIndex % 1000 == 0) {
                	System.out.println(">>>[LOG]: Loaded String " + docIndex + " documents.");
                }           
                docIndex++;
			}
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}
        
        System.out.println(">>>[LOAD]: END");
        
        return docStringsAll;
	}
	
	private InstanceList loadInstanceSubsetList (double dataPercentage, File stopwordFile, boolean isUseTFVector) {
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
             
		InstanceList tflist = new InstanceList(instancePipe); 

		List<Integer> permIndex = new ArrayList<Integer>();
		for (int i = 0; i < docStringsAll.size(); ++i) {
			permIndex.add(i);
		}
		Collections.shuffle(permIndex);
		int subsetNum = (int)(docStringsAll.size() * dataPercentage);
		List<Integer> subIndices = permIndex.subList(0, subsetNum);
		 
		docStrings = new ArrayList<String>();
		docURIs = new ArrayList<String>();
		docCorporaName = new ArrayList<String>();
		docTime = new ArrayList<Integer>();
		
		corporaName2Index = new HashMap<String, Integer>();
        corporaIndex2Name = new HashMap<Integer, String>();
        corporaNumber = 0;
        
        maxTimestampNumber = 0;
        minTimestampNumber = Integer.MAX_VALUE;
        corporaTime2Index = new HashMap<String, Integer>();
        corporaIndex2Time = new HashMap<Integer, String>();
	        
        int docIndex = 0;
		 for (int i = 0; i < subIndices.size(); ++i) {
			 if (docIndex % 1000 == 0) {
             	System.out.println(">>>[LOG]: Mallet parse " + docIndex + " documents.");
             }           
             docIndex++;
			 String text = docStringsAll.get(subIndices.get(i));
			 Instance carrier;
			 carrier = instancePipe.instanceFrom(new Instance (text, 0, null, null));
			 int length = 0;
			 if (isUseTFVector == true) {
				 SparseVector sv = (SparseVector) carrier.getData();
				 int[] index = sv.getIndices();
				 double[] value = sv.getValues();
				 length = Math.min(index.length, value.length);
			 } else {
				 FeatureSequence fs = (FeatureSequence) carrier.getData();
				 int[] indices = fs.toFeatureIndexSequence();
				 length = indices.length;
			 }
			 
			 if (length > 5) {
				 tflist.add (carrier);
				 if (!corporaName2Index.containsKey(docCorporaNameAll.get(subIndices.get(i)))) {
					 corporaName2Index.put(docCorporaNameAll.get(subIndices.get(i)), corporaNumber);
					 corporaIndex2Name.put(corporaNumber, docCorporaNameAll.get(subIndices.get(i)));
					 corporaNumber++;
				 }
				 
				 docStrings.add(docStringsAll.get(subIndices.get(i)));
				 docURIs.add(docURIsAll.get(subIndices.get(i)));
				 docCorporaName.add(docCorporaNameAll.get(subIndices.get(i)));
				 docTime.add(docTimeAll.get(subIndices.get(i)));
				 int year = docTimeAll.get(subIndices.get(i));
				 if (maxTimestampNumber < year) {
					 maxTimestampNumber = year;
				 }
				 if (minTimestampNumber > year) {
					 minTimestampNumber = year;
				 }
			 }
		 }
		 
		 for (int i = minTimestampNumber; i < (maxTimestampNumber + 1); ++i){
			 int year = i;
			 corporaTime2Index.put("" + year, year - minTimestampNumber);
			 corporaIndex2Time.put(year - minTimestampNumber, "" + year);
		 }
		 return tflist;
	}
	
	public void loadFeatureVector(double dataPercentage, File stopwordFile) {
		
		InstanceList tflist = loadInstanceSubsetList (dataPercentage, stopwordFile, true);
        alphabet = tflist.getAlphabet();
        tfVectors = new ArrayList<List<List<OpenIntIntHashMap>>>(); 
        for (int i = 0; i < (maxTimestampNumber - minTimestampNumber + 1); ++i) {
        	tfVectors.add(new ArrayList<List<OpenIntIntHashMap>>());
        	for (int j = 0; j < this.corporaNumber; ++j) {
        		tfVectors.get(i).add(new ArrayList<OpenIntIntHashMap>());
        	}
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
            int time = docTime.get(i);
            tfVectors.get(time - this.minTimestampNumber).get(corporaIndex).add(vector);
//            carrier.unLock();
//            carrier.clearSource();
//            carrier.setData(null);
        }
        System.out.println(">>>[LOG]: Doc number: " + tflist.size());
        System.out.println(">>>[LOG]: Word number: " + totalWord);
        System.out.println(">>>[LOG]: Vocabulary size: " + forwardHash.size());
	}
	
	public void loadFeatureSequence(double dataPercentage, File stopwordFile) {
		
		InstanceList tflist = loadInstanceSubsetList (dataPercentage, stopwordFile, false);
        alphabet = tflist.getAlphabet();
        featureSequences = new ArrayList<List<List<int[]>>>(); 
        for (int i = 0; i < (maxTimestampNumber - minTimestampNumber + 1); ++i) {
        	featureSequences.add(new ArrayList<List<int[]>>());
        	for (int j = 0; j < this.corporaNumber; ++j) {
        		featureSequences.get(i).add(new ArrayList<int[]>());
        	}
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
            int time = docTime.get(i);
            featureSequences.get(time - this.minTimestampNumber).get(corporaIndex).add(indices);
        }
        System.out.println(">>>[LOG]: Doc number: " + tflist.size());
        System.out.println(">>>[LOG]: Word number: " + totalWord);
        System.out.println(">>>[LOG]: Vocabulary size: " + alphabet.size());
	}
	
	public DocTimeCorpora initializeData (
			String inputDirectory, 
			String queryStr, 
			boolean isUseTFVector,
			double dataPercentage, 
			double subSampleRate, 
			File stopwordFile) {
    	loadDataString (inputDirectory, queryStr, dataPercentage);
    	
    	DocTimeCorpora docTimeCorpora = null;
    	
    	if (isUseTFVector) {
    		loadFeatureVector(subSampleRate, stopwordFile);
    		docTimeCorpora = new DocTimeCorpora(alphabet, forwardHash, backwardHash,
        			corporaIndex2Name, corporaName2Index, corporaTime2Index, corporaIndex2Time,
        			null, tfVectors, null, null, "Cobra Bank Data Analysis");
    	} else {
    		loadFeatureSequence(subSampleRate, stopwordFile);
    		docTimeCorpora = new DocTimeCorpora(alphabet, forwardHash, backwardHash,
        			corporaIndex2Name, corporaName2Index, corporaTime2Index, corporaIndex2Time,
        			featureSequences, null, null, null, "Cobra Bank Data Analysis");
    	}
    	return docTimeCorpora;
    }
	
}
