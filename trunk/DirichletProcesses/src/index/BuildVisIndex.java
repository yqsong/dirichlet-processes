package index;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.index.IndexWriter;

/**
 * @author Yangqiu Song
 */
public class BuildVisIndex {

	public void BuildIndex (String filename, String indexDirStr, String typeName, boolean create) {
		File indexFileDir = new File(indexDirStr);
		int docNum = 0;
		int lineID = 0;
		
		IndexWriter luceneWriter = null;
		FileInputStream fin;
		try {
			luceneWriter = new IndexWriter(indexFileDir, new StandardAnalyzer(), create, IndexWriter.MaxFieldLength.UNLIMITED);
			System.out.println("Indexing to directory '" +indexFileDir+ "'...");
			fin = new FileInputStream (filename);
			BufferedReader buf = new BufferedReader(new InputStreamReader(fin));
			String str = buf.readLine();
			while (str != null) {
				docNum++;
				System.out.println("Doc Number: " + docNum);
				
				Document doc = new Document();
				
				String[] tokens = str.split(":");
				if (!tokens[0].trim().equalsIgnoreCase("year")) {
					System.err.println("Error in reading year! " + tokens[0]);
				}
				doc.add(new Field("year", tokens[1].trim(), Field.Store.YES, Field.Index.ANALYZED));
		
				str = buf.readLine();
				tokens = str.split(":");
				if (!tokens[0].trim().equalsIgnoreCase("title")) {
					System.err.println("Error in reading title! " + tokens[0]);
				}
				doc.add(new Field("title", tokens[1].trim(), Field.Store.YES, Field.Index.ANALYZED));

				str = buf.readLine();
				tokens = str.split(":");
				if (!tokens[0].trim().equalsIgnoreCase("link")) {
					System.err.println("Error in reading link! " + tokens[0]);
				}
				String str1 = "";
				for (int i = 1; i < tokens.length; ++i) {
					str1 += tokens[i].trim() + "_";
				}
				doc.add(new Field("link", str, Field.Store.YES, Field.Index.NOT_ANALYZED));
				doc.add(new Field("uri",str1, Field.Store.YES, Field.Index.NOT_ANALYZED));

				str = buf.readLine();
				tokens = str.split(":");
				if (!tokens[0].trim().equalsIgnoreCase("abstract")) {
					System.err.println("Error in reading abstract! " + tokens[0]);
				}
				String abs = "";
				for (int i = 1; i < tokens.length; ++i) {
					abs += tokens[i];
				}
				doc.add(new Field("abstract", abs.trim(), Field.Store.YES, Field.Index.ANALYZED));
				
				str = buf.readLine();
				tokens = str.split(":");
				if (!tokens[0].trim().equalsIgnoreCase("keywords")) {
					System.err.println("Error in reading keywords! " + tokens[0]);
				}
				String keywords = "";
				for (int i = 1; i < tokens.length; ++i) {
					keywords += tokens[i];
				}
				doc.add(new Field("keywords", keywords.trim(), Field.Store.YES, Field.Index.ANALYZED));
				
				str = buf.readLine();
				tokens = str.split(":");
				if (!tokens[0].trim().equalsIgnoreCase("author")) {
					System.err.println("Error in reading author! " + tokens[0]);
				}
				String author = "";
				for (int i = 1; i < tokens.length; ++i) {
					author += tokens[i];
				}
				doc.add(new Field("author", author.trim(), Field.Store.YES, Field.Index.ANALYZED));

				doc.add(new Field("confName", typeName, Field.Store.YES, Field.Index.ANALYZED));
				
				luceneWriter.addDocument(doc);
				
				str = buf.readLine();
				str = buf.readLine();
				lineID += 5;
			}
			buf.close();
		    fin.close();
		    luceneWriter.optimize();
		    luceneWriter.close();
		    System.out.println("Indexing to directory done!" + "\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}
	
	public static void main(String[] args) {
		String textInfoVis = "data/INFOVIS-abstract.txt";
		String textVis =     "data/VIS-abstract.txt";
		String luceneSrc = "data/visinfovis";
		
		String indexDirStr = luceneSrc + "/text_index_2";
		File indexFileDir = new File(indexDirStr);
		if (indexFileDir.exists()) {
			File[] files = indexFileDir.listFiles();
			for (int i = 0; i < files.length; ++i) {
				files[i].delete();
			}
		}
		
		BuildVisIndex biVis = new BuildVisIndex();
		biVis.BuildIndex(textInfoVis, indexDirStr, "infovis", true);
		biVis.BuildIndex(textVis, indexDirStr, "vis", false);
		
	}
	
}
