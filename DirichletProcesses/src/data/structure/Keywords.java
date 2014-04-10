package data.structure;

public class Keywords{

	String[] words;
	int[] freqs;
	public Keywords(String[] words, int[] freqs) {
		// TODO Auto-generated constructor stub
		this.words = words;
		this.freqs = freqs;
	}
	public String[] getWords(){
		return this.words;
	}
	public int[] getFreqs(){
		return this.freqs;
	}
	public int getN(){
		return this.words.length;
	}
	
}