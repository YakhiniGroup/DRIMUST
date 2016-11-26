package motif.search;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import statistics.HGScore;
import statistics.mHG1;

public class MotifSearchValidator {
	
	public static List<String> getAllSequences(int l)
	{
		List<String> all = new ArrayList<String>();
		if(l==0)
		{
			all.add("");
			return all;
		}
		char[] letters = {'A','C','T','G'};
		List<String> almost = getAllSequences(l-1);
		for(String s : almost)
		{
			for(char c : letters)
			{
				all.add(s+c);
			}
		}
		return all;
	}
	
	public static HGScore getEnrichment(List<String> sequences, String m)
	{
		int N = sequences.size();
		int B = 0;
		int[] vec = new int[N];
		int i = 0;
		for(String s : sequences)
		{
			if(s.contains(m))
			{
				vec[i++] = 1;
				B++;
			}
			else
			{
				vec[i++] = 0;
			}
		}
		return mHG1.calculate_mHGT(vec, B);
	}
	
	public static void main(String[] args) throws IOException
	{			
//		FileInputStream fstream = new FileInputStream("C:\\WORK\\thesis\\RBPs\\TargetScan\\scripts\\Galgano_PUM1_targets.txt");//"C:\\WORK\\thesis\\RBPs\\puf\\protein\\yeast_proteome_ranked_puf2.txt"); //"C:\\DRIM\\puf2\\puf2targets.txt"
		FileInputStream fstream = new FileInputStream("C:\\DRIM\\puf2\\puf2targets.txt");
//		FileInputStream fstream = new FileInputStream("C:\\WORK\\thesis\\Meromit\\scripts\\sorted_sequences.txt");
		DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        String strLine;
        int N = 0;
        int totalLength = 0;
        List<String> sequences = new ArrayList<String>();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.isEmpty() || strLine.startsWith(">"))
        	{
        		continue;
        	}
        	N++;
        	{
        		sequences.add(strLine);
        		totalLength += strLine.length();
        	}
        }
        List<String> motifs = getAllSequences(7);
        PrintWriter out = new PrintWriter("enriched_motifs_validator.txt");
        for(String m : motifs)
        {
        	HGScore s = getEnrichment(sequences, m);
        	if(s.score <= 0.0001)
        	{
        		out.println(m + "\tN\t" + s.N + "\tB\t" + s.B + "\tn\t" + s.n + "\tb\t" + s.b + "\tmHG_score\t" + s.score);
        	}
        }
		out.flush();
		out.close();
	}
}
