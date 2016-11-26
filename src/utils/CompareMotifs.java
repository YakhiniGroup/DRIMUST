package utils;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CompareMotifs {
	
	public static int countMatches(String s, String t)
	{
		int matches = 0;
		if(s.length() != t.length())
		{
			System.err.println("s and t must have the same length!");
		}
		for(int i=0; i<s.length(); i++)
		{
			char c1 = s.charAt(i);
			char c2 = t.charAt(i);
			if(c1 == c2)
			{
				matches++;
			}
			else
			{
				if(c1 == 'R')
				{
					if(c2 == 'G' || c2 == 'A') matches++;
				}
				if(c2 == 'R')
				{
					if(c1 == 'G' || c1 == 'A') matches++;
				}
				if(c1 == 'Y')
				{
					if(c2 == 'T' || c2 == 'C') matches++;
				}
				if(c2 == 'Y')
				{
					if(c1 == 'T' || c1 == 'C') matches++;
				}
				if(c1 == 'K')
				{
					if(c2 == 'G' || c2 == 'T') matches++;
				}
				if(c2 == 'K')
				{
					if(c1 == 'G' || c1 == 'T') matches++;
				}
				if(c1 == 'M')
				{
					if(c2 == 'C' || c2 == 'A') matches++;
				}
				if(c2 == 'M')
				{
					if(c1 == 'C' || c1 == 'A') matches++;
				}
				if(c1 == 'S')
				{
					if(c2 == 'C' || c2 == 'G') matches++;
				}
				if(c2 == 'S')
				{
					if(c1 == 'C' || c1 == 'G') matches++;
				}
				if(c1 == 'W')
				{
					if(c2 == 'A' || c2 == 'T') matches++;
				}
				if(c2 == 'W')
				{
					if(c1 == 'A' || c1 == 'T') matches++;
				}
				if(c1 == 'B' && c2 != 'A') matches++;
				if(c2 == 'B' && c1 != 'A') matches++;
				if(c1 == 'D' && c2 != 'C') matches++;
				if(c2 == 'D' && c1 != 'C') matches++;
				if(c1 == 'H' && c2 != 'G') matches++;
				if(c2 == 'H' && c1 != 'G') matches++;
				if(c1 == 'V' && c2 != 'T') matches++;
				if(c2 == 'V' && c1 != 'T') matches++;
			}
		}
		return matches;
	}

	public static int countMaxMatches(String s, String t)
	{
		int maxMatches = 0;
		
		String longer = s, shorter = t;
		if(t.length() > s.length())
		{
			longer = t;
			shorter = s;
		}
		
		for(int i=0; i<longer.length(); i++)
		{
			int matches = 0;
			if(i < shorter.length())
			{
				matches = countMatches(longer.substring(0, i+1), shorter.substring(shorter.length()-i-1));
			}
			else
			{
				matches = countMatches(longer.substring(i-shorter.length()+1, i+1), shorter);
			}
			if(matches > maxMatches)
			{
				maxMatches = matches;
			}
		}
		for(int j=shorter.length()-1; j>0 ; j--)
		{
			int matches = countMatches(longer.substring(longer.length()-j), shorter.substring(0,j));
			if(matches > maxMatches)
			{
				maxMatches = matches;
			}
		}
		
		return maxMatches;
	}
	
	public static Map<String,List<String>> loadKnownMotifs(String fileName) throws IOException
	{
		Map<String,List<String>> motifs = new HashMap<String, List<String>>();
		
		FileInputStream fstream = new FileInputStream(fileName);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		strLine = br.readLine();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.isEmpty())
        	{
        		continue;
        	}
        	String[] columns = strLine.split("\t");
        	String tf = columns[0].toUpperCase();
        	String m1 = columns[1].toUpperCase();
        	String m2 = columns[2].toUpperCase();
        	List<String> l = new ArrayList<String>();
        	l.add(m1); l.add(m2);
        	motifs.put(tf, l);
        }
        return motifs;
	}
	
	public static Map<String,List<String>> loadKnownMotifsEran(String fileName) throws IOException
	{
		Map<String,List<String>> motifs = new HashMap<String, List<String>>();
		
		FileInputStream fstream = new FileInputStream(fileName);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		strLine = br.readLine();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.isEmpty())
        	{
        		continue;
        	}
        	String[] columns = strLine.split("\t");
        	String tf = columns[0].toUpperCase();
        	tf = tf.split("_")[0];
        	String m = columns[1].toUpperCase();
        	if(motifs.containsKey(tf))
        	{
        		motifs.get(tf).add(m);
        	}
        	else
        	{
        		List<String> l = new ArrayList<String>();
            	l.add(m);
            	motifs.put(tf, l);
        	}
        }
        return motifs;
	}
	
	public static Map<String,List<String>> loadOurMotifs(String fileName) throws IOException
	{
		Map<String,List<String>> motifs = new HashMap<String, List<String>>();
		
		FileInputStream fstream = new FileInputStream(fileName);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
//		strLine = br.readLine();
		Map<String,Integer> cnt = new HashMap<String, Integer>();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.isEmpty())
        	{
        		continue;
        	}
        	String[] columns = strLine.split("\t");
        	String tf = columns[0].toUpperCase();
        	if(cnt.containsKey(tf))
        	{
        		cnt.put(tf, cnt.get(tf) + 1);
        	}
        	else
        	{
        		cnt.put(tf,1);
        	}
        	if(cnt.get(tf) > 2) continue;
        	if(columns[3].startsWith("["))
        	{
        		String gaps = columns[3];
        		gaps = gaps.substring(1,gaps.length()-1); //remove the brackets
        		String[] lengths = gaps.split(",");
        		String h1 = columns[1].toUpperCase();
            	String h2 = columns[2].toUpperCase();
//            	System.out.println(columns[0] + "\t" + h1 + "-" + columns[3] + "-" + h2 + "\t" + columns[5] + 
//            			"\t" + columns[7] + "\t" + columns[9] + "\t" + columns[11] + "\t" + columns[15]);
        		for(String len : lengths)
        		{
        			int k = Integer.parseInt(len.trim());
        			String spacer = "";
        			for(int i=0; i<k; i++)
        			{
        				spacer = spacer + "N";
        			}
        			String motif = h1 + spacer + h2;
        			if(motifs.containsKey(tf))
        			{
        				motifs.get(tf).add(motif);
        			}
        			else
        			{
        				List<String> l = new ArrayList<String>();
        	        	l.add(motif);
        	        	motifs.put(tf, l);
        			}
        		}
        	}
        	else
        	{
        		String motif = columns[1];
//        		System.out.println(columns[0] + "\t" + columns[1] + "\t" + columns[5] + 
//            			"\t" + columns[7] + "\t" + columns[9] + "\t" + columns[11] + "\t" + columns[15]);
        		if(motifs.containsKey(tf))
    			{
    				motifs.get(tf).add(motif);
    			}
    			else
    			{
    				List<String> l = new ArrayList<String>();
    	        	l.add(motif);
    	        	motifs.put(tf, l);
    			}
        		
        	}
        }
        return motifs;
	}
	
	public static void printOurBestMotifs(String fileName) throws IOException
	{		
		PrintWriter out = new PrintWriter("our_best.txt");
		FileInputStream fstream = new FileInputStream(fileName);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
//		strLine = br.readLine();
		Map<String,Integer> cnt = new HashMap<String, Integer>();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.isEmpty())
        	{
        		continue;
        	}
        	String[] columns = strLine.split("\t");
        	String tf = columns[0].toUpperCase();
        	if(cnt.containsKey(tf))
        	{
        		cnt.put(tf, cnt.get(tf) + 1);
        	}
        	else
        	{
        		cnt.put(tf,1);
        	}
        	if(cnt.get(tf) > 2) continue;
        	out.println(strLine);
        }
        out.close();
	}
	
	public static void main(String[] argv) throws IOException
	{
//		Map<String,List<String>> known = loadKnownMotifs("C:\\WORK\\thesis\\motif search\\results\\harbison\\known_motifs.txt");
		Map<String,List<String>> known = loadKnownMotifsEran("C:\\WORK\\thesis\\motif search\\results\\harbison\\erans_motifs.txt");
		System.out.println("Eran found motifs for " + known.keySet().size() + " tfs");
		Map<String,List<String>> our = loadOurMotifs("C:\\WORK\\thesis\\motif search\\results\\harbison\\our_motifs.txt");
		printOurBestMotifs("C:\\WORK\\thesis\\motif search\\results\\harbison\\our_motifs.txt");
		for(String tf : our.keySet())
		{
			double best = 0;
			String argBest = "";
			if(known.containsKey(tf))
			{
				for(String k : known.get(tf))
				{
					if(k.isEmpty()) continue;
					for(String o : our.get(tf))
					{
						int matches = countMaxMatches(k, o);
						int minLength = k.length() < o.length() ? k.length() : o.length();
						double rate = (double)matches/(double)minLength;
//						if(tf.equals("ADR1")) System.out.println(k + " " + o + " " + matches + " " + minLength + " " + rate);
						if(rate > best)
						{
							best = rate;
							argBest = k + " " + o;
						}
					}
				}
				int mark = best >= 0.8 ? 1 : 3;
				if(best >= 0.5 && best < 0.8) mark = 2;
				System.out.println(tf + "\t" + best + "\t" + mark + "\t" + argBest);
			}
			else
			{
				System.out.println(tf + "\t" + "0" + "\t0");
			}
		}
	}

}
