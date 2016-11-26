package motif.search;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import clustering.ClusterManager;
import clustering.MotifSet;

import statistics.HGScore;
import utils.ClusterInfo;
import utils.StringOperations;

public class MotifSearchTester {
		
	public static void validate(String motif, List<String> sequences, TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ, boolean dsMode)
	{
		char[] strands;
		if(dsMode)
		{
			char[] a = {'+','-'}; strands = a;
		}
		else
		{
			char[] a = {'+'}; strands = a;
		}
		for(int i : occ.keySet())
		{
			for(char strand : occ.get(i).keySet())
			{
				SortedSet<Integer> ind = occ.get(i).get(strand);
				String seq = sequences.get(i-1);
				if(strand == '-')
				{
					seq = StringOperations.reverseComplement(seq);
				}
				int j = 0;
				int bingo = 0;
				while(seq.substring(j).contains(motif))
				{
					j = seq.indexOf(motif, j);
					if(ind.contains(j))
					{
						bingo++;
					}
					j++;
				}
				if(bingo != ind.size())
				{
					System.err.println(motif + " " + seq + " " + bingo + " " + ind.size() + " " + i + " strand " + strand);
				}
			}
		}
		List<String> seen = new ArrayList<String>();
//		if(motif.equalsIgnoreCase("GCTGAGG")) System.out.println("motif " + occ);
		for(int i=0; i<sequences.size(); i++)
		{
			for(char strand : strands)
			{
				String seq = sequences.get(i);
				if(strand == '-')
				{
					seq = StringOperations.reverseComplement(seq);
				}
				if(seq.contains(motif) && false == occ.containsKey(i+1))
				{
					seen.add(seq);
					System.err.println(motif + " appears in " + seq + " but was not considered " + i + " in " + seq.indexOf(motif) + " strand " + strand);
	//				if(seen.contains(seq)) System.err.println("seen");
				}
			}
		}
	}
	
	public static void outputClusters(String fileName, List<MotifsMapEntry> motiflist,
			List<String> sequences, double threshold, int maxNumMotifs) throws IOException, InterruptedException
	{
		PrintWriter output1 = new PrintWriter(fileName);
		PrintWriter output2 = new PrintWriter(fileName + "_alignment.details.txt");
//		Set<Character> alphabet = new HashSet<Character>();
		SortedSet<MotifSet> clusters = ClusterManager.cluster(asArrayList(motiflist), threshold, maxNumMotifs);
		for(MotifSet set : clusters)
		{
			List<List<String>> mos = ClusterManager.getOccurStringsLong(set, sequences, set.getmHGScore().n);
			if(ClusterManager.chooseAlignerAccToPercentile(mos, 0.8, 0.667) < 0)
			{
				mos = ClusterManager.getOccurStringsShort(set, sequences, set.getmHGScore().n);
			}
			else
			{
				System.out.println("long alignment won " + mos.get(0).get(0));
			}
	        List<List<String>> msa = ClusterManager.alignCluster(mos, output2);
			List<String> motifOutput = new ArrayList<String>();
			Collections.sort(set.getSet());
			String title = "motif\tmHG_score\tcorrected_score\tN\tB\tn\tb\tEnrichment\tOccurrences";
			motifOutput.add(title);
			for(MotifsMapEntry me : set.getSet())
			{
				String motif = me.getMotif();
				HGScore s = me.getmHG();
				double enrichment = (double)((double)s.b/s.n)/(double)((double)s.B/s.N);
				String info = motif + "\t" + s.score + "\t" + s.B*s.score + "\t" + s.N + "\t" + s.B + "\t" + s.n + "\t" + s.b + "\t" + enrichment + "\t";
				for(int j : me.occ.keySet())
				{
					info = info + "(" + (j-1) + ":" + me.occ.get(j).get('+') + ") ";
				}
				motifOutput.add(info);
			}
//			Map<Integer,Map<Character,Integer>> pssm = buildPSSM(alphabet, msa);
//			List<String> pssmPrint = pssmPrint(pssm, alphabet, msa.size());
			//now print all info on cluster
			output1.println("cluster mHG score = " + set.getmHGScore() + " p " + set.getp() + " cluster size = " + set.getSet().size());
			for(List<String> l : msa)
			{
				output1.println(l.get(0));
			}
			output1.println();
			for(String st : motifOutput)
			{
				output1.println(st);
			}
			output1.println();
//			for(String st : pssmPrint)
//			{
//				output.println(st);
//			}
			output1.println();
			output1.println("**********************************************");
		}
		output1.close();
		output2.close();
	}
	
	public static Map<Integer,Map<Character,Integer>> buildPSSM(Set<Character> alphabet, List<String> msa)
	{
		Map<Integer,Map<Character,Integer>> pssm = new HashMap<Integer,Map<Character,Integer>>();
		for(String algn : msa)
		{
			for(int i=0; i<algn.length(); i++)
			{
				char ch = algn.charAt(i);
				if(ch != '-' && false == alphabet.contains(ch))
				{
					alphabet.add(ch);
				}
				if(pssm.containsKey(i))
				{
					if(ch != '-')
					{
						if(pssm.get(i).containsKey(ch))
						{
							pssm.get(i).put(ch, pssm.get(i).get(ch)+1);
						}
						else
						{
							pssm.get(i).put(ch, 1);
						}
					}
				}
				else
				{
					if(ch != '-')
					{
						Map<Character,Integer> m = new HashMap<Character, Integer>();
						m.put(ch, 1);
						pssm.put(i, m);
					}
				}
			}
		}
		return pssm;
	}
	
	public static List<String> pssmPrint(Map<Integer,Map<Character,Integer>> pssm, Set<Character> alphabet, int n)
	{
		List<String> pssmPrint = new ArrayList<String>();
		String title = "\t";
		for(int i=0; i<pssm.keySet().size(); i++)
		{
			title = title + i + "\t";
		}
		pssmPrint.add(title);
		for(char ch : alphabet)
		{
			String strLine = ch + "\t";
			for(int i=0; i<pssm.keySet().size(); i++)
			{
				if(pssm.containsKey(i))
				{
					if(pssm.get(i).containsKey(ch))
					{
						strLine = strLine + (double)pssm.get(i).get(ch)/(double)n + "\t";
					}
					else
					{
						strLine = strLine + "0\t";
					}
				}
				else
				{
					
					System.err.println("very odd...");
				}
			}
			pssmPrint.add(strLine);
		}
		return pssmPrint;
	}
	
	private static List<MotifsMapEntry> filterMotifs(List<MotifsMapEntry> motifs)
	{
		List<MotifsMapEntry> filtered = new ArrayList<MotifsMapEntry>();
		Map<String,MotifsMapEntry> representativeMap = new HashMap<String, MotifsMapEntry>();
		List<List<MotifsMapEntry>> rcPairs = new ArrayList<List<MotifsMapEntry>>();
		for(MotifsMapEntry me : motifs)
		{
			String rcMotif = StringOperations.reverseComplement(me.getMotif());
			if(representativeMap.containsKey(rcMotif))
			{
				List<MotifsMapEntry> pair = new ArrayList<MotifsMapEntry>();
				pair.add(me);
				pair.add(representativeMap.get(rcMotif));
				rcPairs.add(pair);
			}
			else
			{
				representativeMap.put(me.getMotif(), me);
			}
		}
		//out of each pair of reverse complement motifs, we take the member that has more occurrences on the plus strand in the target set (i.e. top n sequences)
		for(List<MotifsMapEntry> pair : rcPairs)
		{
			MotifsMapEntry me1 = pair.get(0);
			TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ = me1.getOcc();
			HGScore score = me1.getmHG();
			Map<Character,Integer> cnt = new HashMap<Character, Integer>();
			cnt.put('+', 0);
			cnt.put('-', 0);
			for(int ind : occ.keySet())
			{
				if(ind >= score.n) break;
				Map<Character,SortedSet<Integer>> m = occ.get(ind);
				for(char strand : m.keySet())
				{
					cnt.put(strand, cnt.get(strand)+m.get(strand).size());
				}
			}
			if(cnt.get('+') >= cnt.get('-'))
			{
				filtered.add(me1);
			}
			else
			{
				filtered.add(pair.get(1));
			}
		}
		return filtered;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException
	{	
		int k1 = 0, k2 = 0;
		boolean dsMode = false;
		float stringency = 0;
		FileInputStream fstream = null;
		PrintWriter out1 = null;
		String fileName2 = "";
		
		if (args.length == 7) 
		{
		    k1 = Integer.parseInt(args[0]);
		    k2 = Integer.parseInt(args[1]);
		    if(args[2].equalsIgnoreCase("ds")) dsMode = true;
		    stringency = Float.parseFloat(args[3]);
		    fstream = new FileInputStream(args[4]);
		    out1 = new PrintWriter(args[5]);
		    fileName2 = args[6];
		}
		else
		{
			System.err.println("You must give seven arguments as input");
	        System.exit(1);
		}
		if(dsMode)
		{
			DSMotifSearchTester.dsSearch(k1, k2, stringency, fstream, out1, fileName2);
			return;
		}
//		FileInputStream fstream = new FileInputStream("C:\\WORK\\thesis\\motif search\\scripts\\pufs\\puf5targets.txt");
		DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        String strLine;
        int N = 0;
        int totalLength = 0;
        List<String> sequences = new ArrayList<String>();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.trim().isEmpty() || strLine.contains("\"\""))
        	{
        		continue;
        	}
        	if(strLine.startsWith(">"))
        	{
        		continue;
        	}
        	N++;
        	sequences.add(strLine);
        	totalLength += strLine.length();
        }
        System.out.println("Testing " + sequences.size()  + " sequences, total length = " + totalLength);
        long before = System.currentTimeMillis();
        
        SearchManager m = new SearchManager(sequences, k1, k2, dsMode);
        m.setThreshold(stringency);
        
        System.out.println("Time to construct the tree: " + ((float)(System.currentTimeMillis()-before)/60000)  + " minutes");
        List<MotifsMapEntry> motifs = Arrays.asList(m.doSearch().values().toArray(new MotifsMapEntry[0]));
        Collections.sort(motifs);
        long after = System.currentTimeMillis();
        System.out.println("Total time: " + ((float)(after-before)/60000)  + " minutes");
		Map<String, MotifsMapEntry> motifInfo = new HashMap<String, MotifsMapEntry>();
		List<MotifsMapEntry> motiflist = new ArrayList<MotifsMapEntry>();
		out1.print("motif\tmHG_score\tcorrected_score\tN\tB\tn\tb\tEnrichment\tOccurrences\n");
        for(MotifsMapEntry me : motifs)
		{
			String motif = me.getMotif();
//			validate(motif, sequences, me.occ, dsMode);
			HGScore s = me.getmHG();
			double enrichment = (double)((double)s.b/s.n)/(double)((double)s.B/s.N);
			out1.print(motif + "\t" + s.score + "\t" + s.B*s.score + "\t" + s.N + "\t" + s.B + "\t" + s.n + "\t" + s.b + "\t" + enrichment + "\t");
			for(int j : me.occ.keySet())
			{
				out1.print("(" + (j-1) + ":" + me.occ.get(j).get('+') + ") ");
			}
			out1.print("\n");
			motifInfo.put(motif, me);
			motiflist.add(me);
		}
		out1.flush();
		out1.close();
		outputClusters(fileName2, motiflist, sequences, stringency, 3);
		long after1 = System.currentTimeMillis();
		System.out.println("Total time of fbfs clustering: " + ((float)(after1-after)/60000)  + " minutes");
	}

	public static List<MotifsMapEntry> asArrayList(List<MotifsMapEntry> l)
	{
		List<MotifsMapEntry> l0 = new ArrayList<MotifsMapEntry>();
		for(MotifsMapEntry e : l)
		{
			l0.add(e);
		}
		return l0;
	}
	
}