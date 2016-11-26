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
import clustering.DSClusterManager;
import clustering.DSMotifSet;
import clustering.MotifSet;

import statistics.HGScore;
import utils.ClusterInfo;
import utils.DSMotif;
import utils.StringOperations;

public class DSMotifSearchTester {
	
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
	
	public static void outputClusters(String fileName, List<DSMotif> motifs,
			List<String> sequences, double threshold, int maxNumMotifs) throws IOException, InterruptedException
	{
		PrintWriter output1 = new PrintWriter(fileName);
		PrintWriter output2 = new PrintWriter(fileName + "_alignment.details.txt");
//		Set<Character> alphabet = new HashSet<Character>();
		SortedSet<DSMotifSet> clusters = DSClusterManager.cluster(motifs, threshold, maxNumMotifs);
		for(DSMotifSet set : clusters)
		{
			List<List<String>> mos = DSClusterManager.getOccurStringsLong(set, sequences, set.getmHGScore().n);
			if(ClusterManager.chooseAlignerAccToPercentile(mos, 0.8, 0.667) < 0)
			{
				mos = DSClusterManager.getOccurStringsShort(set, sequences, set.getmHGScore().n);
			}
			else
			{
				System.out.println("long alignment won " + mos.get(0).get(0));
			}
	        List<List<String>> msa = DSClusterManager.alignCluster(mos, output2);
			List<String> motifOutput = new ArrayList<String>();
			Collections.sort(set.getSet());
			String title = "motif\treverse_complement_motif\tmHG_score\tcorrected_score\tN\tB\tn\tb\tEnrichment\tOccurrences";
			motifOutput.add(title);
			for(DSMotif me : set.getSet())
			{
				HGScore s = me.getScore();
				double enrichment = (double)((double)s.b/s.n)/(double)((double)s.B/s.N);
				String info = me.getMotif1() + "\t" + me.getMotif2() + "\t" + s.score + "\t" + s.B*s.score + "\t" + s.N + "\t" + s.B + "\t" + s.n + "\t" + s.b + "\t" + enrichment + "\t";
				for(int j : me.getOcc().keySet())
				{
					info = info + "(" + (j-1) + ":";
					for(char strand : me.getOcc().get(j).keySet())
					{
						SortedSet<Integer> occs = me.getOcc().get(j).get(strand);
						info = info + "{" + strand + ":" + occs + "}";
					}
					info = info + ") ";
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
	
	public static void dsSearch(int k1, int k2, float stringency, FileInputStream fstream, PrintWriter out1, String fileName) throws IOException, InterruptedException
	{	
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
        
        DSSearchManager m = new DSSearchManager(sequences, k1, k2);
        m.setThreshold(stringency);
        
        System.out.println("Time to construct the tree: " + ((float)(System.currentTimeMillis()-before)/60000)  + " minutes");
        List<DSMotif> motifs = m.doSearch();
		Collections.sort(motifs);
		List<Integer> sequencesLengths = m.getSequencesLengths();
        long after = System.currentTimeMillis();
        System.out.println("Total time: " + ((float)(after-before)/60000)  + " minutes");
		
		out1.print("motif\treverse_complement_motif\tmHG_score\tcorrected_score\tN\tB\tn\tb\tEnrichment\tOccurrences\n");
		
        for(DSMotif me : motifs)
		{
			HGScore s = me.getScore();
			double enrichment = (double)((double)s.b/s.n)/(double)((double)s.B/s.N);
			
			out1.print(me.getMotif1() + "\t" + me.getMotif2() + "\t" + s.score + "\t" + s.B*s.score + "\t" + s.N + "\t" + s.B + "\t" + s.n + "\t" + s.b + "\t" + enrichment + "\t");
			
			for(int j : me.getOcc().keySet())
			{
				out1.print("(" + (j-1) + ":");
				for(char strand : me.getOcc().get(j).keySet())
				{
					SortedSet<Integer> occs = me.getOcc().get(j).get(strand);
					out1.print("{" + strand + ":" + occs + "}");
				}
				out1.print(") ");
			}
			out1.print("\n");
		}
		out1.close();
		
		outputClusters(fileName, motifs, sequences, stringency, 3);
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