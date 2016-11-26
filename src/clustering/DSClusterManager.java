package clustering;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import motif.search.Constants;
import motif.search.MotifsMapEntry;

import utils.DSMotif;
import utils.StringOperations;

public class DSClusterManager extends ClusterManager{
		
	public static DSMotifSet fbfs(List<DSMotif> motifs, double threshold)
	{
		DSMotifSet set = new DSMotifSet(motifs.get(0));
		DSMotif argbest = motifs.get(0);
		motifs.remove(0);
		int cnt = 1;
		while(true)
		{
			boolean flag = false;
			double bestp = set.getp();
			DSMotifSet bestset = set;
			boolean toSwapBest = false;
			for(DSMotif m : motifs)
			{
				boolean toSwap = false;
				double agreementRate1 = calculateMeanClustersMatch(set.getSet(), m.getMotif1());
				double agreementRate2 = calculateMeanClustersMatch(set.getSet(), m.getMotif2());
				if(agreementRate1 < Constants.SIMILARITY && agreementRate2 < Constants.SIMILARITY) continue;
				if(agreementRate2 > agreementRate1)
				{
					toSwap = true;
				}
				DSMotifSet newset = DSMotifSet.addMotif(m, set, threshold);
				double newp = newset.getp();
				if(newp < bestp)
				{
					bestp = newp;
					bestset = newset;
					argbest = m;
					toSwapBest = toSwap;
				}
			}
			if(set != bestset)
			{
				set = bestset;
				if(toSwapBest) argbest.swapMotifs();
				motifs.remove(argbest);
				flag = true;
			}
			//here 'bestset' has to be equal to 'set'
			if(cnt == 5 || false == flag)
			{
				for(DSMotif m : set.getSet())
				{
					DSMotifSet newset = DSMotifSet.removeMotif(m, set, threshold);
					double newp = newset.getp();
					if(newp < bestp)
					{
						bestp = newp;
						bestset = newset;
						argbest = m;
					}
				}
				if(set != bestset)
				{
					set = bestset;
					motifs.add(argbest);
					flag = true;
				}
				//here, too, 'bestset' has to be equal to 'set'
			}
			cnt = (cnt == 5)? 1 : cnt+1;
			if(false == flag) //didn't add or removed any motif
			{
				break;
			}
		}
		return set;        
	}
	
	public static SortedSet<DSMotifSet> cluster(List<DSMotif> motifs, double threshold, int maxNumClusters)
	{
		SortedSet<DSMotifSet> clusters = new TreeSet<DSMotifSet>();
		while(false == motifs.isEmpty() && clusters.size() < maxNumClusters)
		{
			DSMotifSet set = fbfs(motifs, threshold);
			clusters.add(set);
			List<DSMotif> toDelete = new ArrayList<DSMotif>(); 
			for(DSMotif m : motifs)
			{
				double similarity1 = calculateMeanClustersMatch(set.getSet(), m.getMotif1());
				double similarity2 = calculateMeanClustersMatch(set.getSet(), m.getMotif2());
				if(similarity1 >= Constants.SIMILARITY-0.1 || similarity2 >= Constants.SIMILARITY-0.1)
				{
					toDelete.add(m);
				}
			}
			for(DSMotif m : toDelete)
			{
				motifs.remove(m);
			}
		}
		return clusters;        
	}
	
	/**
	 * s and t must be of the same length!
	 */
//	public static int countMatches(String s, String t)
//	{
//		int matches = 0;
//		if(s.length() != t.length())
//		{
//			System.err.println("s and t must have the same length!");
//		}
//		for(int i=0; i<s.length(); i++)
//		{
//			if(s.charAt(i) == t.charAt(i))
//			{
//				matches++;
//			}
//		}
//		return matches;
//	}
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return match b with a, such that b can be modified with additional '-' in its beginning
	 */
//	public static String[] match(String a, String b)
//	{
//		int l = a.length();
//		String c = b;
//		int maxMatches = 0;
//		String argMax = b;
//		for(int i=0; i<l; i++)
//		{
//			int matches = 0;
//			if(a.length() >= c.length())
//			{
//				matches = countMatches(a.substring(0,c.length()), c);
//			}
//			else
//			{
//				matches = countMatches(a, c.substring(0, a.length()));
//			}
//			if(matches > maxMatches)
//			{
//				maxMatches = matches;
//				argMax = c;
//			}
//			c = "-" + c;
//		}
//		String[] result = {argMax, Integer.toString(maxMatches)};
//		return result;
//	}
	
//	public static double countMaxMatches(String s, String t)
//	{
//		int maxMatches = 0;
//		
//		String longer = s, shorter = t;
//		if(t.length() > s.length())
//		{
//			longer = t;
//			shorter = s;
//		}
//		
//		for(int i=0; i<longer.length(); i++)
//		{
//			int matches = 0;
//			if(i < shorter.length())
//			{
//				matches = countMatches(longer.substring(0, i+1), shorter.substring(shorter.length()-i-1));
//			}
//			else
//			{
//				matches = countMatches(longer.substring(i-shorter.length()+1, i+1), shorter);
//			}
//			if(matches > maxMatches)
//			{
//				maxMatches = matches;
//			}
//		}
//		for(int j=shorter.length()-1; j>0 ; j--)
//		{
//			int matches = countMatches(longer.substring(longer.length()-j), shorter.substring(0,j));
//			if(matches > maxMatches)
//			{
//				maxMatches = matches;
//			}
//		}
//		
//		return (double)maxMatches/(double)shorter.length();
//	}
	
	public static double calculateMeanClustersMatch(List<DSMotif> cluster, String s)
	{
		double sumMatches = 0;
		int n = 0;
		for(DSMotif m : cluster)
		{
			String t = m.getMotif1();
			sumMatches += countMaxMatches(s, t);
			n++;
		}
		double mean = sumMatches/(double)n;
		
		return mean;
	}
	
//	public static double calculateMeanMatchOfClusters(List<MotifsMapEntry> c1, List<MotifsMapEntry> c2)
//	{
//		double sumMatches = 0;
//		int n = 0;
//		for(MotifsMapEntry m1 : c1)
//		{
//			for(MotifsMapEntry m2 : c2)
//			{
//				sumMatches += countMaxMatches(m1.getMotif(), m2.getMotif());
//				n++;
//			}
//		}
//		double mean = sumMatches/(double)n;
//		
//		return mean;
//	}
	
//	public static List<String> alignCluster(List<String> cluster)
//	{
//		int maxTotal = 0;
//		List<String> argMax = new ArrayList<String>();
//		for(int i=0; i<cluster.size(); i++)
//		{
//			String m1 = cluster.get(i);
//			int total = 0;
//			List<String> l = new ArrayList<String>();
//			l.add(m1);
//			for(int j=0; j<cluster.size(); j++)
//			{
//				String m2 = cluster.get(j);
//				if(i == j) continue;
//				String[] res = match(m1,m2);
//				total += Integer.parseInt(res[1]);
//				l.add(res[0]);
//			}
//			if(total > maxTotal || cluster.size() == 1)
//			{
//				maxTotal = total;
//				argMax = l;
//			}
//		}
//		//cosmetic phase - needed for presentation with WebLogo
//		int maxL = 0;
//		for(String m : argMax)
//		{
//			if(m.length() > maxL)
//			{
//				maxL = m.length();
//			}
//		}
//		List<String> result = new ArrayList<String>();
//		for(String m : argMax)
//		{
//			while(m.length() < maxL)
//			{
//				m = m + "-";
//			}
//			result.add(m);
//		}
//		return trim(result, maxL);
//	}
	
//	public static List<List<String>> alignCluster(List<List<String>> cluster)
//	{
//		int maxTotal = 0;
//		List<List<String>> argMax = new ArrayList<List<String>>();
//		for(int i=0; i<cluster.size(); i++)
//		{
//			String m1 = cluster.get(i).get(0);
////			String info = cluster.get(i).get(1);
//			int total = 0;
//			List<List<String>> l = new ArrayList<List<String>>();
//			l.add(cluster.get(i)); //adding m1 to l
//			for(int j=0; j<cluster.size(); j++)
//			{
//				String m2 = cluster.get(j).get(0);
//				if(i == j) continue;
//				String[] res = match(m1,m2);
//				total += Integer.parseInt(res[1]);
//				List<String> ll = new ArrayList<String>();
//				ll.add(res[0]); ll.add(cluster.get(j).get(1));
//				l.add(ll);
//			}
//			if(total > maxTotal || cluster.size() == 1)
//			{
//				maxTotal = total;
//				argMax = l;
//			}
//		}
//		//cosmetic phase - needed for presentation with WebLogo
//		int maxL = 0;
//		for(List<String> l : argMax)
//		{
//			String m = l.get(0);
//			if(m.length() > maxL)
//			{
//				maxL = m.length();
//			}
//		}
//		List<List<String>> result = new ArrayList<List<String>>();
//		for(List<String> l : argMax)
//		{
//			String m = l.get(0);
//			while(m.length() < maxL)
//			{
//				m = m + "-";
//			}
//			l.remove(0);
//			l.add(0, m); //will this work?
//			result.add(l);
//		}
//		return trim(result, maxL);
//	}
//	
//	public static List<List<String>> trim(List<List<String>> msa, int l)
//	{
//		List<List<String>> trimmed = new ArrayList<List<String>>();
//		double thresh = (double)msa.size()*0.8;
//		int i = 0;
//		while(i<l)
//		{
//			double cnt = 0;
//			for(List<String> li : msa)
//			{
//				String s = li.get(0);
//				if(s.charAt(i) == '-')
//				{
//					cnt++;
//				}
//			}
//			if(cnt >= thresh)
//			{
//				i++;
//			}
//			else break;
//		}
//		int j = l-1;
//		while(j<l && j>i)
//		{
//			double cnt = 0;
//			for(List<String> li : msa)
//			{
//				String s = li.get(0);
//				if(s.charAt(j) == '-')
//				{
//					cnt++;
//				}
//			}
//			if(cnt >= thresh)
//			{
//				j--;
//			}
//			else break;
//		}
//		for(List<String> li : msa)
//		{
//			String s = li.get(0);
//			List<String> nli = new ArrayList<String>();
//			nli.add(s.substring(i, j+1));
//			nli.add(li.get(1));
//			trimmed.add(nli);
//		}
//		
//		return trimmed;
//	}
	
//	public static List<String> trim(List<String> msa, int l)
//	{
//		List<String> trimmed = new ArrayList<String>();
//		double thresh = (double)msa.size()*0.8;
//		int i = 0;
//		while(i<l)
//		{
//			double cnt = 0;
//			for(String s : msa)
//			{
//				if(s.charAt(i) == '-')
//				{
//					cnt++;
//				}
//			}
//			if(cnt >= thresh)
//			{
//				i++;
//			}
//			else break;
//		}
//		int j = l-1;
//		while(j<l && j>i)
//		{
//			double cnt = 0;
//			for(String s : msa)
//			{
//				if(s.charAt(j) == '-')
//				{
//					cnt++;
//				}
//			}
//			if(cnt >= thresh)
//			{
//				j--;
//			}
//			else break;
//		}
//		for(String s : msa)
//		{
//			trimmed.add((s.substring(i, j+1)));
//		}
//		
//		return trimmed;
//	}
	
	public static List<List<String>> getOccurStringsLong(DSMotifSet set, List<String> sequences, int nstar) 
	{
		Map<Integer,Map<Character,SortedSet<Integer>>> map = new HashMap<Integer,Map<Character,SortedSet<Integer>>>();
		for(DSMotif m : set.getSet())
		{
			int l = m.getMotif1().length();
			TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ = m.getOcc();
			for(int i : occ.keySet())
			{
				if(i > nstar) continue;
				for(char strand : occ.get(i).keySet())
				{
					SortedSet<Integer> s = occ.get(i).get(strand);
					if(strand == '-')
					{
						s = translateOcc(s, sequences.get(i-1).length(), l);
					}
					for(int k : s)
					{
						Map<Character,SortedSet<Integer>> w;
						if(map.containsKey(i))
						{
							w = map.get(i);
						}
						else
						{
							w = new HashMap<Character,SortedSet<Integer>>();
						}
						SortedSet<Integer> x;
						if(w.containsKey(strand))
						{
							x = w.get(strand);
						}
						else
						{
							x = new TreeSet<Integer>();
						}
						for(int j=k; j<k+l; j++)
						{
							x.add(j);
						}
						w.put(strand, x);
						map.put(i, w);
					}
				}
			}
		}
		List<List<String>> mos = new ArrayList<List<String>>();
		for(int i : map.keySet())
		{
			for(char strand : map.get(i).keySet())
			{
				String seq = sequences.get(i-1);
				int len = seq.length();
				if(strand == '-')
				{
					seq = StringOperations.reverseComplement(seq);
				}
				int prev = -1;
				String subseq = "";
				int start = -1;
				for(int j : map.get(i).get(strand))
				{
					if(j == prev+1 || prev == -1)
					{
						subseq = subseq + seq.charAt(j);
					}
					else
					{
						String info;
						if(strand == '+')
						{
							info = new Integer(i-1).toString() + "\t" + strand + "\t" + new Integer(start).toString() + "\t" + new Integer(prev).toString();
						}
						else
						{
							int ind1 = len-1-prev;
							int ind2 = len-1-start;
							info = new Integer(i-1).toString() + "\t" + strand + "\t" + new Integer(ind1).toString() + "\t" + new Integer(ind2).toString();
						}
						List<String> list = new ArrayList<String>();
						list.add(subseq);
						list.add(info);
						mos.add(list);
						subseq = "" + seq.charAt(j);
					}
					if(j == map.get(i).get(strand).last())
					{
						String info;
						if(strand == '+')
						{
							info = new Integer(i-1).toString() + "\t" + strand + "\t" + new Integer(start).toString() + "\t" + new Integer(j).toString();
						}
						else
						{
							int ind1 = len-1-j;
							int ind2 = len-1-start;
							info = new Integer(i-1).toString() + "\t" + strand + "\t" + new Integer(ind1).toString() + "\t" + new Integer(ind2).toString();
						}
						List<String> list = new ArrayList<String>();
						list.add(subseq);
						list.add(info);
						mos.add(list);
					}
					if(prev == -1 || j > prev + 1) //this is the start of an occurrence
					{
						start = j;
					}
					prev = j;
				}
			}
		}
		return mos;		
	}
	
	public static List<List<String>> getOccurStringsShort(DSMotifSet set, List<String> sequences, int nstar) 
	{
		List<MotifsMapEntryLength> mels = new ArrayList<MotifsMapEntryLength>();
		for(DSMotif m : set.getSet())
		{
			MotifsMapEntryLength mel = new MotifsMapEntryLength(m);
			mels.add(mel);
		}
		Collections.sort(mels);
		Collections.reverse(mels);
		Map<Integer,Map<Character,SortedSet<Interval>>> map = new HashMap<Integer, Map<Character,SortedSet<Interval>>>();
		for(MotifsMapEntryLength m : mels)
		{
			int l = m.entry.getMotif().length();
			TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ = m.entry.getOcc();
			for(int i : occ.keySet())
			{
				if(i > nstar) continue;
				for(char strand : occ.get(i).keySet())
				{
					for(int k : occ.get(i).get(strand))
					{
						Map<Character,SortedSet<Interval>> w;
						if(map.containsKey(i))
						{
							w = map.get(i);
						}
						else
						{
							w = new HashMap<Character,SortedSet<Interval>>();
						}
						SortedSet<Interval> x;
						if(w.containsKey(strand))
						{
							x = w.get(strand);
						}
						else
						{
							x = new TreeSet<Interval>();
						}
						boolean flag = true;
						for(Interval in : x)
						{
							if(k >= in.start && k+l-1 <= in.end)
							{
								flag = false;
								break;
							}
							if(k+l-1 < in.start) break;
						}
						if(flag)
						{
							x.add(new Interval(k,k+l-1));
						}
						w.put(strand, x);
						map.put(i, w);
					}
				}
			}
		}
		List<List<String>> mos = new ArrayList<List<String>>();
		for(int i : map.keySet())
		{
			for(char strand : map.get(i).keySet())
			{
				String seq = sequences.get(i-1);
				for(Interval in : map.get(i).get(strand))
				{
					String subseq = seq.substring(in.start, in.end+1);
					if(strand == '-')
					{
						subseq = StringOperations.reverseComplement(subseq);
					}
					String info = new Integer(i-1).toString() + "\t" + strand + "\t" + new Integer(in.start).toString() + "\t" + new Integer(in.end).toString();
					List<String> list = new ArrayList<String>();
					list.add(subseq);
					list.add(info);
					mos.add(list);
				}
			}
		}
		return mos;		
	}
	
	public static SortedSet<Integer> translateOcc(SortedSet<Integer> s, int seqLen, int motifLen)
	{
		SortedSet<Integer> newSet = new TreeSet<Integer>();
		for(int i : s)
		{
			newSet.add(seqLen-motifLen-i);
		}
		return newSet;
	}
	
}
