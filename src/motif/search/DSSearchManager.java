package motif.search;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import statistics.HGScore;
import statistics.mHG;
import suffix.trees.SuffixEdge;
import suffix.trees.SuffixNode;
import suffix.trees.SuffixTree;
import utils.DSMotif;
import utils.StringOperations;

public class DSSearchManager {
	
	List<Integer> lengths = new ArrayList<Integer>();
	SuffixTree tree;
	int k1;
	int k2;
	int N;
	float threshold = (float)0.0001;
	int tested = 0;
//	Map<String, TreeMap<Integer, SortedSet<Integer>>> motifOcc = new HashMap<String, TreeMap<Integer,SortedSet<Integer>>>();
	Set<String> seen = new HashSet<String>();
	List<DSMotif> motifs = new ArrayList<DSMotif>();
	
	public DSSearchManager(List<String> sequences, int k1, int k2)
	{
		this.tree = new SuffixTree(sequences);
		this.lengths = parseLengths(sequences);
		this.k1 = k1;
		this.k2 = k2;
		this.N = sequences.size();
	}
	
	private List<Integer> parseLengths(List<String> sequences)
	{
		List<Integer> lengths = new ArrayList<Integer>();
		for(String s : sequences)
		{
			lengths.add(s.length());
		}
		return lengths;
	}
	
	public void setThreshold(float threshold)
	{
		this.threshold = threshold;
	}
	
	public List<Integer> getSequencesLengths()
	{
		return lengths;
	}
	
	public List<DSMotif> doSearch()
	{
		tested = 0;
		searchKKmers(tree.getRoot(), 0, "");
//		List<DSMotif> motifs = testMotifs();
		System.out.println("tested " + tested + " motifs, out of which " + motifs.size() + " have mHG score better than the threshold");
		return motifs;
	}
	
	public TreeMap<Integer, SortedSet<Integer>> searchKKmers(SuffixNode node, int depth, String currStr)
	{
		TreeMap<Integer, SortedSet<Integer>> indices = new TreeMap<Integer, SortedSet<Integer>>();
		
		if(depth == k2 || (depth >= k1 && depth < k2 && node.isLeaf()))
		{
			indices = extractSequencesIndices(node, null);
			testMotif(currStr, indices);
			return indices;
		}
		if(false == node.isLeaf())
		{
			TreeMap<Integer, SortedSet<Integer>> childIndices;
			for(char c : node.getChildrenCharacters())
			{
				SuffixNode child = node.getChild(c);
				SuffixEdge edge = child.getInEdge();
				int diff = edge.getEndInd()-edge.getStartInd()+1;
				String addition;
				//we can stop in the middle of an edge
				if(depth+diff > k2)
				{
					addition = tree.getSubstring(edge.getStrInd(), edge.getStartInd(), edge.getStartInd()+k2-depth-1);
					childIndices = searchKKmers(child, k2, currStr+addition);
				}
				else
				{
					addition = tree.getSubstring(edge.getStrInd(), edge.getStartInd(), edge.getEndInd());
					childIndices = searchKKmers(child, depth+diff, currStr+addition);
				}
				if(depth >= k1 && depth < k2)
				{
					indices = mergeIndices(indices, childIndices);
				}
			}
			if(depth >= k1 && depth < k2)
			{
				testMotif(currStr, indices);
			}
		}
		return indices;
	}
	
	public static TreeMap<Integer, SortedSet<Integer>> mergeIndices(TreeMap<Integer, SortedSet<Integer>> indices1, TreeMap<Integer, SortedSet<Integer>> indices2)
	{
		TreeMap<Integer, SortedSet<Integer>> indices = new TreeMap<Integer, SortedSet<Integer>>();
		
		Integer[] keys1 = indices1.navigableKeySet().toArray(new Integer[0]);
		Integer[] keys2 = indices2.navigableKeySet().toArray(new Integer[0]);
		for(int i : keys1)
		{
			indices.put(i, new TreeSet<Integer>(indices1.get(i)));
		}
		for(int j : keys2)
		{
			if(indices.containsKey(j))
			{
//				SortedSet<Integer> s = new TreeSet<Integer>(indices2.get(j));
				indices.get(j).addAll(indices2.get(j));
			}
			else
			{
				indices.put(j, new TreeSet<Integer>(indices2.get(j)));
			}
		}
		return indices;
	}
	
//	public void updateMotifOccMap(TreeMap<Integer, SortedSet<Integer>> indices, String motif)
//	{
//		if(motif.charAt(motif.length()-1) == tree.getTerminationSymbol())
//		{
//			return;
//		}
//		TreeMap<Integer, SortedSet<Integer>> indicesCopy = new TreeMap<Integer, SortedSet<Integer>>(indices);
//		if(motifOcc.containsKey(motif))
//		{
//			System.err.println("motif is already present in motifOccMap");
//		}
//		motifOcc.put(motif, indicesCopy);
//	}

//	public List<DSMotif> testMotifs()
//	{
//		List<DSMotif> motifs = new ArrayList<DSMotif>();
//		Set<String> seen = new HashSet<String>();
//		for(String motif : motifOcc.keySet())
//		{
//			if(seen.contains(motif)) continue;
//			System.out.println(motif);
//			TreeMap<Integer, SortedSet<Integer>> occ1 = motifOcc.get(motif);
//			SortedSet<Integer> s = new TreeSet<Integer>(occ1.keySet());
//			String rcMotif = StringOperations.reverseComplement(motif);
//			TreeMap<Integer, SortedSet<Integer>> occ2 = getRCMotifOcc(rcMotif);
//			if(false == motif.equalsIgnoreCase(rcMotif))
//			{
//				s.addAll(occ2.keySet());
//			}
//			seen.add(motif);
//			seen.add(rcMotif);
//			tested++;
//			HGScore score = calculateMotifEnrichment(s);
//			if(score.score*score.B > threshold)
//			{
//				continue;
//			}
//			int n = score.n;
//			int representative = 1;
//			if(occ2.navigableKeySet().headSet(n+1).size() > occ1.navigableKeySet().headSet(n+1).size()) //m2 occurs more than m1 at the top n sequences
//			{
//				representative = 2;
//			}
//			DSMotif dsmotif;
//			if(representative == 1)
//			{
//				TreeMap<Integer,Map<Character,SortedSet<Integer>>> mergedOcc = mergeOcc(occ1, occ2);
//				dsmotif = new DSMotif(motif, rcMotif, score, mergedOcc);
//			}
//			else
//			{
//				TreeMap<Integer,Map<Character,SortedSet<Integer>>> mergedOcc = mergeOcc(occ2, occ1);
//				dsmotif = new DSMotif(rcMotif, motif, score, mergedOcc);
//			}
//			motifs.add(dsmotif);
//		}
//		return motifs;
//	}
	
	public void testMotif(String motif, TreeMap<Integer, SortedSet<Integer>> occ1)
	{
		if(motif.charAt(motif.length()-1) == tree.getTerminationSymbol())
		{
			return;
		}
		if(seen.contains(motif)) return;
		SortedSet<Integer> s = new TreeSet<Integer>(occ1.keySet());
		String rcMotif = StringOperations.reverseComplement(motif);
		TreeMap<Integer, SortedSet<Integer>> occ2 = getRCMotifOcc(rcMotif);
		if(false == motif.equalsIgnoreCase(rcMotif))
		{
			s.addAll(occ2.keySet());
		}
		seen.add(motif);
		seen.add(rcMotif);
		tested++;
		HGScore score = calculateMotifEnrichment(s);
		if(score.score*score.B > threshold)
		{
			return;
		}
		int n = score.n;
		int representative = 1;
		if(occ2.navigableKeySet().headSet(n+1).size() > occ1.navigableKeySet().headSet(n+1).size()) //m2 occurs more than m1 at the top n sequences
		{
			representative = 2;
		}
		DSMotif dsmotif;
		if(representative == 1)
		{
			TreeMap<Integer,Map<Character,SortedSet<Integer>>> mergedOcc = mergeOcc(occ1, occ2);
			dsmotif = new DSMotif(motif, rcMotif, score, mergedOcc);
		}
		else
		{
			TreeMap<Integer,Map<Character,SortedSet<Integer>>> mergedOcc = mergeOcc(occ2, occ1);
			dsmotif = new DSMotif(rcMotif, motif, score, mergedOcc);
		}
		motifs.add(dsmotif);
	}
	
	public TreeMap<Integer, SortedSet<Integer>> getRCMotifOcc(String rcMotif)
	{
//		if(motifOcc.containsKey(rcMotif))
//		{
//			return motifOcc.get(rcMotif);
//		}
		SuffixNode node = tree.walk(rcMotif);
		if(null == node)
		{
			return new TreeMap<Integer, SortedSet<Integer>>();
		}
		TreeMap<Integer, SortedSet<Integer>> occ = extractSequencesIndices(node, null);
		return occ;
	}
	
	public static TreeMap<Integer,Map<Character,SortedSet<Integer>>> mergeOcc(TreeMap<Integer, SortedSet<Integer>> occPlus, TreeMap<Integer, SortedSet<Integer>> occMinus)
	{
		TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ = new TreeMap<Integer, Map<Character,SortedSet<Integer>>>();
		for(int i : occPlus.keySet())
		{
			Map<Character,SortedSet<Integer>> m = new HashMap<Character, SortedSet<Integer>>();
			SortedSet<Integer> s = new TreeSet<Integer>(occPlus.get(i));
			m.put('+', s);
			occ.put(i,m);
		}
		for(int i : occMinus.keySet())
		{
			Map<Character,SortedSet<Integer>> m;
			if(occ.containsKey(i))
			{
				m = occ.get(i);
			}
			else
			{
				m = new HashMap<Character, SortedSet<Integer>>();
			}
			SortedSet<Integer> s = new TreeSet<Integer>(occMinus.get(i));
			m.put('-', s);
			occ.put(i,m);
		}
		return occ;
	}
	
//	public Map<String,HGScore> searchKmers(SuffixNode node, int depth, String currStr, HashMap<String,HGScore> motifsMap)
//	{
//		if(null == motifsMap)
//		{
//			motifsMap = new HashMap<String, HGScore>();
//		}
//		if(depth == k && currStr.charAt(k-1) != tree.getTerminationSymbol())
//		{
//			TreeMap<Integer, Integer> indices = extractSequencesIndices(node, null);
//			//do calculation of mHG
//			HGScore mHG = calculateMotifEnrichment(indices);
//			if(motifsMap.containsKey(currStr))
//			{
//				System.err.println("motif is already present");
//			}
//			if(mHG.score <= threshold)
//			{
//				motifsMap.put(currStr,mHG);
//			}
//			tested++;
//			return motifsMap;
//		}
//		if(false == node.isLeaf())
//		{
//			for(char c : node.getChildrenCharacters())
//			{
//				SuffixNode child = node.getChild(c);
//				SuffixEdge edge = child.getInEdge();
//				int diff = edge.getEndInd()-edge.getStartInd()+1;
//				String addition;
//				//we can stop in the middle of an edge
//				if(depth+diff > k)
//				{
//					addition = tree.getSubstring(edge.getStrInd(), edge.getStartInd(), edge.getStartInd()+k-depth-1);
//					searchKmers(child, k, currStr+addition, motifsMap);
//				}
//				else
//				{
//					addition = tree.getSubstring(edge.getStrInd(), edge.getStartInd(), edge.getEndInd());
//					searchKmers(child, depth+diff, currStr+addition, motifsMap);
//				}	
//			}
//		}
//		return motifsMap;
//	}
	
	public HGScore calculateMotifEnrichment(SortedSet<Integer> ns)
	{
		int B = ns.size();
		int n = 0;
		int b = 1;
		HGScore minHG = new HGScore(N,B,0,0,1);
		Iterator<Integer> iter = ns.iterator();
		int prevn = -1;
		while(iter.hasNext())
		{
			n = iter.next();
			if(n==prevn) System.err.println("double n here");
			HGScore s = mHG.calculate_HGT(N, n, B, b);
			if(s.score < minHG.score)
			{
				minHG = s;
			}
			b++;
		}
		return minHG;
	}
	
	public TreeMap<Integer, SortedSet<Integer>> extractSequencesIndices(SuffixNode node, TreeMap<Integer, SortedSet<Integer>> indices)
	{
		if(null == indices)
		{
			indices = new TreeMap<Integer, SortedSet<Integer>>();
		}
		if(node.isLeaf())
		{
			Set<Integer> containingSequencesInd = node.getAdditionalInfo().getPositions().keySet();
			for(int i : containingSequencesInd)
			{
				if(indices.containsKey(i+1))
				{
					indices.get(i+1).addAll(node.getAdditionalInfo().getPositions().get(i));
				}
				else
				{
					SortedSet<Integer> li = new TreeSet<Integer>();
					li.addAll(node.getAdditionalInfo().getPositions().get(i));
					indices.put(i+1, li);
				}
				//added 1 since the counting of indices starts at 0 and mHG starts from 1
			}
			return indices;
		}
		for(char c : node.getChildrenCharacters())
		{
			SuffixNode child = node.getChild(c);
			indices = extractSequencesIndices(child, indices);
		}
		return indices;
	}
	
}
