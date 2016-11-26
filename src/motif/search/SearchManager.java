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
import utils.StringOperations;

public class SearchManager {
	
	List<Integer> lengths = new ArrayList<Integer>();
	SuffixTree tree;
	int k1;
	int k2;
	int N;
	float threshold = (float)0.0001;
	int tested = 0;
	boolean dsMode = false;
	
	public SearchManager(List<String> sequences, int k1, int k2, boolean dsMode)
	{
		if(dsMode)
		{
			this.dsMode = true;
			this.tree = new SuffixTree(addReverseComplementSeqs(sequences));
		}
		else
		{
			this.tree = new SuffixTree(sequences);
		}
		this.lengths = parseLengths(sequences);
		this.k1 = k1;
		this.k2 = k2;
		this.N = sequences.size();
	}
	
	public SearchManager(List<String> sequences, int k, boolean dsMode)
	{
		if(dsMode)
		{
			this.dsMode = true;
			this.tree = new SuffixTree(addReverseComplementSeqs(sequences));
		}
		else
		{
			this.tree = new SuffixTree(sequences);
		}
		this.lengths = parseLengths(sequences);
		this.k1 = k;
		this.k2 = k;
		this.N = sequences.size();
	}
	
	private List<String> addReverseComplementSeqs(List<String> sequences)
	{
		List<String> list = new ArrayList<String>();
		for(String seq : sequences)
		{
			String s = seq + "#" + StringOperations.reverseComplement(seq);
			list.add(s);
		}
		return list;
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
	
	public Map<String,MotifsMapEntry> doSearch()
	{
		tested = 0;
		Map<String, MotifsMapEntry> motifsMap = new HashMap<String, MotifsMapEntry>();
		searchKKmers(tree.getRoot(), 0, "", motifsMap);
		System.out.println("tested " + tested + " motifs, out of which " + motifsMap.size() + " have mHG score better than the threshold");
//		motifsMap.put("size", new MotifsMapEntry(Integer.toString(tested), null, new TreeMap<Integer,Map<Character,SortedSet<Integer>>>()));
		return motifsMap;
	}
	
	public TreeMap<Integer, SortedSet<Integer>> searchKKmers(SuffixNode node, int depth, String currStr, Map<String, MotifsMapEntry> motifsMap)
	{
		TreeMap<Integer, SortedSet<Integer>> indices = new TreeMap<Integer, SortedSet<Integer>>();
		
		if(depth == k2 || (depth >= k1 && depth < k2 && node.isLeaf()))
		{
			indices = extractSequencesIndices(node, null, currStr);
			testMotif(indices, motifsMap, currStr);
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
					childIndices = searchKKmers(child, k2, currStr+addition, motifsMap);
				}
				else
				{
					addition = tree.getSubstring(edge.getStrInd(), edge.getStartInd(), edge.getEndInd());
					childIndices = searchKKmers(child, depth+diff, currStr+addition, motifsMap);
				}
				if(depth >= k1 && depth < k2)
				{
					indices = mergeIndices(indices, childIndices);
				}
			}
			if(depth >= k1 && depth < k2)
			{
				testMotif(indices, motifsMap, currStr);
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
			indices.put(i, indices1.get(i));
		}
		for(int j : keys2)
		{
			if(indices.containsKey(j))
			{
				indices.get(j).addAll(indices2.get(j));
			}
			else
			{
				indices.put(j, indices2.get(j));
			}
		}
		return indices;
	}

	public void testMotif(TreeMap<Integer, SortedSet<Integer>> indices, Map<String,MotifsMapEntry> motifsMap, String motif)
	{
		if(motif.charAt(motif.length()-1) == tree.getTerminationSymbol() || motif.contains("#"))
		{
			return;
		}
		HGScore mHG = calculateMotifEnrichment(indices);
		if(motifsMap.containsKey(motif))
		{
			System.err.println("motif is already present");
		}
		if(mHG.score*mHG.B <= threshold)
		{
			motifsMap.put(motif, new MotifsMapEntry(motif, mHG, processIndices(indices)));
		}
		tested++;
	}
	
	private TreeMap<Integer,Map<Character,SortedSet<Integer>>> processIndices(TreeMap<Integer, SortedSet<Integer>> indices)
	{
		TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ = new TreeMap<Integer,Map<Character,SortedSet<Integer>>>();
		for(int ind : indices.keySet())
		{
			int l = lengths.get(ind-1);
			SortedSet<Integer> set = indices.get(ind);
			for(int j : set)
			{
				int loc = j;
				char strand = '+';
				if(j >= l)
				{
					loc = j-l-1;
					strand = '-';
				}
				if(occ.containsKey(ind))
				{
					addToMap(occ.get(ind), strand, loc);
				}
				else
				{
					Map<Character,SortedSet<Integer>> m = new HashMap<Character,SortedSet<Integer>>();
					addToMap(m, strand, loc);
					occ.put(ind, m);
				}
				if(strand == '-' && dsMode != true) System.err.println("does not make sense, there must be an error...");
			}
		}
		return occ;
	}
	
	private void addToMap(Map<Character,SortedSet<Integer>> m, char c, int i)
	{
		if(m.containsKey(c))
		{
			m.get(c).add(i);
		}
		else
		{
			SortedSet<Integer> s = new TreeSet<Integer>();
			s.add(i);
			m.put(c, s);
		}
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
	
	public HGScore calculateMotifEnrichment(TreeMap<Integer,SortedSet<Integer>> indices)
	{
		NavigableSet<Integer> ns = indices.navigableKeySet();
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
	
	public TreeMap<Integer, SortedSet<Integer>> extractSequencesIndices(SuffixNode node, TreeMap<Integer, SortedSet<Integer>> indices, String currStr)
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
			indices = extractSequencesIndices(child, indices, currStr);
		}
		return indices;
	}
	
	public static void main(String[] args) throws IOException 
	{
//		List<String> sequences = new ArrayList<String>();
//		sequences.add("BOOKLET");
//		sequences.add("BOOJBARBOO");
//		SearchManager m = new SearchManager(sequences, 3, 4);
//		m.doSearch();
//		System.out.println(mHG.calculate_HGT(854046, 4312, 353510, 1607));
		FileInputStream fstream = new FileInputStream("in1.txt");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
        HGScore mins = new HGScore(0, 0, 0, 0, 1);
        String strLine;
        int B = 0;
        List<Integer> l = new ArrayList<Integer>();
        int i = 1;
		while((strLine = br.readLine()) != null)   
        {
			if(i > 160000) break;
        	if(strLine.equals("1"))
        	{
        		B++;
        		l.add(1);
        	}
        	else
        	{
        		l.add(0);
        	}
        	i++;
        }
//		Collections.reverse(l);
//		Collections.shuffle(l);
		List<Integer> l1 = new ArrayList<Integer>();
		for(int k=1; k<l.size(); k++)
		{
			if(l.get(k-1) == 1) l1.add(k);
		}
		l = l1;
		int b = 1;
		for(int j : l)
		{
			if(j > 50000) break;
			HGScore s = mHG.calculate_HGT(i-1, j, B, b++);
			if(s.score < mins.score)// && s.score >= 0 && s.score <= 1)
			{
				mins = s;
			}
		}
		mins.calcPvalue(10000);
		System.out.println(mins);
	}
}
