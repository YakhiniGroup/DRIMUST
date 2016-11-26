package suffix.trees;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SuffixTreesTester {
	
	public static Map<String,Map<Integer,Integer>> produceSuffixes(List<String> input)
	{
		Map<String,Map<Integer,Integer>> suffixes = new HashMap<String,Map<Integer,Integer>>();
		for(int k=0; k<input.size(); k++)
		{
			String s = input.get(k) + '$';
			for(int i=0; i<s.length(); i++)
			{
				String suffix = s.substring(i);
				if(suffixes.containsKey(suffix))
				{
					suffixes.get(suffix).put(k,i);
				}
				else
				{
					Map<Integer,Integer> m = new HashMap<Integer, Integer>();
					m.put(k, i);
					suffixes.put(suffix, m);
				}
			}
		}
		return suffixes;
	}
	
	public static Map<String,List<String>> getSubstrings(List<String> suffixes, int length)
	{
		Map<String,List<String>> substrings = new HashMap<String, List<String>>();
		for(String s : suffixes)
		{
			String substring = (s.length() <= length) ? s : s.substring(0, length);
			String occ = s.substring(s.indexOf('('));
			String[] occs = occ.split("\\(");
			for(String o : occs)
			{
				if(o.isEmpty()) continue;
				o.replace("(", ""); o.replace(")", ""); o.trim();
				if(substrings.containsKey(substring))
				{
					if(false == substrings.get(substring).contains(o))
					{
						substrings.get(substring).add(o);
					}
				}
				else
				{
					List<String> l = new ArrayList<String>();
					l.add(o);
					substrings.put(substring, l);
				}
			}
		}
		return substrings;
	}
	
	public static void main(String[] args) throws IOException
	{			
		FileInputStream fstream = new FileInputStream("C:\\WORK\\thesis\\motif search\\scripts\\p53_nar.txt");
//		FileInputStream fstream = new FileInputStream("C:\\DRIM\\final\\miR-410.fa");
		DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        String strLine;
        int N = 0;
        List<String> sequences = new ArrayList<String>();
        while((strLine = br.readLine()) != null)   
        {
        	if(strLine.isEmpty() || strLine.startsWith(">"))
        	{
        		continue;
        	}
        	N++;
        	if(N<2000) 
        	{
        		sequences.add(strLine);
        	}
//        	else 
//        	{
//        		break;
//        	}
        }
//        sequences.add("BANANA");
//		sequences.add("ANONA");
//		sequences.add("BEST");
//		sequences.add("HANA");
//        System.out.println(sequences);	
		//produce suffixes naively
//		Map<String,Map<Integer,Integer>> suffixes1 = produceSuffixes(sequences);
//			System.out.println(suffixes1);
//		produce suffixes using suffix trees
		SuffixTree t = new SuffixTree(sequences);
		List<String> suffixes1 = new ArrayList<String>();
		t.postOrder2("",t.getRoot(), suffixes1);
		Map<String, List<String>> m1 = getSubstrings(suffixes1, 14);
//		SuffixTreeTrunc t2 = new SuffixTreeTrunc(sequences, 20);
		SuffixTreeTruncated t2 = new SuffixTreeTruncated(sequences, 20);
		List<String> suffixes2 = new ArrayList<String>();
		t2.postOrder2("",t2.getRoot(), suffixes2);
		Map<String, List<String>> m2 = getSubstrings(suffixes2, 14);
		int cnt1 = 0;
		for(String sbstrg1 : m1.keySet())
		{
			if(m2.containsKey(sbstrg1)) cnt1++;
			else continue;
			for(String occ : m1.get(sbstrg1))
			{
				if(false == m2.get(sbstrg1).contains(occ))
				{
					System.err.println(sbstrg1 + " " + occ + " is not in second tree");
				}
			}
		}
		int cnt2 = 0;
		for(String sbstrg2 : m2.keySet())
		{
			if(m1.containsKey(sbstrg2)) cnt2++;
			else continue;
			for(String occ : m2.get(sbstrg2))
			{
				if(false == m1.get(sbstrg2).contains(occ))
				{
					System.err.println(sbstrg2 + " " + occ + " is not in first tree");
				}
			}
		}
		System.out.println(cnt1 + " " + cnt2 + " " + m1.size() + " " + m2.size());
		System.out.println("success? " + ((cnt1 == cnt2)&&(cnt1 == m1.size())&&(cnt1==m2.size())));
//			System.out.println(suffixes2);
//		int cnt = 0;
//		for(String suffix1 : suffixes1.keySet())
//		{
//			Map<Integer,Integer> positions1 = suffixes1.get(suffix1);
////			System.out.println(suffix1 + " positions:");
////			for(Integer p : positions1.keySet()) 
////			{
////				System.out.println(p + ":" + positions1.get(p));
////			}
//			for(String suffix2 : suffixes2)
//			{
//				String suff = suffix2.substring(0,suffix2.indexOf('('));
//				String rest = suffix2.substring(suffix2.indexOf('('));
//				if(suff.equals(suffix1))
//				{
//					int n = 0;
//					String[] positions2 = rest.split("\t");
//					for(String pos : positions2)
//					{
//						int strInd = Integer.parseInt(pos.substring(pos.indexOf('(')+1,pos.indexOf(':')));
//						int startInd = Integer.parseInt(pos.substring(pos.indexOf(':')+1,pos.indexOf(')')));
////						System.err.println(strInd + " " + startInd);
//						if(positions1.containsKey(strInd) && positions1.get(strInd).intValue() == startInd)
//						{
//							n++;
//						}
//					}
//					if(n == positions1.keySet().size())
//					{
//						cnt++;
//					}
//					else
//					{
//						System.err.println(suffix2);
//					}
//					break;
//				}
//			}
//		}
//		System.out.println("success? " + (cnt == suffixes1.keySet().size()));
    }
}
