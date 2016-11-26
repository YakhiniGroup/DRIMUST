package utils;

import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;

import motif.search.MotifsMapEntry;

import statistics.HGScore;

public class DSMotif extends MotifsMapEntry{
	
	String m2;
	
	public DSMotif(String m1, String m2, HGScore score, TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ) 
	{
		super(m1, score, occ);
		this.m2 = m2;
	}

	public String getMotif1() 
	{
		return super.getMotif();
	}

	public String getMotif2() 
	{
		return m2;
	}

	public void swapMotifs()
	{
		String tmp = super.getMotif();
		super.setMotif(m2);
		m2 = tmp;
		for(int i : super.getOcc().keySet())
		{
			Map<Character,SortedSet<Integer>> m = super.getOcc().get(i);
			Map<Character,SortedSet<Integer>> m0 = new HashMap<Character, SortedSet<Integer>>();
			if(m.containsKey('+'))
			{
				SortedSet<Integer> s1 = m.get('+');
				m0.put('-', s1);
			}
			if(m.containsKey('-'))
			{
				SortedSet<Integer> s2 = m.get('-');
				m0.put('+', s2);
			}
			super.getOcc().put(i,m0);
		}
	}

	public HGScore getScore() 
	{
		return super.getmHG();
	}

	public void setScore(HGScore score) 
	{
		super.setmHG(score);
	}

	public TreeMap<Integer, Map<Character, SortedSet<Integer>>> getOcc() 
	{
		return super.getOcc();
	}

	public void setOcc(TreeMap<Integer, Map<Character, SortedSet<Integer>>> occ) 
	{
		super.setOcc(occ);
	}

}
