package clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import statistics.HGScore;
import statistics.mHG;

import motif.search.MotifsMapEntry;

public class MotifSet implements Comparable<MotifSet>{
	
	List<MotifsMapEntry> set;
	HGScore mHGScore;
	TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ;
	
	public MotifSet()
	{
		this.set = new ArrayList<MotifsMapEntry>();
		this.mHGScore = null;
		this.occ = new TreeMap<Integer, Map<Character,SortedSet<Integer>>>();
	}
	
	public MotifSet(List<MotifsMapEntry> set, HGScore mHGScore, TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ)
	{
		this.set = new ArrayList<MotifsMapEntry>(set);
		this.mHGScore = new HGScore(mHGScore.N, mHGScore.B, mHGScore.n, mHGScore.b, mHGScore.score);
		this.occ = new TreeMap<Integer, Map<Character,SortedSet<Integer>>>();
		for(int i : occ.navigableKeySet())
		{
			Map<Character,SortedSet<Integer>> map = new HashMap<Character, SortedSet<Integer>>();
			for(char strand : occ.get(i).keySet())
			{
				SortedSet<Integer> li = new TreeSet<Integer>(occ.get(i).get(strand));
				map.put(strand, li);
			}
			this.occ.put(i, map);
		}
	}
	
	public MotifSet(MotifsMapEntry m)
	{
		this.set = new ArrayList<MotifsMapEntry>();
		MotifsMapEntry m0 = new MotifsMapEntry(m.getMotif(), m.getmHG(), m.getOcc());
		set.add(m0);
		this.mHGScore = new HGScore(m.getmHG().N, m.getmHG().B, m.getmHG().n, m.getmHG().b, m.getmHG().score);
		this.occ = new TreeMap<Integer, Map<Character,SortedSet<Integer>>>();
		for(int i : m.getOcc().navigableKeySet())
		{
			Map<Character,SortedSet<Integer>> map = new HashMap<Character, SortedSet<Integer>>();;
			for(char strand : m.getOcc().get(i).keySet())
			{
				SortedSet<Integer> li = new TreeSet<Integer>(m.getOcc().get(i).get(strand));
				map.put(strand, li);
			}
			this.occ.put(i, map);
		}
	}

	public List<MotifsMapEntry> getSet() 
	{
		return set;
	}

	public void setSet(List<MotifsMapEntry> set) 
	{
		this.set = set;
	}

	public HGScore getmHGScore() 
	{
		return mHGScore;
	}

	public void setmHGScore(HGScore mHGScore) 
	{
		this.mHGScore = mHGScore;
	}
	
	public double getp() 
	{
		return mHGScore.score*mHGScore.B;
	}

	public TreeMap<Integer, Map<Character,SortedSet<Integer>>> getOcc() 
	{
		return occ;
	}

	public void setOcc(TreeMap<Integer, Map<Character,SortedSet<Integer>>> occ) 
	{
		this.occ = occ;
	}
	
	public int compareTo(MotifSet obj)
	{
		if(this.mHGScore.score*this.mHGScore.B < obj.mHGScore.score*obj.mHGScore.B)
		{
			return -1;
		}
		else if(this.mHGScore.score*this.mHGScore.B == obj.mHGScore.score*obj.mHGScore.B)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
	
	public static MotifSet addMotif(MotifsMapEntry m, MotifSet set, double threshold)
	{
		MotifSet newset = new MotifSet(set.getSet(), set.getmHGScore(), set.getOcc());
		TreeMap<Integer, Map<Character,SortedSet<Integer>>> occ = newset.getOcc();
		newset.getSet().add(m);
		for(int key : m.getOcc().keySet())
		{
			Map<Character,SortedSet<Integer>> map;
			if(occ.containsKey(key))
			{
				map = occ.get(key);
			}
			else
			{
				map = new HashMap<Character, SortedSet<Integer>>();
			}
			for(char strand : m.getOcc().get(key).keySet())
			{
				SortedSet<Integer> s;
				if(map.containsKey(strand))
				{
					s = map.get(strand);
				}
				else
				{
					s = new TreeSet<Integer>();
				}
				s.addAll(m.getOcc().get(key).get(strand));
				map.put(strand, s);
				occ.put(key, map);
			}
		}
		SortedSet<Integer> indices = convertIndices(occ.keySet());
		HGScore score = mHG.calculateMinimalHGIfGood(m.getmHG().N, indices, threshold);
		newset.setmHGScore(score);
		newset.setOcc(occ);
		
		return newset;
	}
	
	public static MotifSet removeMotif(MotifsMapEntry m, MotifSet set, double threshold)
	{
		MotifSet newset = new MotifSet();
		newset.setmHGScore(new HGScore(100, 1, 0, 0, 1));
		for(MotifsMapEntry m0 : set.getSet())
		{
			if(false == m.equals(m0))
			{
				newset = addMotif(m0, newset, threshold);
			}
		}
		
		return newset;
	}
	
	public static SortedSet<Integer> convertIndices(Set<Integer> occ)
	{
		SortedSet<Integer> newOcc = new TreeSet<Integer>(); 
		for(int i : occ)
		{
			newOcc.add(i+1);
		}
		return newOcc;
	}
}