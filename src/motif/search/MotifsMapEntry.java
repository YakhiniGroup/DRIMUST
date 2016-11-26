package motif.search;

import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import statistics.HGScore;

public class MotifsMapEntry implements Comparable<MotifsMapEntry>{
	
	String motif;
	HGScore mHG;
	TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ;
	
	public MotifsMapEntry(String motif, HGScore mHG, TreeMap<Integer,Map<Character,SortedSet<Integer>>> occ)
	{
		this.motif = motif;
		this.mHG = mHG;
		this.occ = new TreeMap<Integer,Map<Character,SortedSet<Integer>>>();
		for(int i : occ.navigableKeySet())
		{
			Map<Character,SortedSet<Integer>> m = occ.get(i);
			for(char strand : m.keySet())
			{
				SortedSet<Integer> li = new TreeSet<Integer>(occ.get(i).get(strand));
				m.put(strand, li);
			}
			this.occ.put(i, m);
		}
	}

	public String getMotif() 
	{
		return motif;
	}

	public void setMotif(String motif) 
	{
		this.motif = motif;
	}

	public HGScore getmHG() 
	{
		return mHG;
	}

	public void setmHG(HGScore mHG) 
	{
		this.mHG = mHG;
	}

	public TreeMap<Integer, Map<Character,SortedSet<Integer>>> getOcc() 
	{
		return occ;
	}

	public void setOcc(TreeMap<Integer, Map<Character,SortedSet<Integer>>> occ) 
	{
		this.occ = occ;
	}
	
	public int compareTo(MotifsMapEntry obj)
	{
		if(this.mHG.score*this.mHG.B < obj.mHG.score*obj.mHG.B)
		{
			return -1;
		}
		else if(this.mHG.score*this.mHG.B == obj.mHG.score*obj.mHG.B)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
}
