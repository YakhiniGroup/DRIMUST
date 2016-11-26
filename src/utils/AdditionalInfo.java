package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AdditionalInfo {
	
	Map<Integer, List<Integer>> positions;
	
	public AdditionalInfo(int startPos, int stringIndex)
	{
		this.positions = new HashMap<Integer, List<Integer>>();
		List<Integer> l = new ArrayList<Integer>();
		l.add(startPos);
		this.positions.put(stringIndex, l);
	}
	
	public AdditionalInfo(Map<Integer, List<Integer>> positions)
	{
		this.positions = positions;
	}
	
	public Map<Integer, List<Integer>> getPositions()
	{
		return positions;
	}
	
	public void addPosition(int startPos, int stringIndex)
	{
		if(positions.containsKey(stringIndex) && false == positions.get(stringIndex).contains(startPos))
		{
			positions.get(stringIndex).add(startPos);
		}
		else if(false == positions.containsKey(stringIndex))
		{
			List<Integer> l = new ArrayList<Integer>();
			l.add(startPos);
			this.positions.put(stringIndex, l);
		}
	}
	
	public String toString()
	{
		String s = "";
		for(int strInd : positions.keySet())
		{
			for(int startPos : positions.get(strInd))
			{
				s = s + "(" + strInd + ":" + startPos + ")\t";
			}
		}
		return s;
	}
}
