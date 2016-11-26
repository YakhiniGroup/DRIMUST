package clustering;

import motif.search.MotifsMapEntry;

public class MotifsMapEntryLength implements Comparable<MotifsMapEntryLength>{
	
	MotifsMapEntry entry;

	public MotifsMapEntryLength(MotifsMapEntry entry) 
	{
		this.entry = entry;
	}

	public int compareTo(MotifsMapEntryLength arg0) 
	{
		if(entry.getMotif().length() < arg0.entry.getMotif().length())
		{
			return -1;
		}
		else if(entry.getMotif().length() == arg0.entry.getMotif().length())
		{
			return 0;
		}
		return 1;
	}
	
}
