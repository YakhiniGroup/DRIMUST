package clustering;

public class Interval implements Comparable<Interval>{
	
	int start;
	int end;
	
	public Interval(int start, int end)
	{
		this.start = start;
		this.end = end;
	}

	public int compareTo(Interval obj) 
	{
		if(start < obj.start)
		{
			return -1;
		}
		else if(start == obj.start)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
}
