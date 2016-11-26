package utils;

import statistics.HGScore;

public class LengthAndScore implements Comparable<LengthAndScore>{
	
	int length;
	HGScore score;
	
	public LengthAndScore (int length, HGScore score) 
	{
		this.length = length;
		this.score = score;
	}

	public int getLength() 
	{
		return length;
	}

	public void setLength(int length) 
	{
		this.length = length;
	}

	public HGScore getScore() 
	{
		return score;
	}

	public void setScore(HGScore score) 
	{
		this.score = score;
	}

	public int compareTo(LengthAndScore arg0) 
	{
		if(this.score.score < arg0.getScore().score)
		{
			return -1;
		}
		else if(this.score.score == arg0.getScore().score)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
}
