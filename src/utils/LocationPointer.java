package utils;

import suffix.trees.SuffixNode;

public class LocationPointer {

	SuffixNode closestParent;
	int steps;
	char firstChar;
	
	public LocationPointer(SuffixNode closestParent, int steps, char firstChar) 
	{
		this.closestParent = closestParent;
		this.steps = steps;
		this.firstChar = firstChar;
	}

	public SuffixNode getClosestParent() 
	{
		return closestParent;
	}

	public int getSteps() 
	{
		return steps;
	}

	public char getFirstChar() 
	{
		return firstChar;
	}

	public void setClosestParent(SuffixNode closestParent) 
	{
		this.closestParent = closestParent;
	}

	public void setSteps(int steps) 
	{
		this.steps = steps;
	}

	public void setFirstChar(char firstChar) 
	{
		this.firstChar = firstChar;
	}
	
	public String toString()
	{
		String s = "closest parent: " + closestParent + " steps: " + steps + " first char : " + firstChar;
		return s;
	}
}
