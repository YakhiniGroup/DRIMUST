package utils;

import suffix.trees.SuffixNode;

public class ExtensionObject {
	
	SuffixNode closestParent;
	int steps;
	char c;
	int rule;
	SuffixNode newNode;
	
	public ExtensionObject(SuffixNode closestParent, int steps, char c, int rule, SuffixNode newNode) 
	{
		this.closestParent = closestParent;
		this.steps = steps;
		this.c = c;
		this.rule = rule;
		this.newNode = newNode;
	}

	public SuffixNode getClosestParent() 
	{
		return closestParent;
	}
	
	public int getSteps() 
	{
		return steps;
	}
	
	public char getChar()
	{
		return c;
	}

	public int getRule() 
	{
		return rule;
	}
	
	public SuffixNode getNewNode()
	{
		return newNode;
	}
	
	public String toString()
	{
		String s = "closest parent is (" + closestParent.getInEdge().getStrInd() + "," + closestParent.getInEdge().getStartInd() + "," + closestParent.getInEdge().getEndInd() + ") " +
				"num steps = " + steps + " rule = " + rule;
		return s;
	}
}
