package suffix.trees;

public class SuffixEdge {

	private int strInd;
	private int startInd;
	private EndPoint endInd;
	
	public SuffixEdge(int strInd, int startInd, int endInd)
	{
		this.strInd = strInd;
		this.startInd = startInd;
		this.endInd = new EndPoint(endInd);
	}
	
	public SuffixEdge(int strInd, int startInd, EndPoint endInd)
	{
		this.strInd = strInd;
		this.startInd = startInd;
		this.endInd = endInd;
	}
	
	public int getStrInd()
	{
		return strInd;
	}
	
	public void setStrInd(int strInd)
	{
		this.strInd = strInd;
	}
	
	public int getStartInd()
	{
		return startInd;
	}
	
	public void setStartInd(int startInd)
	{
		this.startInd = startInd;
	}
	
	public int getEndInd()
	{
		return endInd.getEndPoint();
	}
	
	public void setEndInd(int endInd)
	{
		this.endInd.setEndPoint(endInd);
	}
}
