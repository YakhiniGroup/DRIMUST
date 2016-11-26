package suffix.trees;

public class EndPoint {
	
	private int e;
	
	public EndPoint(int e)
	{
		this.e = e; 
	}
	public int getEndPoint()
	{
		return e;
	}
	public void setEndPoint(int e)
	{
		this.e = e;
	}
	
	public void increment()
	{
		e++;
	}
}
