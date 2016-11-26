package suffix.trees;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import utils.AdditionalInfo;


public class SuffixNode {
	
	private SuffixNode parent;
	private SuffixEdge inEdge;
	private SuffixNode suffixLink;
	private Map<Character,SuffixNode> children;
	private AdditionalInfo info;
	
	public SuffixNode(SuffixNode parent, int strInd, int startInd, EndPoint endInd, SuffixNode suffixLink, Map<Character,SuffixNode> children)
	{
		this.parent = parent;
		this.inEdge = new SuffixEdge(strInd, startInd, endInd);
		this.suffixLink = suffixLink;
		this.children = children;
		this.info = null;
	}
	
	public SuffixNode(SuffixNode parent, int strInd, int startInd, int endInd, SuffixNode suffixLink, Map<Character,SuffixNode> children)
	{
		this.parent = parent;
		this.inEdge = new SuffixEdge(strInd, startInd, endInd);
		this.suffixLink = suffixLink;
		this.children = children;
		this.info = null;
	}
	
	public AdditionalInfo getAdditionalInfo()
	{
		return info;
	}
	
	public void setAdditionalInfo(int startPos, int stringIndex)
	{
		this.info = new AdditionalInfo(startPos, stringIndex);
	}
	
	public void setAdditionalInfo(HashMap<Integer, List<Integer>> positions)
	{
		this.info = new AdditionalInfo(positions);
	}
	
	public void addAdditionalInfo(int startPos, int stringIndex)
	{
		if(null == this.info)
		{
			this.info = new AdditionalInfo(startPos, stringIndex);
		}
		else
		{
			this.info.addPosition(startPos, stringIndex);
		}
	}
	
	public boolean isLeaf()
	{
		return null == children && false == isRoot();  ///?????
	}
	
	public boolean isRoot()
	{
		return null == parent;
	}
	
	public boolean hasSuffixLink()
	{
		return null != suffixLink;
	}
 
	public boolean hasChild(Character c)
	{
		return (null != getChild(c));
	}
	
	public SuffixEdge getInEdge()
	{
		return inEdge;
	}
	
	public SuffixNode getChild(Character c)
	{
		return (null == children ? null : children.get(c));
	}
	
	public void setChild(char c, SuffixNode child)
	{
		if(null == children)
		{
			children = new HashMap<Character, SuffixNode>();
		}
		children.put(c, child);
	}
	
	public void setChildren(Map<Character, SuffixNode> children)
	{
		this.children = children;
	}
	
	public Collection<Character> getChildrenCharacters()
	{
		return children.keySet();
	}
	
	public SuffixNode getParent()
	{
		return parent;
	}
	
	public void setParent(SuffixNode parent)
	{
		this.parent = parent;
	}
	
	public SuffixNode getSuffixLink()
	{
		return suffixLink;
	}
	
	public void setSuffixLink(SuffixNode sv)
	{
		this.suffixLink = sv;
	}
	
	public void setInEdge(SuffixEdge inEdge)
	{
		this.inEdge = inEdge;
	}
	
	public String toString()
	{
		return inEdge.getStrInd() + " : " + inEdge.getStartInd() + ", " + inEdge.getEndInd();
	}
}
