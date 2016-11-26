package suffix.trees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.GroupLayout.SequentialGroup;

import utils.ExtensionObject;
import utils.LocationPointer;
import utils.StringOperations;


public class SuffixTree {
	
	private SuffixNode root;
	private char terminationSymbol = '$';
	private List<String> sequences;
	private EndPoint[] e;
	
	private static final int RULE_1 = 1;
	private static final int RULE_2_CREATE_LEAF = 2;
	private static final int RULE_2_SPLIT_EDGE = 4;
	private static final int RULE_3 = 3;
	
	public SuffixTree(List<String> sequences)
	{
		this.sequences = sequences;
		e = new EndPoint[sequences.size()];
		buildTree();
	}
	
	public char getTerminationSymbol()
	{
		return terminationSymbol;
	}
	
	public void buildTree()
	{
		root = new SuffixNode(null, -1, -1, -1, null, new HashMap<Character,SuffixNode>());
		
		for(int k=0; k<sequences.size(); k++)
		{
			String sequence = sequences.get(k) + terminationSymbol;
			int mk = sequence.length();
			e[k] = new EndPoint(-1);
			ExtensionObject extensionInfo = new ExtensionObject(root, 0, '-', -1, null);
			LocationPointer currentSearchLoc = new LocationPointer(root, 0, '-');
			LocationPointer currentHandleLoc = new LocationPointer(root, 0, '-');
			LocationPointer currentMatchLoc = new LocationPointer(root, 0, '-');
			int jL = 0;
			//handling all suffixes S[j..i] built from S[j..i-1]
			boolean toMatch = true;
			boolean matched = false;
			for(int i=0; i<mk; i++)
			{
				e[k].increment();	
				SuffixNode prevInternalNode = null;
				if(toMatch)
				{
//					System.out.println("string number " + k);
//					System.out.println(currentMatchLoc);
					currentMatchLoc = match(currentMatchLoc, k, i);
				}
				if(null != currentMatchLoc)
				{
					currentSearchLoc = currentMatchLoc;
					matched = true;
					if(i < mk-1) continue;
				}
				else
				{
					toMatch = false;
				}
				for(int j=jL; j<=i; j++)
				{
//					System.out.println("k = " + k + " i = " + i  + " j = "+ j + "#######");
	//				if(j == 2 && i == 2) System.out.println("!!current search loc " + currentSearchLoc);
					if(false == matched)
					{
						currentSearchLoc = find(k, j, i-1, currentHandleLoc, j>jL);		
					}
					matched = false;	
	//				if(j == 3 && i == 5) System.out.println("!!current search loc " + currentSearchLoc);
//					if(j == jL+1)
//					{
	//					System.out.println("doing the weird thing");
	//					extensionInfo.getClosestParent().setSuffixLink(currentSearchLoc.getClosestParent());
//					}
					extensionInfo = appendChar(currentSearchLoc, k, j, i);
//					System.out.println("extension info = " + extensionInfo);
//					if(extensionInfo.getRule() == RULE_3 && i == mk-1 && k > 0 && extensionInfo.getClosestParent().isLeaf())
					if(i == mk-1 && k > 0)
					{
						extensionInfo.getClosestParent().addAdditionalInfo(j, k);
					}
					if(extensionInfo.getRule() == RULE_2_SPLIT_EDGE)
					{
	//					System.out.println(" j = " + j + " i = " + i + " prevInternalNode = " + prevInternalNode + " newNode = " + extensionInfo.getNewNode());
						if(prevInternalNode != null)
						{
	//						System.out.println("setting suffix link!");
							prevInternalNode.setSuffixLink(extensionInfo.getNewNode());
						}
						prevInternalNode = extensionInfo.getNewNode();
					}
					else
					{
						prevInternalNode = null;
					}
					if(j == i)
					{
						extensionInfo.getClosestParent().setSuffixLink(root);
					}
					currentHandleLoc = new LocationPointer(extensionInfo.getClosestParent(), extensionInfo.getSteps(), extensionInfo.getChar());
					
	//				System.out.println("closest parent = " + extensionInfo.getClosestParent());
					if(i < mk-1 && (extensionInfo.getRule() == RULE_3 || j == i))
					{
						jL = j;
						break;
					}				
				}
	
//				System.out.println("i = " + i);
	//			postOrder("",root);
			}
//			System.out.println("k = " + k);
		}
	}
	
	public LocationPointer match(LocationPointer startLoc, int strInd, int i)
	{
//		System.out.println(startLoc + " strInd = " + strInd + " i = " + i);
		SuffixNode startNode = startLoc.getClosestParent(); 
		int steps = startLoc.getSteps();
		char ci = getCharAt(strInd,i);
		if(steps == 0)
		{
			if(startNode.hasChild(ci))
			{
				SuffixNode child = startNode.getChild(ci);
				if(child.getInEdge().getStartInd() < child.getInEdge().getEndInd())
				{
					return new LocationPointer(startNode, 1, ci);
				}
				else
				{
					return new LocationPointer(child, 0, '-');
				}
			}
			else
			{
				return null;
			}
		}
		else
		{
			char edgeChar = startLoc.getFirstChar();
			SuffixEdge inEdge = startNode.getChild(edgeChar).getInEdge();
			char c = getCharAt(inEdge.getStrInd(), inEdge.getStartInd() + steps);
			if(c == ci)
			{
				if(inEdge.getStartInd()+steps < inEdge.getEndInd())
				{
					return new LocationPointer(startNode, steps+1, edgeChar);
				}
				else
				{
					SuffixNode child = startNode.getChild(edgeChar);
					return new LocationPointer(child, 0, '-');
				}
			}
			else
			{
				return null;
			}
		}
	}
	
	public LocationPointer find(int strInd, int start,int end, LocationPointer prevHandleLoc, boolean flag)
	{
		//Last iteration we handled S[j..i-1] which is exactly what we need
		if(false == flag)
		{
			return prevHandleLoc;
		}
		if(start > end)
		{
			return new LocationPointer(root, 0, '-');
		}
//		System.out.println("inside ***find***");
		//Last iteration we handled S[j-1..i] but we need S[j-1..i-1]
		LocationPointer startLoc = stepup(prevHandleLoc);
//		System.out.println("start location: " + startLoc);
		SuffixNode startNode = startLoc.getClosestParent();
		int gamaStrInd = -1; //??
		int gamaStart = -1;
		int gamaEnd = -1;
		if(startLoc.getSteps() > 0)
		{
			char c = startLoc.getFirstChar();
			SuffixEdge edge = startNode.getChild(c).getInEdge();
			gamaStrInd = edge.getStrInd();
			gamaStart = edge.getStartInd();
			gamaEnd = gamaStart + startLoc.getSteps()-1;
		}
		else if(startNode.isLeaf())
		{
			gamaStrInd = startNode.getInEdge().getStrInd();
			gamaStart = startNode.getInEdge().getStartInd();
			gamaEnd = startNode.getInEdge().getEndInd();
			startNode = startNode.getParent();
		}
		SuffixNode sv = null;
		if(startNode.hasSuffixLink())
		{
			sv = startNode.getSuffixLink();
			if(false == sv.isRoot())
			{
//				System.out.println("using suffix links!");
			}
		}
		else if(false == startNode.isRoot())
		{
			startNode = root;
//			System.err.println("the suffix link invariant has been broken for " + startNode);
		}
		if(startNode.isRoot())
		{
			//here startNode = root
			sv = root;
			gamaStrInd = strInd;
			gamaStart = start;
			gamaEnd = end;
		}
//		System.out.println("gamaStart " + gamaStart + " gamaEnd " + gamaEnd);
		return fastscan(sv, gamaStrInd, gamaStart, gamaEnd);
	}
	
	public LocationPointer stepup(LocationPointer loc)
	{
		int steps = loc.getSteps();
		SuffixNode closestParent = loc.getClosestParent();
		char c = loc.getFirstChar();
		if(closestParent.isRoot())
		{
			System.err.println("Should not step up from the root");
		}
		if(steps > 0)
		{
			return new LocationPointer(closestParent, steps-1, c);
		}
		else
		{
			SuffixEdge inEdge = closestParent.getInEdge();
			c = getCharAt(inEdge.getStrInd(), inEdge.getStartInd());
			steps = inEdge.getEndInd() - inEdge.getStartInd();
			closestParent = closestParent.getParent();
			return new LocationPointer(closestParent, steps, c);
		}
	}
	
	public ExtensionObject appendChar(LocationPointer prevLocation, int strInd, int j, int i)
	{
		char toAdd = getCharAt(strInd, i); 
		SuffixNode closestParent = prevLocation.getClosestParent();
		int steps = prevLocation.getSteps();
		char c = prevLocation.getFirstChar();
		if(closestParent.isLeaf())
		{
			return new ExtensionObject(closestParent, 0, '-', RULE_1, null);
		}
		if(steps == 0)
		{
			if(closestParent.hasChild(toAdd))
			{
				SuffixNode child = closestParent.getChild(toAdd);
//				System.out.println("inside appendChar " + child);
				if(child.getInEdge().getEndInd() - child.getInEdge().getStartInd() == 0)
				{
					return new ExtensionObject(child, 0, '-', RULE_3, null);
				}
				else
				{
					return new ExtensionObject(closestParent, 1, toAdd, RULE_3, null);
				}
			}
			else
			{
				if(closestParent.isLeaf())
				{
					System.err.println(closestParent + " doesn't obey any rule");
				}
				SuffixNode leaf = new SuffixNode(closestParent, strInd, i, e[strInd], null, null);
				leaf.setAdditionalInfo(j, strInd);
				closestParent.setChild(toAdd, leaf);
				return new ExtensionObject(leaf, 0, '-', RULE_2_CREATE_LEAF, null);
			}
		}
		else
		{
			SuffixNode child = closestParent.getChild(c);
			int edgeStrInd = child.getInEdge().getStrInd();
			int edgeStartInd = child.getInEdge().getStartInd();
			int edgeEndInd = child.getInEdge().getEndInd();
			if(getCharAt(edgeStrInd, edgeStartInd+steps) == toAdd)
			{
				if(edgeStartInd+steps == edgeEndInd)
				{
					return new ExtensionObject(child, 0, '-', RULE_3, null);
				}
				else
				{
					return new ExtensionObject(closestParent, steps+1, c, RULE_3, null);
				}
			}
			else
			{
				HashMap<Character,SuffixNode> children = new HashMap<Character, SuffixNode>();
				SuffixNode newNode = new SuffixNode(closestParent, edgeStrInd, edgeStartInd, edgeStartInd+steps-1, null, children);
				child.getInEdge().setStartInd(edgeStartInd+steps);
				newNode.setChild(getCharAt(child.getInEdge().getStrInd(), child.getInEdge().getStartInd()), child);
				SuffixNode newLeaf = new SuffixNode(newNode, strInd, i, e[strInd], null, null);
				newLeaf.setAdditionalInfo(j, strInd);
				newNode.setChild(toAdd, newLeaf);
				child.setParent(newNode);
				closestParent.setChild(c, newNode);
				return new ExtensionObject(newLeaf, 0, '-', RULE_2_SPLIT_EDGE, newNode);
			}
		}
	}
	
	/**
	 * 
	 * @param sv
	 * @param gama
	 * @return the node in which the end of gama is located (at the subtree of s(v))
	 * and g as the place of the last character in the edge (counting from 1)
	 */
	public LocationPointer fastscan(SuffixNode sv, int strInd, int gamaStart, int gamaEnd)
	{
		int g = gamaEnd - gamaStart + 1;
		int h = 0;
		if(gamaStart < 0)
		{
			return new LocationPointer(sv, 0, '-');
		}
		while(true)
		{
//			System.out.println("starting node is " + sv + " gamaStart " + gamaStart + " gamaEnd " + gamaEnd);
			SuffixNode child = sv.getChild(getCharAt(strInd, gamaStart + h));
			SuffixEdge inEdge = child.getInEdge();
			int gtag = inEdge.getEndInd() - inEdge.getStartInd() + 1;
			if(gtag < g)
			{
				sv = child;
				g = g - gtag;
				h = h + gtag;
			}
			else if(gtag == g)
			{
				return new LocationPointer(child, 0, '-');
			}
			else
			{
				return new LocationPointer(sv, g, getCharAt(inEdge.getStrInd(), inEdge.getStartInd()));
			}
		}
	}
	
	public SuffixNode walk(String str)
	{
		return dowalk(str, root);
	}
	
	public SuffixNode dowalk(String str, SuffixNode node)
	{
		if(str.isEmpty())
		{
			return node;
		}
		node = node.getChild(str.charAt(0));
		if(null == node)
		{
			return null;
		}
		SuffixEdge edge = node.getInEdge();
		String addition = getSubstring(edge.getStrInd(), edge.getStartInd(), edge.getEndInd());
		if(addition.length() > str.length())
		{
			if(addition.startsWith(str))
			{
				return node;
			}
			return null;
		}
		else if(false == str.startsWith(addition))
		{
			return null;
		}
		int diff = edge.getEndInd()-edge.getStartInd()+1;
		return dowalk(str.substring(diff), node);
	}
	
	public void postOrder(String current, SuffixNode node)
	{
//		if(node.hasSuffixLink()) System.err.println("node " + node + " has suffix link to " + node.getSuffixLink());
		String str;
		if(node.isRoot())
		{
			str = "";
		}
		else
		{
			SuffixEdge inEdge = node.getInEdge();
			str = getSubstring(inEdge.getStrInd(), inEdge.getStartInd(), inEdge.getEndInd());
		}
		if(node.isLeaf())
		{
			System.out.println(current+str+node.getAdditionalInfo());
			return;
		}
		for(Character c : node.getChildrenCharacters())
		{
			postOrder(current+str, node.getChild(c));
		}
	}
	
	public void postOrder2(String current, SuffixNode node, List<String> suffixes)
	{
//		if(node.hasSuffixLink()) System.err.println("node " + node + " has suffix link to " + node.getSuffixLink());
		String str;
		if(node.isRoot())
		{
			str = "";
		}
		else
		{
			SuffixEdge inEdge = node.getInEdge();
			str = getSubstring(inEdge.getStrInd(), inEdge.getStartInd(), inEdge.getEndInd());
		}
		if(node.isLeaf())
		{
			suffixes.add(current+str+node.getAdditionalInfo());
			return;
		}
		for(Character c : node.getChildrenCharacters())
		{
			postOrder2(current+str, node.getChild(c), suffixes);
		}
	}
	
	/**
	 * 
	 * @param ind - the index of the sequence to be parsed
	 * @param start - the starting index of the substring (inclusive)
	 * @param end - the ending index of the substring (inclusive)
	 * @return the substring of the sequence located at 'ind', starting at 'start' and ending at 'end' 
	 */
	public String getSubstring(int strInd, int start, int end)
	{
		if(null != sequences)
		{
			String s = sequences.get(strInd) + terminationSymbol;
//			System.out.println(s + " start = " + start + " end = " + end);
			return s.substring(start, end+1);
		}
		else
		{
			System.err.println("Sequences collection was not initialized");
			return null;
		}
	}
	
	public String getSequence(int strInd)
	{
		if(null != sequences)
		{
			return sequences.get(strInd);
		}
		else
		{
			System.err.println("Sequences collection was not initialized");
			return null;
		}
	}
	
	private char getCharAt(int strInd, int start)
	{
		return getSubstring(strInd, start, start).charAt(0);
	}
	
	public SuffixNode getRoot()
	{
		return root;
	}
	
	public static void main(String[] argv)
	{
//		Map<Character,SuffixNode> children = new HashMap<Character, SuffixNode>();
//		SuffixNode root = new SuffixNode(null, 0, 0, null, children);
//		SuffixNode sv0 = new SuffixNode(root, 0, 1, null, children);
//		SuffixNode sv1 = new SuffixNode(sv0, 2, 3, null, children);
//		SuffixNode sv2 = new SuffixNode(sv1, 4, 6, null, children);
//		SuffixNode sv3 = new SuffixNode(sv2, 7, 10, null, children);
//		SuffixNode sv4 = new SuffixNode(sv3, 11, 11, null, null);
//		children.put('z', sv0);
//		children.put('b', sv1);
//		children.put('d', sv2);
//		children.put('g', sv3);
//		children.put('w', sv4);
//		
//		SuffixTree tree = new SuffixTree("zabcdefghyx");
		
//		ExtensionObject res = tree.walkDown(sv0,"bcdefghy");
//		System.out.println(res.getNode().inEdge.startInd);// + " " + res.getOutEdge() + " " + res.getG());
//		System.out.println(res.getG());
		
		List<String> sequences = new ArrayList<String>();
		sequences.add("BOOKLETBOO");
//		sequences.add("BOOKLET");
		sequences.add("BOOJBARBOO");
		SuffixTree t = new SuffixTree(sequences);
		System.out.println("The tree suffixes are:");
		List<String> suffixes = new ArrayList<String>();
//		t.postOrder2("",t.root, suffixes);
		t.postOrder("",t.root);
		System.out.println(t.walk("KLET"));
//		for(String suffix : suffixes) System.out.println(suffix);
//		for(Character c : t.root.children.keySet())
//		{
//			System.out.println(t.root.getChild('a'));
//		}
//		System.out.println(t.root.children.get('a'));
//		System.out.println(t.root.getChild('a'));
	}
}