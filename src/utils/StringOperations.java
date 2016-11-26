package utils;

public class StringOperations {
	
	public static String reverseComplement(String s)
	{
		String t = "";
		for(int i=s.length()-1; i>=0; i--)
		{
			char c = s.charAt(i);
			if(c == 'A')
			{
				t = t + 'T';
			}
			else if(c == 'C')
			{
				t = t + 'G';
			}
			else if(c == 'G')
			{
				t = t + 'C';
			}
			else if(c == 'T')
			{
				t = t + 'A';
			}
			else
			{
				System.err.println(c + " is not in the DNA alphabet; string = " + s);
			}
		}
		return t;
	}
}
