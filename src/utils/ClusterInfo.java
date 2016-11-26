package utils;

import java.util.List;
import java.util.Set;

import motif.search.MotifsMapEntry;

public class ClusterInfo implements Comparable<ClusterInfo>{
	
	List<String> cluster;
	double bestEnrichment;
	List<MotifsMapEntry> motifs;
	
	public ClusterInfo(List<String> cluster, double bestEnrichment, List<MotifsMapEntry> motifs) 
	{
		this.cluster = cluster;
		this.bestEnrichment = bestEnrichment;
		this.motifs = motifs;
	}

	public List<String> getCluster() 
	{
		return cluster;
	}

	public void setCluster(List<String> cluster) 
	{
		this.cluster = cluster;
	}

	public double getBestEnrichment() 
	{
		return bestEnrichment;
	}

	public void setBestEnrichment(double bestEnrichment) 
	{
		this.bestEnrichment = bestEnrichment;
	}
	
	public List<MotifsMapEntry> getMotifs() 
	{
		return motifs;
	}

	public void setMotifs(List<MotifsMapEntry> motifs) 
	{
		this.motifs = motifs;
	}

	public int compareTo(ClusterInfo arg0) 
	{
		if(this.bestEnrichment < arg0.getBestEnrichment())
		{
			return -1;
		}
		else if(this.bestEnrichment == arg0.getBestEnrichment())
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
}
