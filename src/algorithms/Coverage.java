/*
    File:
        Coverage.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Represents the coverage of the gene entry, i.e. the number of sites
        which are valid in only 1 sequence, in 2 sequences etc.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

import bio.gene.Dataset;
import bio.gene.GeneEntry;


public class Coverage 
{
    // Numbers of sites with coverage 0-10%, 10-20% etc.
    private int[] counts = null;
    private int length   = 0;
    
    private Coverage()
    {
        counts = new int[10];
    }
    
    public int[] getCounts()
    {
        return counts;
    }
    
    public int getTotalLength()
    {
        return length;
    }
    

    public static Coverage calculateCoverage(Dataset ds)
    {
        if(ds==null)
            return null;
        int nGenes = ds.getGenesCount();
        Coverage cov = new Coverage();
        if(nGenes==0)
            return cov;
        for(int i=0;i<nGenes;i++)
        {
            Coverage tmp = calculateCoverage(ds.getGeneEntry(i));
            cov.length += tmp.length;
            for(int n=0;n<10;n++)
                cov.counts[n] += tmp.counts[n];
        }
        return cov;
    }
    
    public static Coverage calculateCoverage(GeneEntry ge)
    {
        if(ge==null)
            return null;
        Coverage cov = new Coverage();
        int nStrains = ge.getStrainsCount();
        if(nStrains==0)
            return cov;
        int l = ge.getStrainEntry(0).getCompleteSequence().length();
        cov.length = l;
        String[] seqs = new String[nStrains];
        for(int i=0;i<nStrains;i++)
            seqs[i] = ge.getStrainEntry(i).getCompleteSequence();
        for(int i=0;i<l;i++)
        {
            SiteComposition sc = calculateSiteComposition(seqs, i);
            int nCount = sc.getValidBasesCount()+sc.getGapsCount();
            float f = (float)nCount/(float)nStrains;
            if(f>=0.9f)
                cov.counts[9]++;
            else if(f>=0.8f)
                cov.counts[8]++;
            else if(f>=0.7f)
                cov.counts[7]++;
            else if(f>=0.6f)
                cov.counts[6]++;
            else if(f>=0.5f)
                cov.counts[5]++;
            else if(f>=0.4f)
                cov.counts[4]++;
            else if(f>=0.3f)
                cov.counts[3]++;
            else if(f>=0.2f)
                cov.counts[2]++;
            else if(f>=0.1f)
                cov.counts[1]++;
            else
                cov.counts[0]++;
        }
        return cov;
    }
    
    private static SiteComposition calculateSiteComposition(String[] seqs, int pos)
    {
        SiteComposition sc = new SiteComposition();
        for(int n=0;n<seqs.length;n++)
            sc.addBase(seqs[n].charAt(pos));
        return sc;
    }
}
