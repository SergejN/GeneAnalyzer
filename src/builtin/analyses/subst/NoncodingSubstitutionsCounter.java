/*
    File:
        NoncodingSubstitutionsCounter.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Counts the substitutions in the non-coding regions and at four-fold
        degenerate sites.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.subst;

import algorithms.SiteComposition;
import java.util.Locale;


public class NoncodingSubstitutionsCounter implements PluginMain.ISubstitutionsCounter
{
    private int nSites = 0;
    private int[] nums = null;          // Number of substitutions.
    private float[] bc = null;          // Base counts: A, C, G, T

    private Locale locale = null;

    public NoncodingSubstitutionsCounter(Locale locale)
    {
        nums = new int[7];
        bc   = new float[4];
        this.locale = locale;
    }

    /**
     *  Adds the substitution. If sco is null, then the site is assumed to be
     *  polymorphic. If sco is not null, the site is assumed to be divergent.
     *
     *  @param scp
     *  @param sco
     */
    public void addSubstitution(SiteComposition scp, SiteComposition sco)
    {
        // If the number of substitutions greater than 1, then the substitution type
        // cannot be figured out.
        int np = scp.getNumberOfPolymorphisms();
        if(np>1)
        {
            nums[6]++;
            return;
        }
        char[] bases = {'A', 'C', 'G', 'T'};
        char base1 = 0;
        char base2 = 0;
        if(sco==null)   // A single population, putative polymorphic site.
        {
            // Increment the number of sites analyzed.
            nSites++;
            // Update the base counts.
            float[] freqs = scp.getBaseFrequencies(false);
            for(int i=0;i<freqs.length;i++)
                bc[i]+=freqs[i];
            // Find the polymorphic bases.
            for(char b:bases)
            {
                if(scp.getBaseCount(b)>0)
                {
                    if(base1==0)
                        base1 = b;
                    base2 = b;
                }
            }
        }
        else    // Two populations, putative divergent site.
        {
            if(np!=0 || sco.getNumberOfPolymorphisms()!=0)
            {
                nums[6]++;
                return;
            }
            // Increment the number of sites analyzed.
            nSites++;
            // Update the base counts.
            float[] freqs = scp.getBaseFrequencies(false);
            for(int i=0;i<freqs.length;i++)
                bc[i]+=freqs[i];
            // Find the divergent bases.
            for(char b:bases)
            {
                if(scp.getBaseCount(b)>0)
                    base1 = b;
                if(sco.getBaseCount(b)>0)
                    base2 = b;
            }
            // If the site is monomorphic, simply return.
            if(base1==base2)
                return;
        }
        // At this point the both bases MUST be different.
        String tmp = new String(new char[]{base1, base2});
        if(tmp.equalsIgnoreCase("AC") || tmp.equalsIgnoreCase("CA"))
            nums[0]++;
        else if(tmp.equalsIgnoreCase("AG") || tmp.equalsIgnoreCase("GA"))
            nums[1]++;
        else if(tmp.equalsIgnoreCase("AT") || tmp.equalsIgnoreCase("TA"))
            nums[2]++;
        else if(tmp.equalsIgnoreCase("CG") || tmp.equalsIgnoreCase("GC"))
            nums[3]++;
        else if(tmp.equalsIgnoreCase("CT") || tmp.equalsIgnoreCase("TC"))
            nums[4]++;
        else if(tmp.equalsIgnoreCase("GT") || tmp.equalsIgnoreCase("TG"))
            nums[5]++;
    }

    public String getBasesCountString()
    {
        return String.format(locale, "%.2f\t%.2f\t%.2f\t%.2f", bc[0], bc[1], bc[2], bc[3]);
    }

    @Override
    public String toString()
    {
        int nTS = nums[1]+nums[4];
        int nTV = nums[0]+nums[2]+nums[3]+nums[5];
        return String.format(locale, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                nSites, nums[0], nums[1], nums[2], nums[3], nums[4], nums[5], nTS, nTV, nums[6]);
    }
}