/*
    File:
        CDSSubstitutionsCounter.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Counts the substitutions in the coding region.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.subst;

import algorithms.CodonComposition;
import algorithms.Path;
import algorithms.SiteComposition;
import bio.gene.dna.Codon;
import bio.gene.dna.ICodonTable;
import java.util.Locale;


public class CDSSubstitutionsCounter implements PluginMain.ISubstitutionsCounter
{
    private int[] nums_syn      = null;
    private int[] nums_nonsyn   = null;
    private float[] bc_syn      = null;
    private float[] bc_nonsyn   = null;
    private float[] sites       = null;
    private Locale locale       = null;
    private ICodonTable ct      = null;

    public CDSSubstitutionsCounter(Locale locale, ICodonTable ct)
    {
        nums_syn    = new int[6];
        nums_nonsyn = new int[6];
        bc_syn      = new float[4];
        bc_nonsyn   = new float[4];
        sites       = new float[2];
        this.ct     = ct;
        this.locale = locale;
    }

    public void addSubstitutions(CodonComposition ccp, CodonComposition cco, boolean bUseTerminal)
    {
        // Update bases count.
        float[] freqs = ccp.getBaseFrequencies(ct, bUseTerminal);
        for(int i=1;i<5;i++)
            bc_syn[i-1]+=freqs[i];
        for(int i=6;i<10;i++)
            bc_nonsyn[i-6]+=freqs[i];
        float[] nsites = ccp.getSitesCounts();
        sites[0]+=nsites[0];
        sites[1]+=nsites[1];
        if(cco==null)   // A single population, putative polymorphic site.
        {
            Path path = ccp.getEvolutionaryPath();
            if(path.getLength()==0)
                return;
            Codon[] codons = path.getPath();
            for(int i=1;i<codons.length;i++)
            {
                String prev = codons[i-1].getSequence();
                String cur  = codons[i].getSequence();
                if(codons[i]==codons[i-1])
                    continue;
                // Since the codons are different, there is an evolutionary path between them.
                Path p = Path.findBestPath(new Codon[]{codons[i], codons[i-1]}, ct, bUseTerminal);
                if(p.getLength()==0)
                    continue;
                // Find the difference between the codons.
                for(int k=0;k<3;k++)
                {
                    char c1 = prev.charAt(k);
                    char c2 = cur.charAt(k);
                    // Decide, whether the difference is syn. or nonsyn.
                    if(c1!=c2)
                    {
                        int type = p.getSubstitutionType(c2, k);
                        boolean bSyn = ( (type & Path.MASK_SYNONYMOUS)>0 );
                        addSubstitution(c1, c2, bSyn);
                    }
                }
            }
        }
        else    // Two populations, putative divergent site.
        {
            Path pp = ccp.getEvolutionaryPath();
            Path po = cco.getEvolutionaryPath();
            if(pp.getLength()==0 || po.getLength()==0)
                return;
            Codon[] cp = pp.getPath();
            Codon[] co = po.getPath();
            MAINLOOP:for(int i=0;i<3;i++)
            {
                SiteComposition scp = ccp.getSiteComposition(i);
                SiteComposition sco = cco.getSiteComposition(i);
                // Check, whether the site has divergent bases.
                int type = SiteComposition.getSiteType(scp, sco);
                if( (type & SiteComposition.ST_DIVERGENT)==0 )
                    continue;
                for(Codon c:cp)
                {
                    char b = c.getSequence().charAt(i);
                    if(scp.getBaseCount(b)>0 && sco.getBaseCount(b)==0)
                    {
                        Path p = Path.findBestPath(c, co, ct, bUseTerminal);
                        type = p.getSubstitutionType(b, i);
                        char otherbase = (char)(type>>24);
                        boolean bSyn = (type & Path.MASK_SYNONYMOUS) > 0;
                        addSubstitution(b, otherbase, bSyn);
                        continue MAINLOOP;
                    }
                }
            }
        }
    }

    public String getBasesCountString()
    {
        return String.format(locale, "%.2f\t%.2f\t%.2f\t%.2f\t\t%.2f\t%.2f\t%.2f\t%.2f",
                bc_syn[0], bc_syn[1], bc_syn[2], bc_syn[3],
                bc_nonsyn[0], bc_nonsyn[1], bc_nonsyn[2], bc_nonsyn[3]);
    }

    @Override
    public String toString()
    {
        int nTS_s = nums_syn[1]+nums_syn[4];
        int nTV_s = nums_syn[0]+nums_syn[2]+nums_syn[3]+nums_syn[5];
        int nTS_n = nums_nonsyn[1]+nums_nonsyn[4];
        int nTV_n = nums_nonsyn[0]+nums_nonsyn[2]+nums_nonsyn[3]+nums_nonsyn[5];
        return String.format(locale, "%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t" +
                                     "%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                             sites[0], nums_syn[0], nums_syn[1], nums_syn[2], nums_syn[3], nums_syn[4], nums_syn[5], nTS_s, nTV_s,
                             sites[1], nums_nonsyn[0], nums_nonsyn[1], nums_nonsyn[2], nums_nonsyn[3], nums_nonsyn[4], nums_nonsyn[5], nTS_n, nTV_n);
    }

    private void addSubstitution(char base1, char base2, boolean bSyn)
    {
        int[] arr = (bSyn) ? nums_syn : nums_nonsyn;
        String tmp = new String(new char[]{base1, base2});
        if(tmp.equalsIgnoreCase("AC") || tmp.equalsIgnoreCase("CA"))
            arr[0]++;
        else if(tmp.equalsIgnoreCase("AG") || tmp.equalsIgnoreCase("GA"))
            arr[1]++;
        else if(tmp.equalsIgnoreCase("AT") || tmp.equalsIgnoreCase("TA"))
            arr[2]++;
        else if(tmp.equalsIgnoreCase("CG") || tmp.equalsIgnoreCase("GC"))
            arr[3]++;
        else if(tmp.equalsIgnoreCase("CT") || tmp.equalsIgnoreCase("TC"))
            arr[4]++;
        else if(tmp.equalsIgnoreCase("GT") || tmp.equalsIgnoreCase("TG"))
            arr[5]++;
    }
}
