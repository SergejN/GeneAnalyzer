/*
    File:
        CodonComposition.java
 *
    Revision:
        1.3.0.2
 *
    Description:
        Represents the composition of a codon in a multiple alignment.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

import bio.gene.dna.Codon;
import bio.gene.dna.ICodonTable;
import java.util.Vector;


public class CodonComposition
{
    private Vector<Codon> codons  = null;
    private ICodonTable ct        = null;

    // Number of gaps at the site.
    private int nGaps             = 0;
    // Number of codons containing a N.
    private int nN                = 0;
    // Number of codons containing a X.
    private int nX                = 0;
    // Total number of codons added.
    private int nTotal            = 0;

    private Path path             = null;

    // Site composition of each codon site.
    private SiteComposition[] comp= null;

    private boolean bUseTerminal  = false;


    /**
     *  Constructs new codon composition instance.
     *
     *  @param ct               codon table
     *  @param bUseTerminal     whether to use terminal codons when generating the path
     */
    public CodonComposition(ICodonTable ct, boolean bUseTerminal)
    {
        codons = new Vector<Codon>();
        this.ct = ct;
        this.bUseTerminal = bUseTerminal;
        comp = new SiteComposition[3];
        for(int i=0;i<3;i++)
            comp[i] = new SiteComposition();
    }

    /**
     *  Adds a single codon. If one of the following is true, the method returns
     *  false and the codon is not added:
     *      - strSequence is null
     *      - the length of strSequence is not 3
     *      - the sequence contains characters other than A,C,G,T,-,N,X
     *
     *  @param strSequence 
     */
    public boolean addCodon(String strSequence)
    {
        if( (strSequence==null) || (!strSequence.toUpperCase().matches("[ACGTNX-]{3}")) )
            return false;
        nTotal++;
        for(int i=0;i<3;i++)
            comp[i].addBase(strSequence.charAt(i));
        // If the codon is a valid codon, add it to the path.
        Codon c = Codon.getCodon(strSequence);
        if(c!=null)
        {
            path = null;
            codons.add(c);
        }
        else
        {
            // Otherwise update the number of gaps, N and X.
            if(strSequence.indexOf('-')>-1)
                nGaps++;
            else
            {
                if(strSequence.indexOf('N')>-1)
                    nN++;
                if(strSequence.indexOf('X')>-1)
                    nX++;
            }
        }
        return true;
    }

    /**
     *  Returns the composition of the specified site, or null, if the site
     *  index is invalid.
     *
     *  @param site
     *  @return
     */
    public SiteComposition getSiteComposition(int site)
    {
        if(site<0 || site>2)
            return null;
        return comp[site];
    }

    /**
     *  Returns the array of valid codons. Added codons which have a gap or N or X
     *  are not contained in this array.
     *
     *  @return
     */
    public Codon[] getValidCodons()
    {
        return (codons.size()>0) ? codons.toArray(new Codon[1]) : new Codon[0];
    }

    /**
     *  Returns the number of valid codons, i.e. the ones without X,-,N.
     *
     *  @return
     */
    public int getValidCodonsCount()
    {
        return codons.size();
    }

    /**
     *  Returns the total number of added codons.
     *
     *  @return
     */
    public int getTotalCodonsCount()
    {
        return nTotal;
    }

    /**
     *  Returns the number of codons which have a gap somewhere in the sequence.
     *
     *  @return
     */
    public int getGapsCount()
    {
        return nGaps;
    }

    /**
     *  Returns the number of codons with N in the sequence.
     *
     *  @return
     */
    public int getNCount()
    {
        return nN;
    }

    /**
     *  Returns the number of codons with X in the sequence.
     *
     *  @return
     */
    public int getXCount()
    {
        return nX;
    }

    /**
     *  Returns the number of synonymous and nonsynonymous sites.
     * 
     *  @return
     */
    public float[] getSitesCounts()
    {
        float[] res = {0.0f, 0.0f};
        for(Codon c:codons)
        {
            float[] tmp =  c.calculateNumberOfSites(ct, bUseTerminal);
            res[0] += tmp[0];
            res[1] += tmp[1];
        }
        res[0] /= codons.size();
        res[1] /= codons.size();
        return res;
    }

    /**
     *  Calculates the frequencies of the bases in the codon.
     *
     *  Algorithm:
     *      For each codon and site within the codon calculate the fraction F
     *      of the substitutions at this site, which lead to a synonymous codon.
     *      E.g. a codon CTT has one synonymous T at the last site, since each
     *      possible substitution CTx leads to a synonymous codon.
     *      Afterwards, add up the numbers for every codon and divide by the
     *      number of codons.
     *
     *  Returns:
     *      An array of 10 elements:
     *          index | value
     *              0 | number of syn. sites
     *              1 | number of syn. A's
     *              2 | number of syn. C's
     *              3 | number of syn. G's
     *              4 | number of syn. T's
     *              5 | number of nonsyn. sites
     *              6 | number of nonsyn. A's
     *              7 | number of nonsyn. C's
     *              8 | number of nonsyn. G's
     *              9 | number of nonsyn. T's
     *
     *  @return
     */
    public float[] getBaseFrequencies(ICodonTable ct, boolean bUseTerminal)
    {
        float[] res = new float[10];
        char[] bases = {'A', 'C', 'G', 'T'};
        for(Codon c:codons)
        {
            float[] sites_tmp = c.calculateNumberOfSites(ct, bUseTerminal);
            res[0] += sites_tmp[0];
            res[5] += sites_tmp[1];
            String strSeq = c.getSequence();
            for(int i=0;i<3;i++)
            {
                int ns = 0;
                int nt = 0;
                char ref = strSeq.charAt(i);
                for(char b:bases)
                {
                    if(b!=ref)
                    {
                        char[] tmp = strSeq.toCharArray();
                        tmp[i] = b;
                        String strCodon = new String(tmp);
                        if(!bUseTerminal && ct.isTerminal(strCodon))
                            continue;
                        nt++;
                        if(ct.areSynonymous(strSeq, strCodon))
                            ns++;
                    }
                }
                // Calculate the frequencies.
                float fs = (float)ns/(float)nt;
                float fn = 1.0f-fs;
                switch(ref)
                {
                    case 'A': res[1]+= fs; res[6]+=fn; break;
                    case 'C': res[2]+= fs; res[7]+=fn; break;
                    case 'G': res[3]+= fs; res[8]+=fn; break;
                    case 'T': res[4]+= fs; res[9]+=fn; break;
                }
            }
        }
        float size = codons.size();
        for(int i=0;i<10;i++)
            res[i] /= size;
        return res;
    }

    /**
     *  Returns the evolutionary path, explaining best the codon diversity at the site.
     *
     *  @return
     */
    public Path getEvolutionaryPath()
    {
        if(path==null)
            path = Path.findBestPath(codons.toArray(new Codon[codons.size()]), ct, bUseTerminal);
        return path;
    }

    /**
     *  Returns the number of polymorphisms. The codons with gaps, N's and X's are ignored.
     *  Array structure:
     *      index   | value
     *          0   | synonymous polymorphisms
     *          1   | nonsynonymous polymorphisms
     *
     *  @return
     */
    public int[] getNumberOfPolymorphisms()
    {
        if(path==null)
            path = Path.findBestPath(codons.toArray(new Codon[codons.size()]), ct, bUseTerminal);
        return path.getPolymorphismsCount();
    }

    /**
     *  Returns the number of singletons. The codons with gaps, N's and X's are ignored.
     *  Array structure:
     *      index   | value
     *          0   | synonymous singletons
     *          1   | nonsynonymous singletons
     *
     *  @return
     */
    public int[] getNumberOfSingletons()
    {
        return getNumberOfSingletons(1.0f);
    }

    /**
     *  Returns the number of singletons. The codons with gaps, N's and X's are ignored.
     *  Array structure:
     *      index   | value
     *          0   | synonymous singletons
     *          1   | nonsynonymous singletons
     *
     *  @param cof
     *  @return
     */
    public int[] getNumberOfSingletons(float cof)
    {
        Codon[] si = getSingletons(cof);
        if(si==null)
            return new int[]{0, 0};
        if(path==null)
            path = Path.findBestPath(codons.toArray(new Codon[codons.size()]), ct, bUseTerminal);
        int[] tmp = {0,0};
        for(Codon c:si)
        {
            if(path.isSubstitutionSynonymous(c))
                tmp[0]++;
            else
                tmp[1]++;
        }
        return tmp;
    }

     /**
     *  Returns an array of singletons. The singletons are searched for
     *  site-wise, not codon-wise, that means, for a codon sample like
     *      ATG
     *      ATC
     *      TTG
     *      TTG
     *  only ATC is a singleton, even though ATG appears only once, because A
     *  appears twice.
     *  cof specifies the cut-off frequency, i.e. singletons have frequency
     *  below it.
     *
     *  @param cof
     *  @param bConstSize 
     *  @return
     */
    public Codon[] getSingletons(float cof)
    {
        if(codons.size()<2)
            return null;
        Vector<Codon> tmp = new Vector<Codon>();
        // Iterate through the sites.
        char[] bases = {'A', 'C', 'G', 'T'};
        for(int i=0;i<3;i++)
        {
            SiteComposition sc = new SiteComposition();
            for(Codon c:codons)
                sc.addBase(c.getSequence().charAt(i));
            // Search for singletons.
            for(char b:bases)
            {
                // If the frequency of a base is below the threshold, add all
                // codons with that base at that site to the singletons list.
                // In the cut-off frequency is not specified, only the real
                // singletons, i.e. bases which appear only once at a certain
                // site, are considered.
                int n = sc.getBaseCount(b);
                if((n>0 && cof<0.5f && (float)n/(float)codons.size()<=cof) || (cof>=0.5f && n==1))
                {
                    for(Codon c:codons)
                    {
                        if(c.getSequence().charAt(i)==b)
                        {
                            if(tmp.indexOf(c)==-1)
                                tmp.add(c);
                        }
                    }
                }
            }
        }
        return (tmp.size()>0) ?  tmp.toArray(new Codon[1]) : null;
    }

    /**
     *  Calculates the number of syn. and nonsyn. transitions
     *  and transversions in the path.
     *  Array structure:
     *      index  | value
     *          0  | number of synonymous transitions
     *          1  | number of syn. transversions
     *          2  | number of nonsyn. transitions
     *          3  | number of nonsyn. transversions
     *
     *  @return
     */
    public int[] getSubstitutionsCount()
    {
        if(path==null)
            path = Path.findBestPath(codons.toArray(new Codon[codons.size()]), ct, bUseTerminal);
        return path.getSubstitutionsCount();
    }
    

    public static final int ST_INVALID              = 0;
    public static final int ST_MONOMORPHIC          = 0x00000001;
    public static final int ST_POLYMORPHIC_FIRST    = 0x00000010;
    public static final int ST_POLYMORPHIC_SECOND   = 0x00000100;
    public static final int ST_DIVERGENT            = 0x00001000;

    /**
     *  Returns the type of the alignment site. The result can be the combination
     *  of the following flags:
     *      ST_INVALID:             if one of the parameters is null. This flag is not
     *                              combined with any other.
     *      ST_MONOMORPHIC:         if the site is monomorphic in both populations.
     *                              This flag is not combined with others either.
     *      ST_POLYMORPHIC_FIRST:   if the site is polymorphic in the first population
     *      ST_POLYMORPHIC_SECOND:  if the site is polymorphic in the second population
     *      ST_DIVERGENT:           if the site is divergent between the populations
     *
     *  @param sc1
     *  @param sc2
     *  @return
     */
    public static int getSiteType(CodonComposition cc1, CodonComposition cc2)
    {
        if(cc1==null || cc1.getValidCodonsCount()==0 || cc2==null || cc2.getValidCodonsCount()==0)
            return ST_INVALID;
        int result = 0;
        // Polymorphisms.
        int[] tmp = cc1.getNumberOfPolymorphisms();
        if(tmp[0]>0 || tmp[1]>0)
            result = result | ST_POLYMORPHIC_FIRST;
        tmp = cc2.getNumberOfPolymorphisms();
        if(tmp[0]>0 || tmp[1]>0)
            result = result | ST_POLYMORPHIC_SECOND;
        // Divergence.
        Codon[] cp1 = cc1.getEvolutionaryPath().getPath();
        Codon[] cp2 = cc2.getEvolutionaryPath().getPath();
        boolean bDivergent = true;
        MAINLOOP: for(Codon c1:cp1)
        {
            for(Codon c2:cp2)
            {
                if(c1==c2)
                {
                    bDivergent = false;
                    break MAINLOOP;
                }
            }
        }
        if(bDivergent)
            result = result | ST_DIVERGENT;
        // If the result is still 0, then the site is monomorphic.
        if(result==0)
            result = ST_MONOMORPHIC;
        return result;
    }

    /**
     *  Merges the two codon composition instances. If the instances use different
     *  codon tables, codon networks or treat terminal codons differently, they
     *  cannot be merged, and the method returns null.
     *
     *  @param other
     *  @return
     */
    public static CodonComposition merge(CodonComposition cc1, CodonComposition cc2)
    {
        if(cc1.bUseTerminal!=cc2.bUseTerminal || cc1.ct!=cc2.ct)
            return null;
        CodonComposition cc = new CodonComposition(cc1.ct, cc1.bUseTerminal);
        cc.codons = new Vector<Codon>();
        cc.codons.addAll(cc1.codons);
        cc.codons.addAll(cc2.codons);
        cc.nGaps = cc1.nGaps+cc2.nGaps;
        cc.nN = cc1.nN+cc2.nN;
        cc.nX = cc1.nX+cc2.nX;
        cc.nTotal = cc1.nTotal+cc2.nTotal;
        cc.path = null;
        return cc;
    }
}
