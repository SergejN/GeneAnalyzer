/*
    File:
        SiteComposition.java
 *
    Revision:
        1.3.0.2
 *
    Description:
        Represents the base composition of a single site in a multiple alignment.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

/**
 *  CAUTION:
 *      DO NOT USE THIS CLASS TO ANALYZE AMBIGUOUS SITES!
 *      It only supports regular bases (A, C, G, T and -) and X and N.
 */
public class SiteComposition
{
    private int[] bc = null;    // Bases counts: A, C, G, T, -, N, X

    
    public SiteComposition()
    {
        bc = new int[7];
    }

    public SiteComposition(char[] bases)
    {
        bc = new int[7];
        for(char c:bases)
            addBase(c);
    }

    /**
     *  Adds a single base. Only supported bases can be added:
     *      - A, C, G, T, and gap ('-')
     *      - X for missing base
     *      - N for any base
     *  If any other base is passed, the method returns false.
     *
     *  @param base
     *  @return
     */
    public boolean addBase(char base)
    {
        switch(Character.toUpperCase(base))
        {
            case 'A': bc[0]++; return true;
            case 'C': bc[1]++; return true;
            case 'G': bc[2]++; return true;
            case 'T': bc[3]++; return true;
            case '-': bc[4]++; return true;
            case 'N': bc[5]++; return true;
            case 'X': bc[6]++; return true;
        }
        return false;
    }

    /**
     *  Returns the total number of VALID bases (A,C,G,T) at this site.
     *
     *  @return
     */
    public int getValidBasesCount()
    {
        return bc[0]+bc[1]+bc[2]+bc[3];
    }

    /**
     *  Returns the total bases count. The gaps are not counted.
     *
     *  @return
     */
    public int getTotalBasesCount()
    {
        return bc[0]+bc[1]+bc[2]+bc[3]+bc[5]+bc[6];
    }

    /**
     *  Returns the number of gaps at this site.
     *
     *  @return
     */
    public int getGapsCount()
    {
        return bc[4];
    }

    /**
     *  Returns the number of times the specified base appear at this site or -1
     *  if the base is not one of the following: A,C,G,T,-,N,X
     *
     *  @param base
     *  @return
     */
    public int getBaseCount(char base)
    {
        switch(Character.toUpperCase(base))
        {
            case 'A': return bc[0];
            case 'C': return bc[1];
            case 'G': return bc[2];
            case 'T': return bc[3];
            case '-': return bc[4];
            case 'N': return bc[5];
            case 'X': return bc[6];
        }
        return -1;
    }

    /**
     *  Returns the frequencies of the single bases at the site.
     *  If bUseAll is true, N and X are considered and their frequencies
     *  are also returned. Otherwise the frequencies are only calculated
     *  using the valid bases: A, C, G and T.
     *
     *  Remarks:
     *      The method returs an array, which has the following structure:
     *          index | value
     *              0 | A
     *              1 | C
     *              2 | G
     *              3 | T
     *       Additional values if bUseAll is true:
     *              4 | N
     *              5 | X
     *
     *  @return
     */
    public float[] getBaseFrequencies(boolean bUseAll)
    {
        if(bUseAll)
        {
            float nTotal = getTotalBasesCount();
            return new float[]{(float)bc[0]/nTotal,     // A
                               (float)bc[1]/nTotal,     // C
                               (float)bc[2]/nTotal,     // G
                               (float)bc[3]/nTotal,     // T
                               (float)bc[5]/nTotal,     // N
                               (float)bc[6]/nTotal};    // X
        }
        else
        {
            float nTotal = getValidBasesCount();
            return new float[]{(float)bc[0]/nTotal,     // A
                               (float)bc[1]/nTotal,     // C
                               (float)bc[2]/nTotal,     // G
                               (float)bc[3]/nTotal};    // T
        }
    }

    /**
     *  Returns the number of polymorphisms at this site.
     *  Gaps, N's and X's are ignored.
     *
     *  @return
     */
    public int getNumberOfPolymorphisms()
    {
        int nTotal = 0;
        for(int i=0;i<4;i++)
            if(bc[i]>0)
                nTotal++;
        return nTotal-1;
    }

    /**
     *  Returns the number of singletons at this site.
     *
     *  @return
     */
    public int getNumberOfSingletons()
    {
        if(getValidBasesCount()<2)
            return 0;
        int nSi = 0;
        for(int i=0;i<4;i++)
            if(bc[i]==1)
                nSi++;
        return nSi;
    }

    /**
     *  Returns the number of bases which have the frequency below the
     *  specified threshold value. The gaps are ignored, but X's and N's are
     *  not. I.e. if the alignment is as follows:
     *      AAA
     *      ATA
     *      A-G
     *      AXN
     *  The frequency of T is 1/3, and the frequency of G is 1/4.
     *
     *
     *  @param cof  cut-off threshold
     *  @return
     */
    public int getNumberOfSingletons(float cof)
    {
        if(cof>=0.5)
            return getNumberOfSingletons();
        int nTotal = getValidBasesCount();
        if(nTotal<2)
            return 0;
        int nSi = 0;
        for(int i=0;i<4;i++)
        {
            float n = (float)bc[i];
            if( (n>0) && (n/(float)nTotal<cof) )
                nSi++;
        }
        return nSi;
    }

    /**
     *  Returns the number of transitions (A<->T and C<->G) at this site.
     *
     *  @return
     */
    public float getNumberOfTransitions()
    {
        int nBases = 0;         // Number of different bases at the site.
        for(int i=0;i<4;i++)
            if(bc[i]>0)
                nBases++;
        switch(nBases)
        {
            case 0:
            case 1: return 0.0f;
            // If there are two different bases, check whether they belong to the same group or not.
            case 2: return ( (bc[0]>0 && bc[2]>0) || (bc[1]>0 && bc[3]>0) ) ? 1.0f : 0.0f;
            // For 3 different bases, there are 2 polymorphic sites. The number of
            // transitions is 1.
            case 3: return 1.0f;
            // For four different bases there are 1.5 transitions and 1.5 transversions.
            case 4: return 1.5f;
            default:
                return 0.0f;
        }
    }

    /**
     *  Returns the number of transversions (A/T<->C/G).
     *
     *  @return
     */
    public float getNumberOfTransversions()
    {
        return getNumberOfPolymorphisms()-getNumberOfTransversions();
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
    public static int getSiteType(SiteComposition sc1, SiteComposition sc2)
    {
        if(sc1==null || sc1.getValidBasesCount()==0 || sc2==null || sc2.getValidBasesCount()==0)
            return ST_INVALID;
        int result = 0;
        // Polymorphisms.
        if(sc1.getNumberOfPolymorphisms()>0)
            result = result | ST_POLYMORPHIC_FIRST;
        if(sc2.getNumberOfPolymorphisms()>0)
            result = result | ST_POLYMORPHIC_SECOND;
        // Divergence.
        if(sc1.getBaseCount('A')*sc2.getBaseCount('A')+
           sc1.getBaseCount('C')*sc2.getBaseCount('C')+
           sc1.getBaseCount('G')*sc2.getBaseCount('G')+
           sc1.getBaseCount('T')*sc2.getBaseCount('T')==0)
            result = result | ST_DIVERGENT;
        // If the result is still 0, then the site is monomorphic.
        if(result==0)
            result = ST_MONOMORPHIC;
        return result;
    }

    /**
     *  Merges the two site composition instances. The bases counts of the both objects are summed up.
     *
     *  @param other
     *  @return
     */
    public static SiteComposition merge(SiteComposition sc1, SiteComposition sc2)
    {
        SiteComposition sc = new SiteComposition();
        sc.bc = new int[sc1.bc.length];
        for(int i=0;i<sc1.bc.length;i++)
            sc.bc[i] = sc1.bc[i]+sc2.bc[i];
        return sc;
    }
}