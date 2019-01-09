/*
    File:
        SitesBlock.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Represents the block of several sites of the alignment. The number of strains
        if the same at all sites of the block.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;


public class SitesBlock
{
    private int nStrains = 0;
    private float pi     = 0.0f;
    private float sites  = 0.0f;
    private int nPoly    = 0;
    private int nSi      = 0;
    private float nTS    = 0;

    private boolean bJC_pi = false;
    private boolean bJC_t  = false;
    private float cof      = 1.0f;
    

    public SitesBlock(int nStrains, boolean bJC_pi, boolean bJC_theta, float cof)
    {
        this.nStrains = nStrains;
        this.bJC_pi = bJC_pi;
        this.bJC_t = bJC_theta;
        this.cof = cof;
    }
    
    /**
     *  Adds the site to the block. If the site has different number of strains
     *  than defined for this block, the method does nothing and returns false.
     * 
     *  @param sc
     *  @param cof  cutoff frequency for counting the singletons
     *  @return
     */
    public boolean addSite(SiteComposition sc)
    {
        if(sc.getValidBasesCount()!=nStrains)
            return false;
        pi += BasicStatistics.calculatePi(sc);
        nPoly += sc.getNumberOfPolymorphisms();
        nSi += sc.getNumberOfSingletons(cof);
        nTS += sc.getNumberOfTransitions();
        sites++;
        return true;
    }

    /**
     *  Returns the number of strains in this block.
     * 
     *  @return
     */
    public int getStrainsCount()
    {
        return nStrains;
    }

    /**
     *  Returns pi of the block.
     *
     *  @return
     */
    public float getPi()
    {
        float _pi = (sites>0) ? pi/sites : 0.0f;
        return (bJC_pi) ? BasicStatistics.correctJC(_pi) : _pi;
    }

    /**
     *  Returns theta of the block.
     *
     *  @return
     */
    public float getTheta()
    {
        float t = BasicStatistics.calculateTheta(nPoly, sites, nStrains);
        return (bJC_t) ? BasicStatistics.correctJC(t) : t;
    }

    /**
     *  Returns the number of sites in the block.
     *
     *  @return
     */
    public float getSitesCount()
    {
        return sites;
    }

    /**
     *  Returns the number of polymorphic sites in the block.
     *
     *  @return
     */
    public int getPolymorphismsCount()
    {
        return nPoly;
    }

    /**
     *  Returns the number of singletons in the block.
     *
     *  @return
     */
    public int getSingletonsCount()
    {
        return nSi;
    }

    /**
     *  Returns the number of transitions in the block. The number of transversions
     *  can be calculated by subtracting this number from the total number of polymorphisms.
     * 
     *  @return
     */
    public float getTransitionsCount()
    {
        return nTS;
    }

    /**
     *  Returns the overall Pi value of all blocks.
     *
     *  @param blocks
     *  @param bJC
     *  @return
     */
    public static float getPi(SitesBlock[] blocks, boolean bJC)
    {
        float pi = 0.0f;
        float n = 0.0f;
        for(SitesBlock sb:blocks)
        {
            pi += sb.pi;
            n += sb.sites;
        }
        return (bJC) ? BasicStatistics.correctJC(pi/n) : pi/n;
    }

    /**
     *  Returns the overall Theta value of all blocks.
     * 
     *  @param blocks
     *  @param bJC 
     *  @return
     */
    public static float getTheta(SitesBlock[] blocks, boolean bJC)
    {
        float t = 0.0f;
        float sites = 0.0f;
        for(SitesBlock sb:blocks)
        {
            float f = BasicStatistics.calculateTheta(sb.nPoly, sb.sites, sb.nStrains);
            if(bJC)
                f = BasicStatistics.correctJC(f);
            t += f*sb.sites;
            sites += sb.sites;
        }
        return t/sites;
    }

    /**
     *  Returns the total number of sites in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static float getSitesCount(SitesBlock[] blocks)
    {
        float n = 0.0f;
        for(SitesBlock sb:blocks)
            n += sb.sites;
        return n;
    }

    /**
     *  Returns the total number of polymorphisms in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static int getPolymorphismsCount(SitesBlock[] blocks)
    {
        int n = 0;
        for(SitesBlock sb:blocks)
            n += sb.nPoly;
        return n;
    }

    /**
     *  Returns the total number of singletons in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static int getSingletonsCount(SitesBlock[] blocks)
    {
        int n = 0;
        for(SitesBlock sb:blocks)
            n += sb.nSi;
        return n;
    }

    /**
     *  Returns the total number of transitions in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static float getTransitionsCount(SitesBlock[] blocks)
    {
        float n = 0.0f;
        for(SitesBlock sb:blocks)
            n += sb.nTS;
        return n;
    }
}
