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

import bio.gene.dna.ICodonTable;


public class CodonsBlock
{
    private int nStrains = 0;
    private float pi_s   = 0.0f;
    private float pi_n   = 0.0f;
    private float sites_s= 0.0f;
    private float sites_n= 0.0f;
    private int nPoly_s  = 0;
    private int nPoly_n  = 0;
    private int nSi_s    = 0;
    private int nSi_n    = 0;
    private int nTS_s    = 0;
    private int nTS_n    = 0;
    private int nTV_s    = 0;
    private int nTV_n    = 0;

    private boolean bJC_pi   = false;
    private boolean bJC_t    = false;
    private float cof        = 1.0f;
    private boolean bUseTerm = false;
    private ICodonTable ct   = null;

    /**
     *  Constructs the codons block.
     *
     *  @param nStrains
     *  @param bJC_pi
     *  @param bJC_theta
     *  @param cof
     *  @param bUseTerm     whether to use terminal codons when calculating the statistics
     *  @param ct
     */
    public CodonsBlock(int nStrains, boolean bJC_pi, boolean bJC_theta, float cof, boolean bUseTerm, ICodonTable ct)
    {
        this.nStrains = nStrains;
        this.bJC_pi = bJC_pi;
        this.bJC_t = bJC_theta;
        this.cof = cof;
        this.bUseTerm = bUseTerm;
        this.ct = ct;
    }

    /**
     *  Adds the site to the block. If the site has different number of strains
     *  than defined for this block, the method does nothing and returns false.
     *
     *  @param sc
     *  @return
     */
    public boolean addCodon(CodonComposition cc)
    {
        if(cc.getValidCodonsCount()!=nStrains)
            return false;
        float[] sites = cc.getSitesCounts();
        float[] tmp = BasicStatistics.calculatePi(cc, ct, bUseTerm);
        pi_s += tmp[0]*sites[0];
        pi_n += tmp[1]*sites[1];
        sites_s += sites[0];
        sites_n += sites[1];
        int[] poly = cc.getNumberOfPolymorphisms();
        nPoly_s += poly[0];
        nPoly_n += poly[1];
        int[] si = cc.getNumberOfSingletons(cof);
        nSi_s += si[0];
        nSi_n += si[1];
        int[] subst = cc.getSubstitutionsCount();
        nTS_s += subst[0];
        nTV_s += subst[1];
        nTS_n += subst[2];
        nTV_n += subst[3];
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
    public float[] getPi()
    {
        float[] res = {0.0f, 0.0f};
        if(sites_s>0)
            res[0] = pi_s/sites_s;
        if(sites_n>0)
            res[1] = pi_n/sites_n;
        if(bJC_pi)
        {
            res[0] = BasicStatistics.correctJC(res[0]);
            res[1] = BasicStatistics.correctJC(res[1]);
        }
        return res;
    }

    /**
     *  Returns theta of the block.
     *
     *  @return
     */
    public float[] getTheta()
    {
        float[] res = {0.0f, 0.0f};
        res[0] = BasicStatistics.calculateTheta(nPoly_s, sites_s, nStrains);
        res[1] = BasicStatistics.calculateTheta(nPoly_n, sites_n, nStrains);
        if(bJC_t)
        {
            res[0] = BasicStatistics.correctJC(res[0]);
            res[1] = BasicStatistics.correctJC(res[1]);
        }
        return res;
    }

    /**
     *  Returns the number of synonymous and nonsynonymous sites in the block.
     *
     *  @return
     */
    public float[] getSitesCount()
    {
        return new float[]{sites_s, sites_n};
    }

    /**
     *  Returns the number of synonymous and nonsynonymous polymorphic sites in the block.
     *
     *  @return
     */
    public int[] getPolymorphismsCount()
    {
        return new int[]{nPoly_s, nPoly_n};
    }

    /**
     *  Returns the number of synonymous and nonsynonymous singletons in the block.
     *
     *  @return
     */
    public int[] getSingletonsCount()
    {
        return new int[]{nSi_s, nSi_n};
    }

    /**
     *  Returns the number of synonymous and nonsynonymous transitions in the block.
     *
     *  @return
     */
    public int[] getTransitionsCount()
    {
        return new int[]{nTS_s, nTS_n};
    }

    /**
     *  Returns the number of synonymous and nonsynonymous transversions in the block.
     *
     *  @return
     */
    public int[] getTransversionsCount()
    {
        return new int[]{nTV_s, nTV_n};
    }

    /**
     *  Returns the overall Pi values of all blocks.
     *
     *  @param blocks
     *  @param bJC
     *  @return
     */
    public static float[] getPi(CodonsBlock[] blocks, boolean bJC)
    {
        float pi_s = 0.0f;
        float pi_n = 0.0f;
        float n_s  = 0.0f;
        float n_n  = 0.0f;
        for(CodonsBlock cb:blocks)
        {
            n_s  += cb.sites_s;
            n_n  += cb.sites_n;
            pi_s += cb.pi_s;
            pi_n += cb.pi_n;
            
        }
        float[] res = {(n_s>0) ? pi_s/n_s : 0.0f, (n_n>0) ? pi_n/n_n : 0.0f};
        if(bJC)
        {
            res[0] = BasicStatistics.correctJC(res[0]);
            res[1] = BasicStatistics.correctJC(res[1]);
        }
        return res;
    }

    /**
     *  Returns the overall Theta values of all blocks.
     *
     *  @param blocks
     *  @param bJC
     *  @return
     */
    public static float[] getTheta(CodonsBlock[] blocks, boolean bJC)
    {
        float t_s = 0.0f;
        float t_n = 0.0f;
        float sites_s = 0.0f;
        float sites_n = 0.0f;
        for(CodonsBlock cb:blocks)
        {
            float f_s = BasicStatistics.calculateTheta(cb.nPoly_s, cb.sites_s, cb.nStrains);
            float f_n = BasicStatistics.calculateTheta(cb.nPoly_n, cb.sites_n, cb.nStrains);
            if(bJC)
            {
                f_s = BasicStatistics.correctJC(f_s);
                f_n = BasicStatistics.correctJC(f_n);
            }
            t_s += f_s*cb.nPoly_s;
            t_n += f_n*cb.nPoly_n;
            sites_s += cb.nPoly_s;
            sites_n += cb.nPoly_n;
        }
        return new float[]{(sites_s>0) ? t_s/sites_s : 0.0f, 
                           (sites_n>0) ? t_n/sites_n : 0.0f};
    }

    /**
     *  Returns the total numbers of sites in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static float[] getSitesCount(CodonsBlock[] blocks)
    {
        float n_s = 0.0f;
        float n_n = 0.0f;
        for(CodonsBlock cb:blocks)
        {
            float[] tmp = cb.getSitesCount();
            n_s += tmp[0];
            n_n += tmp[1];
        }
        return new float[]{n_s, n_n};
    }

    /**
     *  Returns the total numbers of polymorphisms in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static int[] getPolymorphismsCount(CodonsBlock[] blocks)
    {
        int n_s = 0;
        int n_n = 0;
        for(CodonsBlock cb:blocks)
        {
            int[] tmp = cb.getPolymorphismsCount();
            n_s += tmp[0];
            n_n += tmp[1];
        }
        return new int[]{n_s, n_n};
    }

    /**
     *  Returns the total numbers of singletons in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static int[] getSingletonsCount(CodonsBlock[] blocks)
    {
        int n_s = 0;
        int n_n = 0;
        for(CodonsBlock cb:blocks)
        {
            int[] tmp = cb.getSingletonsCount();
            n_s += tmp[0];
            n_n += tmp[1];
        }
        return new int[]{n_s, n_n};
    }

    /**
     *  Returns the total number of transitions in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static int[] getTransitionsCount(CodonsBlock[] blocks)
    {
        int n_s = 0;
        int n_n = 0;
        for(CodonsBlock cb:blocks)
        {
            int[] tmp = cb.getTransitionsCount();
            n_s += tmp[0];
            n_n += tmp[1];
        }
        return new int[]{n_s, n_n};
    }

    /**
     *  Returns the total numbers of transversions in all blocks.
     *
     *  @param blocks
     *  @return
     */
    public static int[] getTransversionsCount(CodonsBlock[] blocks)
    {
        int n_s = 0;
        int n_n = 0;
        for(CodonsBlock cb:blocks)
        {
            int[] tmp = cb.getTransversionsCount();
            n_s += tmp[0];
            n_n += tmp[1];
        }
        return new int[]{n_s, n_n};
    }
}
