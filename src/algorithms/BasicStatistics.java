/*
    File:
        BasicStatistics.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Performs some very basic sequence analyses.
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


public class BasicStatistics 
{
    /**
     *  Corrects the value using Jukes-Cantor correction.
     * 
     *  @param d 
     *  @return 
     */
    public static float correctJC(float d)
    {
        if(d>=0.75f || d<0)
            return d;
        return -0.75f*(float)Math.log(1-4*d/3);
    }
    
    /**
     *  Calculates the Theta value.
     * 
     *  @param nP         number of polymorphic sites
     *  @param sites      length of analyzed sequences/number of sites analyzed
     *  @param nSeqs      number of analyzed sequences
     *  @return 
     */
    public static float calculateTheta(int nP, float sites, int nSeqs)
    {
        // For a single sequence theta cannot be calculated.
        if(nSeqs<2 || nP==0)
            return 0.0f;
        // Sum(1/nSeq-1)
        float a = 0;
        for(float i=1.0f;i<nSeqs;i++)
            a += 1.0f/i;        
        return ((float)nP/sites)/a;
    }
    
    /**
     *  Calculates Tajima's D value.
     * 
     *  @param pi       Pi
     *  @param theta    Theta
     *  @param sites    length of the analyzed sequences/number of sites analyzed
     *  @param nP       number of polymorphic sites
     *  @param nSeqs    number of sequences
     *  @return Tajima's D
     */
    public static float calculateTajD(float pi, float theta, float sites, int nP, int nSeqs)
    {
        // Tajima's D cannot be calculated for less than 4 sequences.
        if(nSeqs<4 || nP==0)
            return Float.NaN;
        // First, some values are pre-calculated for the later use.
        float a1 = 0.0f;
        float a2 = 0.0f;
        for(float i=1.0f;i<nSeqs;i++)
        {
            a1 += 1.0f/i;
            a2 += 1.0f/(i*i);
        }          
        float b1 = (float)(nSeqs+1)/(float)(3*(nSeqs-1));
        float b2 = (float)(2.0f*(nSeqs*nSeqs+nSeqs+3))/(float)(9.0f*nSeqs*(nSeqs-1));
        float c1 = b1 - 1.0f/a1;
        float c2 = b2 - (float)(nSeqs+2)/(nSeqs*a1)+a2/(a1*a1);
        float dc1 = c1/a1;
        float dc2 = c2/(a1*a1+a2);
        // Tajima's D.
        float numerator = sites*(pi-theta); 
        float denominator = (float)Math.sqrt(dc1*nP+dc2*nP*(nP-1));
        return numerator/denominator;
    }
    
    /**
     *  Calculates Tajima's D for multiple sites blocks of non-coding region(s).
     *  Since Tajima's D cannot be calculated for the sample size less than 4, the
     *  blocks with less than 4 strains are ignored.
     *
     *  The formula used is similar to the original Tajima's formula:
     *
     *          sites1*(pi1-theta1) + sites2*(pi2-theta2) + ...
     *  TajD = -------------------------------------------------
     *                      SQRT(V1 + V2 + ...)
     *
     *  whereas V is a variance.
     *
     *  @param blocks
     *  @return
     */
    public static float calculateTajD(SitesBlock[] blocks)
    {
        float numerator = 0.0f;
        float denominator = 0.0f;
        for(SitesBlock sb:blocks)
        {
            int nSeqs = sb.getStrainsCount();
            if(nSeqs<4)
                continue;
            numerator += sb.getSitesCount()*(sb.getPi()-sb.getTheta());
            float a1 = 0.0f;
            float a2 = 0.0f;
            for(float i=1.0f;i<nSeqs;i++)
            {
                a1 += 1.0f/i;
                a2 += 1.0f/(i*i);
            }
            float b1 = (float)(nSeqs+1)/(float)(3*(nSeqs-1));
            float b2 = (float)(2.0f*(nSeqs*nSeqs+nSeqs+3))/(float)(9.0f*nSeqs*(nSeqs-1));
            float c1 = b1 - 1.0f/a1;
            float c2 = b2 - (float)(nSeqs+2)/(nSeqs*a1)+a2/(a1*a1);
            float dc1 = c1/a1;
            float dc2 = c2/(a1*a1+a2);
            int nP = sb.getPolymorphismsCount();
            denominator += dc1*nP+dc2*nP*(nP-1);
        }
        return numerator/(float)Math.sqrt(denominator);
    }

    /**
     *  Calculates Tajima's D for multiple sites blocks of coding region(s).
     *  Since Tajima's D cannot be calculated for the sample size less than 4, the
     *  blocks with less than 4 strains are ignored.
     *
     *  The formula used is similar to the original Tajima's formula:
     *
     *          sites1*(pi1-theta1) + sites2*(pi2-theta2) + ...
     *  TajD = -------------------------------------------------
     *                      SQRT(V1 + V2 + ...)
     *
     *  whereas V is a variance.
     *
     *  @param blocks
     *  @return
     */
    public static float[] calculateTajD(CodonsBlock[] blocks)
    {
        float numerator_s   = 0.0f;
        float numerator_n   = 0.0f;
        float denominator_s = 0.0f;
        float denominator_n = 0.0f;
        for(CodonsBlock cb:blocks)
        {
            int nSeqs = cb.getStrainsCount();
            if(nSeqs<4)
                continue;
            float[] sites = cb.getSitesCount();
            float[] pis   = cb.getPi();
            float[] thetas= cb.getTheta();
            numerator_s += sites[0]*(pis[0]-thetas[0]);
            numerator_n += sites[1]*(pis[1]-thetas[1]);
            float a1 = 0.0f;
            float a2 = 0.0f;
            for(float i=1.0f;i<nSeqs;i++)
            {
                a1 += 1.0f/i;
                a2 += 1.0f/(i*i);
            }
            float b1 = (float)(nSeqs+1)/(float)(3*(nSeqs-1));
            float b2 = (float)(2.0f*(nSeqs*nSeqs+nSeqs+3))/(float)(9.0f*nSeqs*(nSeqs-1));
            float c1 = b1 - 1.0f/a1;
            float c2 = b2 - (float)(nSeqs+2)/(nSeqs*a1)+a2/(a1*a1);
            float dc1 = c1/a1;
            float dc2 = c2/(a1*a1+a2);
            int[] poly = cb.getPolymorphismsCount();
            denominator_s += dc1*poly[0]+dc2*poly[0]*(poly[0]-1);
            denominator_n += dc1*poly[1]+dc2*poly[1]*(poly[1]-1);
        }
        return new float[]{numerator_s/(float)Math.sqrt(denominator_s), numerator_n/(float)Math.sqrt(denominator_n)};
    }

    /**
     *  Calculates the TajD', which is defined as TajD' = (pi-theta)/(pi_min-theta),
     *  whereas pi_min is pi, if all polymorphic sites were singletons.
     * 
     *  @param pi       Pi
     *  @param theta    Theta
     *  @param sites    Number of sites
     *  @param nP       Number of polymorphic sites
     *  @param nSeqs    Number of sequences in the sample
     *  @return
     */
    public static float calculateTajDPrime(float pi, float theta, float sites, int nP, int nSeqs)
    {
        if(nSeqs<2 || nP==0)
            return Float.NaN;
        // First calculate pi_min.
        float pi_min = ((2.0f/(float)nSeqs)*(float)nP)/sites;
        return (pi-theta)/(pi_min-theta);
    }    

    /**
     *  Calculates Tajima's D' for multiple blocks of non-coding region(s).
     *  Since Tajima's D cannot be calculated for the sample size less than 4, the
     *  blocks with less than 4 strains are ignored.
     *
     *  @param blocks
     *  @return
     */
    public static float calculateTajDPrime(SitesBlock[] blocks)
    {
        float numerator = 0.0f;
        float denominator = 0.0f;
        for(SitesBlock sb:blocks)
        {
            int nSeqs = sb.getStrainsCount();
            int nP = sb.getPolymorphismsCount();
            if(nSeqs<4 || nP==0)
                continue;
            numerator += sb.getPi()-sb.getTheta();
            float pi_min = ((2.0f/(float)sb.getStrainsCount())*(float)sb.getPolymorphismsCount())/sb.getSitesCount();
            denominator += pi_min-sb.getTheta();
        }
        return numerator/denominator;
    }

    /**
     *  Calculates Tajima's D' for multiple blocks of coding region(s).
     *  Since Tajima's D cannot be calculated for the sample size less than 4, the
     *  blocks with less than 4 strains are ignored.
     *
     *  @param blocks
     *  @return
     */
    public static float[] calculateTajDPrime(CodonsBlock[] blocks)
    {
        float numerator_s   = 0.0f;
        float numerator_n   = 0.0f;
        float denominator_s = 0.0f;
        float denominator_n = 0.0f;
        for(CodonsBlock cb:blocks)
        {
            int nSeqs = cb.getStrainsCount();
            if(nSeqs<4)
                continue;
            int[] poly = cb.getPolymorphismsCount();            
            float[] pis = cb.getPi();
            float[] thetas = cb.getTheta();
            float[] sites = cb.getSitesCount();
            if(sites[0]>0.0f)
            {
                numerator_s += pis[0]-thetas[0];
                float pi_min_s = ((2.0f/(float)cb.getStrainsCount())*(float)poly[0])/sites[0];
                denominator_s += pi_min_s - thetas[0];
            }
            if(sites[1]>0.0f)
            {
                numerator_n += pis[1]-thetas[1];
                float pi_min_n = ((2.0f/(float)cb.getStrainsCount())*(float)poly[1])/sites[1];
                denominator_n += pi_min_n - thetas[1];
            }
        }
        return new float[]{numerator_s/denominator_s, numerator_n/denominator_n};
    }

    /**
     *  Calculates Pi for non-coding sites.
     * 
     *  Algorithm:
     *  In order to not to compare all sequences pairwise, you can make
     *  use of the information about the nucleotide counts.
     *  The number of pairwise mismatches is, thus, Nm=(N*N-(A*A+C*C+G*G+T*T))/2.
     *  The number of comparisons is Nc=N*(N-1)/2. Pi itself is then calculated
     *  as Nm/Nc.
     * 
     *  @param counts   bases counts
     *  @return Pi
     */
    public static float calculatePi(SiteComposition sc)
    {
        int nTotal = sc.getValidBasesCount();
        if(nTotal<2)
            return 0.0f;
        int as = sc.getBaseCount('A');
        int cs = sc.getBaseCount('C');
        int gs = sc.getBaseCount('G');
        int ts = sc.getBaseCount('T');
        int diffs = nTotal*nTotal - (as*as + cs*cs + gs*gs + ts*ts);
        return (float)diffs/(float)(nTotal*(nTotal - 1));
    }
    
    /**
     *  Calculates Pi for synonymous and nonsynonymous sites using the specified
     *  codons table. The first value in the returned array is Pi_syn and the
     *  second is Pi_nonsyn.
     * 
     *  Algorithm:
     *  Generates all possible pairs of codons and calculates the best evolutionary
     *  path between them using the PathFinder.getBestPath method. Using the best path,
     *  calculates the number of synonymous and nonsynonymous substitutions.
     * 
     *  @param cc
     *  @param ct
     *  @param bTerm    whether or not use terminal codons when generating a path
     *  @return
     */
    public static float[] calculatePi(CodonComposition cc, ICodonTable ct, boolean bTerm)
    {
        Codon[] codons = cc.getValidCodons();
        if(codons.length<2)
            return new float[]{0.0f, 0.0f};
        int ns = 0;         // Syn. substitutions
        int nn = 0;         // Nonsyn. substitutions
        float syn = 0.0f;   // Number of syn. sites
        float non = 0.0f;   // Number of nonsyn. sites
        // Generate the pairs and estimate the number of differences.
        for(int i=0;i<codons.length-1;i++)
        {
            Codon c1 = codons[i];
            float[] nos = c1.calculateNumberOfSites(ct, bTerm);
            syn += nos[0];
            non += nos[1];
            for(int j=i+1;j<codons.length;j++)
            {                
                Codon c2 = codons[j];                
                if(c1!=c2)
                {
                    Path bp = Path.findBestPath(new Codon[]{c1, c2}, ct, bTerm);
                    int[] tmp = (bp!=null) ? bp.getPolymorphismsCount() : new int[]{0,0};
                    ns += tmp[0];
                    nn += tmp[1];
                }
            }
        }
        // Estimate the number of syn. and nonsyn. sites of the last codon.
        float[] nos = codons[codons.length-1].calculateNumberOfSites(ct, bTerm);
        syn += nos[0];
        non += nos[1];
        float[] res = new float[2];
        int nPairs = codons.length*(codons.length-1)/2;
        res[0] = (syn>0.0f) ? ((float)ns/(float)nPairs)/(syn/(float)codons.length) : 0.0f;
        res[1] = (non>0.0f) ? ((float)nn/(float)nPairs)/(non/(float)codons.length) : 0.0f;
        return res;
    }
    
    /**
     *  Calculates the divergence between the population of interest
     *  and the outgroup.
     * 
     *  Algorithm:
     *  You can make use of information about the nucleotide counts.
     *  For each nucleotide X from the outgroup there are exactly N-X
     *  mismatches in the population of interest. The total number of
     *  mismatches between populations Nt is, thus, the sum of mismatches
     *  over all nucleotides. The number of comparisons Nc is simply
     *  the product of sample sizes. K is calculated as Nt/Nc.
     * 
     *  @param countpop
     *  @param countout
     *  @return
     */
    public static float calculateK(SiteComposition pop, SiteComposition out)
    {
        int nTotal1 = pop.getValidBasesCount();
        int nTotal2 = out.getValidBasesCount();
        int[] diffs = new int[]{0, 0, 0, 0};
        char[] bases = {'A', 'C', 'G', 'T'};
        for(int i=0;i<4;i++)
            diffs[i] = out.getBaseCount(bases[i])*(nTotal1-pop.getBaseCount(bases[i]));
        return ((float)(diffs[0]+diffs[1]+diffs[2]+diffs[3]))/(float)(nTotal1*nTotal2);
    }
    
    /**
     *  Calculates K for synonymous and nonsynonymous sites using the specified
     *  codons table. The first value in the returned array is Pi_syn and the
     *  second is Pi_nonsyn.
     * 
     *  Algorithm:
     *  Generates all possible pairs of codons of the population and the outgroup
     *  and calculates the best evolutionary path between them using the 
     *  PathFinder.getBestPath method. Using the best path, calculates the 
     *  number of synonymous and nonsynonymous substitutions.
     * 
     *  Remark:
     *  The number of synonymous and nonsynonymous substitutions of each codon
     *  pair is annotated to both codons. The name of the property is the sequence
     *  of the second codon.
     * 
     *  @param pop
     *  @param out
     *  @param ct
     *  @param bTerm    whether or not use terminal codons when generating a path
     *  @return
     */
    public static float[] calculateK(Codon[] pop, Codon[] out, ICodonTable ct, boolean bTerm)
    {
        float[] res = {0.0f, 0.0f};
        if(pop.length<1 || out.length<1)
            return res;
        int ns = 0;          // Syn. substitutions
        int nn = 0;          // Nonsyn. substitutions
        float syn = 0.0f;    // Number of syn. sites
        float non = 0.0f;    // Number of nonsyn. sites
        boolean flag = true; // Specifies whether or not to estimate the number
                             // of syn. and nonsyn. sites.
        for(Codon c1:pop)
        {
            float[] st = c1.calculateNumberOfSites(ct, bTerm);
            syn += st[0];
            non += st[1];
            for(Codon c2:out)
            {
                if(flag)
                {
                    float[] st2 = c2.calculateNumberOfSites(ct, bTerm);
                    syn += st2[0];
                    non += st2[1];
                }
                // Find the best path between the two codons if they are not equal.
                if(!c1.getSequence().equalsIgnoreCase(c2.getSequence()))
                {
                    Path bp = Path.findBestPath(new Codon[]{c1, c2}, ct, bTerm);
                    int[] tmp = (bp!=null) ? bp.getPolymorphismsCount() : new int[]{0,0};
                    ns += tmp[0];
                    nn += tmp[1];                    
                }
            }
            flag = false;
        }
        int nPairs = pop.length*out.length;
        float n = pop.length+out.length;
        res[0] = (syn>0.0f) ? ((float)ns/(float)nPairs)/(syn/n) : 0.0f;
        res[1] = ((float)nn/(float)nPairs)/(non/n);
        return res;
    }

    /**
     *  Returns the type of substitution:
     *      -1      invalid base, if either base1 or base2 is not a valid base: A,C,G,T
     *       0      if bases are equal
     *       1      if the substitution is a transition: A<->T or C<->G
     *       2      if the substitution is a transvertion: A/T<->C/G
     *       
     *
     *  @param base1
     *  @param base2
     *  @return
     */
    public static int getSubstitutionType(char base1, char base2)
    {
        char b1 = Character.toUpperCase(base1);
        if(b1!='A' && b1!='C' && b1!='G' && b1!='T')
            return -1;
        char b2 = Character.toUpperCase(base2);
        if(b2!='A' && b2!='C' && b2!='G' && b2!='T')
            return -1;
        if(b1==b2)
            return 0;
        if((b1=='A' && b2=='G') || (b1=='G' && b2=='A') || (b1=='C' && b2=='T') || (b1=='T' && b2=='C'))
            return 1;
        else
            return 2;
    }
}
