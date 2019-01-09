/*
    File:
        Codon.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Represents a single codon.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene.dna;

import bio.gene.DataChunk;
import java.util.HashMap;

public class Codon extends DataChunk
{
    private static HashMap<String, Codon> codons = generateCodonsList();

    private final static int NEIGHBORSCOUNT = 9;
    
    private String strSequence  = null;
    private Codon[] neighbors   = null;    


    /**
     *  Creates the Codon instance. The constructor is not accessible for foreign
     *  classes to avoid creation of non-existing codons.
     * 
     *  @param strSequence
     *  @param ct
     */
    private Codon(String strSequence)
    {
        properties = new HashMap<String, Object>();
        this.strSequence = strSequence;
    }

    private void generateNeighbors()
    {
        neighbors = new Codon[NEIGHBORSCOUNT];
        int index = 0;
        char[] bases = new char[]{'A','C','G','T'};
        // For each codon site (0,1,2) generate a codon which only
        // differs in one base at that very site.
        for (int i=0;i<3;i++)
        {
            // Interate through all 4 possible bases.
            for (int j=0;j<4;j++)
            {
                char[] tmp = strSequence.toCharArray();
                tmp[i] = bases[j];
                String strSeq = new String(tmp);
                // If the new codon is not the same as the current one
                // add it to the neighbors list.
                if (!strSeq.equalsIgnoreCase(strSequence))
                {
                    neighbors[index] = getCodon(strSeq);
                    index++;
                }
            }
        }
    }
    
    public String getSequence()
    {
        return strSequence;
    }
    
    /**
     *  Returns the number of neighbors.
     * 
     *  @return
     */
    public int getNeighborsCount()
    {
        return NEIGHBORSCOUNT;
    }
    
    /**
     *  Returns the specified neighbor.
     * 
     *  @param index
     *  @return
     */
    public Codon getNeighbor(int index)
    {
        if(index<0 || index>NEIGHBORSCOUNT-1)
            return null;
        return neighbors[index];
    }

    /**
     *  Calculates the number of syn. and nonsyn. sites in the codon.
     *  The flag bTerm specifies whether or not to exclude the terminal
     *  codons (neighbors) when estimating the number of syn. sites.
     *  Algorithm:
     *      There are no synonymous substitutions at position 2 of any codon.
     *      Thus, this position should not be evaluated in order to reduce
     *      computation costs.
     *        Once the number of synonymous sites of the codon is calculated, the
     *      value is added to the codon property list, in order to not re-calculate
     *      it. Thus, for any number of analyzed codons, the CPU-time consuming
     *      calculations are carried out max. 64 times (by the number of codons).
     *        The calculation itself is described in the paper of Nei: for each
     *        codon position the fraction of synonymous substitutions is calculated.
     *
     *  @param codon
     *  @param ct
     *  @param bTerm
     *  @return number of syn. and nonsyn. sites or null if the site contains a gap.
     */
    public float[] calculateNumberOfSites(ICodonTable ct, boolean bTerm)
    {
        // Check the codon.
        if(strSequence.contains("-"))
            return null;
        // Try to get the annotated value.
        String pref = ct.getName(); // Prefix: the codon table name
        Float sites = (Float)getProperty((bTerm) ? "SynSites_term_"+pref : "SynSites_noterm_"+pref);
        if(sites!=null)
            return new float[]{sites, 3.0f-sites};
        sites = new Float(0.0f);
        // If the value is not annotated, calculate it.
        char[] bases = new char[]{'A','C','G','T'};
        for (int i=0;i<3;i++)
        {
            int nMatches = 0;   // Number of syn. substitutions.
            int nSubst = 0;     // Total number of substitutions.
            for(char base:bases)
            {
                char[] chars = strSequence.toCharArray();
                chars[i] = base;
                String strNewCodon = new String(chars);
                if(strNewCodon.equalsIgnoreCase(strSequence))
                    continue;
                if (!ct.isTerminal(strNewCodon) || bTerm)
                {
                    if(ct.areSynonymous(strNewCodon, strSequence))
                        nMatches++;
                    nSubst++;
                }
            }
            sites  += (float) nMatches / (float)nSubst;
        }
        addProperty((bTerm) ? "SynSites_term_"+pref : "SynSites_noterm_"+pref, sites);
        return new float[]{sites, 3.0f-sites};
    }

    /**
     *  Returns the codon with the specified sequence. If the sequence is invalid
     *  the method returns null.
     *
     *  @param strSequence
     *  @return
     */
    public static Codon getCodon(String strSequence)
    {
        if(!strSequence.matches("[ACGTacgt]{3}"))
            return null;
        Codon c = codons.get(strSequence);
        if(c.neighbors==null)
            c.generateNeighbors();
        return c;
    }

    /**
     *  Generates the codons network and assigns every codon its neighbor.
     *
     *  @return
     */
    private static HashMap<String, Codon> generateCodonsList()
    {
        HashMap<String, Codon> network = new HashMap<String, Codon>();
        String[] seqs = generateCodonSequences();
        for(String strSeq:seqs)
            network.put(strSeq, new Codon(strSeq));
        return network;
    }

    /**
     *  Generates the 64 possible codon sequences.
     *
     *  @return
     */
    public static String[] generateCodonSequences()
    {
        String[] tmp = new String[64];
        int iIndex = 0;
        char[] bases = {'A', 'C', 'G', 'T'};
        // Generate all possible codons.
        int[] idx = {0,0,0};
        while(idx[0]<4)
        {
            tmp[iIndex] = new String(new char[]{bases[idx[0]], bases[idx[1]], bases[idx[2]]});
            // Update indices.
            // 3rd base.
            if(idx[2]<3)
                idx[2]++;
            else
            {
                idx[2]=0;
                idx[1]++;
            }
            // 2nd base.
            if(idx[1]==4)
            {
                idx[1]=0;
                idx[0]++;
            }
            iIndex++;
        } 
        return tmp;
    }
}