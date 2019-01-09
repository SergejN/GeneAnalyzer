/*
    File:
        StrainEntry.java
 *   
    Revision:
        2.2.0.2
 * 
    Description:
        Represents a strain entry of a particular gene.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene;

import java.util.Collections;
import java.util.HashMap;
import java.util.Vector;


public class StrainEntry extends DataChunk
{
    /**
     *  Species.
     */
    private String strSpecies = "";
    
    /**
     *  Strain.
     */
    private String strStrain = "";
    
    /**
     *  Chromosome, e.g. X. Is used to differenciate between genes on
     *  different chromosomes. Thus, different chrosomome for different
     *  species or strain is possible.
     *  
     *  NOTE: do not put chromosomal location here.
     */
    private String strChromosome = "";

    /**
     *  Population(s). Multiple populations are possible.
     */
    private Vector<String> populations = null;

    /**
     *  Vector of GeneRegion objects which represent the sequence
     *  of the entry.
     */
    private Vector<GeneRegion> sequence = null;
    
    
    /**
     *  Constructor.
     * 
     *  @param strSpecies
     *  @param strStrainName
     */
    public StrainEntry(String strSpecies, String strStrainName)
    {
        properties = new HashMap<String, Object>();
        this.strSpecies = (strSpecies!=null) ? strSpecies : "";
        this.strStrain = (strStrainName!=null) ? strStrainName : "";
        sequence = new Vector<GeneRegion>();
        populations = new Vector<String>();
    }    

    /**
     *  Sets species name. If strSpecies is null, the species name is
     *  removed.
     * 
     *  @param strSpecies
     */
    public void setSpeciesName(String strSpecies)
    {
        this.strSpecies = (strSpecies!=null) ? strSpecies : "";
    }
    
    /**
     *  Returns species name.
     * 
     *  @return
     */
    public String getSpeciesName()
    {
        return strSpecies;
    }
    
    /**
     *  Sets the strain name. If strStrainName is null, the strain name
     *  is removed.
     * 
     *  @param strStrainName
     */
    public void setStrainName(String strStrainName)
    {
        this.strStrain = (strStrainName!=null) ? strStrainName : "";
    }
    
    /**
     *  Returns the strain name.
     * 
     *  @return
     */
    public String getStrainName()
    {
        return strStrain;
    }
    
    /**
     *  Sets the chromosome name. If strChromosome is null, the chromosome
     *  name is removed.
     * 
     *  @param strChromosome
     */
    public void setChromosome(String strChromosome)
    {
        this.strChromosome = (strChromosome!=null) ? strChromosome : "";
    }
    
    /**
     *  Returns the chromosome name.
     * 
     *  @return
     */
    public String getChromosome()
    {
        return strChromosome;
    }
    
    /**
     *  Assigns the populations to the strain. If a name is null
     *  it is not added to the list. 
     * 
     *  @param populations
     */
    public void addPopulations(String... populations)
    {
        if (populations != null)
        {
            for (String s : populations)
            {
                if(s!=null)
                    this.populations.add(s);
            }
        }
    }
    
    /**
     *  Returns the array with populations names.
     * 
     *  @return
     */
    public String[] listPopulations()
    {
        if (populations.size() == 0)
        {
            return new String[]{};
        }
        else
        {
            return populations.toArray(new String[1]);
        }
    }
    
    /**
     *  Removes the specified populations from the populations list.
     * 
     *  @param populations
     */
    public void removePopulations(String... populations)
    {
        if (populations != null)
        {
            for (String s : populations)
            {
                this.populations.remove(s);
            }
        }
    }    
    
    /**
     *  Returns true if the strain belongs to the population specified and
     *  false otherwise. If the strPop is null or empty, true is returned for
     *  convenience.
     * 
     *  @param strPop
     *  @return
     */
    public boolean belongsToPopulation(String strPop)
    {
        if ( (strPop == null) || (strPop.isEmpty()) )
        {
            return true;
        }
        for(String s:populations)
        {
            if(s.equalsIgnoreCase(strPop))
                return true;
        }
        return false;
    }
    
    /**
     *  Returns the total regions count.
     * 
     *  @return             the total number of regions
     */
    public int getRegionsCount()
    {
        return sequence.size();
    }
        
    /**
     *  Returns the number of regions of the specified type.
     *  If the parameter strType is null or empty the behavior
     *  is equivalent to that of getRegionsCount without parameters.
     *  Since the type check is performed for each region, this function
     *  is slower than the one without parameters.
     * 
     *  @param strType          region type
     *  @return
     */
    public int getRegionsCount(String strType)
    {
        int n = 0;
        for(GeneRegion gr:sequence)
        {
            if(gr.hasType(strType))
                n++;
        }
        return n;
    }
    
    /**
     *  Returns the iIndex-th region or null if index is invalid.
     * 
     *  @param iIndex           region index
     *  @return
     */
    public GeneRegion getRegion(int iIndex)
    {
        if( (iIndex<sequence.size()) && (iIndex>-1) )
            return sequence.get(iIndex);
        else
            return null;
    }
    
    /**
     *  Returns the gene region the specified site is in. If the site
     *  is less than 1 or greater than the total sequence length, the
     *  method returns null.
     * 
     *  Note:
     *  The sites are 1-based.
     * 
     *  @param iSite
     *  @return
     */
    public GeneRegion getRegionAtSite(int iSite)
    {
        if(iSite<1)
            return null;
        for(GeneRegion reg:sequence)
        {
            if(reg.getStart()<=iSite && reg.getEnd()>=iSite)
                return reg;
        }
        return null;
    }
    
    
    /**
     *  Adds the specified region to the sequence. If the region is null,
     *  the method returns without adding the region. 
     *  NOTE: If an overlapping region already exists 
     *  within the sequence, it is not removed.
     * 
     *  @param strVariant       splice variant
     *  @param region           region to add
     */
    public void addRegion(GeneRegion region)
    {
        if(region==null)
            return;
        sequence.add(region); 
        // Sort the regions.
        Collections.sort(sequence);
    }
    
    /**
     *  Removes the region with the specified index. If the index is
     *  invalid the regions list is untouched.
     * 
     *  @param iIndex
     */
    public void removeRegion(int iIndex)
    {
        if( (iIndex<sequence.size()) && (iIndex>-1) )
            sequence.remove(iIndex);
    }    
    
    /**
     *  Returns the names of the regions within the sequence.
     * 
     *  @return
     */
    public String[] listRegionsNames()
    {
        if(sequence.size()==0)
            return new String[]{};
        String[] names = new String[sequence.size()];
        for(int i=0;i<sequence.size();i++)
            names[i] = sequence.get(i).getType();
        return names;
    }
    
    /**
     *  Returns the complete gene DNA sequence, i.e. the string concatinated
     *  of all region sequences.
     *
     *  NOTE: retrieving the sequence of a not completely initialized strains
     *  is a possible bug source since the retrieved sequence length may differ
     *  from the expected.
     * 
     *  @return
     */
    public String getCompleteSequence()
    {
        StringBuffer sb = new StringBuffer();
        for(GeneRegion reg:sequence)
            sb.append(reg.getSequence());
        return sb.toString();
    }
    
    /**
     *  Returns the coding sequence of the gene, i.e. the sequence which 
     *  is translated into a protein. Thus, the sequence returned does not 
     *  contain any non-exonic sequences.
     *
     *  NOTE: retrieving the sequence of a not completely initialized strains
     *  is a possible bug source since the retrieved sequence length may differ
     *  from the expected.
     *
     *  @return
     */
    public String getCodingSequence()
    {
        StringBuffer sb = new StringBuffer();
        for(GeneRegion reg:sequence)
        {
            if(reg.hasType(GeneRegion.EXON))
                sb.append(reg.getSequence());
        }
        return sb.toString();
    }
        
    /**
     *  Concatinates the sequences of the specified region types and
     *  returns the resulting sequence.
     *
     *  NOTE: retrieving the sequence of a not completely initialized strains
     *  is a possible bug source since the retrieved sequence length may differ
     *  from the expected.
     *
     *  Remarks:
     *  regs must not contain redundand entries, otherwise the resulting sequence
     *  will contain multiple copies of the region sequence.
     * 
     *  @param regs
     *  @return
     */
    public String getSequence(String... regs)
    {
        StringBuffer sb = new StringBuffer();
        for(GeneRegion reg:sequence)
        {
            for(String r:regs)
            {
                if(reg.hasType(r))
                {
                    sb.append(reg.getSequence());
                    break;
                }
            }
        }
        return sb.toString();
    }
}
