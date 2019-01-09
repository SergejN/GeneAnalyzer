/*
    File:
        Dataset.java
 *   
    Revision:
        2.2.0.2
 * 
    Description:
        Represents a complete dataset, including gene entries.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import kernel.QualityChecker;


public class Dataset extends DataChunk
{
    /**
     *  Genes of the data set.
     */
    private Vector<GeneEntry> genes;

    
    /**
     *  Constructor.
     */
    public Dataset()
    {
        properties = new HashMap<String, Object>();
        genes = new Vector<GeneEntry>();
    }
    
    
    /**
     *  Returns the number of genes in the data set.
     * 
     *  @return
     */
    public int getGenesCount()
    {
        return genes.size();
    }
    
    
    /**
     *  Returns the specified gene entry or null if the index is invalid.
     * 
     *  @param iIndex       gene entry index
     *  @return
     */
    public GeneEntry getGeneEntry(int iIndex)
    {
        if( (iIndex<genes.size()) && (iIndex>-1) )
            return genes.get(iIndex);
        else
            return null;
    }
    
    
    /**
     *  Adds the gene entry to the dataset. If a gene entry with the same
     *  common name exists (note, that the alias names are not compared),
     *  the new entry is not added, but its strains are added to the existing
     *  gene entry. If the gene entry is null, the method simply returns.
     * 
     *  @param ge
     */
    public void addGene(GeneEntry ge)
    {
        if(ge==null)
            return;
        // Iterate through the existing gene entries and compare
        // them to the new one.
        for(GeneEntry gene:genes)
        {
            // If the current gene and the new gene are the same object, return.
            if(gene==ge)
                return;
            // If the gene name of the existing and the new gene entries
            // is the same, add the strain entries to existing gene entry.
            if(gene.getCommonName().equalsIgnoreCase(ge.getCommonName()))
            {
                for(int i=0;i<ge.getStrainsCount();i++)
                    gene.addStrain(ge.getStrainEntry(i));
                // Since this is the only way to add a new gene entry to the
                // dataset, the function can return, since the uniqueness of
                // the gene entry is quaranted.
                return;
            }
        }
        // If the gene entry does not exist, add it.
        genes.add(ge); 
    }
    
    
    /**
     *  Removes the specified gene entry from the dataset.
     * 
     *  @param iIndex
     */
    public void removeGene(int iIndex)
    {
        if( (iIndex<genes.size()) && (iIndex>-1) )
            genes.remove(iIndex);
    }
    
    
    /**
     *  Checks whether or not the specified gene is already
     *  present in the data set.
     * 
     *  @param strName
     *  @return
     */
    public boolean hasGene(String strName)
    {
        // Iterate through the genes and compare the common names.
        for(GeneEntry gene:genes)
        {
            if(gene.getCommonName().equalsIgnoreCase(strName))
                return true;
        }
        return false;
    }
    
    
    /**
     *  Returns a non-redundant sorted list of populations present in the dataset.
     * 
     *  @return
     */
    public String[] listPopulations()
    {
        Vector<String> tmp = new Vector<String>();
        for(int i=0;i<genes.size();i++)
        {
            String[] pop = genes.get(i).listPopulations();
            for(String s:pop)
            {
                if(tmp.indexOf(s)==-1)
                    tmp.add(s);
            }
        }
        // Sort.
        String[] pops = tmp.toArray(new String[1]);
        Arrays.sort(pops);
        return pops;
    }
    
    
    /**
     *  Returns a non-redundand list of region names of all genes in the dataset.
     * 
     *  @return
     */
    public String[] listRegionsNames()
    {
        Vector<String> tmp = new Vector<String>();
        for(int i=0;i<genes.size();i++)
        {
            String[] regs = genes.get(i).listRegionsNames();
            for(String s:regs)
            {
                if(tmp.indexOf(s)==-1)
                    tmp.add(s);
            }
        }
        return tmp.toArray(new String[1]);
    }

    
    /**
     *  Merges the dataset with the specified dataset. Redundand gene entries,
     *  strain entries and regions are treated like in methods
     *  addGene, addStrain and addRegion. If the specified dataset is null,
     *  the method simply returns.
     * 
     *  @param other
     */
    public void merge(Dataset other)
    {
        if(other==null)
            return;
        for(int i=0;i<other.getGenesCount();i++)
        {
            addGene(other.getGeneEntry(i));
        }
    }

    /**
     *  Sorts the gene entries of the dataset by name.
     */
    public void sort()
    {
        Collections.sort(genes, new Comparator<GeneEntry>()
                {
                    public int compare(GeneEntry ge1, GeneEntry ge2)
                    {
                        return ge1.getCommonName().compareTo(ge2.getCommonName());
                    }
                });
    }

    /**
     *  Sorts the gene entries of the dataset by their quality level.
     */
    public void sortByQuality()
    {
        Collections.sort(genes, new Comparator<GeneEntry>()
                {
                    public int compare(GeneEntry ge1, GeneEntry ge2)
                    {
                        Integer q1 = (Integer)ge1.getProperty(QualityChecker.QUALITY_LEVEL);
                        Integer q2 = (Integer)ge2.getProperty(QualityChecker.QUALITY_LEVEL);
                        return q1.compareTo(q2);
                    }
                });
    }
}
