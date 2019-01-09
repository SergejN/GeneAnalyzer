/*
    File:
        GeneEntry.java
 *   
    Revision:
        2.2.0.2
 * 
    Description:
        Represents a gene entry. A gene entry holds the StrainEntry
        instances of a particular gene.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene;

import java.util.HashMap;
import java.util.Vector;

public class GeneEntry extends DataChunk
{
    /**
     *  Common and unique gene name, e.g. CG number.
     */
    private String strCommonName = "";
    
    /**
     *  Gene alias, i.e. TX1 for CGxxxx, which is used
     *  to refer to the gene more conviniently.
     */
    private String strAlias = "";
    
    /**
     *  Strains of the gene entry.
     */
    private Vector<StrainEntry> strains = null;
    
    
    /**
     *  Constructor.
     * 
     *  @param strCommonName
     *  @param strAlias
     */
    public GeneEntry(String strCommonName, String strAlias)
    {
        strains = new Vector<StrainEntry>();
        properties = new HashMap<String, Object>();
        this.strCommonName = strCommonName;
        this.strAlias = strAlias;
    }
    
    
    /**
     *  Sets the common name. If strCommonName is null, the common name
     *  is removed.
     * 
     *  @param strCommonName
     */
    public void setCommonName(String strCommonName)
    {
        this.strCommonName = (strCommonName!=null) ? strCommonName : "";
    }
    
    
    /**
     *  Returns the common name.
     * 
     *  @return
     */
    public String getCommonName()
    {
        return strCommonName;
    }

    
    /**
     *  Sets the alias. If strAlias is null, the alias is removed.
     * 
     *  @param strAlias
     */
    public void setAlias(String strAlias)
    {
        this.strAlias = (strAlias!=null) ? strAlias : "";
    }
    
    
    /**
     *  Returns the alias name.
     * 
     *  @return
     */
    public String getAlias()
    {
        return strAlias;
    }    
    
    /**
     *  Adds the strain entry to the gene entry. If a strain entry with
     *  the same species and strain names already exists, the entry is 
     *  not added, but its regions are added to the existing strain entry
     *  instead. If the strain entry is null, the method returns without
     *  adding the entry.
     *  NOTE: the regions to be added are not checked against redundancy. 
     * 
     *  @param se
     */
    public void addStrain(StrainEntry se)
    {
        if(se==null)
            return;
        // Iterate through the existing strain entries and compare
        // them to the new one.
        for(StrainEntry strain:strains)
        {
            // If the current strain and the new strain entry are the same
            // object, return.
            if(strain==se)
                return;
            // If the strain and the species names of the existing and the new
            // strain are the same, add the regions to the existing entry.
            if(strain.getSpeciesName().equalsIgnoreCase(se.getSpeciesName()) &&
               strain.getStrainName().equalsIgnoreCase(se.getStrainName()))
            {                
                for(int i=0;i<se.getRegionsCount();i++)
                    strain.addRegion(se.getRegion(i));
                // Since this is the only way to add a new strain entry to the
                // gene entry, the function can return, since the uniqueness of
                // the strain entry is quaranted.
                return;
            }
        }
        // If the strain entry does not exist, add it.
        strains.add(se);
    }
        
    /**
     *  Removes the specified strain entry.
     * 
     *  @param iIndex
     */
    public void removeStrain(int iIndex)
    {
        if( (iIndex<strains.size()) && (iIndex>-1) )
            strains.remove(iIndex);
    }

    /**
     *  Returns true, if the gene entry contains the strain entry with the specified name.
     *
     *  @param strName
     *  @return
     */
    public boolean hasStrain(String strName)
    {
        // Iterate through the genes and compare the common names.
        for(StrainEntry strain:strains)
        {
            if(strain.getStrainName().equalsIgnoreCase(strName))
                return true;
        }
        return false;
    }
    
    /**
     *  Returns the number of strains in the gene entry.
     * 
     *  @return
     */
    public int getStrainsCount()
    {
        return strains.size();
    }
        
    
    /**
     *  Returns the specified strain entry or null if the index is invalid.
     * 
     *  @param iIndex       strain entry index
     *  @return
     */
    public StrainEntry getStrainEntry(int iIndex)
    {
        if( (iIndex<strains.size()) && (iIndex>-1) )
            return strains.get(iIndex);
        else
            return null;
    }
    
    
    /**
     *  Returns a non-redundand list of population names of the
     *  strains in the gene entry.
     * 
     *  @return
     */
    public String[] listPopulations()
    {
        Vector<String> pops = new Vector<String>();
        for (int i = 0; i < strains.size(); i++)
        {
            String[] _pops = strains.get(i).listPopulations();
            for (String s : _pops)
            {
                if (pops.indexOf(s) == -1)
                {
                    pops.add(s);
                }
            }
        }
        return pops.toArray(new String[1]);
    }

    
    /**
     *  Returns a non-redundand list of region names of all
     *  strains in the gene entry.
     * 
     *  @return
     */
    public String[] listRegionsNames()
    {
        Vector<String> regs = new Vector<String>();
        for (int i=0; i<strains.size(); i++)
        {
            String[] _regs = strains.get(i).listRegionsNames();
            for (String s : _regs)
            {
                if (regs.indexOf(s) == -1)
                {
                    regs.add(s);
                }
            }
        }
        return regs.toArray(new String[1]);
    }
}
