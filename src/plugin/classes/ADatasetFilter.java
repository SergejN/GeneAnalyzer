/*
    File:
        ADatasetFilter.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Abstract basis class which all data set filters should be inherited from.
        A filter iterates through the data set and returns the indices or a new 
        data set of genes, satisfying certain criteria.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin.classes;

import bio.gene.Dataset;
import plugin.PluginType;

public abstract class ADatasetFilter implements IGAPlugin
{
    public PluginType GetType()
    {
        return PluginType.FILTER;
    }    
    
    /**
     *  Returns the indices of the genes in the data set, which satisfy
     *  certain criteria.
     *  
     *  Return value:
     *  The method should return null, if an execution error occurs and
     *  a non-empty array with the first element -1, if the user cancels
     *  the filter.
     *
     *  @param dataset      data set to filter
     *  @param strParams    any additional parameters. Can be empty or null.
     *  @return
     */
    public abstract int[] GetFilteredIndices(Dataset dataset, String strParams);
}
