/*
    File:
        ADatasetAnalyzer.java
 *   
    Revision:
        1.1.0.10
 * 
    Description:
        Abstract basis class which all data set analyzers should be inherited from.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin.classes;

import bio.gene.Dataset;
import kernel.ErrorCode;
import plugin.PluginType;


public abstract class ADatasetAnalyzer implements IGAPlugin
{
    public PluginType GetType()
    {
        return PluginType.ANALYZER;
    }
    

    /**
     *  Analyses the specified data set. 
     *  It is possible to call the analyzer from within a GeneAnalyzer script (*.gas).
     *  Since the point of writing the script is the automatization of analysis
     *  process, the user must provide correct parameters the analyzer might need to
     *  perform correctly. Otherwise the analyzer should show an options window.
     *  The host program does not test the correctness of the syntax!
     * 
     *  @param dataset  data set to analyze
     *  @param params   a string containing the set of parameters in any
     *                  format the analyzer object can understand. Can be empty or null.
     *  @return
     */
    public abstract ErrorCode AnalyzeDataset(Dataset dataset, String params);
}
