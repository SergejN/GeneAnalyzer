/*
    File:
        PluginType.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Possible GeneAnalyzer plugin types.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin;

public enum PluginType 
{
    ANALYZER,       // ADatasetAnalyzer instance
    EXPORTER,       // ADatasetExporter instance
    IMPORTER,       // ADatasetImporter instance
    FILTER          // ADatasetFilter instance   
}
