/*
    File:
        ADatasetExporter.java
 *   
    Revision:
        1.1.0.10
 * 
    Description:
        Abstract basis class which all data set exporters should be inherited from.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin.classes;

import bio.gene.Dataset;
import java.io.File;
import kernel.ErrorCode;
import plugin.PluginType;


public abstract class ADatasetExporter implements IGAPlugin
{
    public PluginType GetType()
    {
        return PluginType.EXPORTER;
    }

    
    /**
     *  Returns the extension (without the '.') of the file, the class can export. 
     *  One exporter class should support only one target file extension. If the class
     *  should export the data into the folder, null should be returned.
     * 
     *  @return
     */
    public abstract String GetFileExtension();
    
    
    /**
     *  Returns the description of the file, the class can export. 
     * 
     *  @return
     */
    public abstract String GetFileDescription();

    
    /**
     *  Exports the dataset into the specified file.
     * 
     *  @param ds           data set to export
     *  @param file         file/directory to save data into
     *  @param strParams    any additional parameters. Can be empty or null.
     *  @return
     */
    public abstract ErrorCode ExportDataset(Dataset ds, File file, String strParams);
}
