/*
    File:
        ADatasetImporter.java
 *   
    Revision:
        1.1.0.10
 * 
    Description:
        Abstract basis class which all data set importers should be inherited from.
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
import plugin.PluginType;


public abstract class ADatasetImporter implements IGAPlugin 
{
    public PluginType GetType()
    {
        return PluginType.IMPORTER;
    }

    
    /**
     *  Returns the extension(s) (without the '.') of the file(s), the class can import. 
     *  Multiple extensions must be separated by '|'. This method should not return null.
     * 
     *  @return
     */
    public abstract String GetFileExtension();
    
    
    /**
     *  Returns the description of the file(s), the class can import. If multiple
     *  file formats are supported, the description should provide a general term
     *  for all of them.
     * 
     *  @return
     */
    public abstract String GetFileDescription();
    
    
    /**
     *  Imports the specified files and combines the data into one big dataset
     *  which is then returned. The files parameter can contain both files and
     *  directories. If a file object represents a file of supported format
     *  it is imported and added to the dataset. If it represents a directory
     *  all supported files from this directory and its underdirectories are
     *  imported.
     * 
     *  @param files        files to import
     *  @param strParams    additional parameters. Can be empty or null.
     *  @return
     */
    public abstract Dataset ImportDataset(File[] files, String strParams);
}
