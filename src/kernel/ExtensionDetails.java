/*
    File:
        ExtensionDetails.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Holds the details about the supported extensions.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel;

import java.io.File;
import javax.swing.filechooser.FileFilter;

public class ExtensionDetails extends FileFilter
{
    private String[] ext   = null;     // Extensions.
    private String strDesc = null;     // Description.
    private String desc    = null;

    private int    index   = -1;
    
    /**
     *  Creates the filter.
     * 
     *  @param strExt
     *  @param strDesc
     *  @param index
     */
    public ExtensionDetails(String strExt, String strDesc, int index)
    {
        this.index = index;
        this.strDesc = strDesc;
        // Check for null pointers.
        if(strExt!=null && !strExt.isEmpty()) // Destination is a file
        {
            this.ext = strExt.split("\\|");
            // Construct description line.
            StringBuffer sb = new StringBuffer();
            sb.append(strDesc+" (");
            for(int i=0;i<ext.length-1;i++)
                sb.append(String.format("*.%s, ", ext[i]));
            sb.append("*."+ext[ext.length-1]+")");
            desc = sb.toString();
        }
        else // Destination is a directory
        {
            this.ext = null;
            desc = strDesc+" (Folder only)";
        }
    }

    @Override
    public boolean accept(File f)
    {
        if (f != null) 
        {
            // Accept directories.
            if (f.isDirectory()) 
            {
                return true;
            }
            // Accept file, if it has the supported extension.
            // If the extension is not specified, the OpenFileDialog
            // is configured to display folders only, so that true is returned
            // in the previous step.
            if(ext!=null)
            {
                String fileName = f.getName();
                int i = fileName.lastIndexOf('.');
                if (i>0 && i<fileName.length()-1)
                {
                    String fe = fileName.substring(i+1);
                    for(String s:ext)
                    {
                        if(fe.equalsIgnoreCase(s))
                            return true;
                    }
                }
            }
        }
        return false;
    }

    @Override
    public String getDescription()
    {
        return desc;
    }
    
    /**
     *  Returns the file extension(s). 
     * 
     *  @return
     */
    public String[] getExtention()
    {
        return ext;
    }

    public int getIndex()
    {
        return index;
    }
}
