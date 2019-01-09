/*
    File:
        SettingsManager.java
 *   
    Revision:
        1.2.0.10
 * 
    Description:
        Manages the application and plugin settings.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */
package kernel;

import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

/*
       XML File Format:

       <settings>
            <host.application>
                <setting.name>value</setting.name>
            </host.application>
            <plugin.name>
                <setting.name>value</setting.name>
            </plugin.name>
       </settings>
 */

public class SettingsManager 
{
    public static final String CODONTABLE   = "CodonTable";
    public static final String CODETYPE     = "CodeType";
    public static final String CT_SIMPLE    = "Simple";
    public static final String CT_EXTENDED  = "Extended";
    public static final String CT_COMPLETE  = "Complete";
    public static final String LOCALE       = "Locale";
    public static final String HEAPSIZE     = "CustomHeapSize";

    // Previous instance.
    private static SettingsManager prev = null;
    private Hashtable<String, String> values = null;

    private SettingsManager(String strFilename)
    {
        // Instance counter to avoid creation of multiple copies.
        values = new Hashtable<String, String>();
        // If the file exists, parse it.
        File f = new File(strFilename);
        if(f.exists())
        {
            SAXBuilder builder = new SAXBuilder();
            try
            {
                Document document = builder.build(f);
                Element rootElement = document.getRootElement();
                List<Element> blocks = rootElement.getChildren();
                for(Element el:blocks)
                {
                    List<Element> vals = el.getChildren();
                    for(Element v:vals)
                    {
                        String strKey = el.getName()+" "+v.getName();
                        values.put(strKey, v.getText());
                    }
                }
            }
            catch(Exception e){/* Ignore exceptions. */}
        }
        prev = this;
    }

    /**
     *  Saves the settings to the file.
     *
     *  @param strFilename      file name
     *  @param bFlush           flag specifying the action:
     *                           - true  the content is flushed to disk
     *                           - false the content is saved and the object
     *                             is destroyed.
     *  @return
     */
    public boolean save(String strFilename, boolean bFlush)
    {
        // Save the settings to disk.
        Element root = new Element("settings");
        String[] keys = values.keySet().toArray(new String[1]);
        if(keys[0]==null)
            return true;
        Arrays.sort(keys);
        Element curr = null;
        for(String strKey:keys)
        {
            // Split the key into block name and setting name.
            String[] tmp = strKey.split(" ");
            if(curr==null)
                curr = new Element(tmp[0]);
            // Current key belongs to the same block.
            if(curr.getName().equalsIgnoreCase(tmp[0]))
            {
                Element el = new Element(tmp[1]);
                el.setText(values.get(strKey));
                curr.addContent(el);
            }
            else    // Otherwise add current element to the file.
            {
                root.addContent(curr);
                curr = new Element(tmp[0]);
                Element el = new Element(tmp[1]);
                el.setText(values.get(strKey));
                curr.addContent(el);
            }
        }
        // Add last element.
        root.addContent(curr);
        // Save XML.
        try
        {
            FileWriter fw = new FileWriter(new File(strFilename));
            Format format = Format.getPrettyFormat();
            new XMLOutputter(format).output(new Document(root), fw);
            fw.close();
        }
        catch(Exception e)
        {
            return false;
        }
        if(bFlush)      // Only flush
            return true;
        else            // Destroy the object
        {
            prev = null;
            return true;
        }
    }

    /**
     *  Retrieves the value for the key strKey from the block strBlock.
     *  If strBlock is null or empty, the block is assumed to be the main
     *  settings block (see Remarks). If strKey is empty or does not exist,
     *  the method returns null.
     *
     *  Remarks:
     *      The plugins should not save their settings in the main settings
     *      block to avoid conflicts with another settings.
     *      Both strBlock and strKey cannot contain ' and ".
     *
     *  @param strBlock
     *  @param strKey
     *  @return
     */
    public String getSetting(String strBlock, String strKey)
    {
        if(strKey==null || strKey.isEmpty())
            return null;
        if(strBlock==null || strBlock.isEmpty())
            strBlock = "host.application";
        return values.get(strBlock.replaceAll(" ", ".")+" "+strKey.replaceAll(" ", "."));
    }

    /**
     *  Adds the specified setting to the settings file.
     *  If strBlock is null or empty, the block is assumed to be the main
     *  settings block (see Remarks).
     *  strKey cannot be null or empty, otherwise the setting is not added.
     *  If strValue is null or empty, the setting is not added.
     *
     *  Remarks:
     *      The plugins should not save their settings in the main settings
     *      block to avoid conflicts with another settings.
     *      Both strBlock and strKey cannot contain ' and ".
     * 
     *  @param strBlock
     *  @param strKey
     *  @param strValue
     */
    public void addSetting(String strBlock, String strKey, String strValue)
    {
        if(strKey==null || strKey.isEmpty())
            return;
        if(strValue==null || strValue.isEmpty())
            return;
        if(strBlock==null || strBlock.isEmpty())
            strBlock = "host.application";
        strKey = strBlock.replaceAll(" ", ".")+" "+strKey.replaceAll(" ", ".");
        values.put(strKey, strValue);
    }

    /**
     *  Removes the specified setting.
     *
     *  @param strBlock
     *  @param strKey
     */
    public void removeSetting(String strBlock, String strKey)
    {
        if(strKey==null || strKey.isEmpty())
            return;
        if(strBlock==null || strBlock.isEmpty())
            strBlock = "host.application";
        strKey = strBlock.replaceAll(" ", ".")+" "+strKey.replaceAll(" ", ".");
        values.remove(strKey);
    }

    /**
     *  Only one instance of SettingsManager is allowed to be created.
     *  If there already is an instance of SettingsManager, the method returns
     *  it. Otherwise a new instance is created.
     *
     *  Remarks:
     *      If the filename is null, empty, or specifies a directory the method
     *      returns null.
     *      If the specified file does not exist, it is created.
     * 
     *  @param strFilename
     *  @return
     */
    public static SettingsManager create(String strFilename)
    {
        // Do not allow multiple copies.
        if(prev!=null)
            return prev;
        // Filename should specify a valid file.
        if(strFilename==null || strFilename.isEmpty())
            return null;
        File f = new File(strFilename);
        if(f.exists() && f.isDirectory())
            return null;
        return new SettingsManager(strFilename);
    }
}
