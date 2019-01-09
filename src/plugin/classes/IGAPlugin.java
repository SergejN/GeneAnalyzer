/*
    File:
        IGAPlugin.java
 *   
    Revision:
        2.2.0.1
 * 
    Description:
        Basic plugin interface, from which all other plugin interfaces
        are inherited. The name of the main class of the plugin, which
        implements this interface should be PluginMain.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package plugin.classes;

import kernel.ErrorCode;
import plugin.AInitData;
import plugin.PluginType;

public interface IGAPlugin 
{
    public static final String PLUGIN_EXTENTION = ".jar";
        
    /**
     *  Initializes the plugin.
     * 
     *  @param initdata     initialization data
     *  @return
     */
    public ErrorCode Initialize(AInitData initdata);
        
    /**
     *  Returns the type of the plugin.
     * 
     *  @return
     */
    public PluginType GetType();
    
    /**
     *  Returns the name to display in the menu.
     * 
     *  @return
     */
    public String GetMenuItemName();
        
    /**
     *  Returns the name of the plugin.
     * 
     *  @return
     */
    public String GetName();
        
    /**
     *  Returns the plugin description.
     * 
     *  @return
     */
    public String GetDescription();

    /**
     *  If the plugin can handle missing data, i.e. sequences with X or N, this method
     *  should return true. If the current setting is to use the genetic code which
     *  allows missing data, the host-application disables the menu items which refer to
     *  the plugins which cannot handle such data.
     *
     *  Note:
     *      Plugins supporting the missing data should also support complete data.
     *
     *  Remarks:
     *      This method should always return true in the classes inherited from <code>ADatasetImporter</code>.
     *
     *  @return
     */
    public boolean SupportsMissingData();

    /**
     *  If the plugin can handle missing data and ambiguous data, i.e. sequences with X, N or any
     *  ambiguity base, this method should return true. If the current setting allows the use of
     *  ambiguity bases, the host-application disables the menu items which refer to
     *  the plugins which cannot handle such data.
     * 
     *  Note:
     *      Plugins supporting the ambiguous data should also support missing data.
     *
     *  Remarks:
     *      This method should always return true in the classes inherited from <code>ADatasetImporter</code>.
     *
     *  @return
     */
    public boolean SupportsAmbiguousData();

    /**
     *  Returns the backbone of the parameter string as follows:
     *      param1=<PARAM1> param2=<PARAM2>
     *  where paramX is the parameter name (paramX can be replaced by any other name)
     *  and <PARAMX> is a place holder of its value.
     *  If the plugin does not need any extra parameters, the string should be
     *  empty, but not null.
     * 
     *  @return
     */
    public String GetParamString();
    
    
    /**
     *  Returns the string description of the last error occured.
     * 
     *  @return
     */
    public String GetLastError();
}
