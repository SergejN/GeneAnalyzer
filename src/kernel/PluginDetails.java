/*
    File:
        PluginDetails.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Holds the details about the plugin.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package kernel;


public class PluginDetails 
{
    public String strMenuItemName = "";
    public String strName         = "";
    public String strDescription  = "";
    public String strParamStr     = "";
    public boolean bSmd           = false;  // Missing data support.
    public boolean bSad           = false;  // Ambiguous data support.
    // This member is always null for ADatasetAnalyzer and ADatasetFilter.
    public ExtensionDetails ed    = null;
}
