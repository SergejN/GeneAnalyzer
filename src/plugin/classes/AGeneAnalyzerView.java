/*
    File:
        AGeneAnalyzerView.java
 *   
    Revision:
        2.2.0.1
 * 
    Description:
        Any class which is intended to be displayed on the views panel of
        GeneAnalyzer must be inheritted from the base class AGeneAnalyzerView.
        The name of the main class of the plugin, which implements this interface
        should be ViewMain.
 
    Note:
        The class extending this abstract class must be fast. Otherwise
        the waiting time is unacceptable. Do not create a view if
        the needed computations are too time consuming.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin.classes;

import bio.gene.GeneEntry;
import javax.swing.JPanel;

public abstract class AGeneAnalyzerView extends JPanel
{
    public static final String PLUGIN_EXTENTION = ".jar";

    /**
     *  Returns the name of the view as it should appear in the tool bar.
     * 
     *  @return
     */
    public abstract String getButtonName();
    
    /**
     *  Returns the hint text which should be displayed when the user
     *  points at the button in the tool bar.
     * 
     *  @return
     */
    public abstract String getButtonHintText();
    
    /**
     *  Updates the view to display the information about the specified
     *  gene entry. The host-application passes null, if the view should
     *  be reset, e.g. after the dataset was unloaded.
     * 
     *  @param ge
     */
    public abstract void displayGeneEntryInfo(GeneEntry ge);
}
