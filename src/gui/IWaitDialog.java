/*
    File:
        IWaitDialog.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        The interface IWaitDialog must be implemented by any wait
        dialog intended to be used with GeneAnalyzer.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui;

public interface IWaitDialog 
{
    public static enum TYPE
    {
        Kernel,
        Import,
        Export,
        Filter,
        Analysis,
        Aligner,
        Custom
    };
    
    /**
     *  Shows the wait dialog.
     * 
     *  Note:
     *  The host-application decides whether or not to display the
     *  wait dialog.
     */
    public void show(TYPE type);
    
    /**
     *  Sets the text to display. If the dialog type is other than TYPE.Custom
     *  the strMainText parameter is ignored. Hint text is displayed under the
     *  main text using different font and color. The strHintText parameter
     *  can be null, if not hint text should be displayed.
     * 
     *  Remarks:
     *  strHintText parameter can be used to display the name of the current
     *  gene etc.
     * 
     *  @param strMainText
     *  @param strHintText
     */
    public void setText(String strMainText, String strHintText);

    /**
     *  Closes the wait dialog.
     */
    public void close();
}
