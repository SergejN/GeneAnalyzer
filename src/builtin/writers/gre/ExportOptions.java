/*
    File:
        ExportOptions.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Holds the export options of Gene Region Exporter. 
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.gre;


public class ExportOptions 
{
    String strType = null;
    String[] pops    = null;
    int maxlen       = Integer.MAX_VALUE;
    int minlen       = 0;
}
