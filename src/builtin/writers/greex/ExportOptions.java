/*
    File:
        ExportOptions.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Holds the export options of extended region exporter. 
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.greex;

public class ExportOptions 
{
    String[] pops    = null;
    String strReg    = null;    
    int nStrains     = Integer.MAX_VALUE;
    String strSites  = null;
    int iMinLen      = 0;
    int iMaxLen      = Integer.MAX_VALUE;
}
