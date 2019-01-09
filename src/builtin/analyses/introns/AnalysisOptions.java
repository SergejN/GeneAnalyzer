/*
    File:
        AnalysisOptions.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Holds the analysis options.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.introns;


public class AnalysisOptions 
{
    public String   strPop    = null;
    public String   strOut    = null;
    public String   strOutput = null;
    public boolean  bNogtag   = false;
    public boolean  bCombine  = false;
    public int[]    sites     = null; 
    public int      maxlen    = Integer.MAX_VALUE;
    public int      minlen    = 0;
    public int      maxstr    = Integer.MAX_VALUE;
    public boolean  bJC_Pi    = false;
    public boolean  bJC_Theta = false;
    public boolean  bJC_K     = false;
    public float    cof       = 1.0f;
    public boolean  bExclAll  = false;
    public boolean  bLenRange = true;

    public boolean  bShowRes  = false;
}
