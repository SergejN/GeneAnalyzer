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

package builtin.analyses.indel;


public class AnalysisOptions 
{
    public String   strPop    = null;
    public String   strOut    = null;
    public String   strDist   = null;
    public String   strOutput = null;
    public String   strType   = null;
    public boolean  bCombine  = false;
    public int      maxlen    = Integer.MAX_VALUE;
    public int      minlen    = 0;
    public int      maxstr    = Integer.MAX_VALUE;
    public boolean  bLenRange = true;

    public boolean  bShowRes  = false;
}
