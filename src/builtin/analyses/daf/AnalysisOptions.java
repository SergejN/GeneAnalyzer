/*
    File:
        AnalysisOptions.java
 *   
    Revision:
        1.0.0.1
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

package builtin.analyses.daf;

public class AnalysisOptions 
{
    public String   strPop    = null;
    public String   strOut    = null;
    public String   strOutput = null;
    public String   strRegion = null;
    public int      iMaxlen   = Integer.MAX_VALUE;
    public int      iMinlen   = 0;
    public int      maxstr    = Integer.MAX_VALUE;
    public boolean  bCombine  = false;
    public boolean  bUseAny   = true;
    public boolean  bNonFfd   = false;
    public int[]    sites     = null;
    public boolean  bSizeConst= false;
    public boolean  bLenRange = true;
    public boolean  bListFreqs= true;

    public boolean bShowRes   = false;
}
