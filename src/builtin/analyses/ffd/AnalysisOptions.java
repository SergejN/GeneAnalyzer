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

package builtin.analyses.ffd;

public class AnalysisOptions 
{
    public String strPop    = null;
    public String strOut    = null;
    public String strOutput = null;
    // Exclusion options.
    public boolean nonffd   = true;
    public boolean nonsyn   = true;
    public boolean gaps     = true;
    public int     maxstr   = Integer.MAX_VALUE;
    // Jukes-Cantor correction.
    public boolean bJC_Pi   = false;
    public boolean bJC_Theta= false;
    public boolean bJC_K    = false;
    // Cut-off frequency.
    public float cof        = 1.0f;
    public boolean bExclAll = false;

    public boolean bShowRes = false;
}
