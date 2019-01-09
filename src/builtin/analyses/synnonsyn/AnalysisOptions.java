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
package builtin.analyses.synnonsyn;

public class AnalysisOptions 
{
    public String strPop    = null;
    public String strOut    = null;
    public String strOutput = null;
    public int     maxstr   = Integer.MAX_VALUE;
    // Jukes-Cantor correction.
    public boolean bJC_pi   = false;
    public boolean bJC_t    = false;
    public boolean bJC_K    = false;
    // Cut-off frequency.
    public float cof        = 1.0f;
    // Terminal codons.
    public boolean bExclTer = true;
    public boolean bUseTer  = false;
    public boolean bExclAll = false;

    public boolean bShowRes = false;
}
