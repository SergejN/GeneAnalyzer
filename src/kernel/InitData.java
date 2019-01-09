/*
    File:
        AInitData.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Plugin initialization data.
        
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package kernel;

import plugin.AInitData;


public class InitData extends AInitData
{
    private Kernel kernel = null;
    
    
    public InitData(Kernel kernel)
    {
        this.kernel = kernel;
    }

    public ErrorCode setCodonTable(String strName)
    {
        return kernel.setCodonTable(strName);
    }
}
