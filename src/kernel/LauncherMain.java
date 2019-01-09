/*
    File:
        LauncherMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Launches the GeneAnalyzer.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel;


public class LauncherMain 
{
    public static void main(String[] args) throws Exception 
    {
        SettingsManager sm = SettingsManager.create("Settings.xml");
        String s = sm.getSetting("", SettingsManager.HEAPSIZE);
        if(s!=null && s.matches("[0-9]+"))
        {
            int nHeapSize = Integer.parseInt(s);
            if(nHeapSize>64)
            {
                Runtime.getRuntime().exec(String.format("java -Xmx%dm -cp GeneAnalyzer.jar kernel.Main", nHeapSize));
            }
            else
            {
                Runtime.getRuntime().exec("java -cp GeneAnalyzer.jar kernel.Main");
            }
        }
        else
        {
            Runtime.getRuntime().exec("java -cp GeneAnalyzer.jar kernel.Main");
        }
    }
}
