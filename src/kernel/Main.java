/*
    File:
        Main.java
 *   
    Revision:
        2.2.0.1
 * 
    Description:
        Application entry point.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel;


import gui.MainForm;
import gui.splash.SplashScreen;
import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.UIManager;


public class Main 
{
    // Host version.
    public final static int HOST_MAJOR_VERSION = 2;
    public final static int HOST_MINOR_VERSION = 4;
    
    public final static String APPTITLE        = "GeneAnalyzer 2.1";

    
    public static void main(String[] args) throws Exception 
    {     
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        JFrame.setDefaultLookAndFeelDecorated(true);
        final SplashScreen splash = new SplashScreen();
        new Thread(new Runnable() {

            public void run()
            {
                splash.setVisible(true);
            }
        }).start();        
        Kernel kernel = new Kernel();
        kernel.initialize();   
        final MainForm mf = new MainForm(kernel);        
        EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    mf.setVisible(true);
                    splash.setVisible(false);
                }
            });
    }    
}
