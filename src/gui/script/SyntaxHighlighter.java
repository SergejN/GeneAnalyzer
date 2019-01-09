/*
    File:
        SyntaxHighlighter.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Highlights the GeneAnalyzer Script syntax.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.script;

import java.awt.Color;
import javax.swing.JTextPane;
import javax.swing.text.DefaultHighlighter;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;
import kernel.script.Instruction;


public class SyntaxHighlighter implements Runnable
{ 
    private Style stInstruction     = null;
    private Style stParameterName   = null;
    private Style stParameterValue  = null;
    private Style stComment         = null;
    private Style stSystem          = null;
    private Style stError           = null;
    
    private JTextPane textPane      = null;    
    private boolean bIsRunning      = true;

    
    public SyntaxHighlighter(JTextPane textPane)
    {
        this.textPane = textPane;
        // Instruction style.
        stInstruction = textPane.addStyle("Instruction", null);
        StyleConstants.setForeground(stInstruction, Color.BLUE);
        StyleConstants.setBold(stInstruction, true);
        // Parameter name style.
        stParameterName = textPane.addStyle("ParamName", null);
        StyleConstants.setForeground(stParameterName, Color.RED);
        // Parameter value style.
        stParameterValue = textPane.addStyle("ParamValue", null);
        StyleConstants.setForeground(stParameterValue, new Color(0, 128, 0));
        // Comment style.
        stComment = textPane.addStyle("Comment", null);
        StyleConstants.setForeground(stComment, new Color(125, 100, 70));
        StyleConstants.setItalic(stComment, true);
        // System style.
        stSystem = textPane.addStyle("System", null);
        StyleConstants.setForeground(stSystem, new Color(0, 0, 90));
        StyleConstants.setBold(stSystem, true);
        // Syntax error style.
        stError = textPane.addStyle("SyntaxError", null);
        StyleConstants.setForeground(stError, Color.GRAY);
        StyleConstants.setItalic(stError, true);
    }
    
    public void run()
    {
        while(bIsRunning)
        {
            String strCode = textPane.getText();
            StyledDocument doc = textPane.getStyledDocument();
            int l = strCode.length();
            int pi = -2;         // Position of the instruction start
            int pc = -2;         // Position of the comment start
            int pv = -2;         // Position of the parameter value start
            int is = -2;
            boolean ignoreSemicolon = false;
            for(int i=0;i<l;i++)
            {
                char c = strCode.charAt(i);
                if(is==-2 && pc==-2)
                    is = i;
                /***************************************************************
                *                       COMMENT
                ***************************************************************/
                // If the char is the newline character after the comment symbol
                // then apply the comment style and reset the position.
                if(c=='\n' && pc>-1)
                {
                    doc.setCharacterAttributes(pc, i-pc, stComment, true);
                    pc = -2;
                }
                // If the character is a slash, check the previous character, and if it
                // is also a slash, update the comment start position.
                else if(c=='/')
                {
                    if(i>0 && strCode.charAt(i-1)=='/')
                    {
                        pc = i-1;
                        is = -2;
                    }
                }
                /***************************************************************
                *                       VALUE
                ***************************************************************/
                else if(c=='"')
                {
                    if(pv>-1)
                    {
                        doc.setCharacterAttributes(pv, i-pv+1, stParameterValue, true);
                        pv = -2;
                    }
                    else
                        pv = i;
                    ignoreSemicolon = !ignoreSemicolon;
                }
                /***************************************************************
                *              PARAMETER NAMES AND INSTRUCTIONS
                ***************************************************************/
                else if(Character.isLetter(c) && pi==-2 && pc==-2 && pv==-2)
                    pi = i;
                else if( (c==' ' || c=='=' || c==';') && pi>-1)
                {
                    String word = strCode.substring(pi, i);
                    if(isInstruction(word))
                        doc.setCharacterAttributes(pi, i-pi, stInstruction, true);
                    else if(isParameterName(word))
                        doc.setCharacterAttributes(pi, i-pi, stParameterName, true);
                    else if(isSystemCall(word))
                        doc.setCharacterAttributes(pi, i-pi, stSystem, true);
                    pi = -2;
                }
                /***************************************************************
                *                     SYNTAX ERRORS
                ***************************************************************/
                if(c==';' && !ignoreSemicolon && is>-1)
                {
                    String strInst = strCode.substring(is, i+1);
                    strInst = strInst.replaceAll("\"\n", "\" ");
                    strInst = strInst.replaceAll("\n", "");
                    if(Instruction.createInstruction(strInst)==null)
                        doc.setCharacterAttributes(is, i-is, stError, true);
                    is = -2;
                }
            }
            try
            {
                Thread.sleep(2000);
            } 
            catch (InterruptedException ex){}
        }
    }

    private boolean isInstruction(String strWord)
    {
        return strWord.equalsIgnoreCase("LOAD") ||
               strWord.equalsIgnoreCase("UNLOAD") ||
               strWord.equalsIgnoreCase("DPGPLOAD") ||
               strWord.equalsIgnoreCase("SAVE") ||
               strWord.equalsIgnoreCase("ANALYZE") ||
               strWord.equalsIgnoreCase("FILTER");
    }
    
    private boolean isSystemCall(String strWord)
    {
        return strWord.equalsIgnoreCase("SYSTEM");
    }
    
    private boolean isParameterName(String strWord)
    {
        return strWord.equalsIgnoreCase("name") ||
               strWord.equalsIgnoreCase("params") ||
               strWord.equalsIgnoreCase("source") ||
               strWord.equalsIgnoreCase("dest");
    }
    
    /***************************************************************************
    *                         PRIVATE CLASSES.                                 *
    ***************************************************************************/
    class KeywordPainter extends DefaultHighlighter.DefaultHighlightPainter 
    {
        public KeywordPainter(Color color) 
        {
            super(color);
        }
    };
}
