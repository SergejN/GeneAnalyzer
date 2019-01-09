/*
    File:
        ScriptInterpreter.java 
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Interprets the GeneAnalyzer scripts.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel.script;

import kernel.*;
import gui.MainForm;
import java.io.IOException;


public class ScriptInterpreter 
{
    private Kernel kernel = null;
    private MainForm mf   = null;
    
    private String strLastErr = null;
    
    /**
     *  Constructs the intepreter.
     * 
     *  @param kernel   kernel to work with
     */
    public ScriptInterpreter(Kernel kernel, MainForm mf)
    {
        this.kernel = kernel;
        this.mf     = mf;
    }
    
    public String getLastErrorString()
    {
        return strLastErr;
    }
    
    /**
     *  Runs the script.
     *  Returns:
     *      Ok                  if succeeds
     *      ExecutionError      if the instruction is malformed
     * 
     *  @param strCode
     *  @return
     */
    public ErrorCode runScript(String strCode)
    {
        // Remove comments.
        strCode = strCode.replaceAll("//[^\n]+\n", "");
        // Remove newlines.
        strCode = strCode.replaceAll("\"\n", "\" ");
        strCode = strCode.replaceAll("\n", "");
        int l = strCode.length();
        StringBuffer sbi = new StringBuffer();
        boolean ignoreSemicolon = false;        // If ignoreSemicolon is true, a semicolon is
                                                // not interpreted as instruction end.
        for(int i=0;i<l;i++)
        {
            char c = strCode.charAt(i);
            sbi.append(c);
            // If a semicolon found, check whether it is within an instruction.
            if(c==';' && !ignoreSemicolon)
            {
                // Remove any spaces at the beginning of the instruction.
                while(sbi.charAt(0)==' ')
                    sbi.deleteCharAt(0);
                Instruction instr = Instruction.createInstruction(sbi.toString());
                if(instr!=null)
                {
                    ErrorCode ec = runInstruction(instr);
                    if(ec==ErrorCode.CancelledByUser)
                    {
                        strLastErr = "The script execution was cancelled by user";
                        return ec;
                    }
                    if(ec==ErrorCode.NoDatasetLoaded)
                    {
                        strLastErr = "No dataset loaded";
                        return ec;
                    }
                    if(ec==ErrorCode.IOError || ec==ErrorCode.ExecutionError)
                    {
                        return ec;
                    }
                    if(ec!=ErrorCode.Ok)
                    {
                        strLastErr = kernel.getLastErrorString();
                        return ec;
                    }
                    mf.updateGeneList();
                }
                else
                {
                    strLastErr = "Unsupported instruction";
                    return ErrorCode.ObjectNotFound;
                }
                // Clear the buffer.
                sbi.setLength(0);
            }
            if(c=='"')
                ignoreSemicolon = !ignoreSemicolon;
        }
        return ErrorCode.Ok;
    }
    
    private ErrorCode runInstruction(Instruction instruction)
    {
        switch(instruction.getType())
        {
            case Load:
                return kernel.importDataset(instruction.getFiles(), 
                                            instruction.getObjectName(), 
                                            instruction.getParameters());
            case DPGPLoad:
                return kernel.importDPGPDataset(instruction.getParameters());
            case Unload:
                kernel.unloadDataset();
                return ErrorCode.Ok;
            case Analyze:
            {
                return kernel.performAnalysis(null,
                                              instruction.getObjectName(), 
                                              instruction.getParameters());
            }
            case Filter:
            {
                return kernel.applyFilter(instruction.getObjectName(), 
                                          instruction.getParameters());
            }
            case Save:
            {
                return kernel.exportDataset(null,
                                            instruction.getObjectName(), 
                                            instruction.getFiles()[0], 
                                            instruction.getParameters());
            }
            case System:
                // Select all instruction
                if(instruction.getCommand().equalsIgnoreCase("Select all"))
                {
                    kernel.selectAll();
                    mf.updateGeneList();
                    return ErrorCode.Ok;
                }
                else if(instruction.getCommand().equalsIgnoreCase("Invert selection"))
                {
                    kernel.invertSelection();
                    mf.updateGeneList();
                    return ErrorCode.Ok;
                }
                else if(instruction.getCommand().equalsIgnoreCase("Sort"))
                {
                    if(instruction.getParameters().equalsIgnoreCase("by name"))
                    {
                        kernel.sort();
                        return ErrorCode.Ok;
                    }
                    else if(instruction.getParameters().equalsIgnoreCase("by quality"))
                    {
                        kernel.sortByQualityLevel();
                        return ErrorCode.Ok;
                    }
                }
                else if(instruction.getCommand().equalsIgnoreCase("Execute"))
                {
                    Runtime rt = Runtime.getRuntime();
                    try
                    {
                        Process proc = rt.exec(instruction.getObjectName());
                        // If the host application should wait.
                        if(instruction.getParameters()==null ||
                           !instruction.getParameters().equalsIgnoreCase("NOWAIT"))
                            proc.waitFor();
                        return ErrorCode.Ok;
                    }
                    catch(IOException e)
                    {
                        strLastErr = "Could not execute the program";
                        return ErrorCode.IOError;
                    }
                    catch(InterruptedException e)
                    {
                        strLastErr = "Execution error";
                        return ErrorCode.ExecutionError;
                    }
                }
                else if(instruction.getCommand().equalsIgnoreCase("Codon table"))
                    return kernel.setCodonTable(instruction.getParameters());
        }
        return ErrorCode.InvalidParameter;
    }
}
