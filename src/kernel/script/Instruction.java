/*
    File:
        Instruction.java 
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Represents a single instruction, supported by GeneAnalyzer.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package kernel.script;

import java.io.File;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Instruction 
{
    private InstructionType type;
    
    /**
     *  Name of the plugin. 
     *  If the instruction does not require a plugin, this field is null.
     */
    private String strObjName = null;
    
    /**
     *  Holds the files the plugin should load/save.
     *  This field is only relevant for instructions of type Load and Save.
     */
    private File[] files = null;
    
    /**
     *  Holds additional instruction parameters. If the instruction does not require
     *  additional parameters, this field is null.
     */
    private String strParams  = null;
    
    /**
     *  System command. This field is only field, if the instruction type is
     *  System.
     */
    private String strCmd     = null;
    
    
    /**
     *  Declare a private constructor, to avoid creating instructions from
     *  outside.
     */
    private Instruction()
    {
    }
    
    /**
     *  Returns the type of the instruction.
     * 
     *  @return
     */
    public InstructionType getType()
    {
        return type;
    }
    
    /**
     *  Returns the plugin name.
     * 
     *  @return
     */
    public String getObjectName()
    {
        return strObjName;
    }
    
    /**
     *  Returns the parameter string.
     * 
     *  @return
     */
    public String getParameters()
    {
        return strParams;
    }   
    
    /**
     *  Returns the files array.
     * 
     *  @return
     */
    public File[] getFiles()
    {
        return files;
    }
    
    /**
     *  Returns the command of the System instruction.
     * 
     *  @return
     */
    public String getCommand()
    {
        return strCmd;
    }
    
    /**
     *  Creates and returns a new Instruction object. If the instruction is not 
     *  recognized, e.g. due to a spelling error like laod instead of load etc, 
     *  the method returns null.
     * 
     *  All currently supported instructions are listed below:
     *      - Load data set.
     *          LOAD name="<NAME>" source="<FILENAME>" params="<PARAMS>";
     *              <NAME>      importer name
     *              <FILENAME>  source: file to load. To load multiple files, repeat the
     *                          source="<FILENAME>" block one time for each file.
     *              <PARAMS>    parameters to pass to the importer
     *      
     *      - Create a dataset from DPGP data (*.VMA files and *.GFF file).
     *          DPGPLOAD params="<PARAMS>";
     *              <PARAMS>    parser parameters
     * 
     *      - Unload data set.
     *          UNLOAD;
     *
     *      - Analyze data set.
     *          ANALYZE name="<NAME>" params="<PARAMS>";
     *              <NAME>      analyzer name
     *              <PARAMS>    parameters to pass to the analyzer
     * 
     *      - Filter the genes.
     *          FILTER name="<NAME>" params="<PARAMS>";
     *              <NAME>      filter name
     *              <PARAMS>    parameters to pass to the filter
     * 
     *      - Save the data set.
     *          SAVE name="<NAME>" dest="<FILENAME>" params="<PARAMS>";
     *              <NAME>      exporter name
     *              <FILENAME>  target file/directory
     *              <PARAMS>    parameters to pass to the exporter
     *      
     *      - Change the codon table.
     *          SYSTEM CODONTABLE "<NAME>";
     *              <NAME> name of the codon table:
     *                  Following names are predefined:
     *                      - Initial: use the codon table GeneAnalyzer was started with.
     *                                 You can change it the program settings.
     *                      - Default: use the default codon table
     *
     *      - Select all genes.
     *          SYSTEM SELECTALL;
     * 
     *      - Sort the dataset.
     *          SYSTEM SORT "<TYPE>";
     *              <TYPE> sort type:
     *                  - By name       sort by name
     *                  - By quality    sort by quality level
     * 
     *      - Invert current selection
     *          SYSTEM INVERTSELECTION;
     * 
     *      - Run external program
     *          SYSTEM EXECUTE <NOWAIT> "<PROGRAM>"
     *              <NOWAIT> if the modifyer NOWAIT is set, the host-application
     *                       waits for the started program to finish otherwise
     *                       it proceeds without waiting.
     *              <PROGRAM> is an arbitrary string to execute
     * 
     *  @param strInstruction
     *  @return
     */
    public static Instruction createInstruction(String strInstruction)
    {
        // Lower case instruction.
        String lci = strInstruction.toLowerCase();
        int fp = 0;
        while(!Character.isLetter(lci.charAt(fp)))
            fp++;
        lci = lci.substring(fp);
        strInstruction = strInstruction.substring(fp);
        // LOAD.
        if(lci.startsWith("load "))
        {
            Matcher m = patLoad.matcher(strInstruction);
            if(m.find())
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.Load;
                inst.strObjName = m.group(1);
                inst.strParams = m.group(3);
                inst.files = parseFilenames(m.group(2));
                inst.strCmd = null;
                return inst;
            }
        }
        // DPGPLOAD.
        else if(lci.startsWith("dpgpload "))
        {
            Matcher m = patDpgpLoad.matcher(strInstruction);
            if(m.find())
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.DPGPLoad;
                inst.strObjName = null;
                inst.strParams = m.group(1);
                inst.files = null;
                inst.strCmd = null;
                return inst;
            }
        }
        // UNLOAD.
        else if(lci.equalsIgnoreCase("unload;"))
        {
            Instruction inst = new Instruction();
            inst.type = InstructionType.Unload;
            inst.strObjName = null;
            inst.strParams = null;
            inst.files = null;
            inst.strCmd = null;
            return inst;
        }
        // ANALYZE.
        else if(lci.startsWith("analyze "))
        {
            Matcher m = patAnalyze.matcher(strInstruction);
            if(m.find())
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.Analyze;
                inst.strObjName = m.group(1);
                inst.strParams = m.group(2);
                inst.files = null;
                inst.strCmd = null;
                return inst;
            }
        }
        // FILTER.
        else if(lci.startsWith("filter "))
        {
            Matcher m = patFilter.matcher(strInstruction);
            if(m.find())
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.Filter;
                inst.strObjName = m.group(1);
                inst.strParams = m.group(2);
                inst.files = null;
                inst.strCmd = null;
                return inst;
            }
        }
        // SAVE.
        else if(lci.startsWith("save "))
        {
            Matcher m = patSave.matcher(strInstruction);
            if(m.find())
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.Save;
                inst.strObjName = m.group(1);
                inst.strParams = m.group(3);
                inst.files = new File[]{new File(m.group(2))};
                inst.strCmd = null;
                return inst;
            }
        }
        // SYSTEM.
        else if(lci.startsWith("system "))
        {
            if(strInstruction.equalsIgnoreCase("system selectall;"))
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.System;
                inst.strObjName = null;
                inst.strParams = null;
                inst.files = null;
                inst.strCmd = "Select all";
                return inst;
            }
            else if(strInstruction.equalsIgnoreCase("system invertselection;"))
            {
                Instruction inst = new Instruction();
                inst.type = InstructionType.System;
                inst.strObjName = null;
                inst.strParams = null;
                inst.files = null;
                inst.strCmd = "Invert selection";
                return inst;
            }
            else if(lci.startsWith("system sort "))
            {
                Pattern p = Pattern.compile("^SYSTEM\\s+SORT\\s+\"([^\"]+)\"",
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
                Matcher m = p.matcher(strInstruction);
                if(m.find())
                {
                    if(m.group(1).equalsIgnoreCase("by name") || 
                       m.group(1).equalsIgnoreCase("by quality"))
                    {
                        Instruction inst = new Instruction();
                        inst.type = InstructionType.System;
                        inst.strObjName = null;
                        inst.strParams = m.group(1);
                        inst.files = null;
                        inst.strCmd = "Sort";
                        return inst;
                    }
                }
            }
            else if(lci.startsWith("system execute "))
            {
                Pattern p = Pattern.compile("^SYSTEM\\s+EXECUTE\\s+(NOWAIT){0,1}\\s*\"([^\"]+)\"",
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
                Matcher m = p.matcher(strInstruction);
                if(m.find())
                {
                    Instruction inst = new Instruction();
                    inst.type = InstructionType.System;                    
                    inst.strObjName = m.group(2);
                    inst.strParams = m.group(1);
                    inst.files = null;
                    inst.strCmd = "Execute";
                    return inst;
                }
            }
            else
            {
                Pattern p = Pattern.compile("^SYSTEM\\s+CODONTABLE\\s+\"([^\"]+)\";", 
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
                Matcher m = p.matcher(strInstruction);
                if(m.find())
                {
                    Instruction inst = new Instruction();
                    inst.type = InstructionType.System;
                    inst.strObjName = null;
                    inst.strParams = m.group(1);
                    inst.files = null;
                    inst.strCmd = "Codon table";
                    return inst;
                }
            }
        }
        return null;
    }
       
    /**
     *  Returns the pattern of the instruction of the specified type.
     * 
     *  @param type
     *  @return
     */
    private static Pattern compileInstructionPatterns(InstructionType type)
    {
        switch(type)
        {
            case Load:
                return Pattern.compile("^load\\s+name\\s*=\\s*\"([^\"]+)\"\\s+(source\\s*=.+)\\s+params\\s*=\\s*\"([^\"]*)\";", 
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);           
            case Save:
                return Pattern.compile("^save\\s+name\\s*=\\s*\"([^\"]+)\"\\s+dest\\s*=\\s*\"([^\"]+)\"\\s+params\\s*=\\s*\"([^\"]*)\";", 
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
            case Filter:
                return Pattern.compile("^filter\\s+name\\s*=\\s*\"([^\"]+)\"\\s+params\\s*=\\s*\"([^\"]*)\";",
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
            case Analyze:
                return Pattern.compile("^analyze\\s+name\\s*=\\s*\"([^\"]+)\"\\s+params\\s*=\\s*\"([^\"]*)\";",
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
            case DPGPLoad:
                return Pattern.compile("^dpgpload\\s+params\\s*=\\s*\"([^\"]+)\";", 
                                            Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
            // Instructions, without the pattern.
            case Unload:
            case System:
                break;
        }
        return null;
    }
    
    private static File[] parseFilenames(String strFiles)
    {
        Pattern p = Pattern.compile("^source\\s*=\\s*\"([^\"]+)\"", Pattern.CASE_INSENSITIVE|Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strFiles);
        Vector<File> tmp = new Vector<File>();
        while(m.find())
            tmp.add(new File(m.group(1)));
        if(tmp.size()==0)
            return null;
        return tmp.toArray(new File[tmp.size()]);
    }
    
    // Patterns of all supported instructions.
    private static Pattern patLoad          = compileInstructionPatterns(InstructionType.Load);
    private static Pattern patDpgpLoad      = compileInstructionPatterns(InstructionType.DPGPLoad);
    private static Pattern patAnalyze       = compileInstructionPatterns(InstructionType.Analyze);
    private static Pattern patFilter        = compileInstructionPatterns(InstructionType.Filter);
    private static Pattern patSave          = compileInstructionPatterns(InstructionType.Save);
}
