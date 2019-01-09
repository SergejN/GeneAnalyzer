/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Parses the data files in John Parsch's FASTA format.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.readers.jpf;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetImporter;


public class PluginMain extends ADatasetImporter
{
    private AInitData initData  = null;
    private String strLastError = "";
    
    @Override
    public String GetFileExtension()
    {
        return "fas";
    }

    
    @Override
    public String GetFileDescription()
    {
        return "Files in John Parsch's FASTA format";
    }

    
    @Override
    public Dataset ImportDataset(File[] files, String strParams)
    {
        initData.wd.show(IWaitDialog.TYPE.Import);
        Dataset ds = new Dataset();
        for(File f:files)
        {
            Dataset tmp = createDataset(f);
            ds.merge(tmp);
        }
        initData.wd.close();
        return ds;
    }

    
    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    
    public String GetMenuItemName()
    {
        return "John Parsch's FASTA";
    }
    

    public String GetName()
    {
        return "John Parsch's FASTA files reader";
    }

    
    public String GetDescription()
    {
        return "Reads the data in John Parsch's FASTA format";
    }

    public boolean SupportsMissingData()
    {
        return true;
    }

    public boolean SupportsAmbiguousData()
    {
        return true;
    }

    public String GetParamString()
    {
        return "";
    }
    
    public String GetLastError()
    {
        return strLastError;
    }
    
    
    private Dataset createDataset(File file)
    {
        // Read the file line by line.
        try
        {
            String strChr = (file.getAbsolutePath().endsWith("X.fas")) ? "X" : "Autosome"; 
            BufferedReader in = new BufferedReader(
                                    new InputStreamReader(
                                        new FileInputStream(file)));
            Dataset ds = new Dataset();
            GeneEntry ge = null;
            
            String strHeader = null;
            StringBuffer seq = new StringBuffer();
            
            String strLine;
            while((strLine = in.readLine())!=null)
            {
                // A new header line candidate is found.
                if(strLine.startsWith(">"))
                {
                    // Create a new StrainEntry instance from the data loaded previously.
                    if(strHeader!=null)
                    {
                        GeneEntry tmp = addStrainEntry(strHeader,seq.toString(), ge, strChr);
                        if(ge!=tmp)
                        {
                            ds.addGene(ge);
                            ge = tmp;
                        }
                    }
                    // Set new header line.
                    strHeader = strLine;
                    seq = new StringBuffer();
                }
                // A line with the DNA sequence.
                else
                {
                    seq.append(strLine);
                }
            }  
            // Add last gene entry and last strain entry.
            if(strHeader!=null)
            {
                ds.addGene(addStrainEntry(strHeader, seq.toString(), ge, strChr));
            }
            return ds;
        }
        // Return null on error and set the error description.
        catch(FileNotFoundException e)
        {
            strLastError = String.format("File %s was not found", file.getAbsolutePath());
            return null;
        }
        catch(IOException e)
        {
            strLastError = String.format("An I/O error occured while reading", file.getAbsolutePath());
            return null;
        }
    }
    
    
    /**
     *  Parses the header line and the sequence and creates a new strain entry.
     *  After that adds this entry to the specified gene entry. If the gene entry
     *  is null, it is created first. If the parsed common gene name is different
     *  than the one of the gene entry, a new gene entry is created.
     *  The method returns the gene entry it added the strain to.
     * 
     *  @param strHeaderLine
     *  @param strSequence
     *  @param ge
     *  @param strChr 
     *  @return 
     */
    private GeneEntry addStrainEntry(String strHeaderLine, String strSequence, GeneEntry ge, String strChr)
    {
        int i = 1; // Last seen region end.
        // Create a new strain entry.
        StrainEntry se = new StrainEntry(null, null);
        // Parse the information.
        // If the header is malformatted, a new strain entry with
        // default (empty) field values is returned. The complete sequence is,
        // however, assigned as "Unnamed region", if it is not an empty string.
        Pattern p = Pattern.compile(">\\s*([A-Za-z\\s]+);\\s+" +    // Species (1)
				    "([A-Za-z_\\d]+);\\s+" +        // Strain  (2)
				    "([A-Z\\d]+);\\s+" +	    // Common name (3)
				    "([A-Z\\d]+);\\s+" + 	    // Internal name (4)
				    "Length=(\\d+);\\s+" +	    // Length (5)
				    "Exons=(\\d+);\\s+" +	    // Exons count (6)
                                    "((\\d+\\.\\.\\d+;\\s*)+)");    // Exon positions (7)
        Matcher m = p.matcher(strHeaderLine);
        if(m.find())
        {
            se.setSpeciesName(m.group(1));
            se.setStrainName(m.group(2));
            // If the parsed common name is different from the one returned by
            // ge.GetCommonName, then create a new gene entry.
            if( (ge==null) || (!ge.getCommonName().equalsIgnoreCase(m.group(3))) )
            {
                ge = new GeneEntry(m.group(3), m.group(4));
            }
            // Add populations.
            if(m.group(2).startsWith("MEL"))
                se.addPopulations("Europe", m.group(1));
            else if(m.group(2).startsWith("ZBMEL"))
                se.addPopulations("Africa", m.group(1));
            else
                // In case of non-melanogaster strains simply use species name
                // as population name.
                se.addPopulations(m.group(1));
            se.setChromosome(strChr);            
            // Create gene regions: in JPF format the only possible region types
            // are intron an exon. If a sequence region is not within the
            // specified exon positions, then it is assumed to be an intron.
            Pattern p_int = Pattern.compile("(\\d+)\\.\\.(\\d+)");
            Matcher m_int = p_int.matcher(m.group(7));
            while(m_int.find())
            {
                int start = Integer.parseInt(m_int.group(1));
                int end   = Integer.parseInt(m_int.group(2));
                // Check whether an intron has to be added before.
                if(start>i)
                {
                    GeneRegion intron = new GeneRegion(GeneRegion.INTRON);
                    intron.setStart(i);
                    intron.setEnd(start-1);
                    // Check whether the substring can be extracted.
                    if(strSequence.length()>start)
                        intron.setSequence(strSequence.substring(i-1, start-1));
                    se.addRegion(intron);                    
                }
                // Add exon.
                GeneRegion exon = new GeneRegion(GeneRegion.EXON);
                exon.setStart(start);
                exon.setEnd(end);
                // Check whether the substring can be extracted.
                if( (strSequence.length()>start) && (end-start>0) )
                    exon.setSequence(strSequence.substring(start-1, end));
                se.addRegion(exon);
                i = end+1;
            }            
        }        
        // Check whether the complete sequence was added.
        if(strSequence.length()>=i)
        {
            GeneRegion unnamed = new GeneRegion(GeneRegion.UNNAMED);
            unnamed.setStart(i);
            unnamed.setEnd(strSequence.length());
            unnamed.setSequence(strSequence.substring(i-1, strSequence.length()));
            se.addRegion(unnamed);
        }
        ge.addStrain(se);
        return ge;
    }
}
