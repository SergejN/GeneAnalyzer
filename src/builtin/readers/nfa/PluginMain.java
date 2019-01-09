/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Parses the data files in GeneAnalyzer native FASTA format.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.readers.nfa;

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
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetImporter;


public class PluginMain extends ADatasetImporter
{
    private AInitData initData = null;
    private String strLastError = "";
    
    @Override
    public String GetFileExtension()
    {
        return "nfa";
    }

    
    @Override
    public String GetFileDescription()
    {
        return "Files in GeneAnalyzer native FASTA format";
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
        return "Native FASTA";
    }
    

    public String GetName()
    {
        return "GeneAnalyzer native FASTA files reader";
    }

    
    public String GetDescription()
    {
        return "Reads the data in GeneAnalyzer native FASTA format";
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
                        GeneEntry tmp = addStrainEntry(strHeader,seq.toString(), ge);
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
                ds.addGene(addStrainEntry(strHeader, seq.toString(), ge));
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
     *  @return 
     */
    private GeneEntry addStrainEntry(String strHeaderLine, String strSequence, GeneEntry ge)
    {
        // Create a new strain entry.
        StrainEntry se = new StrainEntry(null, null);
        // Parse the information.
        // If the header is malformatted, a new strain entry with
        // default (empty) field values is returned. The complete sequence is,
        // however, assigned as "Unnamed region", if it is not an empty string.
        Pattern p = Pattern.compile(">\\s*([^;]+);\\s*" +           // Species (1)
                                    "([^;]+);\\s*" +                // Strain  (2)
                                    "([^;]+);\\s*" +                // CG number (3)
                                    "([^;]*);\\s*" +                // Internal name (4)
                                    "([^;]*);\\s*" +                // Chromosomal location (5)
                                    "Population[s]*=([^;]+);\\s*" + // Population(s) (6)
                                    "((Reg:\\s*[^=]+=\\d+-\\d+;*\\s*)+)");   // Regions (7)
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
            se.setChromosome(m.group(5));
            // Add populations.
            se.addPopulations(m.group(6).split(",\\s*"));
            // Create gene regions.
            Pattern p_int = Pattern.compile("Reg:([^=]+)=(\\d+)-(\\d+)");
            Matcher m_int = p_int.matcher(m.group(7));
            while(m_int.find())
            {
                int start = Integer.parseInt(m_int.group(2));
                int end   = Integer.parseInt(m_int.group(3));
                // Check whether an intron has to be added before.
                GeneRegion region = new GeneRegion(m_int.group(1));
                region.setStart(start);
                region.setEnd(end);
                if( (strSequence.length()>=start) && (end-start>=0) && (strSequence.length()>=end))
                    region.setSequence(strSequence.substring(start-1, end));
                se.addRegion(region);
            }
            // Iterate through the regions to make sure they cover the complete
            // sequence.
            int i = -1; // Last seen end.
            Vector<GeneRegion> tmp = new Vector<GeneRegion>();
            // Sequence gaps before and between the defined regions.
            for(int n=0;n<se.getRegionsCount();n++)
            {
                GeneRegion r = se.getRegion(n);
                if(r.getStart()-i>0 && i>0)
                {
                    GeneRegion un = new GeneRegion(GeneRegion.UNNAMED);
                    un.setStart(i);
                    un.setEnd(r.getStart()-1);
                    if(strSequence.length()>r.getStart())
                        un.setSequence(strSequence.substring(un.getStart()-1, un.getEnd()));
                    tmp.add(un);                    
                }
                i = r.getEnd()+1;
            }
            // Sequence gap at the end of the sequence.
            if(strSequence.length()>=i)
            {
                GeneRegion un = new GeneRegion(GeneRegion.UNNAMED);
                un.setStart(i);
                un.setEnd(strSequence.length());
                un.setSequence(strSequence.substring(i-1, strSequence.length()));
                tmp.add(un);
            }
            // Add all gaps.
            for(GeneRegion r:tmp)
                se.addRegion(r);            
        }
        else
        {
            if(ge==null)
                ge = new GeneEntry("Unknown", "");
            GeneRegion un = new GeneRegion(GeneRegion.UNNAMED);
            un.setStart(1);
            if(strSequence!=null)
                un.setEnd(strSequence.length());
            un.setSequence(strSequence);
            se.addRegion(un);
        }
        ge.addStrain(se);
        return ge;
    }
}
