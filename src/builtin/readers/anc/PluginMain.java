/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Parses the files with ancestral sequence alignment.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.readers.anc;

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
    private AInitData initData = null;
    private String strErr = null;
    
    @Override
    public String GetFileExtension()
    {
        return "reconstr";
    }

    @Override
    public String GetFileDescription()
    {
        return "Ancestral sequence alignment files";
    }

    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Ancestral sequence";
    }

    public String GetName()
    {
        return "Ancestral sequence data importer";
    }

    public String GetDescription()
    {
        return "Imports the data in Peter Andolfatto's ancestral sequence data format";
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
        return strErr;
    }

    @Override
    public Dataset ImportDataset(File[] files, String strParams)
    {
        initData.wd.show(IWaitDialog.TYPE.Import);
        Dataset ds = new Dataset();
        for(File f:files)
            loadData(f, ds);
        initData.wd.close();
        return ds;
    }
    
    private void loadData(File file, Dataset ds)
    {
        try
        {
            BufferedReader in = new BufferedReader(
                                    new InputStreamReader(
                                        new FileInputStream(file)));
            
            GeneEntry ge = null;
            // Intron type is contained in the file name.
            String strType = null;
            Pattern p = Pattern.compile("([0-9A-Za-z]+).+(CG\\d+)");
            Matcher m = p.matcher(file.getName());            
            if(m.find())
            {
                // Look for the gene entry.
                for(int i=0;i<ds.getGenesCount();i++)
                {
                    if(ds.getGeneEntry(i).getCommonName().equalsIgnoreCase(m.group(2))) 
                        ge = ds.getGeneEntry(i);
                }
                if(ge==null)
                {
                    ge = new GeneEntry(m.group(2), "");
                    ds.addGene(ge);
                }
                strType = m.group(1);
                if(strType.equalsIgnoreCase("intron"))
                    strType = GeneRegion.INTRON;
                else if(strType.equalsIgnoreCase("introncat"))
                    strType = GeneRegion.INTRON;
                else if(strType.equalsIgnoreCase("exon"))
                    strType = GeneRegion.EXON;
                else if(strType.equalsIgnoreCase("5UTR"))
                    strType = GeneRegion.UTR5;
                else if(strType.equalsIgnoreCase("3UTR"))
                    strType = GeneRegion.UTR3;
                else if(strType.equalsIgnoreCase("intergenic"))
                    strType = GeneRegion.INTERGENIC;
                else
                    strType = GeneRegion.UNNAMED;
            }
            else
                strType = GeneRegion.UNNAMED;
            
            String strHeader = null;
            StringBuffer seq = new StringBuffer();
            
            String strLine;
            while((strLine = in.readLine())!=null)
            {
                if(strLine.startsWith(">"))
                {
                    if(strHeader!=null)
                    {
                        addStrain(strHeader, ge, seq.toString(), strType);
                    }
                    // Set new header line.
                    strHeader = strLine;
                    seq = new StringBuffer();
                }
                // Sequence data.
                else
                    seq.append(strLine);
            }
            addStrain(strHeader, ge, seq.toString(), strType);
        }
        catch(FileNotFoundException e)
        {
            strErr = String.format("File %s was not found", file.getAbsolutePath());
        }
        catch(IOException e)
        {
            strErr = String.format("An I/O error occured while reading", file.getAbsolutePath());
        }
    }
    
    private void addStrain(String strHeader, GeneEntry ge, String strSeq, String strType)
    {
        Pattern p;
        Matcher m;
        // Extract the strain name.
        if(strHeader.startsWith(">Droso")) 
            p = Pattern.compile(">[^_]+_([^_]+)");
        else
            p = Pattern.compile(">([^_]+)");
        m = p.matcher(strHeader);
        String strStrain = "";
        String strSpec = "";
        if(m.find())
        {
            strStrain = m.group(1);
            if(strStrain.equalsIgnoreCase("ANC"))
                strSpec = "Ancestral";
            else
                strSpec = "Drosophila melanogaster";
        }
        StrainEntry se = null;
        for(int i=0;i<ge.getStrainsCount();i++)
        {
            if(ge.getStrainEntry(i).getStrainName().equalsIgnoreCase(strStrain))
                se = ge.getStrainEntry(i);
        }
        if(se==null)
        {
            se = new StrainEntry(strSpec, strStrain);
            if(strSpec.equalsIgnoreCase("Ancestral"))
                se.addPopulations("Ancestral");
            else
                se.addPopulations("Africa", "Drosophila melanogaster");
            ge.addStrain(se);
        }
        // Create region.
        if(strType.equalsIgnoreCase(GeneRegion.INTRON))
            strSeq = "GT"+strSeq+"AG";
        int iLength = strSeq.length();
        GeneRegion r = new GeneRegion(strType);
        int l = se.getCompleteSequence().length();
        r.setStart(l+1);
        r.setEnd(l+iLength);
        r.setSequence(strSeq);
        se.addRegion(r);
    }
}
