/*
    File:
        PluginMain.java
 *   
    Revision:
        1.2.0.1
 * 
    Description:
        Saves the data into GeneAnalyzer native FASTA format.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.nfa;

import algorithms.SequenceRoutines;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetExporter;


public class PluginMain extends ADatasetExporter
{
    private AInitData initData = null;
    private String strLastErr  = "";
   
    @Override
    public String GetFileExtension()
    {
        return "nfa";
    }

    @Override
    public String GetFileDescription()
    {
        return "GeneAnalyzer native FASTA";
    }

    @Override
    public ErrorCode ExportDataset(Dataset ds, File file, String strParams)
    {
        initData.wd.show(IWaitDialog.TYPE.Export);
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(file));
            int nGenes = ds.getGenesCount();
            for(int i=0;i<nGenes;i++)
            {
                GeneEntry ge = ds.getGeneEntry(i);
                String strCommonName = ge.getCommonName();
                String strInternalName = ge.getAlias();
                int nStrains = ge.getStrainsCount();
                for (int j=0;j<nStrains;j++)
                {
                    StrainEntry se = ge.getStrainEntry(j);
                    StringBuffer regs = new StringBuffer();
                    for (int n=0; n<se.getRegionsCount(); n++)
                    {
                        GeneRegion r = se.getRegion(n);
                        regs.append(String.format("Reg:%s=%d-%d; ",
                                                     r.getType(),
                                                     r.getStart(),
                                                     r.getEnd()));
                    }
                    String strHeader = String.format(">%s; %s; %s; %s; %s; Population=%s; %s",
                                                      se.getSpeciesName(),
                                                      se.getStrainName(),
                                                      strCommonName,
                                                      strInternalName,
                                                      se.getChromosome(),
                                                      getPopulationsString(se),
                                                      regs.toString());
                    String strSequence = SequenceRoutines.getFormattedSequence(se.getCompleteSequence(), 60);
                    out.write(String.format("%s\n%s", strHeader, strSequence));
                }
            }
            out.close();
        }
        catch (IOException e)
        {
            strLastErr = "An I/O error occured while saving the file";
            file.delete();
            initData.wd.close();
            return ErrorCode.IOError;
        }
        initData.wd.close();
        return ErrorCode.Ok;
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
        return "GeneAnalyzer native FASTA exporter";
    }

    public String GetDescription()
    {
        return "Exports the data into GeneAnalyzer native FASTA format";
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
        return strLastErr;
    }
    
    
    /**
     *  Returns the populations of the specified strain entry concatinated
     *  into one comma-separated string.
     * 
     *  @param se
     *  @return
     */
    private String getPopulationsString(StrainEntry se)
    {
        String[] pops = se.listPopulations();
        StringBuffer popbuf = new StringBuffer();
        for(int i=0;i<pops.length;i++)
        {
            popbuf.append(pops[i]);
            if(i<pops.length-1)
                popbuf.append(", ");
        }
        return popbuf.toString();
    }
}
