/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Saves the gene region sequences into a FASTA file.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.gre;

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
        return "grf";
    }

    @Override
    public String GetFileDescription()
    {
        return "Gene region sequence";
    }

    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Gene region sequence";
    }

    public String GetName()
    {
        return "Gene region sequence exporter";
    }

    public String GetDescription()
    {
        return "Exports the gene region sequence into FASTA format";
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
        return "reg='<REG>' pops='<POP1>;<POP2>' minlen='<MIN>' maxlen=<MAX>";
    }
    
    public String GetLastError()
    {
        return strLastErr;
    }
    
    @Override
    public ErrorCode ExportDataset(Dataset ds, File file, String strParams)
    {
        StringBuffer content = new StringBuffer();
        
        OptionsForm of = new OptionsForm();
        ExportOptions opt = of.showOptionsForm(ds.listRegionsNames(), ds.listPopulations());
        // If the user clicked Cancel, cancel the export.
        if(opt==null)
        {
            strLastErr = "Cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Export);
        // Iterate through the data set.
        for(int i=0;i<ds.getGenesCount();i++)
        {
            GeneEntry ge = ds.getGeneEntry(i);
            StrainEntry se = ge.getStrainEntry(0);
            // If the gene entry does not have the specified region, skip it.
            if(se.getRegionsCount(opt.strType)==0)
                continue;
            int counter = 1;      // Indicates the region number.
            // Iterate through the regions.
            for(int n=0;n<se.getRegionsCount();n++)
            {                
                boolean flag = false; // Indicates whether the region was exported.
                GeneRegion reg = se.getRegion(n);                
                // If the region has the desired type and length, add the the 
                // sequences of different strains iteratively.
                if(reg.hasType(opt.strType))
                {
                    int l = reg.getSequence().replaceAll("-", "").length();
                    if(l>opt.maxlen || l<opt.minlen)
                        continue;
                    l = reg.getSequence().length(); // Alignment length.
                    for(int j=0;j<ge.getStrainsCount();j++)
                    {
                        // Check whether the strain belongs to the population.
                        if(belongsToPopulation(ge.getStrainEntry(j), opt.pops))
                        {
                            flag = true;
                            content.append(String.format(">%s; %s; %s; %s %d; Length=%dbp;\n%s\n", 
                                     ge.getStrainEntry(j).getSpeciesName(),
                                     ge.getStrainEntry(j).getStrainName(),
                                     ge.getCommonName(),
                                     opt.strType, counter, l, 
                                     ge.getStrainEntry(j).getRegion(n).getSequence()));
                        }
                    }
                }
                if(flag)
                {
                    counter++;
                    content.append("\n");
                }
            }
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(file));
            out.print(content.toString());
            out.close();
        }
        catch (IOException e)
        {
            strLastErr = "An I/O error occured while saving the file";
            initData.wd.close();
            return ErrorCode.IOError;
        }
        initData.wd.close();
        return ErrorCode.Ok;
    }    
    
    
    /**
     *  Tests whether or not the strain entry belongs to at least one population
     *  from the given list.
     *  
     *  @param se
     *  @param pops
     *  @return
     */
    private boolean belongsToPopulation(StrainEntry se, String[] pops)
    {
        for(String s:pops)
            if(se.belongsToPopulation(s))
                return true;
        return false;
    }
}
