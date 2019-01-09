/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
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

package builtin.writers.nfa2;

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
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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
        return "GeneAnalyzer native FASTA v.2";
    }

    @Override
    public ErrorCode ExportDataset(Dataset ds, File file, String strParams)
    {
        ExportOptions eo = getOptions(strParams, ds.listPopulations());
        if(eo==null)
        {
            strLastErr = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Export);
        try
        {            
            PrintWriter out = new PrintWriter(new FileWriter(file));
            int nGenes = ds.getGenesCount();
            for(int i=0;i<nGenes;i++)
            {
                GeneEntry ge = ds.getGeneEntry(i);
                String strCommonName = ge.getCommonName();
                String strInternalName = (eo.bNoIntNames) ? "" : ge.getAlias();
                StrainEntry[] strains = getStrainEntries(ge, eo);
                if(strains==null)
                    continue;
                int nStrains = strains.length;
                for (int j=0;j<nStrains;j++)
                {
                    StrainEntry se = strains[j];
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
        return "GeneAnalyzer native FASTA exporter v.2";
    }

    public String GetDescription()
    {
        return "Exports the data into GeneAnalyzer native FASTA format using additional options";
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
        return "pops='<POP1>;<POP2>' strlimit='<LIMIT>' nointname='<T/F>'";
    }

    public String GetLastError()
    {
        return strLastErr;
    }

    /**
     *  Returns an array of strains to export.
     *
     *  @param ge
     *  @param eo
     *  @return
     */
    private StrainEntry[] getStrainEntries(GeneEntry ge, ExportOptions eo)
    {
        Vector<StrainEntry> strains = new Vector<StrainEntry>();
        int nStrains = ge.getStrainsCount();
        for(String strPop:eo.pops)
        {
            int n = 0;
            for(int i=0;i<nStrains;i++)
            {
                StrainEntry se = ge.getStrainEntry(i);
                if(se.belongsToPopulation(strPop) && strains.indexOf(se)==-1 && n<eo.nLimit)
                {
                    strains.add(se);
                    n++;
                }
            }
        }
        if(strains.size()==0)
            return null;
        else
            return strains.toArray(new StrainEntry[strains.size()]);
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

    private ExportOptions getOptions(String strParams, String[] pops)
    {
        if(strParams==null || strParams.isEmpty())
            return (new OptionsForm()).getOptions(pops);
        else
        {
            // Try to parse the options.
            Pattern p = Pattern.compile("pops='([^']+)' strlimit='(\\d*)' nointname='(<T/F>)'");
            Matcher m = p.matcher(strParams);
            if(m.find())
            {
                ExportOptions eo = new ExportOptions();
                eo.pops = m.group(1).split("\\s*;\\s*");
                if(!m.group(2).isEmpty())
                    eo.nLimit = Integer.parseInt(m.group(2));
                else
                    eo.nLimit = Integer.MAX_VALUE;
                eo.bNoIntNames = m.group(3).equalsIgnoreCase("T");
                return eo;
            }
            else
                return (new OptionsForm()).getOptions(pops);
        }
    }
}
