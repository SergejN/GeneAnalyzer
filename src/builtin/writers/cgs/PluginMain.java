/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Exports the complete gene sequence in the GeneAnalyzer native FASTA
        format but with all gaps removed.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.cgs;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetExporter;


public class PluginMain extends ADatasetExporter
{
    private AInitData initData = null;
    private String strLastErr = "";

    @Override
    public String GetFileExtension()
    {
        return "cgs";
    }

    @Override
    public String GetFileDescription()
    {
        return "Ungapped FASTA";
    }

    @Override
    public ErrorCode ExportDataset(Dataset ds, File file, String strParams)
    {
        initData.wd.show(IWaitDialog.TYPE.Export);
        StringBuffer content = new StringBuffer();
        for(int i=0;i<ds.getGenesCount();i++)
        {
            GeneEntry ge = ds.getGeneEntry(i);
            String strCommonName = ge.getCommonName();
            String strInternalName = ge.getAlias();
            for (int j = 0; j < ge.getStrainsCount(); j++)
            {
                int offset = 0;
                StringBuffer seq = new StringBuffer();
                StrainEntry se = ge.getStrainEntry(j);
                StringBuffer regs = new StringBuffer();
                for (int n=0; n<se.getRegionsCount(); n++)
                {
                    int oo = offset;        // Old offset.
                    GeneRegion r = se.getRegion(n);
                    String strSeq = r.getSequence().replaceAll("-", "");
                    // Update offset.
                    offset += r.getSequence().length()-strSeq.length();
                    regs.append(String.format("Reg:%s=%d-%d; ",
                                                 r.getType(),
                                                 r.getStart()-oo,
                                                 r.getEnd()-offset));
                    seq.append(strSeq);
                }
                String strHeader = String.format(">%s; %s; %s; %s; %s; Population=%s; %s",
                                                  se.getSpeciesName(),
                                                  se.getStrainName(),
                                                  strCommonName,
                                                  strInternalName,
                                                  se.getChromosome(),
                                                  getPopulationsString(se),
                                                  regs.toString());
                String strSequence = formatSequence(seq.toString());
                content.append(String.format("%s\n%s",
                               strHeader, strSequence));

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

    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Ungapped FASTA";
    }

    public String GetName()
    {
        return "Ungapped FASTA exporter";
    }

    public String GetDescription()
    {
        return "Exports the data into ungapped FASTA format";
    }

    // Only general genetic code is supported: A, C, G, T, -
    public boolean SupportsMissingData()
    {
        return false;
    }

    public boolean SupportsAmbiguousData()
    {
        return false;
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


    /**
     *  Returns the formatted FASTA-like sequence, i.e. with a line break after
     *  every 60bp.
     *
     *  @param strSequence
     *  @return
     */
    private String formatSequence(String strSequence)
    {
        StringBuffer seq = new StringBuffer();
        Pattern p = Pattern.compile("[ACGTacgt-]{1,60}");
        Matcher m = p.matcher(strSequence);
        while (m.find())
        {
            seq.append(m.group());
            seq.append("\n");
        }
        return seq.toString().toLowerCase();
    }
}
