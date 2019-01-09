/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Saves the data into Peter Andolfatto's FASTA format.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.paf;

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
    private AInitData initData      = null;
    private final String EXTENSION  = ".fa";
    private String strLastErr       = "";
   
    @Override
    public String GetFileExtension()
    {
        return null;
    }

    @Override
    public String GetFileDescription()
    {
        return "Peter Andolfatto's FASTA";
    }

    @Override
    public ErrorCode ExportDataset(Dataset ds, File file, String strParams)
    {
        initData.wd.show(IWaitDialog.TYPE.Export);
        for(int i=0;i<ds.getGenesCount();i++)
        {
            StringBuffer content = new StringBuffer();
            GeneEntry ge = ds.getGeneEntry(i);
            String strCommonName = ge.getCommonName();
            // Delete this lines afterwards!!!
            String yak = null;
            String sim = null;
            String sec = null;
            for (int j = 0; j < ge.getStrainsCount(); j++)
            {
                StrainEntry se = ge.getStrainEntry(j);
                if(se.belongsToPopulation("Africa"))
                {
                    // Header line.
                    content.append(String.format(">%s_%s_%s\n", 
                                                 se.getSpeciesName().replaceAll(" ", "."),
                                                 se.getStrainName(),
                                                 strCommonName));
                    content.append(se.getCompleteSequence()+"\n");
                } 
                if(se.belongsToPopulation("Drosophila simulans"))
                {
                    sim = String.format(">%s_%s_%s\n", 
                                                 se.getSpeciesName().replaceAll(" ", "."),
                                                 se.getStrainName(),
                                                 strCommonName);
                    sim += se.getCompleteSequence()+"\n";
                }
                if(se.belongsToPopulation("Drosophila yakuba"))
                {
                    yak = String.format(">%s_%s_%s\n", 
                                                 se.getSpeciesName().replaceAll(" ", "."),
                                                 se.getStrainName(),
                                                 strCommonName);
                    yak += se.getCompleteSequence()+"\n";
                }
                if(se.belongsToPopulation("Drosophila sechellia"))
                {
                    sec = String.format(">%s_%s_%s\n", 
                                                 se.getSpeciesName().replaceAll(" ", "."),
                                                 se.getStrainName(),
                                                 strCommonName);
                    sec += se.getCompleteSequence()+"\n";
                }
            }
            if(sim!=null)
                content.insert(0, sim);
            if(sec!=null)
                content.insert(0, sec);
            content.insert(0, yak);
            // Write regions block.
            content.append("\n;\n");
            for(int n=0;n<ge.getStrainEntry(0).getRegionsCount();n++)
            {
                GeneRegion gr = ge.getStrainEntry(0).getRegion(n);
                String type = gr.getType();
                if(type.equalsIgnoreCase(GeneRegion.EXON))
                    type = "CDS";
                else if(type.equalsIgnoreCase("5'UTR")) 
                    type = "five_UTR";
                content.append(String.format("%s\t%d-%d\n", type, gr.getStart(), gr.getEnd()));
            }
            // Save.
            try
            {
                // Construct filename.
                File output = new File(file.getAbsolutePath()+File.separator+
                                       strCommonName+EXTENSION);
                PrintWriter out = new PrintWriter(new FileWriter(output));
                out.print(content.toString());
                out.close();
            }
            catch (IOException e)
            {
                strLastErr = "An I/O error occured while saving the file";
                new File(file.getAbsolutePath()+File.separator+strCommonName+EXTENSION).delete();
                initData.wd.close();
                return ErrorCode.IOError;
            }
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
        return "Peter Andolfatto's FASTA";
    }

    public String GetName()
    {
        return "Peter Andolfatto's FASTA exporter";
    }

    public String GetDescription()
    {
        return "Exports the data into Peter Andolfatto's FASTA format";
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
}
