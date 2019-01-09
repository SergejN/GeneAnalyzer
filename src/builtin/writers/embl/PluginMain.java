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

package builtin.writers.embl;

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
    private String strLastErr = "";

    @Override
    public String GetFileExtension()
    {
        return "embl";
    }

    @Override
    public String GetFileDescription()
    {
        return "Mega bulk-submission FASTA";
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
            int nStrains = ge.getStrainsCount();
            for (int j=0;j<nStrains;j++)
            {
                StrainEntry se = ge.getStrainEntry(j);
                if(se.getSpeciesName().contains("simulans") || se.getSpeciesName().contains("sechellia") || se.getSpeciesName().contains("yakuba"))
                    continue;
                int offset = 0;
                StringBuffer seq = new StringBuffer();                
                StringBuffer exons  = new StringBuffer();
                int iSeqStart = 0;
                int iCdsStart = -1;
                int nRegs = se.getRegionsCount();
                for (int n=0; n<nRegs; n++)
                {
                    int oo = offset;        // Old offset.
                    GeneRegion r = se.getRegion(n);
                    String strSeq = r.getSequence().replaceAll("-", "");
                    // Update offset.
                    offset += r.getSequence().length()-strSeq.length();
                    seq.append(strSeq);
                    if(r.hasType(GeneRegion.EXON))
                    {
                        if(iSeqStart>0)
                        {
                            exons.append(String.format("%d; %d; ", iSeqStart, r.getEnd()-offset));
                            iSeqStart = 0;
                        }
                        else
                            exons.append(String.format("%d; %d; ", r.getStart()-oo, r.getEnd()-offset));
                        if(iCdsStart==-1)
                            iCdsStart = r.getStart()-oo;
                    }
                    else if(r.hasType(GeneRegion.UTR5))
                    {
                        // If the 5'UTR follows the intergenic region, update the start position of the
                        // first exon.
                        if(n==0 || se.getRegion(n-1).hasType(GeneRegion.INTERGENIC))
                        {
                            iSeqStart = r.getStart()-oo;
                            // If, however, the 5'UTR is followed by an intron, create an additional
                            // exon from the 5'UTR and don't change the exon start.
                            if(n!=nRegs-1 && se.getRegion(n+1).hasType(GeneRegion.INTRON))
                            {
                                iSeqStart = 0;
                                exons.append(String.format("%d; %d; ", r.getStart()-oo, r.getEnd()-offset));
                            }
                        }
                    }
                }
                String strHeader = String.format(">%s; %s; %d; %s",
                                                 se.getStrainName(),
                                                 strCommonName,
                                                 iCdsStart,
                                                 exons.toString());
                String strSequence = SequenceRoutines.getFormattedSequence(seq.toString(), 100);
                content.append(String.format("%s\n%s", strHeader, strSequence));

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
        return "EMBL submission FASTA";
    }

    public String GetName()
    {
        return "EMBL FASTA exporter";
    }

    public String GetDescription()
    {
        return "Exports the data into ungapped FASTA format ready for submission at EMBL";
    }

    // The EMBL exporter only supports general genetic code: A, C, G, T, -.
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
}
