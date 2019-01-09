/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Filters the genes which overlap.

    Note:
        The class requires the property DatasetBuilder.GENOMIC_POSITION to be
        annotated to the gene entry. This property should be an 1x2 array, whereas
        the first value is the start, and the second one the last position in the genome.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.overlap;

import dpgp.DatasetBuilder;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import gui.IWaitDialog;
import java.util.Vector;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetFilter;

public class PluginMain extends ADatasetFilter
{
    private AInitData initData = null;
    private String strErr = "";

    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Overlapping genes";
    }

    public String GetName()
    {
        return "Overlapping genes filter";
    }

    public String GetDescription()
    {
        return "Keeps only the genes which do not overlap with other genes";
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
    public int[] GetFilteredIndices(Dataset dataset, String strParams)
    {
        initData.wd.show(IWaitDialog.TYPE.Filter);
        Vector<Integer> ind = new Vector<Integer>();
        int nGenes = dataset.getGenesCount();
        MAINLOOP: for(int i=0;i<nGenes;i++)
        {
            GeneEntry g1 = dataset.getGeneEntry(i);
            int[] tmp = (int[])g1.getProperty(DatasetBuilder.GENOMIC_POSITION);
            if(tmp==null)
            {
                ind.add(i);
                continue;
            }
            int iStart1 = tmp[0];
            int iEnd1 = tmp[1];
            if(iStart1>iEnd1)
            {
                int k = iEnd1;
                iEnd1 = iStart1;
                iStart1 = k;
            }
            for(int j=0;j<dataset.getGenesCount();j++)
            {
                GeneEntry g2 = dataset.getGeneEntry(j);
                // If the second gene entry is just a splice variant of the same gene as the first gene entry,
                // skip the second gene entry.
                if(isSameGene(g1.getCommonName(), g2.getCommonName()))
                    continue;
                tmp = (int[])g2.getProperty(DatasetBuilder.GENOMIC_POSITION);
                if(tmp==null)
                    continue;
                int iStart2 = tmp[0];
                int iEnd2 = tmp[1];
                if(iStart2>iEnd2)
                {
                    int k = iEnd2;
                    iEnd2 = iStart2;
                    iStart2 = k;
                }
                if( (iStart1<iStart2 && iStart2<iEnd1) || (iStart1<iEnd2 && iEnd2<iEnd1))
                    continue MAINLOOP;
            }
            ind.add(i);
        }
        int[] indices = new int[ind.size()];
        for(int i=0;i<ind.size();i++)
            indices[i] = ind.get(i);
        initData.wd.close();
        return indices;
    }

    private boolean isSameGene(String strCG1, String strCG2)
    {
        return strCG1.substring(0, strCG1.length()-1).equalsIgnoreCase(strCG2.substring(0, strCG2.length()-1));
    }
}
