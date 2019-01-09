/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Filters the genes which have no valid sequence.

    Note:
        When importing DPGP data using GFF files, it can happen, that the specified
        gene lies outside the available alignment. In this case the complete gene sequence
        consists of X's. This filter finds and filters such genes.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.empty;

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
        return "Empty genes";
    }

    public String GetName()
    {
        return "Empty genes filter";
    }

    public String GetDescription()
    {
        return "Filters out the genes with an empty gene sequence";
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
        for(int i=0;i<nGenes;i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            int nEmpty = 0;
            for(int j=0;j<nStrains;j++)
            {
                if(ge.getStrainEntry(j).getCompleteSequence().matches("[XNxn-]+"))
                    nEmpty++;
            }
            // If the number of empty strain sequences if the same as the total
            // number of sequences, skip the gene.
            if(nEmpty==nStrains)
                continue;
            ind.add(i);
        }
        int[] indices = new int[ind.size()];
        for(int i=0;i<ind.size();i++)
            indices[i] = ind.get(i);
        initData.wd.close();
        return indices;
    }
}
