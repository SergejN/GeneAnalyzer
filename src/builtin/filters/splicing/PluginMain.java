/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Keeps only one splice variant of every gene.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.splicing;

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
        return "Splice variants";
    }

    public String GetName()
    {
        return "Splice variants filter";
    }

    public String GetDescription()
    {
        return "Keeps only one splice variant per gene";
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
            GeneEntry ge = dataset.getGeneEntry(i);
            String strName = ge.getCommonName();
            char c = strName.charAt(strName.length()-1);
            // If the last character is a letter, check the splice variants.
            if(Character.isLetter(c))
            {
                // If the splice variant is not variant A, check, whether other
                // splice variants exist.
                if(c!='A')
                {
                    char[] tmp = strName.toCharArray();
                    for(char b='A';b<c;b++)
                    {
                        tmp[strName.length()-1] = b;
                        if(dataset.hasGene(new String(tmp)))
                            continue MAINLOOP;
                    }
                    
                }
            }
            ind.add(i);
        }
        int[] indices = new int[ind.size()];
        for(int i=0;i<ind.size();i++)
            indices[i] = ind.get(i);
        initData.wd.close();
        return indices;
    }
}