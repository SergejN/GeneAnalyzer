/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Filters the genes with specified quality level.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.quality;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import gui.IWaitDialog;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import kernel.QualityChecker;
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
        return "Quality level";
    }

    public String GetName()
    {
        return "Quality level filter";
    }

    public String GetDescription()
    {
        return "Keeps only the genes which have quality level equal or above specified";
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
        return "level='<LEVEL>'";
    }

    public String GetLastError()
    {
        return strErr;
    }

    @Override
    public int[] GetFilteredIndices(Dataset dataset, String strParams)
    {
        int iLevel = getLevel(strParams);
        if(iLevel==-1)
        {
            strErr = "Cancelled by user";
            return new int[]{-1};
        }
        initData.wd.show(IWaitDialog.TYPE.Filter);
        Vector<Integer> ind = new Vector<Integer>();
        for(int i=0;i<dataset.getGenesCount();i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            Integer level = (Integer)ge.getProperty(QualityChecker.QUALITY_LEVEL);
            if(level==null)
                continue;
            if(level<=iLevel)
                ind.add(i);
        }
        int[] indices = new int[ind.size()];
        for(int i=0;i<ind.size();i++)
            indices[i] = ind.get(i);
        initData.wd.close();
        return indices;
    }

    private int getLevel(String strParams)
    {
        if(strParams==null || strParams.isEmpty())
            return (new OptionsDialog()).getLevel();
        else
        {
            Pattern p = Pattern.compile("level='(\\d+)'");
            Matcher m = p.matcher(strParams);
            if(m.find())
                return Integer.parseInt(m.group(1));
            else
                return (new OptionsDialog()).getLevel();
        }
    }
}
