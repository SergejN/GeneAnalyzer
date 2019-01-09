/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Filters out the specified genes.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.name;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import gui.IWaitDialog;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetFilter;


public class PluginMain extends ADatasetFilter
{
    private AInitData initData = null;
    private String strErr = null;
    
    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Gene name";
    }

    public String GetName()
    {
        return "Gene name filter";
    }

    public String GetDescription()
    {
        return "Filters out the genes with the specified names";
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
        return "names='<NAME1>;<NAME2>' action='<E/I>'";
    }

    public String GetLastError()
    {
        return strErr;
    }

    @Override
    public int[] GetFilteredIndices(Dataset dataset, String strParams)
    {
        FilterOptions opt = getOptions(strParams);
        if(opt==null)
        {
            strErr = "Cancelled by user";
            return new int[]{-1};
        }
        initData.wd.show(IWaitDialog.TYPE.Filter);
        // Create regular expressions.
        for(int i=0;i<opt.names.length;i++)
        {
            opt.names[i] = opt.names[i].replace('?', '.');
            opt.names[i] = opt.names[i].replaceAll("\\*", ".*");
        }
        // Iterate through the data set and exclude the specified genes.
        Vector<Integer> ind = new Vector<Integer>();
        MAINLOOP: for(int i=0;i<dataset.getGenesCount();i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            // Check the name. If the name is in the list, treat the index
            // depending on the action to perform.
            for(String s:opt.names)
            {
                if(ge.getCommonName().matches(s))
                {
                    // Exclude.
                    if(opt.bExclude)
                        continue MAINLOOP;
                    else // Include.
                    {
                        ind.add(i);
                        continue MAINLOOP;
                    }
                }
            }
            // If the gene name is unmatched and the specified genes should be
            // INCLUDED, continue MAINLOOP.
            if(!opt.bExclude)
                continue MAINLOOP;
            else // Otherwise add the index to the list.
                ind.add(i);
        }
        int[] indices = new int[ind.size()];
        for(int i=0;i<ind.size();i++)
            indices[i] = ind.get(i);
        initData.wd.close();
        return indices;
    }
    
    private FilterOptions getOptions(String strParams)
    {
        if(strParams==null || strParams.isEmpty())
        {
            return (new OptionsDialog()).getOptions();
        }
        else
        {
            Pattern p = Pattern.compile("names='([^']+)' action='([EeIi])'");
            Matcher m = p.matcher(strParams);
            if(m.find())
            {
                FilterOptions opt = new FilterOptions();
                opt.names = m.group(1).split("\\s*;\\s*");
                opt.bExclude = m.group(2).equalsIgnoreCase("E");
                return opt;
            }
            else
            {
                return (new OptionsDialog()).getOptions();
            }
        }
    }
}
