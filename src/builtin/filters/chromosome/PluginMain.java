/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Filters the genes from the specified chromosome.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.chromosome;

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
    private String strLastErr = "";
    
    
    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Chromosome";
    }

    public String GetName()
    {
        return "Chromosomal location filter";
    }

    public String GetDescription()
    {
        return "Keeps only the genes from the specified chromosome(s)";
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
        return "chr='<CHR1>;<CHR2>'";
    }

    public String GetLastError()
    {
        return strLastErr;
    }

    @Override
    public int[] GetFilteredIndices(Dataset dataset, String strParams)
    {
        String[] chrs = getOptions(strParams);
        if(chrs==null)
        {
            strLastErr = "Cancelled by user";
            return new int[]{-1};
        }
        initData.wd.show(IWaitDialog.TYPE.Filter);
        // Iterate through the dataset.
        Vector<Integer> ind = new Vector<Integer>();        
        for(int i=0;i<dataset.getGenesCount();i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            // Iterate through the data set.
            MAINLOOP: 
            for(int n=0;n<ge.getStrainsCount();n++)
            {
                // If the current gene is on the specified chromosome in 
                // any strain, keep the index and go over to the next gene.
                String strChr = ge.getStrainEntry(n).getChromosome();
                for(int k=0;k<chrs.length;k++)
                {
                    if(strChr.startsWith(chrs[k]))
                    {
                        ind.add(i);
                        break MAINLOOP;
                    }
                }
            }
        }
        int[] indices = new int[ind.size()];
        for(int i=0;i<ind.size();i++)
            indices[i] = ind.get(i);
        initData.wd.close();
        return indices;
    }
    
    /**
     *  Returns the chromosomes to use. If no chromosomes are specified,
     *  null is returned.
     * 
     *  @return
     */
    private String[] getOptions(String strParams)
    {
        // If no parameters are specified.
        if(strParams==null || strParams.isEmpty())
        {
            return (new OptionsForm()).showOptionsForm();
        }
        // Else try to parse the parameter string.
        Pattern p = Pattern.compile("chr='([^']+)'");
        Matcher m = p.matcher(strParams);
        if(m.find())
            return m.group(1).split("\\s*;\\s*");
        else
            return (new OptionsForm()).showOptionsForm();
    }
}
