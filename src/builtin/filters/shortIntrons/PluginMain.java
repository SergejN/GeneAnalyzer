/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Filters the short introns.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.shortIntrons;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
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
        return "Short introns";
    }

    public String GetName()
    {
        return "Short introns filter";
    }

    public String GetDescription()
    {
        return "Keeps only the genes with introns of specified length";
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
        return "max='<MAX>' min='<MIN>' pops='<POP1>;<POP2>'";
    }
    
    public String GetLastError()
    {
        return strLastErr;
    }

    @Override
    public int[] GetFilteredIndices(Dataset dataset, String strParams)
    {
        FilterOptions fo = getOptions(strParams, dataset.listPopulations());
        if(fo==null)
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
            // Since the sequences went through quality control check,
            // the algorithm assumes, that the gene regions are the same
            // in all strain entries, that all regions have a valid sequence
            // and that the gene entry contains at least one strain entry.
            StrainEntry se = ge.getStrainEntry(0);
            for(int n=0;n<se.getRegionsCount();n++)
            {
                // The algorithm adds genes to the list, which have AT LEAST ONE
                // intron satisfying the criteria. Once such an intron is found
                // no further introns are checked.
                GeneRegion gr = se.getRegion(n);
                // Check the region type.
                if(gr.hasType(GeneRegion.INTRON))
                {
                    int lmax = -1;                  // Maximal intron length.
                    int lmin = Integer.MAX_VALUE;   // Minimal intron length.
                    // If random sampling should be used, use se.
                    int l = 0;
                    for(int j=0;j<ge.getStrainsCount();j++)
                    {
                        StrainEntry tmp = ge.getStrainEntry(j);
                        // Check the population.
                        if(!belongsToPopulation(tmp, fo.pops))
                            continue;
                        GeneRegion r = ge.getStrainEntry(j).getRegion(n);
                        l = r.getSequence().replaceAll("-", "").length();
                        if(l>=lmax) lmax = l;
                        if(l<=lmin) lmin = l;
                    }                   
                    // If the intron has proper length, add the index of the
                    // gene entry to the list.
                    if( (lmin>=fo.imin) && (lmax<=fo.imax) && (l>0))
                    {
                        ind.add(i);
                        break;
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
     *  Returns the options for the filtering.
     * 
     *  @param strParams
     *  @param pops
     *  @return
     */
    private FilterOptions getOptions(String strParams, String[] pops)
    {
        if(strParams==null || strParams.isEmpty())
            return (new OptionsForm()).showOptionsForm(pops);
        else
        {
            // Try to parse the options.
            Pattern p = Pattern.compile("max='(\\d*)' min='(\\d*)' pops='([^']+)'");
            Matcher m = p.matcher(strParams);
            if(m.find())
            {
                FilterOptions fo = new FilterOptions();
                fo.imax = (m.group(1).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(1));
                fo.imin = (m.group(2).isEmpty()) ? 0 : Integer.parseInt(m.group(2));
                fo.pops = m.group(3).split("\\s*;\\s*");
                return fo;
            }
            else
                return (new OptionsForm()).showOptionsForm(pops);
        }
    }
    
    /**
     *  Returns true if the strain entry belongs to at least one population
     *  from the list.
     * 
     *  @param se
     *  @param pops
     *  @return
     */
    private boolean belongsToPopulation(StrainEntry se, String[] pops)
    {
        for(String s:pops)
        {
            if(se.belongsToPopulation(s))
                return true;
        }
        return false;
    }
}
