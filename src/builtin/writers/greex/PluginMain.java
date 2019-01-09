/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Extended regions exporter.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.greex;

import algorithms.SequenceRoutines;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetExporter;

public class PluginMain extends ADatasetExporter
{
    private AInitData initData = null;
    private String strErr = null;
    
    @Override
    public String GetFileExtension()
    {
        return null;
    }

    @Override
    public String GetFileDescription()
    {
        return "Gene region";
    }

    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Extended region exporter";
    }

    public String GetName()
    {
        return "Extended Region Exporter";
    }

    public String GetDescription()
    {
        return "Exports the aligned region sequences into separate files";
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
        return "pops='<POP1;POP2>' region='<REG>' strlimit='<LIMIT>' " +
               "range='<POS1>-<POS2>;<POS3>'";
    }

    public String GetLastError()
    {
        return strErr;
    }

    @Override
    public ErrorCode ExportDataset(Dataset ds, File file, String strParams)
    {
        ExportOptions eo = getOptions(ds, strParams);
        if(eo==null)
        {
            strErr = "Cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        // Iterate through the data set.
        initData.wd.show(IWaitDialog.TYPE.Export);
        HashMap<String, Integer> counts = new HashMap<String, Integer>();
        for(int i=0;i<ds.getGenesCount();i++)
        {
            resetHashMap(counts, eo.pops);
            GeneEntry ge = ds.getGeneEntry(i);
            Vector<StrainEntry> strains = new Vector<StrainEntry>();
            for(int n=0;n<ge.getStrainsCount();n++)
            {
                if(proceedStrainEntry(ge.getStrainEntry(n), counts, eo.nStrains))
                    strains.add(ge.getStrainEntry(n)); 
            }
            // Save regions.
            if(strains.size()>0)
                saveRegions(strains, file, eo.strReg, ge.getCommonName(), 
                        eo.iMinLen, eo.iMaxLen, eo.strSites);
        }
        initData.wd.close();
        return ErrorCode.Ok;
    }

    private void resetHashMap(HashMap<String, Integer> hm, String[] pops)
    {
        for(String s:pops)
            hm.put(s, new Integer(0));
    }
    
    private boolean proceedStrainEntry(StrainEntry se, 
                                       HashMap<String, Integer> hm,
                                       int nMax)
    {
        int n = 0;
        String[] pops = se.listPopulations();
        for(String s:pops)
        {
            Integer i = hm.get(s);
            if(i==null)         // The population is not used - skip it.
                continue;
            else
            {
                if(i<nMax)
                {
                    i++;
                    hm.put(s, i);
                    n++;
                }
            }
        }
        // If n>0 then the strain belongs to at least one population and
        // should be exported.
        return n>0;
    }
    
    private void saveRegions(Vector<StrainEntry> strains, 
                             File file, 
                             String strReg,
                             String strGeneName,
                             int iMin,
                             int iMax,
                             String strSites)
    {
        StringBuffer content;
        // Check region type.
        if(strReg.equalsIgnoreCase("CDS")) // export the complete CDS.
        {
            content = new StringBuffer();
            for(StrainEntry se:strains)
            {
                content.append(String.format("%s_%s_%s\n",
                        se.getSpeciesName(), se.getStrainName(), strGeneName));
                content.append(se.getCodingSequence()+"\n");
            }
            String strFilename = file.getAbsolutePath()+File.separator+
                    strGeneName+"_"+strReg+".grs";
            saveToFile(content.toString(), new File(strFilename));
        }
        else // Other region type.
        {
            int n = 0;
            for(int i=0;i<strains.get(0).getRegionsCount();i++)
            {
                if(strains.get(0).getRegion(i).hasType(strReg))
                {
                    n++;
                    // Check sequence length.
                    String[] seqs = new String[strains.size()];
                    for(int k=0;k<seqs.length;k++)
                        seqs[k] = strains.get(k).getRegion(i).getSequence();
                    int[] len = SequenceRoutines.getMinMaxUngappedLength(seqs);
                    if(len[0]<iMin || len[1]>iMax)
                        continue;
                    // Extract the sites.
                    if(strSites!=null && !strSites.isEmpty())
                    {
                        int[] sites = extractSites(strSites);
                        seqs = SequenceRoutines.getUngappedSequences(seqs);
                        seqs = SequenceRoutines.extractAlignedSites(seqs, sites);
                    }
                    content = new StringBuffer();
                    for(int k=0;k<strains.size();k++)
                    {
                        content.append(String.format("%s_%s_%s_%d\n",
                                strains.get(k).getSpeciesName(), strains.get(k).getStrainName(), strGeneName, n));
                        content.append(seqs[k]+"\n");
                    }
                    // Save to file.
                    String strFilename = file.getAbsolutePath()+File.separator+
                        strGeneName+"_"+strReg+Integer.toString(n)+".grs";
                    saveToFile(content.toString(), new File(strFilename));
                }
            }
        }
    }

    private void saveToFile(String strContent, File file)
    {
        System.out.println(file.getName());
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(file));
            out.print(strContent);
            out.close();
        }
        catch (IOException e)
        {
            file.delete();
        }
    }
    
    /**
     *  Returns the export options or null if the user cancels the analysis.
     * 
     *  @param ds
     *  @param strParams
     *  @return
     */
    private ExportOptions getOptions(Dataset ds, String strParams)
    {
        if(strParams==null || strParams.isEmpty())
        {
            String[] regs = ds.listRegionsNames();
            String[] tmp = new String[regs.length+1];
            tmp[0] = "CDS";
            for(int i=1;i<tmp.length;i++)
                tmp[i] = regs[i-1];
            return (new OptionsForm()).showOptionsDialog(ds.listPopulations(), tmp);
        }
        /** Parse the parameters.
         Expected format:
         pops='<POP1;POP2>' region='<REG>' strlimit='<LIMIT>'
         "range='<POS1>-<POS2>;<POS3>'
            - pops:         populations to export
            - region:       region type
            - strlimit:     max. number of strains in each population
            - range:        range. Can be null.
        */
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "region='(.+)'\\s+"+                // 2
                                    "strlimit='(\\d*)'\\s+"+            // 3
                                    "range='([-,;\\d\\s]*)'\\s+",       // 4
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            ExportOptions eo = new ExportOptions();
            eo.pops = m.group(1).split("\\s*;\\s*");
            eo.strReg = m.group(2);
            eo.nStrains = (m.group(3).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(3));
            eo.strSites = m.group(4);
            return eo;
        }
        String[] regs = ds.listRegionsNames();
        String[] tmp = new String[regs.length+1];
        tmp[0] = "CDS";
        for(int i=1;i<tmp.length;i++)
            tmp[i] = regs[i-1];
        return (new OptionsForm()).showOptionsDialog(ds.listPopulations(), tmp);
    }

    /**
     *  Extracts the positions to use from the parameter string.
     *
     *  @param strPositions
     *  @return
     */
    private int[] extractSites(String strSites)
    {
        Vector<Integer> tmp = new Vector<Integer>();
        // Search for pairs.
        Pattern pair = Pattern.compile("(\\d+)-(\\d+)");
        Matcher m = pair.matcher(strSites);
        while(m.find())
        {
            // Make sure, the sites are always positive and the first site
            // is 0.
            int start = Math.max(Integer.parseInt(m.group(1))-1,0);
            int end = Math.max(Integer.parseInt(m.group(2))-1,0);
            // Add range.
            for(int i=Math.min(start, end);i<=Math.max(start,end);i++)
                tmp.add(new Integer(i));
        }
        // Search for single numbers.
        Pattern num = Pattern.compile("[,;\\s]+(\\d+)[,;\\s]*");
        m = num.matcher(strSites);
        while(m.find())
            tmp.add(Math.max(Integer.parseInt(m.group(1))-1,0));
        // Sort positions.
        if(tmp.size()==0)
            return null;
        int[] sites = new int[tmp.size()];
        for(int i=0;i<tmp.size();i++)
            sites[i] = tmp.get(i);
        Arrays.sort(sites);
        return sites;
    }
}