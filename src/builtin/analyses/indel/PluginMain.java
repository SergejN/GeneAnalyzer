/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Performs indels analysis on the region of specified type.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.indel;

import algorithms.SequenceRoutines;
import algorithms.SiteComposition;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetAnalyzer;


public class PluginMain extends ADatasetAnalyzer
{
    private final String HEADER = "Gene\t\tSites\t\tSample size\t#Ins (D)\t#Dels (D)\t#Ins (P)\t#Dels (P)\t\t" +
                                                   "Sample size\t#Ins (D)\t#Dels (D)\t#Ins (P)\t#Dels (P)\t\t#Unk (D)\t#Unk (P)\n";

    private AInitData initData = null;
    private String strErr = "";


    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "InDels";
    }

    public String GetName()
    {
        return "Indels analyzer";
    }

    public String GetDescription()
    {
        return "Performs indels analysis";
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
        return "pop='<POP>' out='<OUT>' dist='<DIST>' regtype='<REG>' maxlen='<MAX>' minlen='<MIN>' "+
               "combine='<T/F>' strlimit='<LIMIT>' lenrange='<T/F>' output='<OUTFILE>'";
    }
    
    public String GetLastError()
    {
        return strErr;
    }

    
    @Override
    public ErrorCode AnalyzeDataset(Dataset dataset, String params)
    {
        AnalysisOptions ao = getOptions(params, dataset.listPopulations(), dataset.listRegionsNames());
        if(ao==null)
        {
            strErr = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Analysis);
        // Outout file content.
        StringBuffer content = new StringBuffer();
        content.append("Analysis type: InDels analysis\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Outgroup: %s\n", ao.strOut));
        content.append(String.format("More distant group: %s\n", ao.strDist));
        content.append(String.format("Region type: %s\n", ao.strType));
        content.append(String.format("\t\t\t\t\t\t%s\t\t\t\t\t\t%s\n", ao.strPop, ao.strOut));
        content.append(HEADER);
        // Iterate through the dataset.
        int nGenes = dataset.getGenesCount();
        for(int i=0;i<nGenes;i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            if(nStrains==0 || ge.getStrainEntry(0).getRegionsCount(ao.strType)==0)
                continue;
            // Iterate through the strains and select the ones which belong either
            // to the population of interest or to the outgroup.
            Vector<StrainEntry> pop = new Vector<StrainEntry>();
            Vector<StrainEntry> out = new Vector<StrainEntry>();
            Vector<StrainEntry> mdg = new Vector<StrainEntry>(); // More distant group.
            for(int n=0;n<nStrains;n++)
            {
                StrainEntry se = ge.getStrainEntry(n);
                if(se.belongsToPopulation(ao.strPop) && pop.size()<ao.maxstr)
                    pop.add(se);
                if(se.belongsToPopulation(ao.strOut) && out.size()<ao.maxstr)
                    out.add(se);
                if(se.belongsToPopulation(ao.strDist) && mdg.size()<ao.maxstr)
                    mdg.add(se);
            }
            // The analysis only makes sence if the number of sequences in all
            // three populations is at least 1.
            if(pop.size()<1 || out.size()<1 || mdg.size()<1)
                continue;
            String[] tmp = analyzeGene(pop, out, mdg, ao);
            if(tmp!=null)
            {
                content.append(ge.getCommonName());
                for(String s:tmp)
                    content.append("\t\t"+s+"\n");
            }
        }
        // Save the data into the file.
        if(ao.strOutput!=null && !ao.strOutput.isEmpty())
        {
            try
            {
                PrintWriter out = new PrintWriter(new FileWriter(new File(ao.strOutput)));
                out.print(content.toString());
                out.close();
            }
            catch (IOException e)
            {
                strErr = "An I/O error occured while saving the file";
                initData.wd.close();
                return ErrorCode.IOError;
            }
        }
        // Display the results if necessary.
        ResultsWindow rw = null;
        if(ao.bShowRes)
        {
            rw = new ResultsWindow();
            rw.setResults(content.toString());
        }
        initData.wd.close();
        if(rw!=null)
            rw.displayResults();
        return ErrorCode.Ok;
    }

    private String[] analyzeGene(Vector<StrainEntry> pop, Vector<StrainEntry> out, 
                                 Vector<StrainEntry> mdg, AnalysisOptions ao)
    {
        Vector<String> res = new Vector<String>();  // Result lines
        String[] popseq = new String[pop.size()];   // Pop. sequences
        String[] outseq = new String[out.size()];   // Out. sequences
        String[] mdgseq = new String[mdg.size()];   // More distant group sequences
        int nRegs = pop.get(0).getRegionsCount();
        for(int i=0;i<nRegs;i++)
        {
            // Check the region type and the length.
            if(pop.get(0).getRegion(i).hasType(ao.strType) && (isRegionLengthOk(pop, out, ao, i)) )
            {
                popseq = generateSequenceSet(popseq, pop, i);
                outseq = generateSequenceSet(outseq, out, i);
                mdgseq = generateSequenceSet(mdgseq, mdg, i);
                // If the regions should be combined, then continue appending
                // the sequences to the array, otherwise run the analysis
                // for one region and clear the array for the next iteration.
                if(ao.bCombine)
                    continue;
                else
                {
                    String tmp = analyzeSequence(popseq, outseq, mdgseq);
                    if(tmp!=null)
                        res.add(tmp);
                    popseq = new String[pop.size()];
                    outseq = new String[out.size()];
                    mdgseq = new String[mdg.size()];
                }
            }
        }
        if(ao.bCombine && popseq[0]!=null)
        {
            String tmp = analyzeSequence(popseq, outseq, mdgseq);
            if(tmp!=null)
                res.add(tmp);
        }
        return (res.size()>0) ? res.toArray(new String[1]) : null;
    }

    /**
     *  Analyzes the sequence. There are three patterns of sequence which can be
     *  combined with one another:
     *
     *      X       X       -
     *      X       -       -
     *
     *  The possible combinations are:
     *
     *      1   2   3           4   5   6           7   8   9
     *      X   X   X           X   X   X           -   -   -
     *      X   X   X           -   -   -           -   -   -
     *
     *      X   X   -           X   X   -           X   X   -
     *      X   -   -           X   -   -           X   -   -
     *
     *      -   D               D   U   I               I   U
     *          O               P       P               O
     *          p               p       p               p
     *
     *  D - deletion
     *  I - insertion
     *  U - unknown
     *  P - population
     *  O - outgroup
     *  d - divergent
     *  p - polymorphic
     *
     *  @param pop
     *  @param out
     *  @param mdg
     *  @return
     */
    private String analyzeSequence(String[] pop, String[] out, String[] mdg)
    {
        boolean b = false;
        for(String s:mdg)
        {
            // If the sequence is not a huge gap, leave the cycle and proceed to analysis.
            if(!s.matches("-+"))
            {
                b = true;
                break;
            }
        }
        if(!b)
            return null;
        int l = pop[0].length();
        int[] vals_p = new int[4];  // Ins(D), Del(D), Ins(P), Del(P)
        int[] vals_o = new int[4];  // Ins(D), Del(D), Ins(P), Del(P)
        int nSites = 0;
        int[] unk  = new int[]{0, 0};   // Unk(D), Unk(P)
        int lsg = -1;   // Last seen gap
        for(int i=0;i<l;i++)
        {
            if(!allPopulationsAvailable(pop, out, mdg, i))
                continue;
            nSites++;
            SiteComposition scp = calculateSiteComposition(pop, i);
            SiteComposition sco = calculateSiteComposition(out, i);
            SiteComposition scd = calculateSiteComposition(mdg, i);
            int nGapsPop = scp.getGapsCount();
            int nGapsOut = sco.getGapsCount();
            int nGapsMdg = scd.getGapsCount();
            int nPop = pop.length-scp.getBaseCount('X');
            int nOut = out.length-sco.getBaseCount('X');
            int nMdg = mdg.length-scd.getBaseCount('X');
            // Case 1: no gaps. Or case 9: only gaps.
            if( (nGapsPop==0 && nGapsOut==0) || (nGapsPop==nPop && nGapsOut==nOut) )
                continue;
            // Check whether new new site has the same pattern as the last one.
            if(lsg>-1 && lsg==i-1)
            {
                if(haveSamePattern(pop, out, mdg, i, lsg))
                {
                    lsg = i;
                    continue;
                }
            }
            lsg = i;
            // Case 2.
            if(nGapsPop==0 && nGapsOut<nOut && nGapsOut>0)
                vals_o[3]++;
            // Case 4.
            else if(nGapsOut==0 && nGapsPop<nPop && nGapsPop>0)
                vals_p[3]++;
            // Case 5.
            else if(nGapsPop>0 && nGapsPop<nPop && nGapsOut>0 && nGapsOut<nOut)
                unk[1]++;
            // Case 6.
            else if(nGapsOut==nOut && nGapsPop>0 && nGapsPop<nPop)
                vals_p[2]++;
            // Case 8.
            else if(nGapsPop==nPop && nGapsOut>0 && nGapsOut<nOut)
                vals_o[2]++;
            // Case 3: second outgroup required.
            else if(nGapsPop==0 && nGapsOut==nOut)
            {
                if(nGapsMdg==0)
                    vals_o[1]++;
                else if(nGapsMdg>0 && nGapsMdg<nMdg)
                    unk[0]++;
                else if(nGapsMdg==nMdg)
                    vals_p[0]++;
            }
            // Case 7: second outgroup required.
            else if(nGapsPop==nPop && nGapsOut==0)
            {
                if(nGapsMdg==0)
                    vals_p[1]++;
                else if(nGapsMdg>0 && nGapsMdg<nMdg)
                    unk[0]++;
                else if(nGapsMdg==nMdg)
                    vals_o[0]++;
            }
        }
        return String.format(initData.locale, "%d\t\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d",
                             nSites, pop.length, vals_p[0], vals_p[1], vals_p[2], vals_p[3],
                             out.length, vals_o[0], vals_o[1], vals_o[2], vals_o[3], unk[0], unk[1]);
    }

    private boolean allPopulationsAvailable(String[] pop, String[] out, String[] mdg, int pos)
    {
        boolean b = false;
        for(String s:pop)
        {
            if(!s.substring(0, pos+1).matches("[-Xx]+") && !s.substring(pos).matches("[-Xx]+"))
            {
                b = true;
                break;
            }
        }
        if(!b)
            return false;
        b = false;
        for(String s:out)
        {
            if(!s.substring(0, pos+1).matches("[-Xx]+") && !s.substring(pos).matches("[-Xx]+"))
            {
                b = true;
                break;
            }
        }
        if(!b)
            return false;
        b = false;
        for(String s:mdg)
        {
            if(!s.substring(0, pos+1).matches("[-Xx]+") && !s.substring(pos).matches("[-Xx]+"))
            {
                b = true;
                break;
            }
        }
        return b;
    }

    /**
     *  Returns true if the old site has the same pattern as the new one.
     *  
     *  @param pop
     *  @param out
     *  @param mdg
     *  @param np
     *  @param op
     *  @return
     */
    private boolean haveSamePattern(String[] pop, String[] out, String[] mdg, int np, int op)
    {
        for(String s:pop)
        {
            char co = s.charAt(op);
            char cn = s.charAt(np);
            if( (co=='-' && cn!='-') || (cn=='-' && co!='-'))
                return false;
        }
        for(String s:out)
        {
            char co = s.charAt(op);
            char cn = s.charAt(np);
            if( (co=='-' && cn!='-') || (cn=='-' && co!='-'))
                return false;
        }
      
        for(String s:mdg)
        {
            char co = s.charAt(op);
            char cn = s.charAt(np);
            if( (co=='-' && cn!='-') || (cn=='-' && co!='-'))
                return false;
        }
        
        return true;
    }

    /**
     *  Extracts the sequence of the specified region and appends it to already extracted region if necessary.
     *
     *  @param old
     *  @param pop
     *  @param iRegion
     *  @return
     */
    private String[] generateSequenceSet(String[] old, Vector<StrainEntry> pop, int iRegion)
    {
        int nStrains = pop.size();
        String[] seqs = new String[nStrains];
        for(int n=0;n<nStrains;n++)
        {
            String seq = pop.get(n).getRegion(iRegion).getSequence();
            seqs[n] = (old[n]==null) ? seq : old[n]+seq;
        }
        return seqs;
    }

    /**
     *  Creates a site composition object of the current site.
     *
     *  @param seqs
     *  @param pos
     *  @return
     */
    private SiteComposition calculateSiteComposition(String[] seqs, int pos)
    {
        SiteComposition sc = new SiteComposition();
        for(int n=0;n<seqs.length;n++)
        {
            char base = seqs[n].charAt(pos);
            sc.addBase(seqs[n].charAt(pos));
        }
        return sc;
    }


    /**
     *  Iterates through both population and outgroup and returns true only if
     *  the minimal and the maxinal length are not beyond the length threshold.
     *
     *  @param pop
     *  @param out
     *  @param ao
     *  @param index
     *  @return
     */
    private boolean isRegionLengthOk(Vector<StrainEntry> pop,
                                     Vector<StrainEntry> out,
                                     AnalysisOptions ao,
                                     int index)
    {        
        int nPopSize = pop.size();
        int nOutSize = out.size();
        String[] seqs = new String[nPopSize+nOutSize];
        for(int i=0;i<nPopSize;i++)
            seqs[i] = pop.get(i).getRegion(index).getSequence();
        for(int i=0;i<nOutSize;i++)
            seqs[nPopSize+i] = out.get(i).getRegion(index).getSequence();
        if(ao.bLenRange)
        {
            int[] tmp = SequenceRoutines.getMinMaxUngappedLength(seqs);
            return ( (tmp[1]<=ao.maxlen) && (tmp[0]>=ao.minlen) );
        }
        else
        {
            seqs = SequenceRoutines.getUngappedSequences(seqs);
            int l = seqs[0].length();
            return ( (l<=ao.maxlen) && (l>=ao.minlen) );
        }
    }
        
    private AnalysisOptions getOptions(String strParams, String[] pops, String[] regs)
    {
        // If the parameters are not specified or empty, show the dialog.
        if(strParams==null || strParams.isEmpty())
        {            
            return (new OptionsDialog(pops, regs)).getOptions();
        }
        // Else parse the parameters.
        /** The parameters line should have the following format:
            pop='<POP>' out='<OUT>' dist='<DIST>' regtype='<REG>' maxlen='<MAX>' minlen='<MIN>'
            combine='<T/F>' strlimit='<LIMIT>' lenrange='<T/F>' output='<OUTFILE>'
            
            - pop:          population of interest
            - out:          outgroup
            - dist:         more distant group to use, if the outgroup is not enough
            - regtype:      type of the region to analyze
            - maxlen:       maximal region length. Use 0 for undefined
            - minlen:       minimal region length. Use 0 for undefined
            - combine:      whether or not to combine the regions
            - strlimit:     maximal number of strains to use
            - lenrange:     whether the minimal and maximal length should fall into the specified range
            - output:       output filename
        */ 
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "out='(.+)'\\s+"+                   // 2
                                    "dist='(.+)'\\s+"+                  // 3
                                    "regtype='(.+)'\\s+"+               // 4
                                    "maxlen='(\\d+)'\\s+"+              // 5
                                    "minlen='(\\d+)'\\s+"+              // 6
                                    "combine='([TF])'\\s+"+             // 7
                                    "strlimit='(\\d*)'\\s+"+            // 8
                                    "lenragnge='([TF])'\\s+"+           // 9
                                    "output='(.+)'$",                   // 10
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strOut = m.group(2);
            ao.strDist = m.group(3);
            ao.strType = m.group(4);
            ao.maxlen = Integer.parseInt(m.group(5));
            if(ao.maxlen==0)
                ao.maxlen = Integer.MAX_VALUE;
            ao.minlen = Integer.parseInt(m.group(6));
            if(ao.maxlen<ao.minlen)
            {
                int tmp = ao.maxlen;
                ao.maxlen = ao.minlen;
                ao.minlen = tmp;
            }
            ao.bCombine  = m.group(7).equalsIgnoreCase("T");
            ao.bLenRange = m.group(9).equalsIgnoreCase("T");
            ao.strOutput = m.group(10);
            // Strains number limit.
            ao.maxstr = (m.group(8).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(8));
            return ao;
        }
        else
            return (new OptionsDialog(pops, regs)).getOptions();
    }
}
