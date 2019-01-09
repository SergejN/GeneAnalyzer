/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Performs short introns analysis. If a gene has multiple introns
        only the short ones are analyzed. 
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.introns;

import algorithms.SitesBlock;
import algorithms.BasicStatistics;
import algorithms.SequenceRoutines;
import algorithms.SiteComposition;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
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
    private final String NAMES  = "\t\t\t%s\t\t\t\t\t\t\t\t\t\t\t%s+%s (%s)\t\t\t\t\t\t\t\t\t%s+%s\n";
    private final String HEADER = "Gene\t\tSample\tSites\tPi\tTheta\tTajD\tTajD'\tP\tSi\tTS\tTV\t"+
                                        "\tSites\tPi\tTheta\tTajD\tTajD'\tP\tSi\tTS\tTV\t"+
                                        "\tK\tP\tTSp\tTVp\tD\tTSd\tTVd\n";

    private AInitData initData = null;
    private String strErr = "";


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
        return "Introns analyzer";
    }

    public String GetDescription()
    {
        return "Performs introns analysis";
    }

    public boolean SupportsMissingData()
    {
        return true;
    }

    public boolean SupportsAmbiguousData()
    {
        return false;
    }

    public String GetParamString()
    {
        return "pop='<POP>' out='<OUT>' maxlen='<MAX>' minlen='<MIN>' nogtag='<T/F>' "+
               "combine='<T/F>' range='<POS1>-<POS2>;<POS3>' strlimit='<LIMIT>' " +
               "si_freq='<FREQ>' jc_pi='<T/F>' jc_t='<T/F>' jc_k='<T/F>' exclAll='<T/F>' lenrange='<T/F>' output='<OUTFILE>'";
    }
    
    public String GetLastError()
    {
        return strErr;
    }

    
    @Override
    public ErrorCode AnalyzeDataset(Dataset dataset, String params)
    {
        AnalysisOptions ao = getOptions(params, dataset.listPopulations());
        if(ao==null)
        {
            strErr = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Analysis);
        // Outout file content.
        StringBuffer content = new StringBuffer();
        content.append("Analysis type: Introns analysis\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Outgroup: %s\n", ao.strOut));
        content.append("Jukes-Cantor corrections:\n");
        content.append(String.format("\tPi: %s\n", ao.bJC_Pi));
        content.append(String.format("\tTheta: %s\n", ao.bJC_Theta));
        content.append(String.format("\tK: %s\n", ao.bJC_K));
        content.append(String.format("Singletons cut-off frequency: %f\n", ao.cof));
        content.append(String.format("Codon table: %s\n", initData.ct.getName()));
        content.append("\nNOTE:\nSample size may not reflect the real sample size used to " +
                            "calculate statistics for each site if the sample contains missing data!\n\n");
        content.append(String.format(NAMES, ao.strPop, ao.strPop, ao.strOut, ao.strPop, ao.strPop, ao.strOut));
        content.append(HEADER);
        // Iterate through the dataset.
        int nGenes = dataset.getGenesCount();
        for(int i=0;i<nGenes;i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            if(nStrains==0 || ge.getStrainEntry(0).getRegionsCount(GeneRegion.INTRON)==0)
                continue;
            // Iterate through the strains and select the ones which belong either
            // to the population of interest or to the outgroup.
            Vector<StrainEntry> pop = new Vector<StrainEntry>();
            Vector<StrainEntry> out = new Vector<StrainEntry>();
            for(int n=0;n<nStrains;n++)
            {
                StrainEntry se = ge.getStrainEntry(n);
                if(se.belongsToPopulation(ao.strPop) && pop.size()<ao.maxstr)
                    pop.add(se);
                if(se.belongsToPopulation(ao.strOut) && out.size()<ao.maxstr)
                    out.add(se);
            }
            // The analysis only makes sence if the number of sequences in the
            // population of interest is at least 2. It is not important here
            // whether the outgroup population is empty or not, since some analyses
            // can be performed even without outgroup.
            if(pop.size()<2)
                continue;
            String[] lines = generateGeneStatistics(pop, out, ao);
            if(lines!=null)
            {
                content.append(ge.getCommonName());
                for(String s:lines)
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

    private String[] generateGeneStatistics(Vector<StrainEntry> pop, Vector<StrainEntry> out, AnalysisOptions ao)
    {
        Vector<String> res = new Vector<String>();  // Result string(s)
        String[] popseq = new String[pop.size()];   // Pop. sequences
        String[] outseq = new String[out.size()];   // Out. sequences
        StrainEntry se = pop.get(0);
        int nRegs = se.getRegionsCount();
        for(int i=0;i<nRegs;i++)
        {
            // Check the region type and the intron length.
            if(se.getRegion(i).hasType(GeneRegion.INTRON) && (isIntronLengthOk(pop, out, ao, i)) )
            {
                // Distinguish between three options:
                //  1. Complete intron sequence
                //  2. Intron sequence w/o GT..AG
                //  3. Specific sites
                // The first two variants are very similar, whereas the last one
                // needs an extra treatment.
                if(ao.sites==null) // Variants 1 and 2.
                {
                    popseq = generateSequenceSet(popseq, pop, i, ao.bNogtag);
                    outseq = generateSequenceSet(outseq, out, i, ao.bNogtag);                                       
                }
                else
                {
                    generateSequenceSet(popseq, outseq, pop, out, i, ao.sites);
                }
                // If the regions should be combined, then continue appending
                // the sequences to the array, otherwise run the analysis
                // for one region and clear the array for the next iteration.
                if(ao.bCombine)
                    continue;
                else
                {
                    res.add(generateSequenceStatistics(popseq, outseq, ao));
                    popseq = new String[pop.size()];
                    outseq = new String[out.size()];                    
                }
            }
        }
        if(ao.bCombine && popseq[0]!=null)
            res.add(generateSequenceStatistics(popseq, outseq, ao));
        return (res.size()>0) ? res.toArray(new String[1]) : null;
    }

    private String generateSequenceStatistics(String[] popseq, String[] outseq, AnalysisOptions ao)
    {
        int l = popseq[0].length();
        SitesBlock[] spb = new SitesBlock[popseq.length]; // Single population
        SitesBlock[] tpb = new SitesBlock[popseq.length]; // Two populations
        for(int i=1;i<=spb.length;i++)
        {
            spb[i-1] = new SitesBlock(i, ao.bJC_Pi, ao.bJC_Theta, ao.cof);
            tpb[i-1] = new SitesBlock(i, ao.bJC_Pi, ao.bJC_Theta, ao.cof);
        }
        float k = 0.0f;
        int nP = 0;
        float tsp = 0.0f;
        int nD = 0;
        float tsd = 0.0f;
        for(int i=0;i<l;i++)
        {
            // For each site find out its composition.
            SiteComposition scp = calculateSiteComposition(popseq, i);
            if(scp==null)
                continue;
            int nvb = scp.getValidBasesCount();
            if(nvb==0 || (nvb<4 && ao.bExclAll) )
                continue;
            spb[nvb-1].addSite(scp);
            SiteComposition sco = calculateSiteComposition(outseq, i);
            if(sco!=null)
            {
                tpb[nvb-1].addSite(scp);
                // Combined site composition.
                SiteComposition scc = SiteComposition.merge(scp, sco);
                float f = BasicStatistics.calculateK(scp, sco);
                // Increase D and TS/TV only if there are at least 2 valid bases at this site.
                if(f>=1.0f && scp.getValidBasesCount()>1)
                {
                    tsd += scc.getNumberOfTransitions();
                    nD++;
                }
                else if(f>0.0f)
                {
                    tsp += scc.getNumberOfTransitions();
                    nP++;
                }
                k += f;                
            }
        }
        k = k/SitesBlock.getSitesCount(tpb);
        if(ao.bJC_K)
            k = BasicStatistics.correctJC(k);
        return String.format(initData.locale, "%d\t%s\t\t%s\t\t%f\t%d\t%.1f\t%.1f\t%d\t%.1f\t%.1f",
                             popseq.length,
                             formatString(spb, ao.bJC_Pi, ao.bJC_Theta),
                             formatString(tpb, ao.bJC_Pi, ao.bJC_Theta),
                             k, nP, tsp, nP-tsp, nD, tsd, nD-tsd);
    }

    /**
     *  Formats the results string.
     *
     *  @param blocks
     *  @param bJC_Pi
     *  @param bJC_Theta
     *  @return
     */
    private String formatString(SitesBlock[] blocks, boolean bJC_Pi, boolean bJC_Theta)
    {
        int n = SitesBlock.getPolymorphismsCount(blocks);
        if(n==0)
        {
            return String.format(initData.locale, "%d\t%f\t%f\t%f\t%f\t%d\t%d\t%.1f\t%.1f",
                                 (int)SitesBlock.getSitesCount(blocks), 0.0f, 0.0f, Float.NaN, Float.NaN, n, 0, 0.0f, 0.0f);
        }
        else
        {
            float f = SitesBlock.getTransitionsCount(blocks);
            return String.format(initData.locale, "%d\t%f\t%f\t%f\t%f\t%d\t%d\t%.1f\t%.1f",
                                 (int)SitesBlock.getSitesCount(blocks), SitesBlock.getPi(blocks, bJC_Pi),
                                 SitesBlock.getTheta(blocks, bJC_Theta), BasicStatistics.calculateTajD(blocks),
                                 BasicStatistics.calculateTajDPrime(blocks), n, SitesBlock.getSingletonsCount(blocks),
                                 f, n-f);
        }
    }

    /**
     *  Extracts the sequence of the specified region and appends it to already extracted region if necessary.
     *
     *  @param old
     *  @param pop
     *  @param iRegion
     *  @param bNogtag
     *  @return
     */
    private String[] generateSequenceSet(String[] old, Vector<StrainEntry> pop, int iRegion, boolean bNogtag)
    {
        int nStrains = pop.size();
        String[] seqs = new String[nStrains];
        for(int n=0;n<nStrains;n++)
        {
            String seq = pop.get(n).getRegion(iRegion).getSequence();
            if(bNogtag)
                seq = seq.substring(2, seq.length()-2);
            seqs[n] = (old[n]==null) ? seq : old[n]+seq;
        }
        return seqs;
    }

    /**
     *  Extracts the sequence of the specified region and extracts only the specified sites from the sequence.
     *  This method skippes the sites, at which there is a gap in any sequence.
     *
     *  @param old_p
     *  @param old_o
     *  @param pop
     *  @param out
     *  @param iRegion
     *  @param sites
     */
    private void generateSequenceSet(String[] old_p, String[] old_o, Vector<StrainEntry> pop, Vector<StrainEntry> out, int iRegion, int[] sites)
    {
        int nPopSize = pop.size();
        int nOutSize = out.size();
        String[] seqs = new String[nPopSize+nOutSize];
        for(int i=0;i<nPopSize;i++)
            seqs[i] = pop.get(i).getRegion(iRegion).getSequence();
        for(int i=0;i<nOutSize;i++)
            seqs[nPopSize+i] = out.get(i).getRegion(iRegion).getSequence();
        // Generate ungapped sequences.
        seqs = SequenceRoutines.getUngappedSequences(seqs);
        // Extract the sites.
        seqs = SequenceRoutines.extractAlignedSites(seqs, sites);
        for(int i=0;i<nPopSize;i++)
            old_p[i] = (old_p[i]==null) ? seqs[i] : old_p[i]+seqs[i];
        for(int i=0;i<nOutSize;i++)
            old_o[i] = (old_o[i]==null) ? seqs[nPopSize+i] : old_o[i]+seqs[nPopSize+i];
    }

    /**
     *  Creates a site composition object of the current site. This method returns
     *  null if there is a gap in any sequence at the specified position.
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
            if(base=='-')
                return null;
            else
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
    private boolean isIntronLengthOk(Vector<StrainEntry> pop,
                                     Vector<StrainEntry> out,
                                     AnalysisOptions ao,
                                     int index)
    {
        // Check the intron location: if the intron is the first or the last
        // annotated region in the strain entry, check its boundaries, and if
        // these are not correct, return false.
        if(index==0 || index==pop.get(0).getRegionsCount()-1)
        {
            // Iterate through the population and count the number of invalid
            // boundaries. If this number exceeds 25% assume the intron to be
            // incomplete and return false.
            //
            // Reason:
            // You can't use a single strain to check the boundaries, since it can have
            // a valid GT or AG pair (which is not a boundary) just by chance. Thus, 25%
            // is a fair value to use.
            int nInvalid = 0;
            int nMax = pop.size()/4;
            for(int i=0;i<pop.size();i++)
            {
                String seq = pop.get(i).getRegion(index).getSequence();
                if(!seq.startsWith("GT") || !seq.endsWith("AG"))
                    nInvalid++;
                if(nInvalid>nMax)
                    return false;
            }
        }
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
        
    private AnalysisOptions getOptions(String strParams, String[] pops)
    {
        // If the parameters are not specified or empty, show the dialog.
        if(strParams==null || strParams.isEmpty())
        {            
            return (new OptionsDialog(pops)).getOptions();
        }
        // Else parse the parameters.
        /** The parameters line should have the following format:
            pop='<POPULATION>' out='<OUTGROUP>' maxlen='<MAX>' minlen='<MIN>' nogtag='<T/F>'
            combine='<T/F>' range='<VALUES>' strlimit='<LIMIT>' si_freq='<FREQ>'
            jc_pi='<T/F>' jc_t='<T/F>' jc_k='<T/F>' exclAll='<T/F>' lenrange='<T/F>' output='<FILENAME>'
            
            - pop:          population of interest
            - out:          outgroup
            - maxlen:       maximal intron length. Use 0 for undefined
            - minlen:       minimal intron length. Use 0 for undefined
            - nogtag:       whether or not to exclude GT..AG
            - combine:      whether or not to combine the regions
            - range:        sites range. If specified, nogtag is ignored!
            - strlimit:     maximal number of strains to use
            - si_freq:      singletons cut-off frequency
            - jc_pi:        whether or not to apply Jukes-Cantor correction to pi
            - jc_t:         whether or not to apply Jukes-Cantor correction to theta
            - jc_k:         whether or not to apply Jukes-Cantor correction to K
            - exclAll:      whether to exclude the blocks with less than 4 strains from all analyses
            - lenrange:     whether the minimal and maximal length should fall into the specified range
            - output:       output filename
        */ 
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "out='(.+)'\\s+"+                   // 2
                                    "maxlen='(\\d*)'\\s+"+              // 3
                                    "minlen='(\\d*)'\\s+"+              // 4
                                    "nogtag='([TF])'\\s+"+              // 5
                                    "combine='([TF])'\\s+"+             // 6
                                    "range='([-,;\\d\\s]*)'\\s+"+       // 7
                                    "strlimit='(\\d*)'\\s+"+            // 8
                                    "si_freq='(\\d*\\.{0,1}\\d*)'\\s+"+ // 9
                                    "jc_pi='([TF])'\\s+"+               // 10
                                    "jc_t='([TF])'\\s+"+                // 11
                                    "jc_k='([TF])'\\s+"+                // 12
                                    "exclAll='([TF])'\\s+"+             // 13
                                    "lenrange='([TF])'\\s+"+            // 14
                                    "output='(.+)'",                    // 15
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strOut = m.group(2);
            if(!m.group(3).isEmpty())
                ao.maxlen = Integer.parseInt(m.group(3));
            else
                ao.maxlen = Integer.MAX_VALUE;
            if(!m.group(4).isEmpty())
                ao.minlen = Integer.parseInt(m.group(4));
            else
                ao.minlen = 0;
            if(ao.maxlen<ao.minlen)
            {
                int tmp = ao.maxlen;
                ao.maxlen = ao.minlen;
                ao.minlen = tmp;
            }
            ao.bNogtag   = m.group(5).equalsIgnoreCase("T");
            ao.bCombine  = m.group(6).equalsIgnoreCase("T");
            ao.bJC_Pi    = m.group(10).equalsIgnoreCase("T");
            ao.bJC_Theta = m.group(11).equalsIgnoreCase("T");
            ao.bJC_K     = m.group(12).equalsIgnoreCase("T");
            ao.bExclAll  = m.group(13).equalsIgnoreCase("T");
            ao.bLenRange = m.group(14).equalsIgnoreCase("T");
            ao.strOutput = m.group(15);
            // If range is specified and valid, set bNogtag to false.
            int[] sites = OptionsDialog.extractSites(m.group(7));
            ao.sites = sites;
            if(sites!=null)
                ao.bNogtag = false;
            // Strains number limit.
            ao.maxstr = (m.group(8).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(8));
            ao.cof = (m.group(9).isEmpty()) ? 1.0f : Float.parseFloat(m.group(9));
            return ao;
        }
        else
            return (new OptionsDialog(pops)).getOptions();
    }
}
