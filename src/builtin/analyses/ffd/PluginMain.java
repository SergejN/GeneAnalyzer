/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Analyses the four-fold degenerate sites.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.ffd;

import algorithms.SitesBlock;
import algorithms.BasicStatistics;
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
    
    private AInitData initData   = null;                                
    private String strLastError = "";
    
    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Four-fold degenerate sites";
    }

    public String GetName()
    {
        return "Four-fold degenerate sites analyzer";
    }

    public String GetDescription()
    {
        return "Performs divergence and polymorphism tests using four-fold degenerate sites.";
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
        return "pop='<POP>' out='<OUT>' nonffd='<T/F>' nonsynffd='<T/F>' gap='<T/F>' "+
               "strlimit='<LIMIT>' si_freq='<FREQ>' jc_pi='<T/F>' jc_t='<T/F>' jc_k='<T/F>' " +
               "exclAll='<T/F>' output='<OUTFILE>'";
    }
    
    public String GetLastError()
    {
        return strLastError;
    }
    
    @Override
    public ErrorCode AnalyzeDataset(Dataset dataset, String params)
    {
        AnalysisOptions ao = getOptions(params, dataset.listPopulations());
        if(ao==null)
        {
            strLastError = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Analysis);
        // Perform analysis.
        StringBuffer content = new StringBuffer();
        content.append("Analysis type: Four-fold degenerate sites\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Outgroup: %s\n", ao.strOut));
        content.append("Jukes-Cantor corrections:\n");
        content.append(String.format("\tPi: %s\n", ao.bJC_Pi));
        content.append(String.format("\tTheta: %s\n", ao.bJC_Theta));
        content.append(String.format("\tK: %s\n", ao.bJC_K));
        content.append(String.format("Singletons cut-off frequency: %f\n", ao.cof));
        content.append(String.format("Codon table: %s\n\n", initData.ct.getName()));
        content.append(String.format(NAMES, ao.strPop, ao.strPop, ao.strOut, ao.strPop, ao.strPop, ao.strOut));
        content.append(HEADER);
        // Iterate through the dataset.
        int nGenes = dataset.getGenesCount();
        for(int i=0;i<nGenes;i++)
        {            
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            if(nStrains==0 || ge.getStrainEntry(0).getRegionsCount(GeneRegion.EXON)==0)
                continue;
            // Iterate through the strains and select the ones which belong either
            // to the population of interest or to the outgroup.
            Vector<String> pop = new Vector<String>();
            Vector<String> out = new Vector<String>();
            for(int n=0;n<nStrains;n++)
            {
                StrainEntry se = ge.getStrainEntry(n);
                if(se.belongsToPopulation(ao.strPop) && pop.size()<ao.maxstr)
                    pop.add(se.getCodingSequence());
                if(se.belongsToPopulation(ao.strOut) && out.size()<ao.maxstr)
                    out.add(se.getCodingSequence());
            }
            // The analysis only makes sence if the number of sequences in the
            // population of interest is at least 2. It is not important here
            // whether the outgroup population is empty or not, since some analyses
            // can be performed even without outgroup.
            if(pop.size()<2)
                continue;
            String strLine = generateGeneStatistics(pop, out, ao);
            if(strLine!=null)
            {
                content.append(ge.getCommonName());
                content.append("\t\t"+strLine+"\n");
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
                strLastError = "An I/O error occured while saving the file";
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

    private String generateGeneStatistics(Vector<String> pop, Vector<String> out, AnalysisOptions ao)
    {
        int length = pop.get(0).length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        // Skip the case, that the ORF is too short to be analyzed at all.
        if(length<3)
            return null;
        // Data blocks.
        SitesBlock[] spb = new SitesBlock[pop.size()]; // Single population
        SitesBlock[] tpb = new SitesBlock[pop.size()]; // Two populations
        for(int n=1;n<=spb.length;n++)
        {
            spb[n-1] = new SitesBlock(n, ao.bJC_Pi, ao.bJC_Theta, ao.cof);
            tpb[n-1] = new SitesBlock(n, ao.bJC_Pi, ao.bJC_Theta, ao.cof);
        }
        // Iterate through the coding sequence.
        float k = 0.0f;
        int nP = 0;
        float tsp = 0.0f;
        int nD = 0;
        float tsd = 0.0f;
        for(int pos=0; pos<length; pos+=3)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos, ao);
            if(scp==null)
                continue;
            int nvb = scp.getValidBasesCount();
            if(nvb==0 || (nvb<4 && ao.bExclAll) )
                continue;
            spb[nvb-1].addSite(scp);
            SiteComposition sco = calculateSiteComposition(out, pos, ao);
            if(sco!=null)
            {
                // If the user selected to not to use nonsyn. FFD sites and the
                // codons in the population of interest and in the outgroup are
                // not synonymous, skip the site.
                if(ao.nonsyn)
                {
                    String strCodon1 = pop.get(0).substring(pos, pos+3);
                    String strCodon2 = out.get(0).substring(pos, pos+3);
                    if(!initData.ct.areSynonymous(strCodon1, strCodon2))
                        continue;
                }
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
                             pop.size(),
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
     *  Creates a site composition object of the current site. This method returns
     *  null if the site should be skipped, depending on the analysis options.
     *  The method also returns null, if there are no FFD codons at the specified site.
     *
     *  @param pop
     *  @param pos
     *  @param ao
     *  @return
     */
    private SiteComposition calculateSiteComposition(Vector<String> pop, int pos, AnalysisOptions ao)
    {
        if(pop.size()==0)
            return null;
        SiteComposition sc = new SiteComposition();
        String strRef = pop.get(0).substring(pos, pos+3);
        boolean b = false; // Flag, specifying, whether any sequence has a FFD codon at the specified site.
        int nStrains = pop.size();
        for(int n=0;n<nStrains;n++)
        {
            String codon = pop.get(n).substring(pos, pos+3);
            // 1. Check whether the sites with at least one gap must be
            // excluded. Otherwise only the codon is skipped.
            if(codon.contains("-"))
            {
                if(ao.gaps)
                    return null;
                else
                    continue;
            }
            if(!codon.matches("XXX"))
            {
                // 2. Check whether the codon is a FFD codon.
                // If it is not and the site must be excluded, then
                // return null, otherwise go to the next iteration, so that
                // the site is not counted.
                if(initData.ct.getFoldFamily(codon)!=4)
                {
                    if(ao.nonffd)
                        return null;
                    else
                        continue;
                }
                // 3. Check whether the codon is syn. to the reference codon.
                // If they are not synonymous, then either reject the whole site
                // and return null or proceed with all sequences.
                if(!initData.ct.areSynonymous(codon, strRef) && ao.nonsyn)
                    return null;
                b = true;
            }
            sc.addBase(codon.charAt(2));
        }
        return (b) ? sc : null;
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
            pop='<POPULATION>' out='<OUTGROUP>' nonffd='<T/F>' nonsynffd='<T/F>' 
            gap='<T/F>' strlimit='<LIMIT>' si_freq='<FREQ>' 
            jc_pi='<T/F>'  jc_t='<T/F>' jc_k='<T/F>' exclAll='<T/F>' output='<FILENAME>'
            
            - pop:          population of interest
            - out:          outgroup
            - nonffd:       whether or not to exclude non ffd sites
            - nonsynffd:    whether or not to exclude nonsyn. ffd sites
            - gap:          whether or not to exclude sites with gaps
            - strlimit:     maximal number of strains to 
            - si_freq:      singletons cut-off frequency
            - jc_pi:        whether or not to apply Jukes-Cantor correction to pi
            - jc_t:         whether or not to apply Jukes-Cantor correction to theta
            - jc_k:         whether or not to apply Jukes-Cantor correction to K
            - exclAll:      whether to exclude the blocks with less than 4 strains from all analyses
            - output:       output filename
        */ 
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "out='(.+)'\\s+"+                   // 2
                                    "nonffd='([TF])'\\s+"+              // 3
                                    "nonsynffd='([TF])'\\s+"+           // 4
                                    "gap='([TF])'\\s+"+                 // 5
                                    "strlimit='(\\d*)'\\s+"+            // 6
                                    "si_freq='(\\d*\\.{0,1}\\d*)'\\s+"+ // 7
                                    "jc_pi='([TF])'\\s+"+               // 8
                                    "jc_t='([TF])'\\s+"+                // 9
                                    "jc_k='([TF])'\\s+"+                // 10
                                    "exclAll='([TF])'\\s+"+             // 11
                                    "output='(.+)'$",                   // 12
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strOut = m.group(2);
            ao.nonffd = m.group(3).equalsIgnoreCase("T");
            ao.nonsyn = m.group(4).equalsIgnoreCase("T");
            ao.gaps   = m.group(5).equalsIgnoreCase("T");
            ao.maxstr = (m.group(6).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(6));
            ao.cof = (m.group(7).isEmpty()) ? 1.0f : Float.parseFloat(m.group(7));
            ao.bJC_Pi = m.group(8).equalsIgnoreCase("T");
            ao.bJC_Theta  = m.group(9).equalsIgnoreCase("T");
            ao.bJC_K  = m.group(10).equalsIgnoreCase("T");
            ao.bExclAll  = m.group(11).equalsIgnoreCase("T");
            ao.strOutput = m.group(12);
            return ao;
        }
        else
        {
            return (new OptionsDialog(pops)).getOptions();
        }
    }
}
