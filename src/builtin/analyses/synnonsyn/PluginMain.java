/*
    File:
        PluginMain.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Analyses the synonymous and nonsynonymous sites.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.synnonsyn;

import algorithms.BasicStatistics;
import algorithms.CodonComposition;
import algorithms.CodonsBlock;
import algorithms.Path;
import algorithms.SiteComposition;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import bio.gene.dna.Codon;
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
    private final String NAMES  = "\t\t\t%s\t\t\t"+
                                  "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s+%s (%s)" +
                                  "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s+%s\n";
    private final String HEADER = "Gene\t\tSample\tSyn.Sites\tPi(s)\tTheta(s)\tTajD(s)\tTajD'(s)\tP(s)\tSi(s)\tTS(s)\tTV(s)\t\t"+
                                                  "Nonsyn.Sites\tPi(n)\tTheta(n)\tTajD(n)\tTajD'(n)\tP(n)\tSi(n)\tTS(n)\tTV(n)\t\t" +
                                                  "Syn.Sites\tPi(s)\tTheta(s)\tTajD(s)\tTajD'(s)\tP(s)\tSi(s)\tTS(s)\tTV(s)\t\t"+
                                                  "Nonsyn.Sites\tPi(n)\tTheta(n)\tTajD(n)\tTajD'(n)\tP(n)\tSi(n)\tTS(n)\tTV(n)\t\t" +
                                                  "Syn.Sites\tK(s)\tP(s)\tTSp(s)\tTVp(s)\tD(s)\tTSd(s)\tTVd(s)\t\t"+
                                                  "Nonsyn.Sites\tK(n)\tP(n)\tTSp(n)\tTVp(n)\tD(n)\tTSd(n)\tTVd(n)\n";
    
    private String strLastErr = "";
    private AInitData initData = null;
    
    public ErrorCode Initialize(AInitData initdata)
    {
        this.initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Synonymous and nonsynonymous sites";
    }

    public String GetName()
    {
        return "Synonymous and nonsynonymous sites analyzer";
    }

    public String GetDescription()
    {
        return "Analyzes the synonymous and nonsynonymous sites within the CDS";
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
        return "pop='<POP>' out='<OUT>' strlimit='<LIMIT>' si_freq='<FREQ>' " +
               "jc_pi='<T/F>' jc_t='<T/F>' jc_k='<T/F>' use_term='<T/F>' " +
               "excl_term='<T/F>' exclAll='<T/F>' output='<OUTFILE>'";
    }

    public String GetLastError()
    {
        return strLastErr;
    }

    @Override
    public ErrorCode AnalyzeDataset(Dataset dataset, String params)
    {
        AnalysisOptions ao = getOptions(params, dataset.listPopulations());
        if(ao==null)
        {
            strLastErr = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Analysis);
        // Perform analysis.
        StringBuffer content = new StringBuffer();
        content.append("Analysis type: Synonymous and nonsynonymous sites\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Outgroup: %s\n", ao.strOut));
        content.append("Jukes-Cantor corrections:\n");
        content.append(String.format("\tPi: %s\n", ao.bJC_pi));
        content.append(String.format("\tTheta: %s\n", ao.bJC_t));
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
                strLastErr = "An I/O error occured while saving the file";
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

    /**
     *  Generates the statistics line.
     *
     *  @param pop
     *  @param out
     *  @param ao
     *  @return
     */
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
        CodonsBlock[] spb = new CodonsBlock[pop.size()]; // Single population
        CodonsBlock[] tpb = new CodonsBlock[pop.size()]; // Two populations
        for(int n=1;n<=spb.length;n++)
        {
            spb[n-1] = new CodonsBlock(n, ao.bJC_pi, ao.bJC_t, ao.cof, ao.bUseTer, initData.ct);
            tpb[n-1] = new CodonsBlock(n, ao.bJC_pi, ao.bJC_t, ao.cof, ao.bUseTer, initData.ct);
        }
        // Iterate through the coding sequence.
        float[] ks = {0.0f, 0.0f};
        int[] poly = {0, 0};
        int[] tsps = {0, 0};
        int[] tvps = {0, 0};
        int[] Ds   = {0, 0};
        int[] tsds = {0, 0};
        int[] tvds = {0, 0};
        float[] sites = {0.0f, 0.0f};
        for(int pos=0; pos<length; pos+=3)
        {
            CodonComposition ccp = calculateCodonComposition(pop, pos, ao.bUseTer, ao.bExclTer);
            if(ccp==null)
                continue;
            int nvb = ccp.getValidCodonsCount();
            if(nvb==0 || (nvb<4 && ao.bExclAll) )
                continue;
            spb[nvb-1].addCodon(ccp);
            CodonComposition cco = calculateCodonComposition(out, pos, ao.bUseTer, ao.bExclTer);
            if(cco!=null)
            {
                tpb[nvb-1].addCodon(ccp);
                float[] tmp = BasicStatistics.calculateK(ccp.getValidCodons(), cco.getValidCodons(), initData.ct, ao.bUseTer);
                // Combined codon composition.
                CodonComposition ccc = CodonComposition.merge(ccp, cco);
                float[] cs = ccc.getSitesCounts(); // Combined sites count.
                sites[0] += cs[0];
                sites[1] += cs[1];                
                ks[0] += tmp[0]*cs[0];
                ks[1] += tmp[1]*cs[1];
                if(isSiteDivergent(cco, ccp) && ccp.getValidCodonsCount()>1)
                {
                    int[] nn = countSubstitutions(ccp, cco, ao.bUseTer);
                    Ds[0] += nn[0];
                    tsds[0] += nn[1];
                    tvds[0] += nn[2];
                    Ds[1] += nn[3];
                    tsds[1] += nn[4];
                    tvds[1] += nn[5];
                }
                // Statistics: polymorphic sites.
                int[] p_p = ccp.getNumberOfPolymorphisms();
                int[] p_o = cco.getNumberOfPolymorphisms();
                poly[0] += p_p[0]+p_o[0];
                poly[1] += p_p[1]+p_o[1];
                // TS and TV.
                int[] subst_p = ccp.getSubstitutionsCount();
                int[] subst_o = cco.getSubstitutionsCount();
                tsps[0] += subst_p[0]+subst_o[0];
                tvps[0] += subst_p[1]+subst_o[1];
                tsps[1] += subst_p[2]+subst_o[2];
                tvps[1] += subst_p[3]+subst_o[3];
            }
        }
        ks[0] = ks[0]/sites[0];
        ks[1] = ks[1]/sites[1];
        if(ao.bJC_K)
        {
            ks[0] = BasicStatistics.correctJC(ks[0]);
            ks[1] = BasicStatistics.correctJC(ks[1]);
        }
        String[] twopop = new String[2];
        for(int i=0;i<2;i++)
        {
            twopop[i] = String.format(initData.locale, "%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d",
                    sites[i], ks[i], poly[i], tsps[i], tvps[i], Ds[i], tsds[i], tvds[i]);
        }
        return String.format(initData.locale, "%d\t%s\t\t%s\t\t%s\t\t%s",
                             pop.size(),
                             formatResultsString(spb, ao.bJC_pi, ao.bJC_t),
                             formatResultsString(tpb, ao.bJC_pi, ao.bJC_t),
                             twopop[0], twopop[1]);
    }

    /**
     *  Formats the results string of one population.
     *
     *  @param blocks
     *  @param bJC_Pi
     *  @param bJC_Theta
     *  @return
     */
    private String formatResultsString(CodonsBlock[] blocks, boolean bJC_Pi, boolean bJC_Theta)
    {
        float[] sites = CodonsBlock.getSitesCount(blocks);
        float[] pis   = CodonsBlock.getPi(blocks, bJC_Pi);
        float[] thetas= CodonsBlock.getTheta(blocks, bJC_Theta);
        float[] tds   = BasicStatistics.calculateTajD(blocks);
        float[] tdsp  = BasicStatistics.calculateTajDPrime(blocks);
        int[]   poly  = CodonsBlock.getPolymorphismsCount(blocks);
        int[]   si    = CodonsBlock.getSingletonsCount(blocks);
        int[]   tss   = CodonsBlock.getTransitionsCount(blocks);
        int[]   tvs   = CodonsBlock.getTransversionsCount(blocks);
        String[] res  = new String[2];
        for(int i=0;i<2;i++)
        {
            res[i] = String.format(initData.locale, "%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d",
                    sites[i], pis[i], thetas[i], tds[i], tdsp[i], poly[i], si[i], tss[i], tvs[i]);
        }
        return String.format("%s\t\t%s", res[0], res[1]);
    }

    /**
     *  Returns the true if the site is divergent and false otherwise.
     * 
     *  @param ccp
     *  @param cco
     *  @return
     */
    private boolean isSiteDivergent(CodonComposition cc1, CodonComposition cc2)
    {
        Path pp = cc1.getEvolutionaryPath();
        Path po = cc2.getEvolutionaryPath();
        if(pp.getLength()==0 || po.getLength()==0)
            return false;
        Codon[] codons = pp.getPath();
        for(Codon c:codons)
        {
            // If the evolutionary path of the outgroup contains the current codon
            // then the sites are not divergent.
            if(po.indexOf(c)>-1)
                return false;
        }
        return true;
    }

    /**
     *  Returns the number of substitutions:
     *      index | value
     *          0 | D syn
     *          1 | TS syn
     *          2 | TV syn
     *          3 | D nonsyn
     *          4 | TS nonsyn
     *          5 | TV nonsyn
     * @param ccp
     * @param cco
     * @param bUseTerm 
     * @return
     */
    private int[] countSubstitutions(CodonComposition ccp, CodonComposition cco, boolean bUseTerm)
    {
        int[] res = {0, 0, 0, 0, 0, 0};
        Codon[] codons_p = ccp.getEvolutionaryPath().getPath();
        Codon[] codons_o = cco.getEvolutionaryPath().getPath();
        boolean bIsTS = false;
        MAINLOOP:for(int i=0;i<3;i++)
        {
            SiteComposition scp = ccp.getSiteComposition(i);
            // Check, whether the site has divergent bases.
            SiteComposition sco = cco.getSiteComposition(i);
            if(scp.getBaseCount('A')*sco.getBaseCount('A')+
               scp.getBaseCount('C')*sco.getBaseCount('C')+
               scp.getBaseCount('G')*sco.getBaseCount('G')+
               scp.getBaseCount('T')*sco.getBaseCount('T')!=0)
               continue MAINLOOP;
            for(Codon c:codons_p)
            {
                char b = c.getSequence().charAt(i);
                if(scp.getBaseCount(b)>0 && sco.getBaseCount(b)==0)
                {
                    Path p = findBestPath(c, codons_o, bUseTerm);
                    int type = p.getSubstitutionType(b, i);
                    if( (type & Path.MASK_TRANSITION) > 0)
                        bIsTS = true;
                    if( (type & Path.MASK_SYNONYMOUS) > 0)
                    {
                        res[0]++;
                        if( (type & Path.MASK_TRANSITION) > 0)
                            res[1]++;
                        else if( (type & Path.MASK_TRANSVERSION) > 0)
                            res[2]++;
                        continue MAINLOOP;
                    }
                }
            }
            // Since the site is divergent, there MUST be a substitution. At this
            // point it is obvious, that the substitution is nonsynonymous, otherwise
            // this point would not be reached.
            res[3]++;
            if(bIsTS)
                res[4]++;
            else
                res[5]++;
        }
        return res;
    }

    /**
     *  Finds the best path from the divergent codon to any of the codons
     *  of another population.
     *
     *  @param d
     *  @param c
     *  @param bUseTerminal
     *  @return
     */
    private Path findBestPath(Codon d, Codon[] codons, boolean bUseTerminal)
    {
        int nn = Integer.MAX_VALUE;
        Path bp = null;
        for(Codon c:codons)
        {
            Path p = Path.findBestPath(new Codon[]{d, c}, initData.ct, bUseTerminal);
            // There are two general rules, how the best path is determined:
            //      - the new path as many nonsyn. substitutions as the old one, but is shorter
            //      - the new path has less nonsyn. substitutions than the old one
            int tmp = p.getPolymorphismsCount()[1];
            if(bp==null || (tmp==nn && p.getLength()<bp.getLength()) || tmp<nn)
            {
                bp = p;
                nn = tmp;
            }
        }
        return bp;
    }

    /**
     *  Calculates the codons composition of the specified site.
     * 
     *  @param pop
     *  @param pos
     *  @param bUseTerm     whether to use the terminal codons when generating the path
     *  @param bExclTerm    whether to exclude the last terminal codon from the analysis: i.e. the last codon of the CDS will
     *                      not be used if it is a terminal codon
     *  @return
     */
    private CodonComposition calculateCodonComposition(Vector<String> pop, int pos, boolean bUseTerm, boolean bExclTerm)
    {
        if(pop.size()==0)
            return null;
        boolean bIsLast = pos==pop.get(0).length()-3;
        CodonComposition cc = new CodonComposition(initData.ct, bUseTerm);
        int nStrains = pop.size();
        for(int n=0;n<nStrains;n++)
        {
            String strCodon = pop.get(n).substring(pos, pos+3);
            if(strCodon.contains("-") || (bIsLast && bExclTerm && initData.ct.isTerminal(strCodon)))
                return null;
            cc.addCodon(strCodon);
        }
        return cc;
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
             pop='<POP>' out='<OUT>' strlimit='<LIMIT>' si_freq='<FREQ>' 
             jc_pi='<T/F>' jc_t='<T/F>' jc_k='<T/F>' use_term='<T/F>' excl_term='<T/F>'
             exclAll='<T/F>' output='<OUTFILE>'
            
            - pop:          population of interest
            - out:          outgroup
            - strlimit:     maximal number of strains to 
            - si_freq:      singletons cut-off frequency
            - jc_pi:        whether or not to apply Jukes-Cantor correction to pi
            - jc_t:         whether or not to apply Jukes-Cantor correction to theta
            - jc_k:         whether or not to apply Jukes-Cantor correction to K
            - use_term:     whether or not to use terminal codons when generating a path
            - excl_term:    whether or not to exclude terminal codons from analysis
            - exclAll:      whether to exclude the blocks with less than 4 strains from all analyses
            - output:       output filename
        */ 
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "out='(.+)'\\s+"+                   // 2
                                    "strlimit='(\\d*)'\\s+"+            // 3
                                    "si_freq='(\\d*\\.{0,1}\\d*)'\\s+"+ // 4
                                    "jc_pi='([TF])'\\s+"+               // 5
                                    "jc_t='([TF])'\\s+"+                // 6
                                    "jc_k='([TF])'\\s+"+                // 7
                                    "use_term='([TF])'\\s+"+            // 8
                                    "excl_term='([TF])'\\s+"+           // 9
                                    "exclAll='([TF])'\\s+"+             // 10
                                    "output='(.+)'$",                   // 11
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strOut = m.group(2);
            ao.maxstr = (m.group(3).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(3));
            ao.cof = (m.group(4).isEmpty()) ? 1.0f : Float.parseFloat(m.group(4));
            ao.bJC_pi = m.group(5).equalsIgnoreCase("T");
            ao.bJC_t  = m.group(6).equalsIgnoreCase("T");
            ao.bJC_K  = m.group(7).equalsIgnoreCase("T");
            ao.bUseTer= m.group(8).equalsIgnoreCase("T");
            ao.bExclTer = m.group(9).equalsIgnoreCase("T");
            ao.bExclAll = m.group(10).equalsIgnoreCase("T");
            ao.strOutput = m.group(11);
            return ao;
        }
        else
        {
            return (new OptionsDialog(pops)).getOptions();
        }
    }
}
