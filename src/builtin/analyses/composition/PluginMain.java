/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Calculates the number of the bases in the entire dataset using
        the specified region(s).
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.composition;

import algorithms.CodonComposition;
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
    private AInitData initData = null;
    private String strLastErr  = null;

    public ErrorCode Initialize(AInitData initdata)
    {
        initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Bases composition";
    }

    public String GetName()
    {
        return "Bases composition analyzer";
    }

    public String GetDescription()
    {
        return "Calculates the number of the bases in the entire dataset using" +
               " the specified region(s)";
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
        return "pop='<POP>' type='<NAME>' maxlen='<MAXLEN>' minlen='<MINLEN>' " +
               "range='<POS1>-<POS2>;<POS3>' any='<T/F>' exclnonffd='<T/F>' nogtag='<T/F>' strlimit='<LIMIT>' " +
               "constsize='<T/F>' useterm='<T/F>' exclterm='<T/F>' lenrange='<T/F>' output='<OUTFILE>'";
    }

    public String GetLastError()
    {
        return strLastErr;
    }

    @Override
    public ErrorCode AnalyzeDataset(Dataset dataset, String params)
    {
        AnalysisOptions ao = getOptions(params, dataset.listPopulations(),
                                            dataset.listRegionsNames());
        if(ao==null)
        {
            strLastErr = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Analysis);
        // Perform analysis.
        StringBuffer content = new StringBuffer();
        content.append("Analysis type: Bases composition\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Region type: %s\n", ao.strRegion));
        content.append(String.format("Range: %d - %d\n", ao.iMinlen, ao.iMaxlen));
        content.append(String.format("Codon table: %s\n\n", initData.ct.getName()));
        boolean bCDS = ao.strRegion.equalsIgnoreCase("CDS");
        // Write the header.
        if(bCDS)
            content.append("Syn.sites\tA\tC\tG\tT\t\tNonsyn.sites\tA\tC\tG\tT\n");
        else
            content.append("Sites\tA\tC\tG\tT\tN\tX\n");
        float[] nums = (bCDS) ? new float[10] : new float[7];
        for(int i=0;i<dataset.getGenesCount();i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            if(nStrains==0)
                continue;
            Vector<StrainEntry> pop = new Vector<StrainEntry>();
            for(int n=0;n<nStrains;n++)
            {
                StrainEntry se = ge.getStrainEntry(n);
                if(se.belongsToPopulation(ao.strPop) && pop.size()<ao.maxstr)
                    pop.add(se);
            }
            if(pop.size()<1)
                continue;
            float[] tmp = analyzeGene(pop, ao);
            if(tmp!=null)
            {
                for(int n=0;n<tmp.length;n++)
                    nums[n] += tmp[n];
            }
        }
        if(bCDS)
        {
            content.append(String.format(initData.locale, "%f\t%.2f\t%.2f\t%.2f\t%.2f\t\t%f\t%.2f\t%.2f\t%.2f\t%.2f",
                     nums[0], nums[1], nums[2], nums[3], nums[4],
                     nums[5], nums[6], nums[7], nums[8], nums[9]));
        }
        else
        {
            content.append(String.format(initData.locale, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
                    (int)nums[0], nums[1], nums[2], nums[3], nums[4], nums[5], nums[6]));
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
            rw.setResults(content.toString(), bCDS);
        }
        initData.wd.close();
        if(rw!=null)
            rw.displayResults();
        return ErrorCode.Ok;
    }

    private float[] analyzeGene(Vector<StrainEntry> pop, AnalysisOptions ao)
    {
        String[] popseq = new String[pop.size()];   // Pop. sequences
        // If the type is FFD, extract the complete coding sequence.
        if(ao.strRegion.equalsIgnoreCase("CDS") || ao.strRegion.equalsIgnoreCase("FFD"))
        {
            for(int i=0;i<pop.size();i++)
            {
                String strCds = pop.get(i).getCodingSequence();
                if(strCds==null || strCds.isEmpty())
                    return null;
                popseq[i] = strCds;
            }
            if(ao.strRegion.equalsIgnoreCase("CDS"))
                return countBasesCDS(popseq, ao);
            else
                return countBasesFFD(popseq, ao);
        }
        else // Non-coding regions.
        {
            int nRegs = pop.get(0).getRegionsCount();
            for(int i=0;i<nRegs;i++)
            {
                // Check the region type and the length.
                if(pop.get(0).getRegion(i).hasType(ao.strRegion) && (isRegionLengthOk(pop, ao, i)) )
                {
                    // Distinguish between two options:
                    //  1. Complete region sequence
                    //  2. Specific sites
                    if(ao.sites==null) // Option 1.
                        popseq = generateSequenceSet(popseq, pop, i, ao.bNoGtag);
                    else    // Option 2.
                        generateSequenceSet(popseq, pop, i, ao.sites);
                }
            }
            // If no region(s) found, return null.
            if(popseq[0]==null)
                return null;
            return countBasesNoncoding(popseq, ao);
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
    private String[] generateSequenceSet(String[] old, Vector<StrainEntry> pop, int iRegion, boolean bNoGtag)
    {
        int nStrains = pop.size();
        String[] seqs = new String[nStrains];
        for(int n=0;n<nStrains;n++)
        {
            String seq = pop.get(n).getRegion(iRegion).getSequence();
            if(bNoGtag)
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
    private void generateSequenceSet(String[] old_p, Vector<StrainEntry> pop, int iRegion, int[] sites)
    {
        int nPopSize = pop.size();
        String[] seqs = new String[nPopSize];
        for(int i=0;i<nPopSize;i++)
            seqs[i] = pop.get(i).getRegion(iRegion).getSequence();
        // Generate ungapped sequences.
        seqs = SequenceRoutines.getUngappedSequences(seqs);
        // Extract the sites.
        seqs = SequenceRoutines.extractAlignedSites(seqs, sites);
        for(int i=0;i<nPopSize;i++)
            old_p[i] = (old_p[i]==null) ? seqs[i] : old_p[i]+seqs[i];
    }

    /**
     *  Iterates through the population returns true only if the minimal and the maxinal
     *  length are not beyond the length threshold.
     *
     *  @param pop
     *  @param ao
     *  @param index
     *  @return
     */
    private boolean isRegionLengthOk(Vector<StrainEntry> pop, AnalysisOptions ao, int index)
    {
        int nPopSize = pop.size();
        String[] seqs = new String[nPopSize];
        for(int i=0;i<nPopSize;i++)
            seqs[i] = pop.get(i).getRegion(index).getSequence();
        if(ao.bLenRange)
        {
            int[] tmp = SequenceRoutines.getMinMaxUngappedLength(seqs);
            return ( (tmp[1]<=ao.iMaxlen) && (tmp[0]>=ao.iMinlen) );
        }
        else
        {
            seqs = SequenceRoutines.getUngappedSequences(seqs);
            int l = seqs[0].length();
            return ( (l<=ao.iMaxlen) && (l>=ao.iMinlen) );
        }
    }

    /*****************************************************************************************
    *                                  ANALYSES                                              *
    *****************************************************************************************/
    private float[] countBasesCDS(String[] pop, AnalysisOptions ao)
    {
        int length = pop[0].length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        if(length<3)
            return null;
        float[] res = new float[10];
        for(int pos=0;pos<length;pos+=3)
        {
            CodonComposition ccp = calculateCodonComposition(pop, pos, ao.bUseTerm, ao.bExclTerm);
            if(ccp==null)
                continue;
            float[] tmp = ccp.getBaseFrequencies(initData.ct, ao.bUseTerm);
            for(int i=0;i<10;i++)
                res[i]+=tmp[i];
        }
        return res;
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
    private CodonComposition calculateCodonComposition(String[] pop, int pos, boolean bUseTerm, boolean bExclTerm)
    {
        boolean bIsLast = pos==pop[0].length()-3;
        CodonComposition cc = new CodonComposition(initData.ct, bUseTerm);
        int nStrains = pop.length;
        for(int n=0;n<nStrains;n++)
        {
            String strCodon = pop[n].substring(pos, pos+3);
            if(strCodon.contains("-") || (bIsLast && bExclTerm && initData.ct.isTerminal(strCodon)))
                return null;
            cc.addCodon(strCodon);
        }
        return cc;
    }

    private float[] countBasesFFD(String[] pop, AnalysisOptions ao)
    {
        int length = pop[0].length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        if(length<3)
            return null;
        float[] res = (ao.bSizeConst) ? new float[7] : new float[5];
        for(int pos=0;pos<length;pos+=3)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos, ao);
            if(scp==null)
                continue;
            res[0]++;
            float[] tmp = scp.getBaseFrequencies(ao.bSizeConst);
            for(int i=0;i<tmp.length;i++)
                res[i+1] += tmp[i];
        }
        return res;
    }

    private SiteComposition calculateSiteComposition(String[] pop, int pos, AnalysisOptions ao)
    {
        SiteComposition sc = new SiteComposition();
        String strRef = pop[0].substring(pos, pos+3);
        boolean b = false; // Flag, specifying, whether any sequence has a FFD codon at the specified site.
        for(int n=0;n<pop.length;n++)
        {
            String codon = pop[n].substring(pos, pos+3);
            // 1. If there is a gap in any strain return null.
            if(codon.contains("-"))
                return null;
            // 2. Check whether the codon is a FFD codon.
            // If it is not and the option "Ignore strains with non-FFD sites" is not selected,
            // return null to exclude the entire site. Otherwise just skip the single strain.
            if(initData.ct.getFoldFamily(codon)!=4)
            {
                if(ao.bNonFfd)
                    continue;
                else
                    return null;
            }
            // 3. Check whether the codon is syn. to the reference codon.
            // If they are not synonymous, then either reject the whole site
            // and return null or proceed with all sequences.
            if(!initData.ct.areSynonymous(codon, strRef) && !ao.bUseAny)
                return null;
            b = true;
            sc.addBase(codon.charAt(2));
        }
        return (b) ? sc : null;
    }

    private float[] countBasesNoncoding(String[] pop, AnalysisOptions ao)
    {
        int length = pop[0].length();
        float[] res = (ao.bSizeConst) ? new float[7] : new float[5];
        for(int pos=0;pos<length;pos++)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos);
            if(scp==null)
                continue;
            res[0]++;
            float[] tmp = scp.getBaseFrequencies(ao.bSizeConst);
            for(int i=0;i<tmp.length;i++)
                res[i+1] += tmp[i];
        }
        return res;
    }

    private SiteComposition calculateSiteComposition(String[] pop, int pos)
    {
        SiteComposition sc = new SiteComposition();
        for(String s:pop)
        {
            char base = s.charAt(pos);
            if(base=='-')
                return null;
            sc.addBase(base);
        }
        return sc;
    }

    /**
     *  Returns the analysis options.
     *
     *  @param strParams
     *  @param pops
     *  @param regs
     *  @return
     */
    private AnalysisOptions getOptions(String strParams, String[] pops, String[] regs)
    {
        // If the parameters are not specified or empty, show the dialog.
        if(strParams==null || strParams.isEmpty())
        {
            return (new OptionsDialog(pops, regs)).getOptions();
        }
        // If the parameter string is not empty, parse it.
        /** The parameters line should have the following format:
            pop='<POP>' type='<NAME>' maxlen='<MAXLEN>' minlen='<MINLEN>'
            range='<POS1-POS2; POS3>' any='<T/F>' exclnonffd='<T/F>' strlimit='<LIMIT>' nogtag='<T/F>'
            combine='<T/F>' constsize='<T/F>' useterm='<T/F>' exclterm='<T/F>'  lenrange='<T/F>' output='<OUTFILE>'

            - pop:          population of interest
            - type:         region type
            - maxlen:       maximal intron length. Use 0 for undefined (see Remarks)
            - minlen:       minimal intron length. Use 0 for undefined (see Remarks)
            - range:        specific sites to use. See Remarks section for details.
            - any:          whether or not to use any (i.e. both syn. and nonsyn.) FFD sites.
            - exclnonffd:   whether to exclude the strains with non-FFD codons.
            - nogtag:       whether to exclude the first GT and the last AG from the intron sequence
            - strlimit:     maximal number of strains to use
            - constsize:    whether to assume the population size to be constant when calculating the frequencies
            - useterm:      whether to use terminal codons when estimating the number of sites
            - exclterm:     whether to exclude the terminal codon of the gene from the analysis
            - lenrange:     whether the minimal and maximal length should fall into the specified range
            - output:       output filename

            Remarks:
            If type is CDS or FFD, the parameters minlen, maxlen and range are ignored.
        */
        // Parse the parameter string.
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "type='(.+)'\\s+" +                 // 2
                                    "maxlen='(\\d+)'\\s+"+              // 3
                                    "minlen='(\\d+)'\\s+"+              // 4
                                    "range='([-,;\\d\\s]*)'\\s+"+       // 5
                                    "any='([TF])'\\s+"+                 // 6
                                    "exclnonffd='([TF])'\\s+"+          // 7
                                    "nogtag='([TF])'\\s+"+              // 8
                                    "strlimit='(\\d*)'\\s+"+            // 9
                                    "constsize='([TF])'\\s+"+           // 10
                                    "useterm='([TF])'\\s+"+             // 11
                                    "exclterm='([TF])'\\s+"+            // 12
                                    "lenrange='([TF])'\\s+"+            // 13
                                    "output='(.+)'$",                   // 14
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strRegion = m.group(2);
            ao.iMaxlen = Integer.parseInt(m.group(3));
            if(ao.iMaxlen==0)
                ao.iMaxlen = Integer.MAX_VALUE;
            ao.iMinlen = Integer.parseInt(m.group(4));
            if(ao.iMaxlen<ao.iMinlen)
            {
                int tmp = ao.iMaxlen;
                ao.iMaxlen = ao.iMinlen;
                ao.iMinlen = tmp;
            }
            ao.maxstr = (m.group(9).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(9));
            ao.bSizeConst = m.group(10).equalsIgnoreCase("T");
            ao.bUseTerm = m.group(11).equalsIgnoreCase("T");
            ao.bExclTerm = m.group(12).equalsIgnoreCase("T");
            ao.bLenRange = m.group(12).equalsIgnoreCase("T");
            ao.strOutput = m.group(14);
            // Check whether the region type is CDS or FFD.
            if(ao.strRegion.equalsIgnoreCase("CDS") || ao.strRegion.equalsIgnoreCase("FFD"))
            {
                ao.iMaxlen = Integer.MAX_VALUE;
                ao.iMinlen = 0;
            }
            // Range.
            ao.sites = OptionsDialog.extractSites(m.group(5));
            ao.bUseAny = m.group(6).equalsIgnoreCase("T");
            ao.bNonFfd = m.group(7).equalsIgnoreCase("T");
            ao.bNoGtag = (ao.sites==null) ? m.group(8).equalsIgnoreCase("T") : false;
            return ao;
        }
        else
        {
            return (new OptionsDialog(pops, regs)).getOptions();
        }
    }
}
