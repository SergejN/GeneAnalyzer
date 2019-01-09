/*
    File:
        PluginMain.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Calculates the number of substitution pairs in the entire dataset using
        the specified region(s):
            - A<->C
            - A<->G
            - A<->T
            - C<->G
            - C<->T
            - G<->T
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.subst;

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
    public static interface ISubstitutionsCounter
    {
        public String getBasesCountString();
    };

    
    private AInitData initData = null;
    private String strLastErr  = null;

    public ErrorCode Initialize(AInitData initdata)
    {
        initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Substitutions count";
    }

    public String GetName()
    {
        return "Substitutions count analyzer";
    }

    public String GetDescription()
    {
        return "Calculates the number of substitution pairs in the entire dataset using"+
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
        return "pop='<POP>' out='<OUT>' type='<NAME>' maxlen='<MAXLEN>' minlen='<MINLEN>' " +
               "range='<POS1>-<POS2>;<POS3>' any='<T/F>' exclnonffd='<T/F>' nogtag='<T/F>' strlimit='<LIMIT>' " +
               "useterm='<T/F>' exclterm='<T/F>' lenrange='<T/F>' output='<OUTFILE>'";
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
        content.append("Analysis type: Substitution pairs\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Outgroup: %s\n", ao.strPop));
        content.append(String.format("Region type: %s\n", ao.strRegion));
        content.append(String.format("Range: %d - %d\n", ao.iMinlen, ao.iMaxlen));
        content.append(String.format("Codon table: %s\n\n", initData.ct.getName()));
        boolean bCds = ao.strRegion.equalsIgnoreCase("CDS");
        // Write the header.
        if(bCds)
            content.append("Type\t\tA<->C\tA<->G\tA<->T\tC<->G\tC<->T\tG<->T\tTransitions\tTransversions" +
                           "\t\tA<->C\tA<->G\tA<->T\tC<->G\tC<->T\tG<->T\tTransitions\tTransversions\n");
        else
            content.append("Type\t\tSites\tA<->C\tA<->G\tA<->T\tC<->G\tC<->T\tG<->T\tTransitions\tTransversions\tUnknown\n");
        ISubstitutionsCounter[] subst = new ISubstitutionsCounter[3];
        for(int i=0;i<3;i++)
            subst[i] = (bCds) ? new CDSSubstitutionsCounter(initData.locale, initData.ct)
                              : new NoncodingSubstitutionsCounter(initData.locale);
        for(int i=0;i<dataset.getGenesCount();i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            if(nStrains==0)
                continue;
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
            if(pop.size()<1 || out.size()<1)
                continue;
            analyzeGene(pop, out, subst, ao);
        }
        content.append(String.format(initData.locale, "Polymorphic sites (population)\t\t%s\n", subst[0].toString()));
        content.append(String.format(initData.locale, "Polymorphic sites (population+outgroup)\t\t%s\n", subst[1].toString()));
        content.append(String.format(initData.locale, "Divergent sites\t\t%s\n", subst[2].toString()));
        content.append("\n\nBase counts:\t\t");
        content.append((bCds) ? "Syn. A\tSyn. C\tSyn. G\tSyn. T\t\tNonsyn. A\tNonsyn. C\tNonsyn. G\tNonsyn. T\n" : "A\tC\tG\tT\n");
        content.append(String.format(initData.locale, "Polymorphic sites (population)\t\t%s\n", subst[0].getBasesCountString()));
        content.append(String.format(initData.locale, "Polymorphic sites (population+outgroup)\t\t%s\n", subst[1].getBasesCountString()));
        content.append(String.format(initData.locale, "Divergent sites\t\t%s\n", subst[2].getBasesCountString()));
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
            rw.setResults(content.toString(), bCds);
        }
        initData.wd.close();
        if(rw!=null)
            rw.displayResults();
        return ErrorCode.Ok;
    }

    /**
     *  Calculates the number of substitutions in a gene.
     *
     *  @param pop
     *  @param out
     *  @param subst
     *  @param ao
     *  @return
     */
    private void analyzeGene(Vector<StrainEntry> pop, Vector<StrainEntry> out, ISubstitutionsCounter[] subst, AnalysisOptions ao)
    {
        String[] popseq = new String[pop.size()];   // Pop. sequences
        String[] outseq = new String[out.size()];   // Out. sequences
        // If the type is CDS or FFD, extract the complete coding sequence.
        if(ao.strRegion.equalsIgnoreCase("CDS") || ao.strRegion.equalsIgnoreCase("FFD"))
        {
            for(int i=0;i<pop.size();i++)
            {
                String strCds = pop.get(i).getCodingSequence();
                if(strCds==null || strCds.isEmpty())
                    return;
                popseq[i] = strCds;
            }
            for(int i=0;i<out.size();i++)
            {
                String strCds = out.get(i).getCodingSequence();
                if(strCds==null || strCds.isEmpty())
                    return;
                outseq[i] = strCds;
            }
            if(ao.strRegion.equalsIgnoreCase("CDS"))
                countSubstitutionsCDS(popseq, outseq, subst, ao);
            else
                countSubstitutionsFFD(popseq, outseq, subst, ao);
        }
        else // Non-coding regions.
        {
            int nRegs = pop.get(0).getRegionsCount();
            for(int i=0;i<nRegs;i++)
            {
                // Check the region type and the intron length.
                if(pop.get(0).getRegion(i).hasType(ao.strRegion) && (isRegionLengthOk(pop, out, ao, i)) )
                {
                    // Distinguish between two options:
                    //  1. Complete region sequence
                    //  2. Specific sites
                    if(ao.sites==null) // Option 1.
                    {
                        popseq = generateSequenceSet(popseq, pop, i, ao.bNoGtag);
                        outseq = generateSequenceSet(outseq, out, i, ao.bNoGtag);
                    }
                    else    // Option 2.
                        generateSequenceSet(popseq, outseq, pop, out, i, ao.sites);
                }
            }
            // If no region(s) found, return.
            if(popseq[0]==null)
                return;
            countSubstitutionsNoncoding(popseq, outseq, subst, ao);
        }
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
     *  Iterates through both population and outgroup and returns true only if
     *  the minimal and the maxinal length are not beyond the length threshold.
     *
     *  @param pop
     *  @param out
     *  @param ao
     *  @param index
     *  @return
     */
    private boolean isRegionLengthOk(Vector<StrainEntry> pop, Vector<StrainEntry> out, AnalysisOptions ao, int index)
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
    *                                  CODING SEQUENCE                                       *
    *****************************************************************************************/
    private void countSubstitutionsCDS(String[] pop, String[] out, ISubstitutionsCounter[] subst, AnalysisOptions ao)
    {
        int length = pop[0].length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        if(length<3)
            return;
        for(int pos=0;pos<length;pos+=3)
        {
            CodonComposition ccp = calculateCodonComposition(pop, pos, ao.bUseTerm, ao.bExclTerm);
            if(ccp==null)
                continue;
            // Calculate the number of substitutions for the population of interest only.
            ((CDSSubstitutionsCounter)subst[0]).addSubstitutions(ccp, null, ao.bUseTerm);
            CodonComposition cco = calculateCodonComposition(out, pos, ao.bUseTerm, ao.bExclTerm);
            if(cco!=null)
            {
                int type = CodonComposition.getSiteType(ccp, cco);
                // If the site is monomorphic between the populations add it to both counters.
                if(type==CodonComposition.ST_MONOMORPHIC)
                {
                    ((CDSSubstitutionsCounter)subst[1]).addSubstitutions(ccp, null, ao.bUseTerm);
                    ((CDSSubstitutionsCounter)subst[2]).addSubstitutions(ccp, cco, ao.bUseTerm);
                }
                // If the site is divergent, add it to the counter of divergent sites.
                else if( (type & SiteComposition.ST_DIVERGENT)>0 )
                {
                    ((CDSSubstitutionsCounter)subst[2]).addSubstitutions(ccp, cco, ao.bUseTerm);
                }
                // If the site is polymorphic in the population, but not divergent, add it to the
                // first counter only.
                else if( (type & SiteComposition.ST_POLYMORPHIC_FIRST)>0 )
                {
                    ((CDSSubstitutionsCounter)subst[1]).addSubstitutions(ccp, null, ao.bUseTerm);
                }
            }
        }
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

    /*****************************************************************************************
    *                            FOUR-FOLD DEGENERATE SITES                                  *
    *****************************************************************************************/
    private void countSubstitutionsFFD(String[] pop, String[] out, ISubstitutionsCounter[] subst, AnalysisOptions ao)
    {
        int length = pop[0].length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        if(length<3)
            return;
        for(int pos=0;pos<length;pos+=3)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos, ao);
            if(scp==null)
                continue;
            // Calculate the number of substitutions for the population of interest only.
            ((NoncodingSubstitutionsCounter)subst[0]).addSubstitution(scp, null);
            SiteComposition sco = calculateSiteComposition(out, pos, ao);
            if(sco!=null)
            {
                // If the user selected to not to use nonsyn. FFD sites and the
                // codons in the population of interest and in the outgroup are
                // not synonymous, skip the site.
                if(!ao.bUseAny)
                {
                    String strCodon1 = pop[0].substring(pos, pos+3);
                    String strCodon2 = out[0].substring(pos, pos+3);
                    if(!initData.ct.areSynonymous(strCodon1, strCodon2))
                        continue;
                }
                int type = SiteComposition.getSiteType(scp, sco);
                // If the site is monomorphic between the populations add it to both counters.
                if(type==SiteComposition.ST_MONOMORPHIC)
                {
                    ((NoncodingSubstitutionsCounter)subst[1]).addSubstitution(scp, null);
                    ((NoncodingSubstitutionsCounter)subst[2]).addSubstitution(scp, sco);
                }
                // If the site is divergent, add it to the counter of divergent sites.
                else if( (type & SiteComposition.ST_DIVERGENT)>0 )
                {
                    ((NoncodingSubstitutionsCounter)subst[2]).addSubstitution(scp, sco);
                }
                // If the site is polymorphic in the population, but not divergent, add it to the
                // first counter only.
                else if( (type & SiteComposition.ST_POLYMORPHIC_FIRST)>0 )
                {
                    ((NoncodingSubstitutionsCounter)subst[1]).addSubstitution(scp, null);
                }
            }
        }
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
            if(!ao.bUseAny && !initData.ct.areSynonymous(codon, strRef))
                return null;
            b = true;
            sc.addBase(codon.charAt(2));
        }
        return (b) ? sc : null;
    }

    /*****************************************************************************************
    *                                  NONCODING REGIONS                                     *
    *****************************************************************************************/
    private void countSubstitutionsNoncoding(String[] pop, String[] out, ISubstitutionsCounter[] subst, AnalysisOptions ao)
    {
        int length = pop[0].length();
        for(int pos=0;pos<length;pos++)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos);
            if(scp==null)
                continue;
            // Calculate the number of substitutions for the population of interest only.
            ((NoncodingSubstitutionsCounter)subst[0]).addSubstitution(scp, null);
            SiteComposition sco = calculateSiteComposition(out, pos);
            if(sco!=null)
            {
                int type = SiteComposition.getSiteType(scp, sco);
                // If the site is monomorphic between the populations add it to both counters.
                if(type==SiteComposition.ST_MONOMORPHIC)
                {
                    ((NoncodingSubstitutionsCounter)subst[1]).addSubstitution(scp, null);
                    ((NoncodingSubstitutionsCounter)subst[2]).addSubstitution(scp, sco);
                }
                // If the site is divergent, add it to the counter of divergent sites.
                else if( (type & SiteComposition.ST_DIVERGENT)>0 )
                {
                    ((NoncodingSubstitutionsCounter)subst[2]).addSubstitution(scp, sco);
                }
                // If the site is polymorphic in the population, but not divergent, add it to the
                // first counter only.
                else if( (type & SiteComposition.ST_POLYMORPHIC_FIRST)>0 )
                {
                    ((NoncodingSubstitutionsCounter)subst[1]).addSubstitution(scp, null);
                }
            }
        }
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
            pop='<POP>' out='<OUT>' type='<NAME>' maxlen='<MAXLEN>' minlen='<MINLEN>'
            range='<POS1-POS2; POS3>' any='<T/F>' exclnonffd='<T/F>' nogtag='<T/F>' strlimit='<LIMIT>'
            useterm='<T/F>' exclterm='<T/F>' lenrange='<T/F>' output='<OUTFILE>'

            - pop:          population of interest
            - out:          outgroup
            - type:         region type
            - maxlen:       maximal intron length. Use 0 for undefined (see Remarks)
            - minlen:       minimal intron length. Use 0 for undefined (see Remarks)
            - range:        specific sites to use. See Remarks section for details.
            - any:          whether or not to use any (i.e. both syn. and nonsyn.) FFD sites.
            - exclnonffd:   whether to exclude the strains with non-FFD codons.
            - nogtag:       whether to exclude the first GT and the last AG from the intron sequence
            - strlimit:     maximal number of strains to use                      
            - useterm:      whether to use terminal codons when estimating the number of sites
            - exclterm:     whether to exclude the terminal codon of the gene from the analysis
            - lenrange:     whether the minimal and maximal length should fall into the specified range
            - output:       output filename

            Remarks:
            If type is CDS or FFD, the parameters minlen, maxlen and range are ignored.
        */
        // Parse the parameter string.
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "out='(.+)'\\s+"+                   // 2
                                    "type='(.+)'\\s+" +                 // 3
                                    "maxlen='(\\d+)'\\s+"+              // 4
                                    "minlen='(\\d+)'\\s+"+              // 5
                                    "range='([-,;\\d\\s]*)'\\s+"+       // 6
                                    "any='([TF])'\\s+"+                 // 7
                                    "exclnonffd='([TF])'\\s+"+          // 8
                                    "nogtag='([TF])'\\s+"+              // 9
                                    "strlimit='(\\d*)'\\s+"+            // 10
                                    "useterm='([TF])'\\s+"+             // 11
                                    "exclterm='([TF])'\\s+"+            // 12
                                    "lenragnge='([TF])'\\s+"+           // 13
                                    "output='(.+)'$",                   // 14
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strOut = m.group(2);
            ao.strRegion = m.group(3);
            ao.iMaxlen = Integer.parseInt(m.group(4));
            if(ao.iMaxlen==0)
                ao.iMaxlen = Integer.MAX_VALUE;
            ao.iMinlen = Integer.parseInt(m.group(5));
            if(ao.iMaxlen<ao.iMinlen)
            {
                int tmp = ao.iMaxlen;
                ao.iMaxlen = ao.iMinlen;
                ao.iMinlen = tmp;
            }
            ao.maxstr = (m.group(10).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(10));
            ao.strOutput = m.group(14);
            // Check whether the region type is FFD.
            if(ao.strRegion.equalsIgnoreCase("FFD"))
            {
                ao.iMaxlen = Integer.MAX_VALUE;
                ao.iMinlen = 0;
            }
            // Range.
            ao.sites = OptionsDialog.extractSites(m.group(6));
            ao.bUseAny = m.group(7).equalsIgnoreCase("T");
            ao.bNonFfd = m.group(8).equalsIgnoreCase("T");
            ao.bUseTerm = m.group(11).equalsIgnoreCase("T");
            ao.bExclTerm = m.group(12).equalsIgnoreCase("T");
            ao.bLenRange = m.group(13).equalsIgnoreCase("T");
            ao.bNoGtag = (ao.sites==null) ? m.group(9).equalsIgnoreCase("T") : false;
            return ao;
        }
        else
        {
            return (new OptionsDialog(pops, regs)).getOptions();
        }
    }
}