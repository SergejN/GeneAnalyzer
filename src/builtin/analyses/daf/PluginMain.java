/*
    File:
        PluginMain.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Estimates the frequencies of derived mutation in the specified region.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.daf;

import algorithms.SequenceRoutines;
import algorithms.SiteComposition;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.StrainEntry;
import bio.gene.dna.ICodonTable.TYPE;
import gui.IWaitDialog;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetAnalyzer;


public class PluginMain extends ADatasetAnalyzer
{
    /**
     *  Represents the derived mutation in the coding region.
     */
    private class CDSSite
    {
        char anc;
        char derived;
        int  nCount;
    };

    private String strErr = null;
    private Vector<Float> freqs = null;
    private int nMutations = 0;
    private AInitData initData = null;
    

    public ErrorCode Initialize(AInitData initdata)
    {
        initData = initdata;
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Derived allele frequency";
    }

    public String GetName()
    {
        return "Derived allele frequency analyzer";
    }

    public String GetDescription()
    {
        return "Calculates the frequencies of derived alleles";
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
               "range='<POS1>-<POS2>;<POS3>' any='<T/F>' exclnonffd='<T/F>' strlimit='<LIMIT>' " +
               "combine='<T/F>' freqs='<T/F>' constsize='<T/F>' lenrange='<T/F>' output='<OUTFILE>'";
    }

    public String GetLastError()
    {
        return strErr;
    }

    @Override
    public ErrorCode AnalyzeDataset(Dataset dataset, String params)
    {
        AnalysisOptions ao = getOptions(params, dataset.listPopulations(), 
                                            dataset.listRegionsNames());
        if(ao==null)
        {
            strErr = "Analysis cancelled by user";
            return ErrorCode.CancelledByUser;
        }
        initData.wd.show(IWaitDialog.TYPE.Analysis);
        nMutations = 0;
        // Perform analysis.
        StringBuffer content = new StringBuffer();
        content.append("Analysis type: Derived allele frequency\n");
        content.append(String.format("Population of interest: %s\n", ao.strPop));
        content.append(String.format("Outgroup: %s\n", ao.strOut));
        content.append(String.format("Region type: %s\n", ao.strRegion));
        content.append(String.format("Range: %d - %d\n", ao.iMinlen, ao.iMaxlen));
        content.append(String.format("Codon table: %s\n\n", initData.ct.getName()));
        // Write the header.
        content.append("Gene\t\tSamp.size\tSites\tP\tDAF\t\tSubstitutions\n");
        // Single frequencies.
        freqs = new Vector<Float>();
        // Iterate through the data set and analyze the genes.
        for(int i=0;i<dataset.getGenesCount();i++)
        {
            GeneEntry ge = dataset.getGeneEntry(i);
            int nStrains = ge.getStrainsCount();
            if(nStrains==0)
                continue;
            // Iterate through the gene entry to select the relevant strains.
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
            // If the number of strains in the population of interest is less
            // than two or the number of strains in the outgroup is less than
            // one, the analysis of derived allele frequencies cannot be performed.
            if(pop.size()<2 || out.size()<1)
                continue;
            // Once, the vectors of data are ready, analysis can be performed.
            // Add results to the file content.
            String[] lines = analyzeGene(pop, out, ao);
            if(lines!=null)
            {
                content.append(ge.getCommonName());
                for(String s:lines)
                    content.append("\t\t"+s+"\n");
            }           
        }
        // Add single frequencies.
        if(ao.bListFreqs)
        {
            content.append("\n\nSingle frequencies");
            for(Float f:freqs)
                content.append(String.format(initData.locale, "\n%.3f", f));
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
            rw.setResults(content.toString(), nMutations);
        }
        initData.wd.close();
        if(rw!=null)
            rw.displayResults();
        return ErrorCode.Ok;
    }

    private String[] analyzeGene(Vector<StrainEntry> pop, Vector<StrainEntry> out, AnalysisOptions ao)
    {
        Vector<String> res = new Vector<String>(); // Result lines
        String[] popseq = new String[pop.size()];   // Pop. sequences
        String[] outseq = new String[out.size()];   // Out. sequences
        // If the type is CDS or FFD, extract the complete coding sequence.
        if(ao.strRegion.equalsIgnoreCase("CDS") || ao.strRegion.equalsIgnoreCase("FFD"))
        {
            for(int i=0;i<pop.size();i++)
            {
                String strCds = pop.get(i).getCodingSequence();
                if(strCds==null || strCds.isEmpty())
                    return null;
                popseq[i] = strCds;
            }
            for(int i=0;i<out.size();i++)
            {
                String strCds = out.get(i).getCodingSequence();
                if(strCds==null || strCds.isEmpty())
                    return null;
                outseq[i] = strCds;
            }
            if(ao.strRegion.equalsIgnoreCase("CDS"))
            {
                String s = calculateStatisticsSynNonsyn(popseq, outseq, ao);
                if(s!=null)
                    res.add(s);
            }
            else
            {
                String s = calculateStatisticsFFD(popseq, outseq, ao);
                if(s!=null)
                    res.add(s);
            }
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
                        popseq = generateSequenceSet(popseq, pop, i);
                        outseq = generateSequenceSet(outseq, out, i);
                    }
                    else    // Option 2.
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
                        res.add(calculateStatisticsNoncoding(popseq, outseq, ao));
                        popseq = new String[pop.size()];
                        outseq = new String[out.size()];
                    }
                }
            }
            if(ao.bCombine && popseq[0]!=null)
                res.add(calculateStatisticsNoncoding(popseq, outseq, ao));
        }
        return (res.size()!=0) ? res.toArray(new String[1]) : null;
    }

    /**********************************************************************************************
    *                           SYNONYMOUS/NONSYNONYMOUS SITES                                    *
    **********************************************************************************************/
    /**
     *  Calculates the derived allele frequencies using synonymous/nonsyn. sites.
     *
     *  @param pop
     *  @param out
     *  @param ao 
     *  @return
     */
    private String calculateStatisticsSynNonsyn(String[] pop, String[] out, AnalysisOptions ao)
    {
        int length = pop[0].length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        if(length<3)
            return null;
        int nSites = 0;
        int nDerived = 0;
        float freq = 0.0f;
        StringBuffer sb = new StringBuffer();
        for(int pos=0;pos<length;pos+=3)
        {
            char[] scp = calculateAAComposition(pop, pos);
            if(scp==null)
                continue;
            char[] sco = calculateAAComposition(out, pos);
            if(sco==null)
                continue;
            nSites++;
            CDSSite[] derived = findDerivedMutations(scp, sco);
            // Divergent or monomorphic site.
            if(derived==null)
                continue;
            else // Polymorphic site.
            {
                for(CDSSite site:derived)
                {
                    float f = (ao.bSizeConst) ? (float)site.nCount/(float)pop.length
                                              : (float)site.nCount/(float)scp.length;
                    freqs.add(f);
                    sb.append(String.format(initData.locale, "%.3f\t(%c->%c)\t",
                                            f, site.anc, site.derived));
                    freq += f;
                    nDerived++;
                }
            }
        }
        // Add average frequency.
        if(nDerived>0)
            freq = freq/(float)nDerived;
        else
            freq = Float.NaN;
        if(nDerived>nMutations)
            nMutations = nDerived;
        sb.insert(0, String.format(initData.locale, "%d\t%d\t%d\t%.3f\t\t",
                                   pop.length, nSites, nDerived, freq));
        return sb.toString();
    }

    /**
     *  Returns the amino-acid composition of the specified site. The site should
     *  be the first site of a codon.
     *  If there is a gap in any sequence in the codon specified by its first site,
     *  the method returns null. If the codon sequence contains X or N, the codon
     *  is skipped, and, thus, the length of the resulting array can be less
     *  than the length of the sequences array.
     *  If there is no valid codon at the specified site, the method also returns null.
     * 
     *  @param seqs
     *  @param pos
     *  @return
     */
    private char[] calculateAAComposition(String[] seqs, int pos)
    {
        char[] aas = new char[seqs.length];
        int n = 0;
        for(int i=0;i<aas.length;i++)
        {
            String strCodon = seqs[i].substring(pos, pos+3);
            if(strCodon.indexOf('-')>-1)
                return null;
            if(strCodon.indexOf('X')>-1 || strCodon.indexOf('N')>-1)
                continue;
            aas[n] = initData.ct.getAminoAcid(strCodon, TYPE.OneLetterCode).charAt(0);
            n++;
        }
        if(n==0)
            return null;
        return Arrays.copyOf(aas, n);
    }


    /**
     *  Finds the derived mutations and counts how often the derived mutation appears
     *  in the population.
     *
     *  If the site is divergent or monomorphic, the method returns null.
     *
     *  @param scp
     *  @param sco
     *  @return
     */
    private CDSSite[] findDerivedMutations(char[] scp, char[] sco)
    {
        char anc  = 0;   // Currently best ancestral amino acid.
        int  nAnc = 0;   // Currently best ancestral amino acid count:
                         // # of matches in the population of interest.
        for(char c:sco)
        {
            int n = 0;
            if(c==anc)
                continue;
            // For each amino acid from the outgroup count the number of times it matches
            // an amino acid from the population of interest.
            for(char b:scp)
                if(c==b)
                    n++;
            if(n>nAnc)
            {
                nAnc=n;
                anc=c;
            }
        }
        // If there is no ancestral amino acid, then the site is divergent.
        if(anc==0)
            return null;
        // Otherwise calculate the frequency of the derived mutations.
        HashMap<Character, Integer> tmp = new HashMap<Character, Integer>();
        for(char b:scp)
        {
            if(b==anc)
                continue;
            Integer i = tmp.get(b);
            if(i==null)
                i = new Integer(1);
            else
                i++;
            tmp.put(b, i);
        }
        if(tmp.size()==0)
            return null;
        CDSSite[] sites = new CDSSite[tmp.size()];
        int i = 0;
        Iterator<Character> iterator = tmp.keySet().iterator();
        while(iterator.hasNext())
        {
            Character aa = iterator.next();
            CDSSite site = new CDSSite();
            site.anc = anc;
            site.derived = aa;
            site.nCount = tmp.get(aa);
            sites[i] = site;
            i++;
        }
        return sites;
    }

    /**********************************************************************************************
    *                           FOUR-FOLD DEGENERATE SITES                                        *
    **********************************************************************************************/
    /**
     *  Calculates the derived allele frequencies of four-fold degenerate sites.
     *
     *  @param pop
     *  @param out
     *  @param ao
     *  @return
     */
    private String calculateStatisticsFFD(String[] pop, String[] out, AnalysisOptions ao)
    {
        int length = pop[0].length();
        // For the case that not the complete gene sequence is present
        // use only the part of it, so that the ORF is intact.
        if(length%3!=0)
            length = (length/3)*3;
        if(length<3)
            return null;
        int nSites = 0;
        int nDerived = 0;
        float freq = 0.0f;
        StringBuffer sb = new StringBuffer();
        for(int pos=0;pos<length;pos+=3)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos, ao);
            // If the site is divergent/monomorphic or if there are less than 2
            // valid bases at the site, then it cannot be analyzed.
            if(scp==null || scp.getValidBasesCount()<2)
                continue;
            SiteComposition sco = calculateSiteComposition(out, pos, ao);
            if(sco==null)
                continue;
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
            nSites++;
            CDSSite[] derived = findDerivedMutations(scp, sco);
            // Divergent or monomorphic site.
            if(derived==null)
                continue;
            else // Polymorphic site.
            {
                for(CDSSite site:derived)
                {
                    float f = (ao.bSizeConst) ? (float)site.nCount/(float)pop.length
                                              : (float)site.nCount/(float)scp.getTotalBasesCount();
                    freqs.add(f);
                    sb.append(String.format(initData.locale, "%.3f\t(%c->%c)\t",
                                            f, site.anc, site.derived));
                    freq += f;
                    nDerived++;
                }
            }
        }
        // Add average frequency.
        if(nDerived>0)
            freq = freq/(float)nDerived;
        else
            freq = Float.NaN;
        if(nDerived>nMutations)
            nMutations = nDerived;
        sb.insert(0, String.format(initData.locale, "%d\t%d\t%d\t%.3f\t\t",
                                   pop.length, nSites, nDerived, freq));
        return sb.toString();
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

    /**********************************************************************************************
    *                                  NON-CODING REGIONS                                         *
    **********************************************************************************************/
    /**
     *  Calculates the derived allele frequencies in introns or other noncoding regions.
     *
     *  @param pop
     *  @param out
     *  @param ao
     *  @return
     */
    private String calculateStatisticsNoncoding(String[] pop, String[] out, AnalysisOptions ao)
    {
        int length = pop[0].length();
        int nSites = 0;
        int nDerived = 0;
        float freq = 0.0f;
        StringBuffer sb = new StringBuffer();
        for(int pos=0;pos<length;pos++)
        {
            SiteComposition scp = calculateSiteComposition(pop, pos);
            // If the site is divergent/monomorphic or if there are less than 2
            // valid bases at the site, then it cannot be analyzed.
            if(scp==null || scp.getValidBasesCount()<2)
                continue;
            SiteComposition sco = calculateSiteComposition(out, pos);
            if(sco==null)
                continue;
            nSites++;
            CDSSite[] derived = findDerivedMutations(scp, sco);
            // Divergent or monomorphic site.
            if(derived==null)
                continue;
            else // Polymorphic site.
            {
                for(CDSSite site:derived)
                {
                    float f = (ao.bSizeConst) ? (float)site.nCount/(float)pop.length
                                              : (float)site.nCount/(float)scp.getTotalBasesCount();
                    freqs.add(f);
                    sb.append(String.format(initData.locale, "%.3f\t(%c->%c)\t",
                                            f, site.anc, site.derived));
                    freq += f;
                    nDerived++;
                }
            }
        }
        // Add average frequency.
        if(nDerived>0)
            freq = freq/(float)nDerived;
        else
            freq = Float.NaN;
        if(nDerived>nMutations)
            nMutations = nDerived;
        sb.insert(0, String.format(initData.locale, "%d\t%d\t%d\t%.3f\t\t",
                                   pop.length, nSites, nDerived, freq));
        return sb.toString();
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

    /**
     *  Extracts the sequence of the specified region and appends it to already extracted region if necessary.
     *
     *  @param old
     *  @param pop
     *  @param iRegion
     *  @param bNogtag
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

    /*********************************************************************************************/

    /**
     *  Finds the derived mutations and counts how often the derived mutation appears
     *  in the population.
     *
     *  If the site is divergent or monomorphic, the method returns null.
     *
     *  @param scp
     *  @param sco
     *  @return
     */
    private CDSSite[] findDerivedMutations(SiteComposition scp, SiteComposition sco)
    {
        // If the site is monomorphic, return null.
        int nP = scp.getNumberOfPolymorphisms();
        if(nP==0)
            return null;
        char anc  = 0;   // Currently best ancestral base.
        int  nAnc = 0;   // Currently best ancestral base count:
                         // # of matches in the population of interest.
        char[] bases = {'A', 'C', 'G', 'T'};
        for(char c:bases)
        {
            if(sco.getBaseCount(c)>0)
            {
                int n = scp.getBaseCount(c);
                if(n>nAnc)
                {
                    nAnc = n;
                    anc = c;
                }
            }
        }
        // If there is no ancestral base return null.
        if(anc==0)
            return null;
        CDSSite[] sites = new CDSSite[nP];
        int i = 0;
        for(char c:bases)
        {
            if(c==anc)
                continue;
            int n = scp.getBaseCount(c);
            if(n>0)
            {
                CDSSite site = new CDSSite();
                site.anc = anc;
                site.derived = c;
                site.nCount = n;
                sites[i] = site;
                i++;
            }
        }
        return sites;
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
            range='<POS1-POS2; POS3>' any='<T/F>' exclnonffd='<T/F>' strlimit='<LIMIT>' combine='<T/F>' freqs='<T/F>'
            constsize='<T/F>' lenrange='<T/F>' output='<OUTFILE>'
            
            - pop:          population of interest
            - out:          outgroup
            - type:         region type
            - maxlen:       maximal intron length. Use 0 for undefined (see Remarks)
            - minlen:       minimal intron length. Use 0 for undefined (see Remarks)
            - range:        specific sites to use. See Remarks section for details.
            - any:          whether or not to use any (i.e. both syn. and nonsyn.) FFD sites.
            - exclnonffd:   whether to exclude the strains with non-FFD codons.
            - strlimit:     maximal number of strains to use
            - combine:      whether or not to combine multiple regions (see Remarks)
            - freqs:        whether to list single frequencies
            - constsize:    whether to assume the population size to be constant when calculating the frequencies
            - lenrange:     whether the minimal and maximal length should fall into the specified range
            - output:       output filename
          
            Remarks:
            If type is CDS or FFD, the parameters minlen, maxlen, combine and 
            range are ignored.
        */ 
        // Parse the parameter string.
        Pattern p = Pattern.compile("pop='(.+)'\\s+"+                   // 1
                                    "out='(.+)'\\s+"+                   // 2
                                    "type='(.+)'\\s+" +                 // 3         
                                    "maxlen='(\\d*)'\\s+"+              // 4
                                    "minlen='(\\d*)'\\s+"+              // 5
                                    "range='([-,;\\d\\s]*)'\\s+"+       // 6
                                    "any='([TF])'\\s+"+                 // 7
                                    "exclnonffd='([TF])'\\s+"+          // 8
                                    "strlimit='(\\d*)'\\s+"+            // 9
                                    "combine='([TF])'\\s+"+             // 10
                                    "freqs='([TF])'\\s+"+               // 11
                                    "constsize='([TF])'\\s+"+           // 12
                                    "lenrange='([TF])'\\s+"+            // 13
                                    "output='(.+)'$",                   // 14
                                    Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            AnalysisOptions ao = new AnalysisOptions();
            ao.strPop = m.group(1);
            ao.strOut = m.group(2);
            ao.strRegion = m.group(3);
            if(!m.group(4).isEmpty())
                ao.iMaxlen = Integer.parseInt(m.group(4));
            else
                ao.iMaxlen = Integer.MAX_VALUE;
            if(!m.group(5).isEmpty())
                ao.iMinlen = Integer.parseInt(m.group(5));
            else
                ao.iMinlen = 0;
            if(ao.iMaxlen<ao.iMinlen)
            {
                int tmp = ao.iMaxlen;
                ao.iMaxlen = ao.iMinlen;
                ao.iMinlen = tmp;
            }
            ao.maxstr = (m.group(9).isEmpty()) ? Integer.MAX_VALUE : Integer.parseInt(m.group(9));
            ao.bCombine  = m.group(10).equalsIgnoreCase("T");
            ao.bListFreqs = m.group(11).equalsIgnoreCase("T");
            ao.bSizeConst = m.group(12).equalsIgnoreCase("T");
            ao.bLenRange = m.group(13).equalsIgnoreCase("T");
            ao.strOutput = m.group(14);
            // Check whether the region type is CDS or FFD.
            if(ao.strRegion.equalsIgnoreCase("CDS") || ao.strRegion.equalsIgnoreCase("FFD"))
            {
                ao.iMaxlen = Integer.MAX_VALUE;
                ao.iMinlen = 0;
                ao.bCombine = false;
            }
            // Range.
            ao.sites = OptionsDialog.extractSites(m.group(6));
            ao.bUseAny = m.group(7).equalsIgnoreCase("T");
            ao.bNonFfd = m.group(8).equalsIgnoreCase("T");
            return ao;
        }
        else
        {
            return (new OptionsDialog(pops, regs)).getOptions();
        }
    }
}