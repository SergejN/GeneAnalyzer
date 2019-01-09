/*
    File:
        DatasetFixer.java
 *
    Revision:
        1.3.0.1
 *
    Description:
        Fixes some annotation or sequence mistakes in the gene entries of the
        dataset.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package kernel.fix;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.util.Vector;
import kernel.InitData;

public class DatasetFixer
{
    public class RemovedRegion
    {
        public int iStart     = 0;
        public int iEnd       = 0;
        public String[] seqs  = null;
    };


    public static final String DELETED_BASES    = "Deleted bases";
    public static final String MASKED_CODONS    = "Masked codons";

    private InitData initData = null;

    
    public DatasetFixer(InitData initData)
    {
        this.initData = initData;
    }

    /**
     *  Masks the premature terminal codons in the coding sequence with three
     *  consequent gaps ("---") and returns the reference to the dataset.
     *
     *  @param ds
     *  @return
     */
    public Dataset maskPrematureTerminalCodons(Dataset ds)
    {
        int nGenes = ds.getGenesCount();
        initData.wd.show(IWaitDialog.TYPE.Kernel);
        for(int i=0;i<nGenes;i++)
            maskPrematureTerminalCodons(ds.getGeneEntry(i));
        initData.wd.close();
        return ds;
    }

    /**
     *  Masks the premature terminal codons in the coding sequence with three
     *  consequent gaps ("---") and returns the reference to the gene entry.
     *
     *  @param ge
     *  @return
     */
    private GeneEntry maskPrematureTerminalCodons(GeneEntry ge)
    {
        int nStrains = ge.getStrainsCount();
        for(int i=0;i<nStrains;i++)
        {
            StrainEntry se = ge.getStrainEntry(i);
            int nRegs = se.getRegionsCount();
            int iFrame = 0;
            int iCdsLength = se.getCodingSequence().length();
            int iCheckedLength = 0;
            int ile = 0;
            for(int n=0;n<nRegs;n++)
            {
                GeneRegion r = se.getRegion(n);
                if(r.hasType(GeneRegion.EXON))
                {                    
                    String strSeq = r.getSequence();
                    int l = strSeq.length();
                    iCheckedLength += l;
                    int pos = (iFrame==0) ? 0 : (3-iFrame);
                    // If the frame is not 0, check whether the last bases of the last exon and
                    // the beginning of the current exon build a PTC.
                    if(iFrame>0)
                    {
                        String tmp = se.getRegion(ile).getSequence();
                        String strCodon = tmp.substring(tmp.length()-iFrame);
                        strCodon += strSeq.substring(0, pos);
                        if(initData.ct.isTerminal(strCodon))
                        {
                            for(int k=0;k<pos;k++)
                                r.replaceBase(k, '-');
                            GeneRegion prev = se.getRegion(ile);
                            int pl = prev.getSequence().length();
                            for(int k=1;k<=iFrame;k++)
                                se.getRegion(ile).replaceBase(pl-k, '-');
                        }
                    }
                    while(pos<=l-3)
                    {
                        String strCodon = strSeq.substring(pos, pos+3);
                        // If the codon is terminal, check, whether it is in the last
                        if(initData.ct.isTerminal(strCodon))
                        {
                            if(pos!=l-3 || iCheckedLength!=iCdsLength)
                            {
                                r.replaceBase(pos, '-');
                                r.replaceBase(pos+1, '-');
                                r.replaceBase(pos+2, '-');
                            }
                        }
                        pos+=3;
                    }
                    iFrame = (iFrame+l%3)%3;
                    ile = n;
                }
            }
        }
        return ge;
    }

    /**
     *  Removes ORF disrupting insertions.
     *  E.g. assume the alignment of a gene entry looks as follows:
     *      ATGAAA-TAA
     *      ATGAAA-TAA
     *      ATGAAAATAA
     *  The last sequence obviously disrupts the ORF by introducing an additional
     *  base A. This will be removed.
     *
     *  Remarks:
     *      Only the insertions in the coding sequence will be removed. The sequence
     *      of all non-coding regions is untouched.
     *      A RemovedRegion instance is added to the gene region of the reference
     *      strain entry.
     *
     *  @param ds
     *  @return 
     */
    public Dataset removeORFDisruptingInsertions(Dataset ds)
    {
        boolean b = false;                          // Shows whether or not the wait dialog is displayed.
        int nGenes = ds.getGenesCount();
        DatasetFixerRemoveDlg.Options opt = null;
        initData.wd.show(IWaitDialog.TYPE.Kernel);
        b = true;
        for(int i=0;i<nGenes;i++)
        {
            GeneEntry ge = ds.getGeneEntry(i);
            if(opt==null)
            {
                initData.wd.close();
                b = false;
                opt = new DatasetFixerRemoveDlg().getFixOptions(ge);
            }
            if(opt==null)
                return ds;
            // Find the reference strain.
            int iIndex = getReferenceStrainIndex(opt, ge, initData.wd);
            if(iIndex==-1)
                return ds;
            if(!b)
            {
                initData.wd.show(IWaitDialog.TYPE.Kernel);
                b = true;
            }
            removeORFDisruptingInsertions(ge, iIndex);
            if(!opt.bKeep)
                opt = null;
        }
        initData.wd.close();
        return ds;
    }

    /**
     *  Removes ORF disrupting insertions.
     *  E.g. assume the alignment of a gene entry looks as follows:
     *      ATGAAA-TAA
     *      ATGAAA-TAA
     *      ATGAAAATAA
     *  The last sequence obviously disrupts the ORF by introducing an additional
     *  base A. This will be removed.
     *
     *  Note:
     *      Only the insertions in the coding sequence will be removed. The sequence
     *      of all non-coding regions is untouched. If the number of consequent gaps
     *      id a multiple of 3 and the gaps are in-frame, they are not removed.
     *
     *  Algorithm:
     *  Iterate through the CDS and look for blocks of gaps in the reference strain.
     *  Once such a block was found, check the following:
     *      - if the start of the block is the first base in a codon and the end of
     *        the block is the last base in the same or in other codon - keep the block,
     *        since it does not disrupt the ORF of the reference strain.
     *      - if the block starts in one codon and ends in another codon, delete all
     *        codons which contain bases of the block
     *      - if the block starts and ends in the same codon, remove the bases of the
     *        block only
     *
     *  Example 1:
     *      xxx ATG xxx
     *      xxx --- xxx
     *    These theree gaps do not disrupt the ORF of the strain -> sites not removed
     *
     *  Example 2:
     *      xxx A-- xxx     xxx -ATG xxx    xxx AA- --- --A xxx
     *      xxx -T- xxx     xxx AATG xxx    xxx XXX XXX XXX xxx
     *      xxx --G xxx     xxx AATG xxx    xxx XXX XXX XXX xxx
     *    In these cases the sites with gaps disrupt the ORF and, thus, are removed.
     * 
     *  @param ds
     *  @return
     */
    private GeneEntry removeORFDisruptingInsertions(GeneEntry ge, int iRefStrain)
    {
        int nStrains = ge.getStrainsCount();
        if(nStrains==0)
            return ge;
        StrainEntry se = ge.getStrainEntry(0);        
        int iFrame = 0;     // Current frame.        
        // Iterate through the gene regions to find the exons.
        Vector<Integer> sites = new Vector<Integer>();
        int nRegs = se.getRegionsCount();
        for(int i=0;i<nRegs;i++)
        {
            if(se.getRegion(i).hasType(GeneRegion.EXON))
            {
                int nRemoved = 0;   // Number of already removed bases.
                String[] seqs = new String[nStrains];
                for(int n=0;n<nStrains;n++)
                    seqs[n] = ge.getStrainEntry(n).getRegion(i).getSequence();
                // Iterate through the sequences and search for potential ORF distuptions.
                int pos = 0;
                int iLength = seqs[0].length();                
                while(pos<iLength)
                {
                    // Add the gaps to the block of gaps.
                    if(seqs[iRefStrain].charAt(pos)=='-')
                        sites.add(pos-nRemoved);
                    // Once the end of the block is reached, remove the ORF disrupting insertions.
                    else if(sites.size()>0)
                    {
                        nRemoved += removeBases(ge, i, sites, iFrame);
                        sites.clear();
                    }
                    pos++;
                }
        if(ge.getCommonName().startsWith("CG14341"))
        {
            int j = 0;
            ge.getStrainEntry(0).getRegion(i).getSequence();
        }
                // If the sequence ends with a block of gaps, remove those, too.
                if(sites.size()>0)
                {
                    nRemoved += removeBases(ge, i, sites, iFrame);
                    sites.clear();
                }
                // Update the frame.
                iFrame = (iFrame+seqs[0].length()%3)%3;
                // Update positions.
                for(int n=0;n<nStrains;n++)
                {
                    for(int k=i+1;k<nRegs;k++)
                    {
                        GeneRegion reg = ge.getStrainEntry(n).getRegion(k);
                        reg.setStart(reg.getStart()-nRemoved);
                        reg.setEnd(reg.getEnd()-nRemoved);
                    }
                }
            }
        }
        return ge;
    }

    private int removeBases(GeneEntry ge, int iRegion, Vector<Integer> sites, int iFrame)
    {
        // If the block covers one or more codons completely, keep the whole block.
        int iStart = sites.get(0);
        int iStartFrame = (iFrame+iStart%3)%3;
        int iEnd = sites.get(sites.size()-1);
        int iEndFrame = (iFrame+iEnd%3)%3;
        if(iStartFrame==0 && iEndFrame==2)
            return 0;
        // Otherwise remove the specified sites.
        int nStrains = ge.getStrainsCount();
        int nSites = sites.size();
        int nRemoved = 0;
        for(int i=0;i<nStrains;i++)
        {
            GeneRegion reg = ge.getStrainEntry(i).getRegion(iRegion);
            nRemoved = reg.removeBases(iStart, nSites);
        }        
        return nRemoved;
    }

    /**
     *  Iterates though the strains and returns the index of the reference strain within
     *  the gene entry.
     *  If the strain with the specified name is not found, the options dialog is displayed.
     *  If the user then clicks CANCEL, the method returns -1.
     *
     *  @param opt
     *  @param ge
     *  @return
     */
    private int getReferenceStrainIndex(DatasetFixerRemoveDlg.Options opt, GeneEntry ge, IWaitDialog wd)
    {
        int nStrains = ge.getStrainsCount();
        for(int n=0;n<nStrains;n++)
        {
            if(ge.getStrainEntry(n).getStrainName().equalsIgnoreCase(opt.strRefStrain))
                return n;
        }
        wd.close();
        DatasetFixerRemoveDlg.Options tmp = new DatasetFixerRemoveDlg().getFixOptions(ge);
        if(tmp==null)
            return -1;
        for(int n=0;n<nStrains;n++)
        {
            if(ge.getStrainEntry(n).getStrainName().equalsIgnoreCase(tmp.strRefStrain))
            {
                if(tmp.bKeep)
                {
                    opt.strRefStrain = tmp.strRefStrain;
                    opt.bKeep = tmp.bKeep;
                }
                wd.show(IWaitDialog.TYPE.Kernel);
                return n;
            }
        }        
        return -1;
    }
}
