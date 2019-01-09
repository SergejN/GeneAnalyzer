/*
    File:
        GeneFusioner.java
 *   
    Revision:
        1.4.0.1
 * 
    Description:
        Fuses the genes of the dataset.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms.alignment;

import algorithms.SequenceBuffer;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.util.Vector;
import plugin.classes.IAligner;


public class GeneFusioner 
{
    private String strLastError = "";
    private IWaitDialog wd      = null;


    public GeneFusioner(IWaitDialog wd)
    {
        this.wd = wd;
    }


    /**
     *  Returns the string description of the last error occured.
     * 
     *  @return
     */
    public String getLastError()
    {
        return strLastError;
    }

    /**
     *  The agorithm first aligns all genes with each other to identify the best
     *  reciprocal hits. These hits are now displayed to allow the user to select
     *  which ones keep.
     *  If the dataset is null or empty, the method returns the dataset passed.
     *
     *  @param ds
     *  @param aligner
     *  @return
     */
    public Dataset fuseGenes(Dataset ds, IAligner[] aligners)
    {
        if(ds==null || ds.getGenesCount()==0)
        {
            strLastError = "The dataset is null or empty";
            return ds;
        }
        AlignerSelectionDialog.Options opt = (new AlignerSelectionDialog()).getOptions(aligners);
        if(opt==null)
        {
            strLastError = "Cancelled by user";
            return ds;
        }
        wd.show(IWaitDialog.TYPE.Aligner);
        IAligner aligner = opt.aligner;
        int nGenes = ds.getGenesCount();
        int[] pairs = new int[nGenes];
        Alignment[] alignments = new Alignment[nGenes];
        for(int i=0;i<nGenes;i++)
        {
            GeneEntry ge1 = ds.getGeneEntry(i);
            pairs[i] = -1;
            alignments[i] = null;
            for(int j=0;j<nGenes;j++)
            {
                if(i==j)
                    continue;
                GeneEntry ge2 = ds.getGeneEntry(j);
                if(opt.bUseSim && areSimilar(ge1.getCommonName(), ge2.getCommonName()))
                {
                    Alignment ali = constructAlignment(aligner, ge1, ge2, "");
                    if(ali==null)
                        continue;
                    if(alignments[i]==null || ali.getScore()>alignments[i].getScore())
                    {
                        alignments[i] = ali;
                        pairs[i] = j;
                    }
                }
            }
        }
        // Process alignments.
        Dataset dataset = new Dataset();
        for(int i=0;i<nGenes;i++)
        {
            // The best match for the gene i.
            int bm = pairs[i];
            // Reciprocal best hit.
            if(bm!=-1 && pairs[bm]==i)
            {
                String strCommonName = ds.getGeneEntry(i).getCommonName()+"_"+ds.getGeneEntry(bm).getCommonName();
                GeneEntry ge = new GeneEntry(strCommonName, "");
                int nStrains = alignments[i].getStrainsCount();
                for(int n=0;n<nStrains;n++)
                    ge.addStrain(alignments[i].getStrainEntry(n));
                dataset.addGene(ge);
                pairs[i] = -1;
            }
        }
        wd.close();
        return dataset;
    }

    /**
     *  Aligns two genes and returns the alignment of the strain entries of both genes.
     *  The order of the strain entries in the alignment is the same as if the strain
     *  entries of both genes, beginning with the first strain entry of ge1, were added
     *  one by one to the alignment.
     *
     *  @param ge1
     *  @param ge2
     *  @param strParams
     *  @return
     */
    private Alignment constructAlignment(IAligner aligner, GeneEntry ge1, GeneEntry ge2, String strParams)
    {
        if(ge1.getStrainsCount()==0 || ge2.getStrainsCount()==0)
        {
            strLastError = "One of the genes does not have any strains";
            return null;
        }
        // Find the best initial pairwise alignment.
        Alignment ba = null;
        int[] indices = {-1, -1};   // Indices of the strains which align best.
        StrainEntry[] entries = new StrainEntry[2];
        int nStrains1 = ge1.getStrainsCount();
        int nStrains2 = ge2.getStrainsCount();
        for(int i=0;i<nStrains1;i++)
        {
            entries[0] = new StrainEntry("Strain 1", "");
            GeneRegion reg1 = new GeneRegion(GeneRegion.UNNAMED);
            reg1.setStart(1);
            String strSeq1 = ge1.getStrainEntry(i).getCompleteSequence();
            reg1.setEnd(strSeq1.length());
            reg1.setSequence(strSeq1);
            entries[0].addRegion(reg1);
            // If the entry sequence is not valid, skip it, since the alignment does not make any sence.
            if(strSeq1.matches("[NXnx-]+"))
                continue;
            for(int j=0;j<nStrains2;j++)
            {
                entries[1] = new StrainEntry("Strain 2", "");
                GeneRegion reg2 = new GeneRegion(GeneRegion.UNNAMED);
                reg2.setStart(1);
                String strSeq2 = ge2.getStrainEntry(j).getCompleteSequence();
                reg2.setEnd(strSeq2.length());
                reg2.setSequence(strSeq2);
                entries[1].addRegion(reg2);
                // If the entry sequence is not valid, skip it, since the alignment does not make any sence.
                if(entries[1].getCompleteSequence().matches("[NXnx-]+"))
                    continue;
                Alignment ali = aligner.alignStrains(entries, strParams);                
                if(ba==null || ali.getScore()>ba.getScore())
                {
                    ba = ali;
                    indices[0] = i;
                    indices[1] = j;
                }
            }
        }
        // If no alignment was found, return null.
        if(ba==null)
        {
            strLastError = "One of the genes consists of invalid sequences only";
            return null;
        }
        // Align the sequences.
        entries[0] = ge1.getStrainEntry(indices[0]);
        entries[1] = ge2.getStrainEntry(indices[1]);
        Alignment ali = aligner.alignStrains(entries, strParams);
        // Otherwise, extend the alignments to all strains.
        Vector<StrainEntry> tmp = new Vector<StrainEntry>();
        StrainEntry[] ent = extendAlignment(ge1, ali.getStrainEntry(0), indices[0]);
        for(StrainEntry se:ent)
            tmp.add(se);
        ent = extendAlignment(ge2, ali.getStrainEntry(1), indices[1]);
        for(StrainEntry se:ent)
            tmp.add(se);
        return new Alignment(ba.getScore(), tmp.toArray(new StrainEntry[tmp.size()]), "Pairwise gene alignment", "");
    }

    /**
     *  Extends the alignment to all strains of a gene entry.
     *
     *  @param ge
     *  @param ae
     *  @param iIndex
     *  @return
     */
    private StrainEntry[] extendAlignment(GeneEntry ge, StrainEntry ae, int iIndex)
    {
        String as = ae.getCompleteSequence();
        int length = as.length();
        int nStrains = ge.getStrainsCount();
        // Copy old sequences into an array and create a new sequences buffer.
        String[] old = new String[nStrains];
        SequenceBuffer[] seqs = new SequenceBuffer[nStrains];
        for(int i=0;i<nStrains;i++)
        {
            old[i] = ge.getStrainEntry(i).getCompleteSequence();
            seqs[i] = new SequenceBuffer(length);
        }
        // Iterate through the newly aligned sequence and compare it to the old sequence.
        // The only difference, which can appear, is a gap in the NEW sequence, due to new
        // alignment. In this case introduce a gap into all old sequences. Otherwise copy the
        // old bases into the buffer.
        int po = 0; // Position in the OLD sequence
        int pn = 0; // Position in the NEW sequence
        // Length of the old sequence.
        int ol = ge.getStrainEntry(0).getCompleteSequence().length();
        while(pn<length)
        {
            if(po>=ol || old[iIndex].charAt(po)!=as.charAt(pn))
            {
                for(int i=0;i<nStrains;i++)
                    seqs[i].appendBase('-');
                pn++;
            }
            // The old sequence and the new one are equal at this position.
            else
            {
                for(int i=0;i<nStrains;i++)
                    seqs[i].appendBase(old[i].charAt(po));
                po++;
                pn++;
            }
        }
        // Create new strain entries.
        StrainEntry[] entries = new StrainEntry[nStrains];
        int nRegions = ae.getRegionsCount();
        for(int n=0;n<nStrains;n++)
        {
            StrainEntry orig = ge.getStrainEntry(n);
            entries[n] = new StrainEntry(orig.getSpeciesName(), orig.getStrainName());
            entries[n].setChromosome(orig.getChromosome());
            entries[n].addPopulations(orig.listPopulations());
            for(int i=0;i<nRegions;i++)
            {
                GeneRegion template = ae.getRegion(i);
                GeneRegion reg = new GeneRegion(template.getType());
                reg.setStart(template.getStart());
                reg.setEnd(template.getEnd());
                reg.setSequence(seqs[n].substring(template.getStart()-1, template.getEnd()));
                entries[n].addRegion(reg);
            }
        }
        return entries;
    }

    /**
     *  Returns true if at least 75% of the shorter string matches with the longer one.
     *
     *  @param strName1
     *  @param strName2
     *  @return
     */
    private boolean areSimilar(String strName1, String strName2)
    {
        int l = Math.min(strName1.length(), strName2.length());
        int i = 0;
        while(i<l && strName1.charAt(i)==strName2.charAt(i))
            i++;
        return (float)i>=0.75f*l;
    }
}