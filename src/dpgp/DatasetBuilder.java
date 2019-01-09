/*
    File:
        DatasetBuilder.java
 *
    Revision:
        1.0.0.0
 *
    Description:
        Builds the dataset from the alignment by exctracting the specified regions.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package dpgp;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.StrainEntry;
import java.io.Serializable;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;


public class DatasetBuilder
{
    public static class Fragment implements Serializable
    {
        public String strType;
        public int iStart;
        public int iEnd;
        public boolean bLeading;

        @Override
        public Fragment clone()
        {
            Fragment f = new Fragment();
            f.strType = this.strType;
            f.iStart = this.iStart;
            f.iEnd = this.iEnd;
            f.bLeading = this.bLeading;
            return f;
        }
    };

    public static final String GENOMIC_POSITION = "Genomic position";

    private String strLastErr = null;

    /**
     *  Builds a dataset using the the specified multiple alignment and the HashMap
     *  of fragments.
     *  The keys of the HashMap must be the gene names and the values are the vectors
     *  of fragments.
     *  The method creates one GeneEntry per key using the fragments from the corresponding
     *  fragments vector.
     *
     *  @param alignment
     *  @param fragments
     *  @param strChromosome
     *  @return
     */
    public Dataset build(DPGPAlignment alignment, HashMap<String, Vector<Fragment>> fragments, String strChromosome)
    {
        if(alignment==null || fragments==null || fragments.size()==0 || strChromosome==null)
        {
            strLastErr = "Invalid parameter";
            return null;
        }
        Dataset dataset = new Dataset();
        Iterator it = fragments.keySet().iterator();
        while(it.hasNext())
        {
            String strGene = (String)it.next();
            Vector<Fragment> frags = fragments.get(strGene);
            Dataset tmp = build(alignment, strGene, frags, strChromosome);
            if(tmp==null)
                return null;
            else
                dataset.merge(tmp);
        }
        return dataset;
    }

    /**
     *  Builds a dataset consisting of a single GeneEntry using the provided multiple
     *  alignment and the vector of fragments.
     *
     *  @param alignment
     *  @param strGene 
     *  @param fragments
     *  @param strChromosome
     * @return
     */
    public Dataset build(DPGPAlignment alignment, String strGene, Vector<Fragment> fragments, String strChromosome)
    {
        if(alignment==null || strGene==null || strGene.isEmpty() || fragments==null || fragments.isEmpty() || strChromosome==null)
        {
            strLastErr = "Invalid parameter";
            return null;
        }
        // Sort the fragments.
        Collections.sort(fragments, new Comparator<Fragment>() 
        {
            public int compare(Fragment o1, Fragment o2)
            {
                if(o1.bLeading)
                    return (o1.iStart<o2.iStart ? -1 : (o1.iStart==o2.iStart ? 0 : 1));
                else
                    return (o1.iStart>o2.iStart ? -1 : (o1.iStart==o2.iStart ? 0 : 1));
            }
        });
        // Create regions templates.
        SequenceFragment[] templates = new SequenceFragment[fragments.size()];
        boolean b = fragments.get(0).bLeading;
        int iStart = (b) ? fragments.get(0).iStart : fragments.get(0).iEnd;
        for(int i=0;i<fragments.size();i++)
        {
            Fragment f = fragments.get(i);
            SequenceFragment sf = (b) ? new SequenceFragment(f.strType, alignment, 0, f.iStart, f.iEnd)
                                      : new SequenceFragment(f.strType, alignment, 0, f.iEnd, f.iStart);
            sf.setStart(Math.abs(sf.getGenomicStart()-iStart)+1);
            sf.setEnd(0);
            templates[i] = sf;
        }
        for(int i=0;i<templates.length;i++)
        {
            SequenceFragment f = templates[i];
            f.setEnd(f.getStart()+f.getSequenceLength()-1);
            for(int j=i+1;j<templates.length;j++)
            {
                SequenceFragment f2 = templates[j];
                int iOffset = 0;
                // If the fragments overlap.
                if((b) ? f2.getGenomicStart()<f.getGenomicEnd() : f2.getGenomicStart()>f.getGenomicEnd())
                    iOffset = f.getGapsCount(f.getGenomicStart(), f2.getGenomicStart());
                else // Otherwise move the region to the right.
                    iOffset = f.getGapsCount(f.getGenomicStart(), f.getGenomicEnd());
                f2.setStart(f2.getStart()+iOffset);
            }
        }
        Dataset ds = new Dataset();
        GeneEntry ge = build(alignment, strGene, templates, strChromosome);
        ds.addGene(ge);
        return ds;
    }

    private GeneEntry build(DPGPAlignment alignment, String strGene, SequenceFragment[] templates, String strChromosome)
    {
        GeneEntry ge = new GeneEntry(strGene, "");
        String[] specs = alignment.getSpeciesNames();
        String[] strains = alignment.getStrainsNames();
        int nStrains = alignment.getStrainsCount();
        for(int i=0;i<nStrains;i++)
        {
            StrainEntry se = new StrainEntry(specs[i], strains[i]);
            se.setChromosome(strChromosome);
            se.addPopulations(specs[i]);
            for(SequenceFragment sf:templates)
            {
                SequenceFragment frag = new SequenceFragment(sf.getType(), alignment, i, sf.getGenomicStart(), sf.getGenomicEnd());
                frag.setStart(sf.getStart());
                frag.setEnd(sf.getEnd());
                se.addRegion(frag);
            }
            ge.addStrain(se);
        }
        if(ge.getStrainsCount()>0)
        {
            StrainEntry se = ge.getStrainEntry(0);
            if(se.getRegionsCount()>0)
            {   
                SequenceFragment f = (SequenceFragment)se.getRegion(0);
                SequenceFragment l = (SequenceFragment)se.getRegion(se.getRegionsCount()-1);
                ge.addProperty(GENOMIC_POSITION, new int[]{f.getGenomicStart(), l.getGenomicEnd()});
            }
        }
        return ge;
    }

    /**
     *  Returns the string description of the last error occured.
     *
     *  @return
     */
    public String getLastErrorString()
    {
        return strLastErr;
    }
}
