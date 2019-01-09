/*
    File:
        QualityChecker.java
 *   
    Revision:
        1.3.0.1
 * 
    Description:
        Performs a low-level sequence quality check. The performed tests are:
            - properties: i.e. all class fields are filled
                          properly and completely
            - regions: i.e. all strains of a gene entry have the same
                       regions with the same boundaries
            - ORF: i.e. the coding sequences of the genes have a valid ORF
            - PTC: i.e. no coding sequence has a premature terminal codon
            - Introns: i.e. introns have valid intron boundaries (GT..AG)
            - Gene boundaries: i.e. coding sequences have valid start and stop codons
 
        Output:
            The results of the quality check are annotated on gene entry (summary)
            and on the faulty strain entry (detailed description). The annotation
            name is defined in PROPERTY_NAME and the type is String. If the
            strain/gene entry is correct, getProperty(PROPERTY_NAME) returns null,
            i.e. nothing is annotated.
 
            Every gene entry is assigned a quality value (see Remarks for details):
                0       no errors
                1       invalid ORF
                2       premature terminal codon
                3       more remarkable deficits, might still be usable
                4       more severe deficits, generally not usable
                5       fatal errors, cannot be used
 
        Remarks:
            QUALITY LEVEL 0 means, that the gene entry does not have any annotation
            deficits and can be analyzed without any restrictions.
            QUALITY LEVEL 1 means, that the gene has some not very important
            annotation errors, like missing start/stop codon, intron boundaries etc.
            The analysis of the gene is still possible in most cases.
            QUALITY LEVEL 2 means, that the coding region does not have a valid
            ORF. LEVEL 2 is only assigned to strains, which do not end in exon or
            intron, i.e. there is at least one more gene region after the last
            exon or intron - for example, 3'UTR.
            QUALITY LEVEL 3 means, that the strain has a premature terminal codon.
            The analyses are still possible, however, a PTC is a good indication
            that the sequence might be incorrect.
            QUALITY LEVEL 4 means, that the gene has overlapping regions, or that
            the start/end of the regions are incorrect, or that the sequence length
            is not end-start+1. The entry should not be exported but it may be analyzed
            if this information is not required.
            QUALITY LEVEL 5 means, that the gene entry has severe annotation
            problems, like missing sequence of a gene region, no regions
            in the strain entry or no strain entries in the gene entry. The
            further analysis of such genes is highly discouraging.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;


public class QualityChecker 
{
    /**
     *  Name of the property holding the description of the annotation
     *  deficits.
     */
    public static final String QUALITY_DESCRIPTION = "QualityCheckResult";
    
    /**
     *  Name of the property holding the integer value, specifying the
     *  level of data quality.
     */
    public static final String QUALITY_LEVEL = "QualityLevel";
    
    /**
     *  Quality levels.
     */
    public static final int QUALITY_LEVEL_0 = 0;
    public static final int QUALITY_LEVEL_1 = 1;
    public static final int QUALITY_LEVEL_2 = 2;
    public static final int QUALITY_LEVEL_3 = 3;
    public static final int QUALITY_LEVEL_4 = 4;
    public static final int QUALITY_LEVEL_5 = 5;

    private InitData initData = null;

    public QualityChecker(InitData initData)
    {
        this.initData = initData;
    }

    /**
     *  Validates the complete dataset. It is reasonable to call
     *  validateDataset at the very beginning, when a dataset is
     *  loaded. However, if you edit one gene entry, it is more
     *  efficient validate that single entry.
     * 
     *  @param ds
     */
    public void validateDataset(Dataset ds)
    {
        for(int i=0;i<ds.getGenesCount();i++)
            validateGeneEntry(ds.getGeneEntry(i));
    }    
    
    /**
     *  Validates a single gene entry.
     *
     *  For a gene entry following properties are tested:
     *      - common name
     *      - strains
     *
     *  @param ge
     *  @return     quality level of the gene entry
     */
    public int validateGeneEntry(GeneEntry ge)
    {
        int quality = QUALITY_LEVEL_0;
        StringBuffer annot = new StringBuffer();
        // Common name.
        if(ge.getCommonName().isEmpty())
        {
            annot.append("Missing species name;");
            quality = QUALITY_LEVEL_1;
        }
        // Strains.
        if(ge.getStrainsCount()==0)
        {
            annot.append("No strains specified;");
            quality = QUALITY_LEVEL_5;
        }
        else // Validate strains.
        {
            for(int i=0;i<ge.getStrainsCount();i++)
            {
                // Strain entry quality.
                StrainEntry se = ge.getStrainEntry(i);
                quality = Math.max(quality, validateStrainEntry(se));
                // All strain entries in the gene entry must have identical regions
                // i.e. same type and same start and end positions.
                if(i>0)
                {
                    StrainEntry prev = ge.getStrainEntry(i-1);
                    // If the strains have unequal number of regions, set quality
                    // level to QUALITY_LEVEL_5.
                    if(se.getRegionsCount()!=prev.getRegionsCount())
                    {
                        quality = QUALITY_LEVEL_5;
                        annot.append(String.format("Strains %s and %s have different number of regions;",
                                prev.getStrainName(), se.getStrainName()));
                    }
                    else // Compare the regions
                    {
                        for(int n=0;n<se.getRegionsCount();n++)
                        {
                            GeneRegion rc = se.getRegion(n);
                            GeneRegion rp = prev.getRegion(n);
                            if(!rc.getType().equalsIgnoreCase(rp.getType()) ||  // Different type
                                rc.getStart()!=rp.getStart() ||                 // Different start positions
                                rc.getEnd()!=rp.getEnd())
                            {
                                quality = QUALITY_LEVEL_5;
                                annot.append(String.format("Region #%d differs in %s and %s;",
                                        n+1, prev.getStrainName(), se.getStrainName()));
                            }
                        }
                    }
                }
            }
        }
        // Assign the quality.
        ge.addProperty(QUALITY_LEVEL, quality);
        ge.addProperty(QUALITY_DESCRIPTION, annot.toString());
        return quality;
    }    
    
    /**
     *  Validates a single strain entry.
     *
     *  For a strain entry following properties are tested:
     *      - species name
     *      - strain name
     *      - chromosome
     *      - populations
     *      - sequence
     *
     *  @param se
     *  @return     quality level of the strain entry
     */
    public int validateStrainEntry(StrainEntry se)
    {
        int quality = QUALITY_LEVEL_0;
        StringBuffer annot = new StringBuffer();
        // Species name.
        if(se.getSpeciesName().isEmpty())
        {
            annot.append("Missing species name;");
            quality = QUALITY_LEVEL_1;
        }
        // Strain name.
        if(se.getStrainName().isEmpty())
        {
            annot.append("Missing strain name;");
            quality = QUALITY_LEVEL_1;
        }
        // Chromosome.
        if(se.getChromosome().isEmpty())
        {
            annot.append("Missing chromosomal location;");
            quality = QUALITY_LEVEL_1;
        }
        // Populations.
        if(se.listPopulations().length==0)
        {
            annot.append("No populations assigned;");
            quality = QUALITY_LEVEL_4;
        }
        // Sequence.
        if(se.getRegionsCount()==0)
        {
            annot.append("No gene regions specified;");
            quality = QUALITY_LEVEL_5;
        }
        else    // ORF, PTC, overlapping regions and gaps.
        {
            // Coding sequence quality.
            String strSeq = se.getCodingSequence();
            // Intact ORF.
            if(strSeq.length()%3!=0)
            {
                // If the last region of the strain entry is neither an exon nor
                // an intron, then the sequence does not have a valid ORF.
                GeneRegion lr = se.getRegion(se.getRegionsCount()-1);
                if(!lr.hasType(GeneRegion.EXON) && !lr.hasType(GeneRegion.INTRON))
                {
                    annot.append(String.format("Invalid ORF (%d bp);", strSeq.length()));
                    quality = Math.max(quality, QUALITY_LEVEL_2);
                }
                else
                {
                    annot.append(String.format("Incomplete ORF (%d bp);", strSeq.length()));
                    quality = Math.max(quality, QUALITY_LEVEL_1);
                }
            }
            // PTC.
            int pos = (strSeq.length()>2) ? 0 : 3;
            boolean bStop = false;  // Flag specifying whether a proper terminal codon was found.
            int length = (strSeq.length()/3)*3;
            while(pos<length)
            {
                String strCodon = strSeq.substring(pos, pos+3);
                // If the codon is terminal ...
                if(initData.ct.isTerminal(strCodon))
                {
                    // ... and the codon is not the last codon in the sequence -> PTC
                    if(pos<strSeq.length()-3)
                    {
                        annot.append(String.format("A premature terminal codon found at position %d of the coding sequence;", pos+1));
                        quality = Math.max(quality, QUALITY_LEVEL_3);
                    }
                    else    // ... otherwise the codon is a proper terminal codon
                        bStop = true;
                }
                pos += 3;
            }
            // Proper terminal codon.
            if(!bStop)
            {
                annot.append("Stop codon missing;");
                quality = Math.max(quality, QUALITY_LEVEL_1);
            }
            // Start codon.
            if(length>2)
            {
                if(!initData.ct.isStartCodon(strSeq.substring(0,3)))
                {
                    annot.append("Start codon missing;");
                    quality = Math.max(quality, QUALITY_LEVEL_1);
                }
            }
            // Overall sequence quality.
            if(se.getRegion(0).getStart()>1)
            {
                annot.append(String.format("The first region (%s) does not begin at sequence start",
                        se.getRegion(0).getType()));
                quality = Math.max(quality, QUALITY_LEVEL_4);
            }
            int ile = se.getRegion(0).getEnd(); // last end position.
            quality = Math.max(quality, validateGeneRegion(se.getRegion(0)));
            for(int i=1;i<se.getRegionsCount();i++)
            {
                GeneRegion reg = se.getRegion(i);
                quality = Math.max(quality, validateGeneRegion(reg));
                // Gap.
                if(reg.getStart()>ile+1)
                {
                    annot.append(String.format("Gap between regions %d (%s) and %d (%s);",
                            i, se.getRegion(i-1).getType(), i+1, se.getRegion(i).getType()));
                    quality = Math.max(quality, QUALITY_LEVEL_4);
                }
                // Overlap.
                if(reg.getStart()<=ile)
                {
                    annot.append(String.format("Overlap between regions %d (%s) and %d (%s);",
                            i, se.getRegion(i-1).getType(), i+1, se.getRegion(i).getType()));
                    quality = Math.max(quality, QUALITY_LEVEL_4);
                }
                ile = reg.getEnd();
            }
        }
        // Assign the quality.
        se.addProperty(QUALITY_LEVEL, quality);
        se.addProperty(QUALITY_DESCRIPTION, annot.toString());
        return quality;
    }
    
    /**
     *  Validates a single gene region.
     *
     *  For a region following properties are tested:
     *      - region type
     *      - sequence: the sequence must be available and have the length end-start+1
     *      - start>0
     *      - end>0 and end>start
     * 
     *  @param gr
     *  @return     quality level of the gene region
     */
    public int validateGeneRegion(GeneRegion gr)
    {
        int quality = QUALITY_LEVEL_0;
        StringBuffer annot = new StringBuffer();
        // Region type.
        if(gr.getType().isEmpty())
        {
            annot.append("Region type missing;");
            quality = QUALITY_LEVEL_4;
        }
        // Region start/end.
        if(gr.getStart()<1 || gr.getEnd()<gr.getStart())
        {
            annot.append("Invalid start/end position;");
            quality = QUALITY_LEVEL_4;
        } 
        // Sequence.
        if(gr.getSequence().isEmpty())
        {
            annot.append("Missing sequence;");
            quality = QUALITY_LEVEL_5;
        }  
        // If the code type is not set, initialize it to SettingsManager.CT_SIMPLE.
        if(initData.sm.getSetting("", SettingsManager.CODETYPE)==null)
            initData.sm.addSetting("", SettingsManager.CODETYPE, SettingsManager.CT_SIMPLE);
        // Check for invalid chracters.
        else if(initData.sm.getSetting("", SettingsManager.CODETYPE).equalsIgnoreCase(SettingsManager.CT_SIMPLE) &&
                !gr.getSequence().matches("[ACGTacgt-]+"))
        {
            annot.append("Sequence contains invalid characters;");
            quality = QUALITY_LEVEL_5;
        }
        else if(initData.sm.getSetting("", SettingsManager.CODETYPE).equalsIgnoreCase(SettingsManager.CT_EXTENDED) &&
                !gr.getSequence().matches("[ACGTNXacgtnx-]+"))
        {
            annot.append("Sequence contains invalid characters;");
            quality = QUALITY_LEVEL_5;
        }
        else if(initData.sm.getSetting("", SettingsManager.CODETYPE).equalsIgnoreCase(SettingsManager.CT_COMPLETE) &&
                !gr.getSequence().matches("[ACGTSWRYKMBVHDNUXacgtswrykmbvhdnux-]+"))
        {
            annot.append("Sequence contains invalid characters;");
            quality = QUALITY_LEVEL_5;
        }
        // Sequence length.
        int iLength = gr.getEnd()-gr.getStart()+1;
        if(gr.getSequence().length()!=iLength)
        {
            annot.append("Actual sequence length differs from the one given by the start/end;");
            quality = Math.max(quality, QUALITY_LEVEL_4);
        }
        // Assign the quality.
        gr.addProperty(QUALITY_LEVEL, quality);
        gr.addProperty(QUALITY_DESCRIPTION, annot.toString());
        return quality;
    }
}
