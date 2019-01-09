/*
    File:
        GFFParser.java
 *
    Revision:
        1.0.0.0
 *
    Description:
        Parses GFF (General Feature Format) files.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package dpgp;

import bio.gene.GeneRegion;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class GFFParser
{
    private final String CDS = "GFFRealCDS";

    private String strErr = null;
    private String strWarn = null;

    /**
     *  Returns the text description of the last error or null if no error occured.
     *
     *  @return
     */
    public String getLastErrorString()
    {
        return strErr;
    }

    /**
     *  Returns the list of warnings if any.
     *
     *  @return
     */
    public String getWarnings()
    {
        return strWarn;
    }

    /**
     *  Parses the provided GFF file and returns the HashMap of gene fragments.
     *
     *  Remarks:
     *  The start position of the fragment MUST ALWAYS be less than or equal to
     *  the end position, whatever strand the gene is on. The strand should only
     *  be specified by the strand character (+ or -).
     *
     *  @param gffFile
     *  @param strChromosome 
     *  @return
     */
    public HashMap<String, Vector<DatasetBuilder.Fragment>> createFragments(File gffFile, String strChromosome)
    {
        // Check whether the file exists.
        if(!gffFile.exists())
        {
            strErr = String.format("%s does not exist.", gffFile.getAbsolutePath());
            return null;
        }
        strWarn = "";
        // Create gene entries.
        Pattern p = Pattern.compile("([^\\s]+)\\s+" +                  // 1st column: chromosome
                                    "[^\\s]+\\s+" +                    // 2nd column
                                    "([^\\s]+)\\s+" +                  // 3rd column: region
                                    "(\\d+)\\s+(\\d+)\\s+" +           // 4th and 5th columns: start and end
                                    "[^\\s]+\\s+" +                    // 6th column
                                    "([+-])\\s+" +                     // 7th column: strand
                                    "[^\\s]+\\s+" +                    // 8th column
                                    "ID=[^;]*(CG\\d+).*Parent=([^;]+)",   // 9th column: gene name
                                    Pattern.CASE_INSENSITIVE);
        // Hash map of precursor gene regions, which only specify the boundaries
        // and the region type.
        HashMap<String, Vector<DatasetBuilder.Fragment>> frags = new HashMap<String, Vector<DatasetBuilder.Fragment>>();
        try
        {            
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(gffFile)));
            String strLine;
            while((strLine = in.readLine())!=null)
            {
                Matcher m = p.matcher(strLine);
                if(m.find())
                {
                    boolean bLeading = m.group(5).charAt(0)=='+';
                    // Check the chromosome, since GFF can contain entries for different chromosomes.
                    if(!m.group(1).equalsIgnoreCase(strChromosome))
                        continue;
                    String strType = m.group(2);
                    int iStart = Integer.parseInt(m.group(3));
                    int iEnd = Integer.parseInt(m.group(4));
                    // Found out the type.
                    if(strType.equalsIgnoreCase("five_prime_UTR"))
                        strType = GeneRegion.UTR5;
                    else if(strType.equalsIgnoreCase("exon"))
                        strType = GeneRegion.EXON;
                    else if(strType.equalsIgnoreCase("intron"))
                        strType = GeneRegion.INTRON;
                    else if(strType.equalsIgnoreCase("three_prime_UTR"))
                        strType = GeneRegion.UTR3;
                    else if(strType.equalsIgnoreCase("CDS"))
                    {
                        strType = CDS;
                        if(bLeading)
                            iEnd += 3;
                        else
                            iStart -=3;
                    }
                    else                    // Unsupported region
                        continue;
                    // Create new gene region.
                    DatasetBuilder.Fragment fr = new DatasetBuilder.Fragment();
                    fr.strType = strType;
                    fr.iStart = Math.min(iStart, iEnd);
                    fr.iEnd = Math.max(iStart, iEnd);
                    fr.bLeading = bLeading;
                    // Splice variants list.
                    String[] sv = m.group(7).split("\\s*,\\s*");
                    Arrays.sort(sv);
                    // Iterate through the list of parents and add the region to each parent gene.
                    for(String s:sv)
                    {
                        String strParent = s.replaceAll("-R", "");
                        Vector<DatasetBuilder.Fragment> tmp = frags.get(strParent);
                        if(tmp==null)
                            tmp = new Vector<DatasetBuilder.Fragment>();
                        tmp.add(fr.clone());
                        frags.put(strParent, tmp);
                    }
                }
            }
        }
        catch(Exception e)
        {
            strErr = "An error occured while parsing the GFF file";
            return null;
        }
        processFragments(frags);
        return frags;
    }

    /**
     *  Iterates through the genes and adjusts the start position of the first
     *  exon if there is a 5' UTR, and the end position of the last exon, if there
     *  if a 3'UTR.
     *
     *  Remarks:
     *  The start and end of the coding sequence is also adjusted using the entry
     *  CDS of the original GFF file.
     *
     *  @param frags
     */
    private void processFragments(HashMap<String, Vector<DatasetBuilder.Fragment>> frags)
    {
        Iterator it = frags.keySet().iterator();
        while(it.hasNext())
        {
            String strGene = (String)it.next();
            Vector<DatasetBuilder.Fragment> fragments = frags.get(strGene);
            Vector<DatasetBuilder.Fragment> toRemove = new Vector<DatasetBuilder.Fragment>();
            int nFrags = fragments.size();            
            // Iterate through the regions and find overlapping exons and UTR's.
            DatasetBuilder.Fragment cds = null;
            for(int i=0;i<nFrags;i++)
            {
                DatasetBuilder.Fragment f1 = fragments.get(i);
                if(f1.strType.equalsIgnoreCase(GeneRegion.UTR5) || f1.strType.equalsIgnoreCase(GeneRegion.UTR3))
                {
                    for(int j=0;j<nFrags;j++)
                    {
                        DatasetBuilder.Fragment f2 = fragments.get(j);
                        if(f2.strType.equalsIgnoreCase(GeneRegion.EXON))
                        {
                            if(f2.iStart==f1.iStart)
                                f2.iStart=f1.iEnd+1;
                            if(f2.iEnd==f1.iEnd)
                                f2.iEnd=f1.iStart-1;
                            if(f2.iEnd<f2.iStart)
                                toRemove.add(f2);
                        }
                    }
                }
                else if(f1.strType.equalsIgnoreCase(CDS))
                {
                    cds = f1;
                    toRemove.add(f1);
                }
            }
            nFrags = toRemove.size();
            for(int i=0;i<nFrags;i++)
                fragments.remove(toRemove.get(i));
            // Adjust positions.
            if(cds==null)
                continue;
            // Find the indices of the left-most and the right-most exons.
            DatasetBuilder.Fragment lme = null;
            DatasetBuilder.Fragment rme = null;
            nFrags = fragments.size();
            for(int i=0;i<nFrags;i++)
            {
                DatasetBuilder.Fragment f = fragments.get(i);
                if(f.strType.equalsIgnoreCase(GeneRegion.EXON))
                {
                    if(lme==null || f.iStart<lme.iStart)
                        lme = f;
                    if(rme==null || f.iEnd>rme.iEnd)
                        rme = f;
                }
            }
            // The mRNA fragment BEFORE the left-most exon.
            if(cds.iStart>lme.iStart)
            {
                DatasetBuilder.Fragment f = new DatasetBuilder.Fragment();
                f.bLeading = lme.bLeading;
                f.iStart = lme.iStart;
                f.iEnd = cds.iStart-1;
                f.strType = GeneRegion.mRNA;
                lme.iStart = cds.iStart;
                fragments.add(f);
            }
            else if(cds.iStart<lme.iStart)
                strWarn += String.format("Gene %s: annotated CDS starts BEFORE the first exon;", strGene);
            // The mRNA fragment AFTER the right-most exon.
            if(cds.iEnd<rme.iEnd)
            {
                DatasetBuilder.Fragment f = new DatasetBuilder.Fragment();
                f.bLeading = rme.bLeading;
                f.iStart = cds.iEnd+1;
                f.iEnd = rme.iEnd;
                f.strType = GeneRegion.mRNA;
                rme.iEnd = cds.iEnd;
                fragments.add(f);
            }
            else if(cds.iEnd>rme.iEnd)
                strWarn += String.format("Gene %s: annotated CDS ends AFTER the last exon;", strGene);
        }
    }
}
