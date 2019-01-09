/*
    File:
        VMAParser.java
 *
    Revision:
        1.0.0.0
 *
    Description:
        Parses Vertical Multiple DPGPAlignment.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package dpgp;

import algorithms.SequenceBuffer;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class VMAParser
{
    public static class RawAlingment
    {
        public int iStart = 0;
        public int iEnd   = 0;
        public SequenceBuffer smel = null;
        public SequenceBuffer sseq = null;
    };

    private String strErrMsg    = "";

    /**
     *  Parses the VMA files and returns the resulting multiple alignment.
     *  Sites with quality below iThreshold are replaces by 'N' in other species.
     *
     *  Remarks:
     *      filenames, species and strains arrays must have the same length and cannot be null.
     *
     *  VMA format:
     *      Column      | Value
     *          1       | Position in D. mel. genome
     *          2       | D. mel. base
     *          3       | D. sim. base
     *          4       | Quality score
     *          ...     | Several irrelevant columns
     *      All columns are TAB-separated.
     *
     *
     *  @param files
     *  @param iThreshold
     *  @param species
     *  @param strains
     *  @return
     */
    public DPGPAlignment parseFiles(File[] files, int iThreshold, String[] species, String[] strains)
    {
        if(files==null || species==null || strains==null ||
           files.length!=species.length || species.length!=strains.length)
        {
            strErrMsg = "Arrays have unequal length";
            return null;
        }
        DPGPAlignment ali = new DPGPAlignment();
        for(int i=0;i<files.length;i++)
        {
            try
            {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
                RawAlingment ra = new RawAlingment();                
                String strLine;               
                ra.smel = new SequenceBuffer(5000000);      // Melanogaster sequence
                ra.sseq = new SequenceBuffer(5000000);      // Other sequence.
                // Read the data.
                int cp = 0;                                 // Current position with respect to the mel. genome.
                while((strLine = in.readLine())!=null)
                {
                    char base;    // D. sim. base.
                    // Read the line until the first blank.
                    int pos = 0;
                    while(strLine.charAt(pos)!=' ' && strLine.charAt(pos)!='\t')
                        pos++;
                    // If the part is a position, check whether there is a gap between genome
                    // positions.
                    if(strLine.charAt(0)!='I')
                    {
                        int ne = Integer.parseInt(strLine.substring(0, pos)); // New end position.
                        if(ra.iStart==0)                // Set alignment start.
                            ra.iStart = ne;
                        ra.iEnd = ne;                   // Update alignment end.
                        if(cp==0)
                            cp=ne;
                        else
                        {
                            // If there is missing data between the previously parsed site and the
                            // current site, insert missing data character.
                            if(ne-1!=cp)
                            {
                                int tmp = ne-cp;
                                ra.smel.appendBase(DPGPAlignment.MISSING, tmp-1);
                                ra.sseq.appendBase(DPGPAlignment.MISSING, tmp-1);
                            }
                        }
                        cp=ne;
                    }
                    // Skip the blanks until the next base.
                    while(strLine.charAt(pos)==' ' || strLine.charAt(pos)=='\t')
                        pos++;
                    // The current position is the D.mel. base - add it anyway.
                    ra.smel.appendBase(strLine.charAt(pos));
                    pos++;
                    // Skip the blanks until the next base.
                    while(strLine.charAt(pos)==' ' || strLine.charAt(pos)=='\t')
                        pos++;
                    // The current position is the base of D. sim.
                    base = strLine.charAt(pos);
                    pos++;
                    // Skip the blanks until the next column - Q-score.
                    while(strLine.charAt(pos)==' ' || strLine.charAt(pos)=='\t')
                        pos++;
                    int pss = pos; // Position: score start.
                    int l = strLine.length();
                    while(pos<l && strLine.charAt(pos)!=' ' && strLine.charAt(pos)!='\t')
                        pos++;
                    // Parse the score.
                    int iScore = Integer.parseInt(strLine.substring(pss, pos));
                    ra.sseq.appendBase((iScore>=iThreshold) ? base : DPGPAlignment.UNKNOWN);
                }
                ali.addAlignment(species[i], strains[i], ra);
            }
            catch(Exception e)
            {
                strErrMsg = "An error occured while parsing the files";
                return null;
            }
        }
        return ali;
    }

    /**
     *  Returns the text description of the last error occured.
     *
     *  @return
     */
    public String getLastErrorString()
    {
        return strErrMsg;
    }
}
