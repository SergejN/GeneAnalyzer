/*
    File:
        DPGPAlignment.java
 *
    Revision:
        1.2.0.0
 *
    Description:
        Represents the multiple sequence alignment derived from multiple VMA files.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package dpgp;

import algorithms.SequenceBuffer;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Vector;
import kernel.ErrorCode;


public class DPGPAlignment implements Serializable
{
    private static final long serialVersionUID = 1407200916451L;

    /*
        Following notation should be used in the alignment:
            - alignment gap:    there is a gap in some sequences but not in all
            - missing data:     no data is available for this site in any sequence
            - unknown:          the data at this site is poor, but it is not a gap and isn't missing

        Any class which uses this class should use these constants instead of common representation (e.g.
        '-' for a gap), otherwise the alignment might be incorrect.
    */
    public static final char GAP        = '-';
    public static final char MISSING    = '.';
    public static final char UNKNOWN    = 'N';
    public static final char NODATA     = 'x';


    // The first sequence in the vector is a melanogaster sequence.
    private SequenceBuffer[] seqs       = null;
    private Vector<String> species      = null;
    private Vector<String> strains      = null;
    private int iStart                  = 0;        // Left-most aligned position in the melanogaster genome
    private int iEnd                    = 0;        // Right-most aligned position in the melanogaster genom
    private Vector<Integer> index       = null;     // Position indices: position of every 1000th base of D.mel. genome
                                                    // in the alignment.

    public DPGPAlignment()
    {
        species = new Vector<String>();
        species.add("Drosophila melanogaster");
        strains = new Vector<String>();
        strains.add("D.mel");
        index = new Vector<Integer>(10000);
    }

    /**
     *  Returns the left-most and the right-most aligned positions with respect to
     *  the melanogaster genome.
     *
     *  @return
     */
    public int[] getBoundaries()
    {
        return new int[]{iStart, iEnd};
    }

    /**
     *  Extracts the specified region. iStart and iEnd are positions relative
     *  to the melanogaster genome. If one of the parameters is non-positive
     *  the method returns null. If iStart is smaller than iEnd, the method
     *  assumes, that the gene lies at the complementary strand and returns
     *  the reverse complement of the actual alignment sequences.
     *
     *  Remarks:
     *      If iStart is less than the start of the alignment, the beginning of
     *      the alignment is padded with gaps ('-');
     *      If iEnd is greater than the end of the alignment, the end of
     *      the alignment is padded with gaps ('-');
     *
     *  @param iStart
     *  @param iEnd
     *  @return
     */
    public SequenceBuffer[] extractRegion(int iStart, int iEnd)
    {
        if(iStart<0 || iEnd<0)
            return null;
        SequenceBuffer[] buf = new SequenceBuffer[seqs.length];
        // Case 1: region lies completely outside of the alignment.
        if( (iStart<this.iStart && iEnd<this.iStart) || (iStart>this.iEnd && iEnd>this.iEnd) )
        {
            int nBases = iEnd-iStart+1;
            for(int i=0;i<seqs.length;i++)
            {
                buf[i] = new SequenceBuffer(nBases);
                buf[i].appendBase(NODATA, nBases);
                return buf;
            }
        }
        for(int i=0;i<seqs.length;i++)
            buf[i] = new SequenceBuffer(Math.abs(iEnd-iStart)+1);
        int ls = Math.min(iStart, iEnd);    // Local start.
        int le = Math.max(iStart, iEnd);    // Local end.
        // If the beginning of the region does not overlap with the alignment
        //  => fill with NODATA.
        if(ls<this.iStart)
        {
            for(int i=0;i<seqs.length;i++)
                buf[i].appendBase(NODATA, this.iStart-ls);
            ls = 0;
        }
        else
            ls = genToAli(ls);
        // Copy overlapping part.
        if(le>this.iEnd)
            le = seqs[0].length()-1;
        else
            le = genToAli(le);
        for(int i=0;i<seqs.length;i++)
            buf[i].appendBases(seqs[i].subsequence(ls, le), 0, le-ls);
        // If the end of the region does not overlap with the alignment
        //  => fill with NODATA.
        le = Math.max(iStart, iEnd);
        if(le>this.iEnd)
            for(int i=0;i<seqs.length;i++)
                buf[i].appendBase(NODATA, le-this.iEnd);
        // Check, whether the sequence must be reverse complemented.
        if(iStart>iEnd)
            for(int i=0;i<seqs.length;i++)
                buf[i].reverseComplement();
        for(int i=0;i<seqs.length;i++)
            buf[i].replaceAll(MISSING, NODATA);
        return buf;
    }

    /**
     *  Extracts and returns the sequence of the single strain.
     *  If you intend to extract the same region of multiple strains,
     *  it is better to use extractRegion() returning the sequences array.
     *  If iEnd is less than iStart the method returns the reversed complement
     *  of the region sequence.
     *
     *  Remarks:
     *  If one of the following is true, the method returns null.
     *      - iStrainIndex is an invalid index
     *      - iStart or iEnd is negative
     *
     *  - If the region lies completely outside of the alignment, the resulting
     *  sequence buffer is filled with NODATA character ('x'):
     *
     *      |--------|   |----------------------------|
     *        region     |----------------------------|
     *                   |----------------------------|
     *                   |----------------------------|
     *                              alignment
     *
     *     |----------------------------|   |--------|
     *     |----------------------------|     region
     *     |----------------------------|
     *     |----------------------------|
     *                alignment
     *  
     *  - If only part of the region overlaps with the alignment, the non-overlapping
     *    part is filled with NODATA character ('x').
     *
     *  @param iStart
     *  @param iEnd
     *  @param iStrainIndex
     *  @return
     */
    public SequenceBuffer extractRegion(int iStart, int iEnd, int iStrainIndex)
    {
        if(iStrainIndex<0 || iStrainIndex>=seqs.length || iStart<0 || iEnd<0)
            return null;
        SequenceBuffer sb = null;
        // Case 1: region lies completely outside of the alignment.
        if( (iStart<this.iStart && iEnd<this.iStart) || (iStart>this.iEnd && iEnd>this.iEnd) )
        {
            int nBases = Math.abs(iEnd-iStart)+1;
            sb = new SequenceBuffer(nBases);
            sb.appendBase(NODATA, nBases);
            return sb;
        }
        sb = new SequenceBuffer(Math.abs(iEnd-iStart)+1);
        int ls = Math.min(iStart, iEnd);    // Local start.
        int le = Math.max(iStart, iEnd);    // Local end.
        // If the beginning of the region does not overlap with the alignment
        //  => fill with NODATA.        
        if(ls<this.iStart)
        {
            sb.appendBase(NODATA, this.iStart-ls);
            ls = 0;
        }
        else
            ls = genToAli(ls);
        // Copy overlapping part.
        if(le>this.iEnd)
            le = seqs[0].length()-1;
        else
            le = genToAli(le);
        sb.appendBases(seqs[iStrainIndex].subsequence(ls, le), 0, le-ls);
        // If the end of the region does not overlap with the alignment
        //  => fill with NODATA.
        le = Math.max(iStart, iEnd);
        if(le>this.iEnd)
            sb.appendBase(NODATA, le-this.iEnd);
        // Check, whether the sequence must be reverse complemented.
        if(iStart>iEnd)
            sb.reverseComplement();
        sb.replaceAll(MISSING, NODATA);
        return sb;
    }

    /**
     *  Returns the list of species names in the same order as the alignments were
     *  added with Drosophila melanogaster at the first position.
     *
     *  @return
     */
    public String[] getSpeciesNames()
    {
        return species.toArray(new String[1]);
    }

    /**
     *  Returns the list of strain names in the same order as the alignments were
     *  added with D.mel. at the first position.
     *
     *  @return
     */
    public String[] getStrainsNames()
    {
        return strains.toArray(new String[1]);
    }

    /**
     *  Returns the total number of strains in the alignment.
     *
     *  @return
     */
    public int getStrainsCount()
    {
        return seqs.length;
    }

    /**
     *  Saves the alignment to the binary file, so that it can be loaded and
     *  used later.
     *
     *  @param strFilename
     *  @return
     */
    public ErrorCode saveToBinaryFile(String strFilename)
    {
        try
        {
            ObjectOutputStream oout = new ObjectOutputStream (new FileOutputStream(strFilename));
            oout.writeObject(this);
            return ErrorCode.Ok;
        }
        catch (Exception e)
        {
            return ErrorCode.IOError;
        }
    }

    /**
     *  Loads the alignment object from the saved file and saves its data into
     *  the alignment object.
     *
     *  @param strFilename
     *  @param alignment
     *  @return
     */
    public static ErrorCode loadFromBinaryFile(String strFilename, DPGPAlignment alignment)
    {
        File f = new File(strFilename);
        if(!f.exists())
            return ErrorCode.FileDoesNotExist;
        try
        {
            ObjectInputStream oin = new ObjectInputStream(new FileInputStream(strFilename));
            DPGPAlignment ali = (DPGPAlignment)oin.readObject();
            alignment.iEnd = ali.iEnd;
            alignment.iStart = ali.iStart;
            alignment.seqs = ali.seqs;
            alignment.species = ali.species;
            alignment.strains = ali.strains;
            alignment.index   = ali.index;
            return ErrorCode.Ok;
        }
        catch(Exception e)
        {
            return ErrorCode.IOError;
        }
    }

    /**
     *  Saves the complete alignment into text file with specified number of characters per line.
     *  If the alignment is empty, the method returns ErrorCode.InvalidParameter.
     *
     *  @param strFilename
     *  @return
     */
    public ErrorCode saveToTextFile(String strFilename, int nBlockSize)
    {
        if(seqs==null)
            return ErrorCode.InvalidParameter;
        try
        {
            String strMissing = "\\"+Character.toString(MISSING);
            String strNodata = Character.toString(NODATA);
            // Construct filename.
            File output = new File(strFilename);
            PrintWriter out = new PrintWriter(new FileWriter(output));
            int pos = 0;
            int maxpos = seqs[0].length();
            int sp = iStart;                // Start position.
            while(pos<maxpos)
            {
                int np = sp; // New position.
                for(int i=0;i<seqs.length;i++)
                {
                    String strSeq = seqs[i].substring(pos, Math.min(pos+nBlockSize, maxpos));
                    // Calculate the genome position.
                    if(i==0)
                    {
                        for(int l=0;l<strSeq.length();l++)
                        {
                            char c = strSeq.charAt(l);
                            if(c!=GAP)
                                np++;
                        }
                    }
                    String strLine = String.format("%20s   %s",
                            strains.get(i)+"\t"+Integer.toString(sp),
                            strSeq.replaceAll(strMissing, strNodata));
                    out.println(strLine);
                }
                out.println();
                pos += nBlockSize;
                sp = np;
            }
            out.close();
            return ErrorCode.Ok;
        }
        catch (IOException e)
        {
            return ErrorCode.IOError;
        }
    }

    /**
     *  Adds new pairwise alignment to the multiple alignment and re-aligns all sequences.
     *
     *  @param strSpecies   species name
     *  @param strStrain    strain name
     *  @param ali          2 lines alignment
     */
    public void addAlignment(String strSpecies, String strStrain, VMAParser.RawAlingment ali)
    {
        strains.add(strStrain);
        species.add(strSpecies);
        // Add the first sequence.
        if(seqs==null)
        {
            seqs = new SequenceBuffer[2];
            seqs[0] = ali.smel.clone();
            seqs[1] = ali.sseq.clone();
            iStart = ali.iStart;
            iEnd = ali.iEnd;
        }
        else    // There is at least one pairwise alignment.
        {
            // Estimate the number of gaps in the existing multiple aligment.
            int nGapOld = seqs[0].length() - (iEnd-iStart+1);
            // Estimate the number of gaps in the new alignment.
            int nGapNew = ali.smel.length() - (ali.iEnd-ali.iStart+1);
            // The upper boundary of the length of the resulting multiple alignment
            // is the distance from left-most to right-most position plus the total
            // number of gaps in both old and new alignment, since a gap in either alignment
            // can yield the gap in another one, thus, increasing the total length.
            int iLength = (Math.max(iEnd, ali.iEnd)-Math.min(iStart, ali.iStart)+1)+nGapOld+nGapNew;
            SequenceBuffer[] tmp = new SequenceBuffer[seqs.length+1];
            for(int i=0;i<tmp.length;i++)
                tmp[i] = new SequenceBuffer(iLength);
            /*****************************************************************
            *                       Re-align the sequences                   *
            *****************************************************************/
            // Two alignments, new and old:
            //          |----------|
            //                  |-----------|
            int start = Math.min(iStart, ali.iStart);
            int pn = start - ali.iStart;
            int po = start - iStart;
            /***************************************************************
            *               DPGPAlignment beginning.                           *
            ***************************************************************/
            // If the new pairwise alignment starts before the old one, add bases
            // from the mel. sequence of the new alignment to the combined alignment.
            while(po<0 && pn<ali.smel.length())
            {
                tmp[0].appendBase(ali.smel.baseAt(pn));
                for(int i=1;i<tmp.length-1;i++)
                   tmp[i].appendBase(GAP);
                tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                if(ali.smel.baseAt(pn)!=GAP)
                    po++;
                pn++;
            }
            // If the new alignment starts beyond the old start, copy the bases
            // from the old alignment and put gaps into the last sequence.
            while(pn<0 && po<seqs[0].length())
            {
                for(int i=0;i<tmp.length-1;i++)
                   tmp[i].appendBase(seqs[i].baseAt(po));
                tmp[tmp.length-1].appendBase(GAP);
                if(seqs[0].baseAt(po)!=GAP)
                    pn++;
                po++;
            }
            // If the old and the new alignment do not have common sites:
            //      |-------|
            //                |---------|
            // then either pn or po will still be negative after the previous
            // while loop. In this case fill the missing part of the sequence
            // with gaps.
            // However, this SHOULD NORMALLY NOT be the case, so it's done
            // extra and not within the previous while loops.
            if(po<0)
            {
                for(int i=0; i<tmp.length;i++)
                    tmp[i].appendBase(MISSING, -po);
                po = 0;
            }
            if(pn<0)
            {
                for(int i=0; i<tmp.length;i++)
                    tmp[i].appendBase(MISSING, -pn);
                pn = 0;
            }
            /*****************************************************************
            *               Common part.                                     *
            *                                                                *
            *   There are 3 possible patterns of a site:                     *
            *       -  BB  both mel. and another sequence(s) have a base     *
            *       -  -B  gap in mel. sequence, and no gap in other seqs.   *
            *       -  ..  complete site is missing                          *
            *                                                                *
            *   Thus, there are 9 possible combinations of site patterns in  *
            *   the old and new alignment.                                   *
            *****************************************************************/
            while(pn<ali.smel.length() && po<seqs[0].length())
            {
                char mo = seqs[0].baseAt(po);   // Mel. base in old alignment.
                char mn = ali.smel.baseAt(pn);  // Mel. base in new alignment.
                // ----------------------------------------------------- //
                // Pattern #1: BB - BB
                if(mo>MISSING && mn>MISSING && mo==mn)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(seqs[i].baseAt(po));
                    tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                    po++;
                    pn++;
                }
                // Pattern #2: BB - ..
                else if(mo>MISSING && mn==MISSING)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(seqs[i].baseAt(po));
                    tmp[tmp.length-1].appendBase(MISSING);
                    po++;
                    pn++;
                }
                // Pattern #3: BB - -B
                else if(mo>MISSING && mn==GAP)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(GAP);
                    tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                    pn++;
                }
                // ----------------------------------------------------- //
                // Pattern #4: -B - BB
                else if(mo==GAP && mn>MISSING)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(seqs[i].baseAt(po));
                    tmp[tmp.length-1].appendBase(GAP);
                    po++;
                }
                // Pattern #5: -B - ..
                else if(mo==GAP && mn==MISSING)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(seqs[i].baseAt(po));
                    tmp[tmp.length-1].appendBase(MISSING);
                    po++;
                }
                // Pattern #6: -B - -B
                else if(mo==GAP && mn==GAP)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(seqs[i].baseAt(po));
                    tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                    po++;
                    pn++;
                }
                // ----------------------------------------------------- //
                // Pattern #7: .. - BB
                else if(mo==MISSING && mn>MISSING)
                {
                    tmp[0].appendBase(mn);
                    for(int i=1;i<tmp.length-1;i++)
                        tmp[i].appendBase(MISSING);
                    tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                    po++;
                    pn++;
                }
                // Pattern #8: .. - -B
                else if(mo==MISSING && mn==GAP)
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(GAP);
                    tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                    pn++;
                }
                // Pattern #9: .. - ..
                else if(mo==MISSING && mn==MISSING)
                {
                    for(int i=0;i<tmp.length;i++)
                        tmp[i].appendBase(MISSING);
                    po++;
                    pn++;
                }
                else // mn!=mo and the data is not missing ==> something wrong with the alignments.
                {
                    for(int i=0;i<tmp.length-1;i++)
                        tmp[i].appendBase(Character.toLowerCase(seqs[i].baseAt(po)));
                    tmp[tmp.length-1].appendBase(Character.toLowerCase(ali.sseq.baseAt(pn)));
                    pn++;
                    po++;
                }
            }
            /***************************************************************
            *               DPGPAlignment end.                                 *
            ***************************************************************/
            while(po<seqs[0].length())
            {
                for(int i=0;i<tmp.length-1;i++)
                   tmp[i].appendBase(seqs[i].baseAt(po));
                tmp[tmp.length-1].appendBase(GAP);
                pn++;
                po++;
            }
            // If the new alignment starts beyond the old start, copy the bases
            // from the old alignment and put gaps into the last sequence.
            while(pn<ali.sseq.length())
            {
                tmp[0].appendBase(ali.smel.baseAt(pn));
                for(int i=1;i<tmp.length-1;i++)
                   tmp[i].appendBase(GAP);
                tmp[tmp.length-1].appendBase(ali.sseq.baseAt(pn));
                pn++;
                po++;
            }
            seqs = tmp;
            iStart = Math.min(iStart, ali.iStart);
            iEnd = Math.max(iEnd, ali.iEnd);
        }
        createIndex();
    }

    /**
     *  Creates an index vector. Each value in the vector corresponds to a
     *  1000-bp step in the melanogaster sequence, i.e. the first value is
     *  the position of the 1000th site of D.mel. sequence in the sequence
     *  buffer (this must not be position 1000, since there might be gaps in
     *  between), the second value is the 2000th site and so on.
     */
    private void createIndex()
    {
        index.clear();
        index.add(0);
        int n=0;
        for(int i=0;i<seqs[0].length();i++)
        {
            if(seqs[0].baseAt(i)!=GAP)
                n++;
            if(n==1000)
            {
                index.add(i);
                n=0;
            }
        }
    }

    /**
     *  Converts the genomic position to the alignment position
     *  and returns the converted value.
     *
     *  @param iPosition
     *  @return
     */
    private int genToAli(int iPosition)
    {
        if(iPosition<iStart || iPosition>iEnd)
            return iPosition;
        // Find the 1000-bp block the position is in.
        int ind = (iPosition-this.iStart)/1000;
        int i = index.get(ind);          // Position in the alignment to start searching from.
        int pos = this.iStart+ind*1000-1;  // Position in the mel. genome the search starts at
        if(ind==0)
            pos++;
        while(pos<iPosition)
        {
            i++;
            if(seqs[0].baseAt(i)!=GAP)
                pos++;            
        }
        return i;
    }
}