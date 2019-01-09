/*
    File:
        Alignment.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Represents the sequence alignment.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms.alignment;

import bio.gene.StrainEntry;


public class Alignment
{
    private StrainEntry[] strains = null;
    private float score           = 0.0f;
    private String strOptions     = null;
    private String strType        = null;


    /**
     *  Creates a new alignment. The strains are backed by the initial
     *  strains array, so that changes in the initial array would also
     *  affect the alignment. strType specifies the type of the aligner used.
     *  strOptions is a string containing the information about the options used 
     *  by the aligner.
     * 
     *  @param score
     *  @param strains
     *  @param strType 
     *  @param strOptions
     */
    public Alignment(float score, StrainEntry[] strains, String strType, String strOptions)
    {
        if(strains==null)
            throw new NullPointerException("Strains array cannot be null");
        this.strains = strains;
        this.score = score;
        this.strType = strType;
        this.strOptions = strOptions;
    }

    @Override
    public String toString()
    {
        return toString(0);
    }

    /**
     *  Returns the string representation of the alignment. The blockSize specifies
     *  the number of bases after which a new-line should be inserted.
     *  
     *  Note:
     *  A non-positive value of blockSize means no newlines are required.
     *
     *  @param blockSize
     *  @return
     */
    public String toString(int blockSize)
    {
        // If strains is an empty array, print the error message.
        if(strains.length==0)
            return "The alignment is empty";       
        int length = strains[0].getCompleteSequence().length();
        if(blockSize<1)
            blockSize = length;
        StringBuffer ali = new StringBuffer();        
        // Sequences.
        int maxlen = 0; // Maximal strain name length.
        String[] seqs = new String[strains.length];
        for(int i=0;i<strains.length;i++)
        {
            seqs[i] = strains[i].getCompleteSequence();
            if(strains[i].getStrainName().length()>maxlen)
                maxlen = strains[i].getStrainName().length();
        }
        // Match string.
        String ms = generateMatchString(seqs);
        // Maximal number of digits in the number of sites.
        int nDigits = 1+(int)(Math.floor(Math.log10(length)));
        int pos = 0;
        while(pos<length)
        {
            if(strains.length!=2)
            {
                for(int i=0;i<seqs.length;i++)
                {
                    String strSeq = seqs[i].substring(pos, Math.min(pos+blockSize, length));
                    String strLine = String.format("%-"+Integer.toString(maxlen)+"s   %-"+Integer.toString(nDigits)+"s %s\n",
                            strains[i].getStrainName(), Integer.toString(pos+1), strSeq);
                    ali.append(strLine);
                }
                if(strains.length>2)
                {
                    String strSeq = ms.substring(pos, Math.min(pos+blockSize, length));
                    ali.append(String.format("%"+Integer.toString(maxlen+nDigits+3)+"s %s\n", "", strSeq));
                }
            }
            else
            {
                String strSeq = seqs[0].substring(pos, Math.min(pos+blockSize, length));
                String strLine = String.format("%-"+Integer.toString(maxlen)+"s   %-"+Integer.toString(nDigits)+"s %s\n",
                                               strains[0].getStrainName(), Integer.toString(pos+1), strSeq);
                ali.append(strLine);
                strSeq = ms.substring(pos, Math.min(pos+blockSize, length));
                ali.append(String.format("%"+Integer.toString(maxlen+nDigits+3)+"s %s\n", "", strSeq));
                strSeq = seqs[1].substring(pos, Math.min(pos+blockSize, length));
                strLine = String.format("%-"+Integer.toString(maxlen)+"s   %-"+Integer.toString(nDigits)+"s %s\n",
                                        strains[1].getStrainName(), Integer.toString(pos+1), strSeq);
                ali.append(strLine);
            }
            ali.append("\n\n");
            pos += blockSize;
        }
        return ali.toString();
    }

    /**
     *  Returns the number of aligned sequences.
     *
     *  @return
     */
    public int getStrainsCount()
    {
        return strains.length;
    }

    /**
     *  Returns the specified strain entry or null if index is an invalid index.
     *
     *  @param index
     *  @return
     */
    public StrainEntry getStrainEntry(int index)
    {
        if(index<0 || index>=strains.length)
            return null;
        return strains[index];
    }
    
    /**
     *  Returns the alignment score.
     * 
     *  @return
     */
    public float getScore()
    {
        return score;
    }
    
    /**
     *  Returns the options the alignment was created with.
     * 
     *  @return
     */
    public String getOptions()
    {
        return strOptions;
    }
    
    /**
     *  Returns the alignment type.
     * 
     *  @return
     */
    public String getType()
    {
        return strType;
    }
    
    /**
     *  Returns the number of gaps in the alignment.
     * 
     *  @return
     */
    public int getGapsCount()
    {
        String[] seqs = new String[strains.length];
        for(int i=0;i<strains.length;i++)
            seqs[i] = strains[i].getCompleteSequence();
        int length = seqs[0].length();
        int nGaps = 0;
        MAINLOOP:for(int i=0;i<length;i++)
        {
            for(int n=0;n<seqs.length;n++)
            {
                if(seqs[n].charAt(i)=='-')
                {
                    nGaps++;
                    continue MAINLOOP;
                }
            }
        }
        return nGaps;
    }
    
    /**
     *  Returns the number of matches in the alignment.
     * 
     *  @return
     */
    public int getMatchesCount()
    {
        String[] seqs = new String[strains.length];
        for(int i=0;i<strains.length;i++)
            seqs[i] = strains[i].getCompleteSequence();
        int length = seqs[0].length();
        int nMatches = 0;
        MAINLOOP:for(int i=0;i<length;i++)
        {
            char ref = seqs[0].charAt(i);
            for(int n=1;n<seqs.length;n++)
            {
                if(seqs[n].charAt(i)!=ref)
                    continue MAINLOOP;
            }
            nMatches++;
        }
        return nMatches;
    }
    
    /**
     *  Returns the alignment length.
     * 
     *  @return
     */
    public int getLength()
    {
        if(strains.length>0)
            return strains[0].getCompleteSequence().length();
        else
            return 0;
    }

    
    private String generateMatchString(String[] seqs)
    {
        int length = seqs[0].length();
        char[] matches = new char[length];      // Match string.
        char mc = (seqs.length==2) ? '|' : '*'; // Match char.
        MAINLOOP:for(int i=0;i<length;i++)
        {
            char ref = seqs[0].charAt(i);
            if(ref=='-')
            {
                matches[i] = ' ';
                continue MAINLOOP;
            }
            for(int n=1;n<strains.length;n++)
            {
                // Mismatch.
                char c = seqs[n].charAt(i);
                if(c!=ref || c=='-')
                {
                    matches[i] = ' ';
                    continue MAINLOOP;
                }
            }
            // Match.
            matches[i] = mc;
        }
        return new String(matches);
    }
}
