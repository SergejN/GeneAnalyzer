/*
    File:
        SequenceRoutines.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Basic routines to work with DNA sequences.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

public class SequenceRoutines 
{
    /**
     *  Iterates through the sample site-wise and removes the sites with
     *  a gap in any sequence from the sample.
     * 
     *  @param sequences    sequences to upgap. The sequences must be of the
     *                      same length.
     *  @return             the array of ungapped sequences of the equal length.
     * 
     *  Note:
     *      If the sequences are not of equal length, the method can
     *      throw an IndexOutOfBounds exception.
     *      If sequences is null, the method throws a NullPointerException.
     *      If sequences is an empty array, the method returns an array of
     *      empty strings.
     */
    public static String[] getUngappedSequences(String[] sequences)
    {
        SequenceBuffer[] seqs = new SequenceBuffer[sequences.length];
        // Initialize the sequences.
        for(int i=0;i<seqs.length;i++)
            seqs[i] = new SequenceBuffer();
        int iLength = sequences[0].length();
        LENGTHLOOP: for(int i=0;i<iLength;i++)
        {
            for(int j=0;j<sequences.length;j++)
            {
                if(sequences[j].charAt(i)=='-')
                    continue LENGTHLOOP;
            }
            // If no gap found, add the site.
            for(int j=0;j<seqs.length;j++)
                seqs[j].appendBase(sequences[j].charAt(i));
        }
        // Result.
        String[] res = new String[seqs.length];
        for(int i=0;i<res.length;i++)
            res[i] = seqs[i].toString();
        return res;
    }
    
    /**
     *  Given aligned ungapped sequences, the method returns aligned specified
     *  sites. If sites is null or an empty array, the method returns the
     *  initial sequences.
     *  If ungappedSequences is null or an empty array, the initial sequences
     *  are returned.
     *  
     *  @param ungappedSeqs
     *  @param sites
     *  @return
     */
    public static String[] extractAlignedSites(String[] ungappedSeqs, int[] sites)
    {
        if(sites==null || sites.length==0 || ungappedSeqs==null || ungappedSeqs.length==0)
            return ungappedSeqs;
        SequenceBuffer[] seqs = new SequenceBuffer[ungappedSeqs.length];
        for(int i=0;i<seqs.length;i++)
            seqs[i] = new SequenceBuffer();
        // Iterate through the sites.
        for(int i:sites)
        {
            // If the position is within the sequence, extract it.
            if(i>0 && i<ungappedSeqs[0].length())
            {
                for(int n=0;n<seqs.length;n++)
                {
                    seqs[n].appendBase(ungappedSeqs[n].charAt(i));
                }
            }
        }
        String[] result = new String[seqs.length];
        for(int i=0;i<seqs.length;i++)
            result[i] = seqs[i].toString();
        return result;
    }
    
    /**
     *  Returns the minimal and the maximal sequence length with all
     *  gaps removed.
     * 
     *  @param seqs
     *  @return
     */
    public static int[] getMinMaxUngappedLength(String[] seqs)
    {
        int iMin = Integer.MAX_VALUE;
        int iMax = 0;
        for(String s:seqs)
        {
            int i = s.replaceAll("-", "").length();
            if(i<iMin)
                iMin = i;
            if(i>iMax)
                iMax = i;
        }
        return new int[]{iMin, iMax};
    }

    /**
     *  Formats the sequence by inserting a new line every nBlockLength characters.
     *
     *  @param strSequence
     *  @param nBlockLength
     *  @return
     */
    public static String getFormattedSequence(String strSequence, int nBlockLength)
    {
        StringBuffer sb = new StringBuffer();
        int pos = 0;
        while(pos<strSequence.length())
        {
            String tmp = strSequence.substring(pos, Math.min(pos+nBlockLength, strSequence.length()));
            sb.append(tmp);
            sb.append("\n");
            pos += nBlockLength;
        }
        return sb.toString();
    }

    /**
     *  Returns the string consisting of the strRepeat strings concatinated nCount times.
     *
     *  @param strRepeat
     *  @param nCount
     *  @return
     */
    public static String generateRepetetiveSequence(String strRepeat, int nCount)
    {
        if(strRepeat==null || strRepeat.isEmpty() || nCount<1)
            return "";
        int nLength = strRepeat.length();
        char[] buf = new char[nCount*nLength];
        char[] pattern = strRepeat.toCharArray();
        for(int i=0;i<nCount;i++)
            System.arraycopy(pattern, 0, buf, i*nLength, nLength);
        return new String(buf);
    }
}
