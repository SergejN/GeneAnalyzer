/*
    File:
        SequenceBuffer.java
 *
    Revision:
        1.1.0.0
 *
    Description:
        A mutable sequence of characters representing the gene sequence.
        A sequence buffer is like a StringBuffer, but it is not thread-safe,
        and, thus, should not be used by multiple threads. It is, however, faster
        than the StringBuffer, due to removed thread-safety.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

import java.io.Serializable;
import java.util.Arrays;

public class SequenceBuffer implements Serializable
{
    private static final long serialVersionUID  = 1407200916451L;
    private final static int DEFAULT_CAPACITY   = 500;
    
    private int nBases = 0;         // Number of bases stored in the buffer.
    private byte[] bases = null;    // Array of bases.

    /**
     *  Creates a sequence buffer with the default initial number of bases (500).
     */
    public SequenceBuffer()
    {
        this(DEFAULT_CAPACITY);
    }

    /**
     *  Creates a sequence buffer with the specified initial number of bases.
     *
     *  @param nBases
     */
    public SequenceBuffer(int nBases)
    {
        bases = new byte[nBases];
    }

    /**
     *  Creates a sequence buffer containing the specified string.
     *
     *  @param strString
     */
    public SequenceBuffer(String strString)
    {
        bases = strString.getBytes();
        nBases = strString.length();
    }

    /**
     *  Resizes the buffer to the desired size. All previously saved bases
     *  are transfered to the new buffer. If the new size is equal or less
     *  than the current size, the method simply returns.
     * 
     *  @param nBases
     */
    public void ensureCapacity(int nBases)
    {
        if(nBases<=bases.length)
            return;
        bases = Arrays.copyOf(bases, nBases);
    }

    /**
     *  Clears the buffer.
     */
    public void clearBuffer()
    {
        nBases = 0;
    }

    /**
     *  Appends a new base to the end of the sequence.
     *
     *  Remarks:
     *      The base is automatically converted to upper-case if possible.
     *
     *  @param base
     */
    public void appendBase(char base)
    {
        if(nBases==bases.length)
            bases = Arrays.copyOf(bases, 2*bases.length);
        bases[nBases] = (byte)base;
        nBases++;
    }

    /**
     *  Appends the new base to the end of the sequence nCount times.
     *  This method is faster, than calling appendBase(base) nCount times.
     *
     *  Remarks:
     *      The base is automatically converted to upper-case if possible.
     * 
     *  @param base
     *  @param nCount
     */
    public void appendBase(char base, int nCount)
    {
        // Calculate new length.
        int inl = nBases+nCount;
        if(inl>bases.length)
            bases = Arrays.copyOf(bases, 2*inl);
        Arrays.fill(bases, nBases, nBases+nCount, (byte)base);
        nBases += nCount;
    }

    /**
     *  Appends the specified part of another sequence buffer including iStart
     *  and iEnd positions to the sequence buffer. Throws exceptions if the
     *  start and end positions are out of bounds of if src is null.
     *
     *  @param src
     *  @param iStart
     *  @param iEnd
     */
    public void appendBases(SequenceBuffer src, int iStart, int iEnd)
    {
        int nCount = iEnd-iStart+1;
        int inl = nBases+nCount;
        if(inl>bases.length)
            bases = Arrays.copyOf(bases, 2*inl);
        System.arraycopy(src.bases, iStart, bases, nBases, nCount);
        nBases += nCount;
    }

    /**
     *  Returns the base at the specified position.
     *
     *  Throws IndexOutOfBoundsException if iSite is less than 0 or greater
     *  than or equal to the sequence length.
     *
     *  @param iSite
     *  @return
     */
    public char baseAt(int iSite)
    {
        if(iSite<0 || iSite>=nBases)
            throw new IndexOutOfBoundsException();
        else
            return (char)bases[iSite];
    }

    /**
     *  Sets the base at the specified position.
     *
     *  Throws IndexOutOfBoundsException if iSite is less than 0 or greater
     *  than or equal to the sequence length.
     *
     *  Remarks:
     *      The base is automatically converted to upper-case if possible.
     *
     *  @param iSite
     *  @return
     */
    public void setBaseAt(int iSite, char base)
    {
        if(iSite<0 || iSite>=nBases)
            throw new IndexOutOfBoundsException();
        else
            bases[iSite] = (byte)base;
    }

    /**
     *  Removes the specified site.
     *
     *  Throws IndexOutOfBoundsException if iSite is less than 0 or greater
     *  than or equal to the sequence length.
     *
     *  @param iSite
     */
    public void removeBase(int iSite)
    {
        if(iSite<0 || iSite>=nBases)
            throw new IndexOutOfBoundsException();
        else
        {
            byte[] tmp = new byte[bases.length];
            if(iSite>0)
                System.arraycopy(bases, 0, tmp, 0, iSite);
            System.arraycopy(bases, iSite+1, tmp, iSite, nBases-iSite-1);
            nBases--;
            bases = tmp;
        }
    }

    /**
     *  Removes the block of bases beginning at the position iSite and returns the
     *  number of removed bases.
     *
     *  Throws IndexOutOfBoundsException if iStart is less than 0 or greater
     *  than or equal to the sequence length. If iStart+nCount exceeds the sequence
     *  length, the method effectively creates a subsequence from 0 to iStart-1.
     *  If nCount is non-positive, the method does nothing and returns 0.
     *
     *  @param iStart
     *  @param nCount
     *  @return
     */
    public int removeBases(int iStart, int nCount)
    {
        if(nCount<1)
            return 0;
        if(iStart<0 || iStart>=nBases)
            throw new IndexOutOfBoundsException();
        if(iStart+nCount>nBases)
            nCount = nBases-iStart;
        byte[] tmp = new byte[bases.length];
        if(iStart>0)
            System.arraycopy(bases, 0, tmp, 0, iStart);        
        System.arraycopy(bases, iStart+nCount, tmp, iStart, nBases-iStart-nCount);
        nBases = nBases-nCount;
        bases = tmp;
        return nCount;
    }

    /**
     *  Inserts the base at the specified site. If the site is invalid the method 
     *  does nothing. The base is converted into uppercase.
     *
     *  Throws IndexOutOfBoundsException if iStart is less than 0 or greater
     *  than the sequence length.
     *
     *  @param iSite
     *  @param base
     */
    public void insertBase(int iSite, char base)
    {
        if(iSite<0 || iSite>nBases)
            throw new IndexOutOfBoundsException();
        byte[] tmp = null;
        if(nBases+1>=bases.length)
            tmp = new byte[bases.length*2];
        else
            tmp = new byte[bases.length];
        if(iSite>0)
            System.arraycopy(bases, 0, tmp, 0, iSite);
        tmp[iSite] = (byte)Character.toUpperCase(base);
        System.arraycopy(bases, iSite, tmp, iSite+1, nBases-iSite);
        nBases++;
        bases = tmp;
    }

    /**
     *  Returns the sequence length.
     *
     *  @return
     */
    public int length()
    {
        return nBases;
    }

    /**
     *  Reverses the sequence IN PLACE and returns the reference to this object.
     *
     *  @return
     */
    public SequenceBuffer reverse()
    {
        int lim = (nBases+1)/2;
        for(int i=0;i<lim;i++)
        {
            byte c = bases[i];
            bases[i] = bases[nBases-1-i];
            bases[nBases-1-i] = c;
        }
        return this;
    }

    /**
     *  Replaces each base in the sequence by its complement IN PLACE and returns
     *  the reference to this object.
     *
     *  Remarks:
     *      Only "real" bases (A, C, G and T) can be complemented.
     *
     * @return
     */
    public SequenceBuffer complement()
    {
        for(int i=0;i<nBases;i++)
        {
            switch(bases[i])
            {
                case 'A': bases[i] = 'T'; break;
                case 'C': bases[i] = 'G'; break;
                case 'G': bases[i] = 'C'; break;
                case 'T': bases[i] = 'A'; break;
                case 'a': bases[i] = 't'; break;
                case 'c': bases[i] = 'g'; break;
                case 'g': bases[i] = 'c'; break;
                case 't': bases[i] = 'a'; break;
            }
        }
        return this;
    }

    /**
     *  Reverses the sequence and replaces each base by its complement IN PLACE
     *  and returns the reference to this object. The effect of calling this method
     *  is the same as of calling complement() and then reverse().
     * 
     *  @return
     */
    public SequenceBuffer reverseComplement()
    {
        int lim = (nBases+1)/2;
        for(int i=0;i<lim;i++)
        {
            byte c1;
            switch(bases[i])
            {
                case 'A': c1 = 'T'; break;
                case 'C': c1 = 'G'; break;
                case 'G': c1 = 'C'; break;
                case 'T': c1 = 'A'; break;
                case 'a': c1 = 't'; break;
                case 'c': c1 = 'g'; break;
                case 'g': c1 = 'c'; break;
                case 't': c1 = 'a'; break;
                default: c1 = bases[i]; break;
            }
            byte c2;
            switch(bases[nBases-1-i])
            {
                case 'A': c2 = 'T'; break;
                case 'C': c2 = 'G'; break;
                case 'G': c2 = 'C'; break;
                case 'T': c2 = 'A'; break;
                case 'a': c2 = 't'; break;
                case 'c': c2 = 'g'; break;
                case 'g': c2 = 'c'; break;
                case 't': c2 = 'a'; break;
                default: c2 = bases[nBases-1-i]; break;
            }
            bases[i] = c2;
            bases[nBases-1-i] = c1;
        }
        return this;
    }

    @Override
    /**
     *  Returns the string representation of the sequence.
     */
    public String toString()
    {
        return new String(bases, 0, nBases);
    }

    /**
     *  Returns the specified substring.
     *
     *  @param beginIndex
     *  @param endIndex
     *  @return
     */
    public String substring(int beginIndex, int endIndex)
    {
        return new String(bases, beginIndex, endIndex-beginIndex);
    }

    /**
     *  Returns the sequence buffer containing the specified part of the
     *  original sequence.
     *
     *  @param beginIndex
     *  @param endIndex
     *  @return
     */
    public SequenceBuffer subsequence(int beginIndex, int endIndex)
    {
        SequenceBuffer sb = new SequenceBuffer(endIndex-beginIndex+1);
        sb.appendBases(this, beginIndex, endIndex);
        return sb;
    }

    public SequenceBuffer replaceAll(char orig, char rep)
    {
        for(int i=0;i<nBases;i++)
        {
            if(bases[i]==orig)
                bases[i]=(byte)rep;
        }
        return this;
    }

    /**
     *  Returns the number of times the character c appears in the sequence.
     * 
     *  @param c
     *  @return
     */
    public int getCharacterCount(char c)
    {
        return getCharacterCount(c, 0, nBases-1);
    }
    
    /**
     *  Returns the number of times the charater c appears in the sequence
     *  between positions from and to.
     * 
     *  @param c
     *  @param from
     *  @param to
     *  @return
     */
    public int getCharacterCount(char c, int from, int to)
    {
        if(from<0 || to<0)
            return 0;
        if(to>=nBases)
            to = nBases-1;
        int n = 0;
        for(int i=from;i<=to;i++)
        {
            if(bases[i]==(byte)c)
                n++;
        }
        return n;
    }

    /**
     *  Returns true if the sequence buffer is empty and false otherwise.
     *
     *  @return
     */
    public boolean isEmpty()
    {
        return nBases==0;
    }
    
    @Override
    /**
     *  Clones the sequence buffer.
     */
    public SequenceBuffer clone()
    {
        SequenceBuffer sb = new SequenceBuffer(nBases);
        sb.bases = Arrays.copyOf(bases, nBases);
        sb.nBases = nBases;
        return sb;
    }
}
