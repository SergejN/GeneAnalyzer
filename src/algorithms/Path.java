/*
    File:
        Path.java
 *
    Revision:
        1.3.0.1
 *
    Description:
        Represents one evolutionary path from one codon to another.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

import bio.gene.dna.Codon;
import bio.gene.dna.ICodonTable;
import java.util.Vector;


public class Path
{
    public static final int MASK_SYNONYMOUS         = 0x00000001;
    public static final int MASK_NONSYNONYMOUS      = 0x00000010;
    public static final int MASK_TRANSITION         = 0x00000100;
    public static final int MASK_TRANSVERSION       = 0x00001000;

    
    private ICodonTable ct      = null;

    private Vector<Codon> path  = null;
    // Flags, specifying whether the codon is observed.
    private Vector<Boolean> flg = null;

    
    public Path(ICodonTable ct)
    {
        path = new Vector<Codon>();
        flg  = new Vector<Boolean>();
        this.ct = ct;
    }

    /**
     *  Tests whether the path contains a codon and returns
     *  its index in the path. If the codon is not present
     *  the method returns -1;
     *
     *  @param c
     *  @return
     */
    public int indexOf(Codon c)
    {
        return path.indexOf(c);
    }

    /**
     *  Adds the specified codon to the path. The flag bObserved
     *  specifies whether or not the codon is observed in the sample.
     *  If this flag is false the codon is assumed to be a transient
     *  codon. When calculating the number of polymorphic sites in the
     *  path only the observed codons are considered.
     *
     *  If codon is null, the method returns false, and the codon is not added.
     *
     *  @param codon
     *  @param bObserved
     *  @return
     */
    public boolean addCodon(Codon codon, boolean bObserved)
    {
        if(codon==null)
            return false;
        path.add(codon);
        flg.add(bObserved);
        return true;
    }

    /**
     *  Checks whether the path is contains all specified codons.
     *
     *  @param codons
     *  @return
     */
    public boolean containsAll(Codon[] codons)
    {
        for(Codon c:codons)
        {
            if(path.indexOf(c)==-1)
                return false;
        }
        return true;
    }

    /**
     *  Returns the length of the path.
     *
     *  @return
     */
    public int getLength()
    {
        return path.size();
    }

    /**
     *  Returns the codons the path consists of. The order of the codons
     *  is preserved in the array.
     *
     *  @return
     */
    public Codon[] getPath()
    {
        return path.toArray(new Codon[1]);
    }

    /**
     *  Returns the mask specifying which codons are observed and which are
     *  only transient. The length of the returned array is the same as the
     *  length of the array returned by the getPath method.
     *
     *  @return
     */
    public Boolean[] getObservedMask()
    {
        return flg.toArray(new Boolean[1]);
    }

    /**
     *  Returns the number of synonymous and nonsynonymous substitutions in the
     *  path.
     *  Array structure:
     *      index  | value
     *          0  | number of synonymous substitutions
     *          1  | number of nonsynonymous substitutions
     *
     *  @return
     */
    public int[] getPolymorphismsCount()
    {
        if(path.size()==0)
            return new int[]{0,0};
        // First, count the number of syn. changes.
        int nSyn = 0;
        for(int i=1;i<path.size();i++)
        {
            Codon prev = path.get(i-1);
            Codon cur  = path.get(i);
            if(ct.areSynonymous(prev.getSequence(), cur.getSequence()))
                nSyn++;
        }
        // There are one less substitutions in the path than there are codons.
        // The number of nonsyn. substitutions is the total number of substitutions
        // minus the number of syn. substitutions.
        int nNonsyn = path.size()-1-nSyn;
        return new int[]{nSyn, nNonsyn};
    }

    /**
     *  Calculates the number of syn. and nonsyn. transitions
     *  and transversions in the path.
     *  Array structure:
     *      index  | value
     *          0  | number of synonymous transitions
     *          1  | number of syn. transversions
     *          2  | number of nonsyn. transitions
     *          3  | number of nonsyn. transversions
     *
     *  @return
     */
    public int[] getSubstitutionsCount()
    {
        int[] res = new int[4];
        // Iterate through the path.
        int offset = 0;
        for(int i=1;i<path.size();i++)
        {
            String prev = path.get(i-1).getSequence();
            String cur  = path.get(i).getSequence();
            offset = (ct.areSynonymous(prev, cur)) ? 0 : 2;
            // Find the difference between the codons.
            SEARCHLOOP:for(int k=0;k<3;k++)
            {
                switch(BasicStatistics.getSubstitutionType(prev.charAt(k), cur.charAt(k)))
                {
                    // Invalid base - impossible
                    case -1:
                    // Bases are the same
                    case 0:
                        break;
                    // Transition.
                    case 1:
                        res[offset]++;
                        break SEARCHLOOP;
                    case 2:
                        res[offset+1]++;
                        break SEARCHLOOP;
                }
            }
        }
        return res;
    }

    /**
     *  Returns true only if the specified codon is not null,
     *  present in the path and contains a synonymous mutation with
     *  respect to the path.
     *
     *  Remarks:
     *      If possible, the codon is compared to both its successor and
     *      its ancestor. If any of them is synonymous to the specified
     *      codon, the method returns true.
     *
     *  @param c
     *  @return
     */
    public boolean isSubstitutionSynonymous(Codon c)
    {
        if(c==null || path.size()<2)
            return false;
        int idx = path.indexOf(c);
        if(idx==-1)
            return false;
        // The codon is the first codon in the path.
        if(idx==0)
        {
            Codon next = path.get(idx+1);
            return ct.areSynonymous(c.getSequence(), next.getSequence());
        }
        // The codon is the last codon in the path.
        else if(idx==path.size()-1)
        {
            Codon prev = path.get(idx-1);
            return ct.areSynonymous(c.getSequence(), prev.getSequence());
        }
        // The codon is somewhere in-between.
        else
        {
            // If the previous codon is syn. there is no need to check
            // the next one.
            Codon prev = path.get(idx-1);
            if(ct.areSynonymous(c.getSequence(), prev.getSequence()))
                return true;
            else
            {
                // If the previous codon is nonsyn. check whether the next
                // codon is syn.
                Codon next = path.get(idx+1);
                return ct.areSynonymous(c.getSequence(), next.getSequence());
            }
        }
    }

    /**
     *  Returns the type of the amino acid substitution caused by the specified
     *  base at the specified site.
     *
     *  Remarks:
     *      If possible, the algorithm finds both the substitution leading to
     *      the specified base and the one leading away from it. If any of those
     *      is synonymous, the method returns 1.
     *
     *  Returns: 
     *      0       the base is invalid
     *      -1      the specified position is monomorphic in the path
     *      -2      base does not appear in any codon of the path at the specified site
     *    If the base is ok, the following flags can be set:
     *      MASK_SYNONYMOUS       the substitution synonymous
     *      MASK_NONSYNONYMOUS    the substitution is nonsynonymous
     *      MASK_TRANSITION       the base substitution is a transition
     *      MASK_TRANSVERSION     the base substitution is a transversion
     *
     *  The high byte of the result contains the base, used to determine the type.
     *      
     *  @param base
     *  @param iSite
     *  @return
     */
    public int getSubstitutionType(char base, int iSite)
    {
        int nCodons = path.size();
        if(nCodons<2)
            return -1;
        base = Character.toUpperCase(base);
        char otherbase = 0;
        // Invalid base.
        if(base!='A' && base!='C' && base!='G' && base!='T')
            return 0;
        int type = 0; // Substitution type: transition or transversion.
        // First and last codon with the specified base.
        int first = -1;
        int last  = -1;        
        for(int i=0;i<nCodons;i++)
        {
            if(path.get(i).getSequence().charAt(iSite)==base)
            {
                if(first==-1)
                    first = i;
                last = i;
            }
        }
        // Base does not appear in the path.
        if(first==-1)
            return -2;
        // Specified site is monomorphic.
        if(first==0 && last==nCodons-1)
            return -1;
        // If the first codon with the specified base is the first codon in the path,
        // compare the last codon with the specified base with the codon following it.
        if(first==0)
        {
            String c1 = path.get(last).getSequence();
            String c2 = path.get(last+1).getSequence();
            otherbase = c2.charAt(iSite);
            type = BasicStatistics.getSubstitutionType(c1.charAt(iSite), otherbase);
            if(ct.areSynonymous(c1, c2))
            {
                int result = MASK_SYNONYMOUS;                
                if(type==1) // Transition.
                    return result | MASK_TRANSITION | ((int)otherbase<<24);
                else if(type==2)
                    return result | MASK_TRANSVERSION | ((int)otherbase<<24);
            }
        }
        // If the last codon with the specified base is the last codon in the path,
        // compare the first codon with the specified base with the previous codon.
        else if(last==nCodons-1)
        {
            String c1 = path.get(last).getSequence();
            String c2 = path.get(last+1).getSequence();
            otherbase = c2.charAt(iSite);
            type = BasicStatistics.getSubstitutionType(c1.charAt(iSite), otherbase);
            if(ct.areSynonymous(c1, c2))
            {
                int result = MASK_SYNONYMOUS;
                if(type==1) // Transition.
                    return result | MASK_TRANSITION | ((int)otherbase<<24);
                else if(type==2)
                    return result | MASK_TRANSVERSION | ((int)otherbase<<24);
            }
        }
        // Otherwise, if both, first and last codons, are in the middle of the
        // path, compare both to their neighbors.
        else
        {
            String f = path.get(first).getSequence();
            String l = path.get(last).getSequence();
            String other = path.get(first-1).getSequence();
            otherbase = other.charAt(iSite);
            type = BasicStatistics.getSubstitutionType(f.charAt(iSite), otherbase);
            if(ct.areSynonymous(f, other))
            {
                int result = MASK_SYNONYMOUS;
                if(type==1) // Transition.
                    return result | MASK_TRANSITION | ((int)otherbase<<24);
                else if(type==2)
                    return result | MASK_TRANSVERSION | ((int)otherbase<<24);
            }
            other = path.get(last+1).getSequence();
            otherbase = other.charAt(iSite);
            type = BasicStatistics.getSubstitutionType(l.charAt(iSite), otherbase);
            if(ct.areSynonymous(l, other))
            {
                int result = MASK_SYNONYMOUS;
                if(type==1) // Transition.
                    return result | MASK_TRANSITION | ((int)otherbase<<24);
                else if(type==2)
                    return result | MASK_TRANSVERSION | ((int)otherbase<<24);
            }
        }
        if(type==1)
            return MASK_NONSYNONYMOUS | MASK_TRANSITION | ((int)otherbase<<24);
        else if(type==2)
            return MASK_NONSYNONYMOUS | MASK_TRANSVERSION | ((int)otherbase<<24);
        else
            return 0;
    }

    /**
     *  Returns true if the path contains terminal codons.
     * 
     *  @return
     */
    public boolean containsTerminalCodons()
    {
        for(Codon c:path)
        {
            if(ct.isTerminal(c.getSequence()))
                return true;
        }
        return false;
    }

    @Override
    /**
     *  Returns the string representation of the path.
     *  If the codon is observable, the arrow pointing to is is '->'. If
     *  the codon is transient, the arrow is '=>'.
     */
    public String toString()
    {
        if(path.size()==0)
            return "";
        StringBuffer sb = new StringBuffer();
        // The first codon is always observable.
        sb.append(path.get(0).getSequence());
        // Iterate through the rest of the codons.
        int n = path.size();
        for(int i=1;i<n;i++)
        {
            sb.append((flg.get(i)) ? "->" : "=>");
            sb.append(path.get(i).getSequence());
        }
        return sb.toString();
    }

    @Override
    public Path clone()
    {
        Path p = new Path(ct);
        p.flg = (Vector<Boolean>)flg.clone();
        p.path = (Vector<Codon>)path.clone();
        return p;
    }

    /**
     *  Finds all possible evolutionary paths which can explain the observed
     *  codons pattern. Only the shortest possible paths are considered,
     *  so that there are at most as many steps as there are differences
     *  between the codons. If the codon array contains redundand codons
     *  the resulting path array contains redundand paths.
     *
     *  Note:
     *      The generated paths may contain terminal codons.
     *
     *  @param codons
     *  @param ct
     *  @return
     */
    public static Path[] generateAllPaths(Codon[] codons, ICodonTable ct)
    {
        // There can be at most as many mutations as there are differences
        // between the codons.
        int nSteps = calculateMutationsCount(codons);
        // Every codon can be the ancestral one. Iterate through the array.
        // If there are site-wise substitutions which appear in several different
        // codons, then such a substitution is only counted once, meaning that the
        // path will not be found, e.g.:
        //      CGG
        //      AAG
        //      CGA
        //      CAA
        // In this case iteratively increase the number of allowed mutations until
        // the path is found.
        // In the number of steps exceeds 64, then the path cannot be found for
        // some reason.
        Vector<Path> paths = new Vector<Path>();
        while(paths.size()==0 && nSteps<65)
        {
            for(Codon c:codons)
                generatePaths(c, new Path(ct), paths, codons, nSteps);
            nSteps++;
        }
        return paths.toArray(new Path[1]);
    }

    /**
     *  Finds the best evolutionary path which can explain the observed
     *  codons pattern. If two or more paths have the same number
     *  of nonsynonymous substitutions one of them is selected at random.
     *
     *  Remarks:
     *  The method returns null if one of the following is true:
     *      - codons is null or empty
     *      - ct is null
     *
     *  @param codons
     *  @param ct
     *  @param bUseTerminal whether or not to use paths with terminal codons
     *  @return
     */
    public static Path findBestPath(Codon[] codons, ICodonTable ct, boolean bUseTerminal)
    {
        if(codons==null || codons.length==0 || ct==null)
            return null;
        Path[] paths = generateAllPaths(codons, ct);
        Path bp = null;
        int i = Integer.MAX_VALUE;
        for(Path p:paths)
        {
            if(!bUseTerminal && p.containsTerminalCodons())
                continue;
            if(p.getPolymorphismsCount()[1]<i)
            {
                bp = p;
                i = p.getPolymorphismsCount()[1];
            }
        }
        return (bp!=null) ? bp : new Path(ct);
    }

    /**
     *  Finds the best evolutionary path from the specified codon to any of the
     *  codons from the codons array.
     *
     *  Remarks:
     *  The method returns null if one of the following is true:
     *      - codon is null
     *      - codons is null or empty
     *      - ct is null
     *
     *  @param codon
     *  @param codons
     *  @param ct
     *  @param bUseTerminal
     *  @return
     */
    public static Path findBestPath(Codon codon, Codon[] codons, ICodonTable ct, boolean bUseTerminal)
    {
        if(codon==null || codons==null || codons.length==0 || ct==null)
            return null;
        int nn = Integer.MAX_VALUE;
        Path bp = null;
        for(Codon c:codons)
        {
            Path p = Path.findBestPath(new Codon[]{codon, c}, ct, bUseTerminal);
            // There are two general rules, how the best path is determined:
            //      - the new path as many nonsyn. substitutions as the old one, but is shorter
            //      - the new path has less nonsyn. substitutions than the old one
            int tmp = p.getPolymorphismsCount()[1];
            if(bp==null || (tmp==nn && p.getLength()<bp.getLength()) || tmp<nn)
            {
                bp = p;
                nn = tmp;
            }
        }
        return bp;
    }


    /**
     *  Calculates the number of differences between the codons.
     *  This number corresponds to the required number of mutations
     *  which have to occur in order to generate all specified codons.
     *
     *  @param codons
     *  @return
     */
    private static int calculateMutationsCount(Codon[] codons)
    {
        int nMutations = 0;
        for(int j=0;j<3;j++)
        {
            SiteComposition sc = new SiteComposition();
            for(Codon c:codons)
                sc.addBase(c.getSequence().charAt(j));
            nMutations += sc.getNumberOfPolymorphisms();
        }
        return nMutations;
    }

    /**
     *  Searches for the complete evolutionary paths between the specified
     *  codons recursively.
     *
     *  @param ancestor     ancestral codon
     *  @param current      current path
     *  @param paths        already found complete paths
     *  @param codons       codons the path must include
     *  @param nSteps       maximal number of mutations/steps
     */
    private static void generatePaths(Codon ancestor, Path current, Vector<Path> paths, final Codon[] codons, int nSteps)
    {
        // Add the ancestral codon to the path. If the ancestral codon
        // is in the array of codons, then it is observable.
        boolean bObserved = false;
        for(Codon c:codons)
        {
            if(c==ancestor)
            {
                bObserved = true;
                break;
            }
        }
        current.addCodon(ancestor, bObserved);
        // If the maximal number of steps is reached, validate the current path
        // and add it to the current path.
        if(nSteps==0)
        {
            if(current.containsAll(codons))
                paths.add(current);
        }
        else
        {
            // Iterate through all neighbors of the ancestral codon and add them
            // to the path, if they are not already in.
            for(int i=0;i<ancestor.getNeighborsCount();i++)
            {
                Codon c = ancestor.getNeighbor(i);
                if(current.indexOf(c)==-1)
                    generatePaths(c, current.clone(), paths, codons, nSteps-1);
            }
        }
    }
}
