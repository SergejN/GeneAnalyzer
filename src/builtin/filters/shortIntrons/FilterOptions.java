/*
    File:
        FilterOptions.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Allows the user to select the min. and max. intron length and the
        population(s) to check the length for.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.filters.shortIntrons;


public class FilterOptions 
{
    String[] pops = null;
    int      imax = Integer.MAX_VALUE;
    int      imin = 0;
}
