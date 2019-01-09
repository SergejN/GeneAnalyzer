/*
    File:
        IAligner.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        The interface IAligner must be implemented by every aligner.
        The name of the main class of the plugin, which implements this interface
        should be AlignerMain.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin.classes;

import algorithms.alignment.Alignment;
import bio.gene.StrainEntry;
import kernel.ErrorCode;
import plugin.AInitData;

public interface IAligner
{
    public static final String PLUGIN_EXTENTION = ".jar";

    public static final String AMBIGUOUS_REGION = "Ambiguous";


    /**
     *  Initializes the aligner.
     *
     *  @param initData
     *  @return
     */
    public ErrorCode initialize(AInitData initData);

    /**
     *  Returns the aligner/algorithm name.
     * 
     *  @return
     */
    public String getName();

    /**
     *  Returns the description of the aligner.
     *
     *  @return
     */
    public String getDescription();

    /**
     *  Returns the number of sequences which can be aligned by the aligner.
     *  All aligners must support at least 2 sequences.
     * 
     *  @return
     */
    public int getSupportedSequenceNumber();

    /**
     *  Returns the template of the parameter string.
     *
     *  @return
     */
    public String getParamString();

    /**
     *  Returns the description of the last error occured.
     *
     *  @return
     */
    public String getLastError();
    
    /**
     *  Aligns the complete sequences of multiple strain entries and then splits the
     *  alignment into regions using the strain entries as reference. The method
     *  returns the Alignment object, containing the aligned strain entries.
     *  If entries is null or empty, the method should return null.
     *  The size of the entries array should be the same as returned by the method
     *  getSupportedSequenceNumber. The class implementing the interface is, however,
     *  free to decide, how to handle the case, if the size of the entries array
     *  is not equal to the number returned by getSupportedSequenceNumber. The method
     *  should also try to annotate the regions to the new alignment. The regions, in which
     *  the region type cannot be determined due to annotation differences in the initial
     *  strain entries, should be annotated as AMBIGUOUS.
     *
     *  If this method is called for the first time, the class should save the settings
     *  for the later use. However, if strParams is null in a subsequent call, the
     *  options dialog should be shown again, if any. If strParams is empty but not null
     *  in the subsequent call, the class should use the previously saved parameters. However,
     *  if the strParam is not empty, the class should use the parameters specified in the
     *  parameter string if possible.
     *
     *  @param entries
     *  @param strParams
     *  @return
     */
    public Alignment alignStrains(StrainEntry[] entries, String strParams);
}