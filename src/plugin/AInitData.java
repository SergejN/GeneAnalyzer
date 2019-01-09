/*
    File:
        AInitData.java
 *
    Revision:
        2.3.0.1
 *
    Description:
        Plugin initialization data.
        Further versions of GeneAnalyzer might use some descendants of this
        base class. It is a good practice to first try to cast the object to
        the desired class, if you expect to get a descendant of the base class.

 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin;


import bio.gene.dna.ICodonTable;
import gui.IWaitDialog;
import java.util.Locale;
import kernel.ErrorCode;
import kernel.SettingsManager;
import plugin.classes.IAligner;

public abstract class AInitData
{
    public SettingsManager          sm          = null;
    public Locale                   locale      = null;
    public ICodonTable              ct          = null;
    // Wait dialog. The decision whether or not to display the wait dialog
    // is made by the host-application. E.g. the wait dialog is not displayed
    // when running the script or during initialization. However, this field
    // is quaranteed to not to be null.
    public IWaitDialog              wd          = null;
    public String[]                 codontables = null;
    public IAligner[]               aligners    = null;


    /**
     *  Sets the codon table.
     *
     *  @param strName      specifies the name of the table to use or can be
     *                      one of the following constants:
     *                          Initial - sets the codon table to the initial value
     *                                    which can be changed in Tools->Options
     *                          Default - sets the codon table to the universal
     *                                    codon table
     *  @return
     */
    public abstract ErrorCode setCodonTable(String strName);
}
