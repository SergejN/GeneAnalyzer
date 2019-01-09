/*
    File:
        EditorPanel.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Scrollable editor panel.

 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package gui.editor;

import bio.gene.GeneEntry;
import java.awt.Dimension;
import java.awt.Rectangle;
import javax.swing.JPanel;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;

public class EditorPanel extends JPanel implements Scrollable
{
    // Show modes.
    public static final int SM_ALL      = RegionRenderer.SM_ALL;       // Show all
    public static final int SM_NO_IDENT = RegionRenderer.SM_NO_IDENT;  // Mask identical
    
    // Highlight modes.
    public static final int HM_REGIONS  = RegionRenderer.HM_REGIONS;   // Highlight regions
    
    private static final int OFFSET_TOP  = 10;
    
    private GeneEntry ge = null;

    /**
     *  Draws the alignment of the sequences of the gene entry and highlights
     *  different geen regions.
     *
     *  @param ge
     */
    public void displayAnnotation(GeneEntry ge)
    {
        if(ge==null || ge.getStrainsCount()==0 || ge.getStrainEntry(0).getRegionsCount()==0)
            return;
        this.ge = ge;
        update();
    }

    /**
     *  Updates the editor.
     */
    public void update()
    {
        this.removeAll();
        this.setLayout(null);
        int nRegs = ge.getStrainEntry(0).getRegionsCount();
        int x = 0;
        int h = 0;
        for(int i=0;i<nRegs;i++)
        {
            RegionRenderer rr = new RegionRenderer(ge, i);
            if(h==0)
                h=rr.getHeight();
            rr.setBounds(x, OFFSET_TOP, rr.getWidth(), h);
            x += rr.getWidth()-1;
            this.add(rr);
        }
        this.setPreferredSize(new Dimension(x+3*RegionRenderer.CELLWIDTH, OFFSET_TOP+h));
    }
    
    /**
     *  Sets the display mode.
     *  showMode can be one of the following:
     *      SM_ALL          if all bases should be displayed
     *      SM_NO_IDENT     if identical bases should be masked by '.'
     *   
     *  @param showMode
     */
    public void setDisplayMode(int showMode)
    {
        int nComp = getComponentCount();
        for(int i=0;i<nComp;i++)
            ((RegionRenderer)getComponent(i)).setDisplayMode(showMode);
        invalidate();
    }
    
    /**
     *  Sets the highlight mask.
     * 
     *  highlightMask can be the combination of the following flags (see Remarks):
     *       HM_REGIONS     if regions should be highlighted
     * 
     *  Remarks:
     *      Currently only the HM_REGIONS flag is supported, but new flags might
     *      be added in the later versions.
     * 
     *  @param highlightMask
     */
    public void setHighlightMask(int highlightMask)
    {
        int nComp = getComponentCount();
        for(int i=0;i<nComp;i++)
            ((RegionRenderer)getComponent(i)).setHighlightMask(highlightMask);
        invalidate();
    }
    
    public Dimension getPreferredScrollableViewportSize()
    {
        return getPreferredSize();
    }

    public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction)
    {
        if(orientation==SwingConstants.VERTICAL)
            return RegionRenderer.CELLHEIGHT;
        else
            return 3*RegionRenderer.CELLWIDTH;
    }

    public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction)
    {
        if(orientation==SwingConstants.VERTICAL)
            return RegionRenderer.CELLHEIGHT;
        else
            return 40*RegionRenderer.CELLWIDTH;
    }

    public boolean getScrollableTracksViewportWidth()
    {
        return false;
    }

    public boolean getScrollableTracksViewportHeight()
    {
        return false;
    }

}
