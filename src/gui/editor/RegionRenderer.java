/*
    File:
        RegionRenderer.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Renders a single gene region in the editor. 

 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.editor;

import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import javax.swing.JComponent;
import plugin.classes.IAligner;


public class RegionRenderer extends JComponent
{  
    // Show modes.
    public static final int SM_ALL      = 0;    // Show all
    public static final int SM_NO_IDENT = 1;    // Mask identical
    
    // Highlight modes.
    public static final int HM_REGIONS  = 0x00000001;    // Highlight regions
    
    public static final int CELLWIDTH      = 10;
    public static final int CELLHEIGHT     = 15;
    
    private final int ALPHA         = 75;
    private final Color EXON_EVEN   = new Color(0,255,0,ALPHA);
    private final Color EXON_ODD    = new Color(0,128,0,ALPHA);

    private final int COORDHEIGHT   = 20;
    private final int TICK_STEP_GENE= 10;
    private final int TICK_STEP_CDS = 9;
            
    private int displayMode   = SM_ALL;
    private int highlightMode = HM_REGIONS;

    private BufferedImage img = null;
    private GeneEntry ge      = null;
    private int iRegion       = 0;

    private int iRef          = 0;
    
    public RegionRenderer(GeneEntry ge, int iRegion)
    {        
        this.ge = ge;
        this.iRegion = iRegion;
        // Find the reference sequence.
        int nStrains = ge.getStrainsCount();
        for(int i=0;i<nStrains;i++)
        {
            String s = ge.getStrainEntry(i).getRegion(iRegion).getSequence();
            if(!s.matches("[NXnx-]+"))
            {
                iRef = i;
                break;
            }
        }
        paintRegion();
    }
    
    @Override
    public int getWidth()
    {
        return img.getWidth();
    }
    
    @Override
    public int getHeight()
    {
        return img.getHeight();
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
        this.displayMode = showMode;   
        paintRegion();
        repaint();
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
        this.highlightMode = highlightMask;
        paintRegion();
        repaint();
    }
    
    @Override
    public void paint(Graphics g)
    {
        g.drawImage(img, 0, 0, null);
    }
    
    private void paintRegion()
    {
        // Create image.
        int width = ge.getStrainEntry(0).getRegion(iRegion).getSequence().length()*CELLWIDTH;
        int height = ge.getStrainsCount()*CELLHEIGHT+2*COORDHEIGHT;
        img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = (Graphics2D)img.getGraphics();
        int nStrains = ge.getStrainsCount();
        // Draw background.
        GeneRegion reg = ge.getStrainEntry(0).getRegion(iRegion);
        String strType = reg.getType();        
        if(!strType.equalsIgnoreCase(GeneRegion.EXON))
        {
            if( (highlightMode & HM_REGIONS) >0 )
            {
                Color[] tmp = getColors(strType);
                GradientPaint gp = new GradientPaint(0, 0, tmp[1], 0, img.getHeight(), tmp[0]);
                g.setPaint(gp);
                g.fillRect(0, 0, img.getWidth(), img.getHeight());            
            }
        }
        else
        {
            // Calculate the frame of the first base of the exon.
            int iCdsLength = 0;
            StrainEntry se = ge.getStrainEntry(0);
            for(int i=0;i<iRegion;i++)
            {
                GeneRegion r = se.getRegion(i);
                if(r.hasType(GeneRegion.EXON))
                    iCdsLength += r.getSequence().length();
            }
            // Highlight the exon if necessary.
            if( (highlightMode & HM_REGIONS) >0 )
            {
                int w = 3*CELLWIDTH;
                int h = img.getHeight();
                BufferedImage patch = new BufferedImage(2*w, h, BufferedImage.TYPE_INT_ARGB);
                Graphics2D gi = (Graphics2D)patch.getGraphics();
                GradientPaint gp = new GradientPaint(0, 0, Color.WHITE, 0, img.getHeight(), EXON_ODD);
                gi.setPaint(gp);
                gi.fillRect(0, 0, w, h); 
                gp = new GradientPaint(0, 0, Color.WHITE, 0, img.getHeight(), EXON_EVEN);
                gi.setPaint(gp);
                gi.fillRect(w, 0, w, h);                 
                int iPos = iCdsLength%6;   // Position in the patch to paint first.
                // First, draw incomplete 6-base-blocks.
                if(iPos>0)
                {
                    g.drawImage(patch, 0, 0, (6-iPos)*CELLWIDTH, h,
                                       iPos*CELLWIDTH, 0, 2*w, h, null);
                }
                int x = ((6-iPos)%6)*CELLWIDTH;
                int nBases = se.getRegion(iRegion).getSequence().length();//-iPos;
                int nInterations = (int)Math.ceil((double)nBases/6.0);
                w = 2*w;
                for(int i=0;i<nInterations;i++)
                {
                    g.drawImage(patch, x, 0, x+w, h, 0, 0, w, h, null);
                    x += w;
                }
            }
            // Draw CDS ticks.
            g.setFont(new Font("Garamond", Font.PLAIN, 12));
            g.setColor(new Color(75,75,75));
            int l = reg.getSequence().length();
            for(int i=1;i<=l;i++)
            {
                if((iCdsLength+i)%TICK_STEP_CDS==0 || (i==1 && iCdsLength==0))
                {
                    String s = Integer.toString(iCdsLength+i);
                    int sw = g.getFontMetrics().bytesWidth(s.getBytes(), 0, s.length());
                    int cx = (i-1)*CELLWIDTH;
                    int cy = img.getHeight();
                    if((cx+(CELLWIDTH-sw)/2>=0) && (cx+(CELLWIDTH-sw)/2+sw<=img.getWidth()))
                    {
                        g.drawString(s, cx+(CELLWIDTH-sw)/2, cy-5);
                        g.drawString("'", cx+CELLWIDTH/2, cy-10);
                    }
                }
            }
        }
        // Draw frame.
        g.setColor(new Color(200,200,200));
        g.drawRect(0, 0, img.getWidth()-1, img.getHeight()-1);
        // Draw sequence.
        g.setFont(new Font("Courier New", Font.PLAIN, 16));
        g.setColor(Color.BLACK);
        if(displayMode==SM_ALL)
        {
            for(int i=0;i<nStrains;i++)
                g.drawString(ge.getStrainEntry(i).getRegion(iRegion).getSequence(), 1, (i+1)*CELLHEIGHT+COORDHEIGHT);
        }
        else
        {
            String strRef = ge.getStrainEntry(iRef).getRegion(iRegion).getSequence();
            for(int i=0;i<nStrains;i++)
            {
                if(i==iRef)
                    g.drawString(strRef, 1, (i+1)*CELLHEIGHT+COORDHEIGHT);
                else
                {
                    char[] tmp = ge.getStrainEntry(i).getRegion(iRegion).getSequence().toCharArray();
                    for(int j=0;j<tmp.length;j++)
                    {
                        if(tmp[j]==strRef.charAt(j) && tmp[j]!='-')
                            tmp[j] = '.';
                    }
                    g.drawString(new String(tmp), 1, (i+1)*CELLHEIGHT+COORDHEIGHT);
                }
            }
        }
        // Draw coordinates.
        g.setFont(new Font("Garamond", Font.PLAIN, 12));
        g.setColor(new Color(75,75,75));
        int cx = 0;
        for(int i=reg.getStart();i<=reg.getEnd();i++)
        {
            if(i%TICK_STEP_GENE==0 || i==1)
            {
                String s = Integer.toString(i);
                int sw = g.getFontMetrics().bytesWidth(s.getBytes(), 0, s.length());
                int cy = COORDHEIGHT;
                if((cx+(CELLWIDTH-sw)/2>=0) && (cx+(CELLWIDTH-sw)/2+sw<=img.getWidth()))
                {
                    g.drawString(s, cx+(CELLWIDTH-sw)/2, cy-5);
                    g.drawString("'", cx+CELLWIDTH/2, cy+8);
                }
            }
            cx += CELLWIDTH;
        }
    }
    
    /**
     *  Returns the start and end color of the gradient for each region type.
     *
     *  @param strType
     *  @return
     */
    private Color[] getColors(String strType)
    {
        if(strType.equalsIgnoreCase(GeneRegion.INTRON))
            return new Color[]{new Color(255,0,0,ALPHA), Color.WHITE};
        else if(strType.equalsIgnoreCase(GeneRegion.UTR5))
            return new Color[]{new Color(255,200,0,ALPHA), Color.WHITE};
        else if(strType.equalsIgnoreCase(GeneRegion.UTR3))
            return new Color[]{new Color(255,255,0,ALPHA), Color.WHITE};
        else if(strType.equalsIgnoreCase(GeneRegion.INTERGENIC))
            return new Color[]{new Color(128,128,128,ALPHA), Color.WHITE};
        else if(strType.equalsIgnoreCase(GeneRegion.mRNA))
            return new Color[]{new Color(220,30,225,ALPHA), Color.WHITE};
        else if(strType.equalsIgnoreCase(GeneRegion.UNNAMED))
            return new Color[]{new Color(225,225,225,ALPHA), Color.WHITE};
        else if(strType.equalsIgnoreCase(IAligner.AMBIGUOUS_REGION))
            return new Color[]{new Color(245,110,10,ALPHA), Color.WHITE};
        else
            return new Color[]{new Color(70,75,200,ALPHA), Color.WHITE};
    }
}
