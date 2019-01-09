/*
    File:
        AnnotationRenderer.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        In the GeneAnnotationView AnnotationRenderer renders the relative location
        of the gene regions by representing them by the rectangles. If the regions
        overlap, the overlapping region is painted under the region it overlaps with.

        |----||----||----------|
            |------|

            ^ overlap

 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.views;

import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import javax.swing.JPanel;
import plugin.classes.IAligner;

public class AnnotationRenderer extends JPanel
{
    private static final int OFFSET = 15;
    private static final int RECT_HEIGHT = 30;

    private final Color UTR5        = new Color(255,150,0);
    private final Color UTR3        = new Color(255,255,0);
    private final Color mRNA        = new Color(220,30,225);
    private final Color EXON        = new Color(0,255,0);
    private final Color INTRON      = new Color(255,0,0);
    private final Color INTERGENIC  = new Color(128,128,128);
    private final Color UNNAMED     = new Color(225,225,225);
    private final Color OTHER       = new Color(55,55,230);
    private final Color AMBIGUOUS   = new Color(245,110,10);

    private BufferedImage bim = null;

    public AnnotationRenderer()
    {
        super();        
    }

    @Override
    public void paint(Graphics g)
    {
        Color tmp = g.getColor();
        g.setColor(Color.LIGHT_GRAY);
        g.clearRect(0, 0, getParent().getWidth(), getParent().getHeight());
        g.drawImage(bim, OFFSET, OFFSET, null);
        g.setColor(tmp);
    }

    /**
     *  Resets the renderer and clears the image.
     */
    public void reset()
    {
        bim = null;
        repaint();
    }

    public void paintRegions(StrainEntry se)
    {
        int w = getParent().getWidth()-2*OFFSET;
        int h = getParent().getHeight()-2*OFFSET;
        if(w<1 || h<1 || se==null)
            return;
        bim = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = (Graphics2D)bim.getGraphics();
        int y = paintLegend(g);
        int[] xs = new int[se.getRegionsCount()];       // X-coordinates of the regions.
        int[] ws = new int[se.getRegionsCount()];       // Widths of the regions.
        int[] ls = new int[se.getRegionsCount()];       // Levels of the regions.
        byte[] rlm = new byte[se.getRegionsCount()];    // Region level mask.
        int iTotalLength = se.getCompleteSequence().length();
        for(int i=0;i<se.getRegionsCount();i++)
        {
            GeneRegion reg = se.getRegion(i);
            Color c = getColor(reg.getType());
            int iLength = reg.getEnd()-reg.getStart()+1;
            // Calculate region rectangle width.
            ws[i] = (int)((float)iLength/(float)iTotalLength*(bim.getWidth()-1));
            // Count overlaps.
            // Since the regions are sorted already, only consider the previously painted regions.
            Arrays.fill(rlm, (byte)0);
            for(int j=i-1;j>=0;j--)
            {
                // If the previous region overlaps with the current one.
                if(se.getRegion(j).getEnd()>=reg.getStart())
                {
                    if(se.getRegion(j).getStart()==reg.getStart())
                        xs[i] = xs[j];
                    rlm[ls[j]]=1;
                }
                else if(se.getRegion(j).getEnd()+1==reg.getStart())
                    xs[i] = xs[j]+ws[j];
            }
            for(int k=0;k<rlm.length;k++)
            {
                if(rlm[k]==0)
                {
                    ls[i]=k;
                    break;
                }
            }
            // Calculate the X-coordinate of the region rectangle.
            if(xs[i]==0)
                xs[i] = (i==0) ? 1 : (int)((float)reg.getStart()/(float)iTotalLength*(bim.getWidth()-1));
            paintRegion(g, xs[i], y+ls[i]*RECT_HEIGHT, ws[i], RECT_HEIGHT, c);
        }
        repaint();
    }

    /**
     *  Paints the region.
     *
     *  @param g
     *  @param x
     *  @param y
     *  @param w
     *  @param h
     *  @param c
     */
    private void paintRegion(Graphics2D g, int x, int y, int w, int h, Color c)
    {
        Color tmp = g.getColor();
        Color from = Color.WHITE;
        GradientPaint gp = new GradientPaint(x, y, from, x, y+h, c);
        g.setPaint(gp);
        if(x+w>=bim.getWidth()-2)
            w = bim.getWidth()-x-2;
        g.fillRect(x, y, w, h);
        g.setColor(Color.BLACK);
        g.drawRect(x, y, w, h);
        g.setColor(tmp);
    }

    /**
     *  Paints the legend and returns the height of the header image.
     *
     *  @param g
     *  @return
     */
    private int paintLegend(Graphics2D g)
    {
        Color tmp = g.getColor();
        int h = 10;
        int w = 20;
        g.setColor(Color.BLACK);
        String[] regs = {GeneRegion.UTR5, GeneRegion.EXON, GeneRegion.INTRON,
                         GeneRegion.UTR3, GeneRegion.INTERGENIC, GeneRegion.UNNAMED,
                         GeneRegion.mRNA, "Other"};
        int y = OFFSET;
        for(String s:regs)
        {
            paintRegion(g, 1, y, w, h, getColor(s));
            g.drawString(s, w+5, y+h);
            y += h+5;
        }
        g.setColor(tmp);
        return y+10;
    }

    /**
     *  Returns the color of each region type.
     *
     *  @param strType
     *  @return
     */
    private Color getColor(String strType)
    {
        if(strType.equalsIgnoreCase(GeneRegion.EXON))
            return EXON;
        else if(strType.equalsIgnoreCase(GeneRegion.INTRON))
            return INTRON;
        else if(strType.equalsIgnoreCase(GeneRegion.UTR5))
            return UTR5;
        else if(strType.equalsIgnoreCase(GeneRegion.UTR3))
            return UTR3;
        else if(strType.equalsIgnoreCase(GeneRegion.INTERGENIC))
            return INTERGENIC;
        else if(strType.equalsIgnoreCase(GeneRegion.UNNAMED))
            return UNNAMED;
        else if(strType.equalsIgnoreCase(GeneRegion.mRNA))
            return mRNA;
        else if(strType.equalsIgnoreCase(IAligner.AMBIGUOUS_REGION))
            return AMBIGUOUS;
        else
            return OTHER;
    }
}
