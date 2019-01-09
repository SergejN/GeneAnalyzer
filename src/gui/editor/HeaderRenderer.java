/*
    File:
        HeaderRenderer.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Renders the header containing the names of the strains.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.editor;

import bio.gene.GeneEntry;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import javax.swing.JComponent;

public class HeaderRenderer extends JComponent
{
    private final int CELLHEIGHT    = 15;
    private final int OFFSET_LEFT   = 10;
    private final int COORDHEIGHT   = 20;
    private final int IMGWIDTH      = OFFSET_LEFT+10*10;
    private final int RADIUS        = 5;

    private BufferedImage img       = null;

    
    public HeaderRenderer(GeneEntry ge)
    {
        int nStrains = ge.getStrainsCount();
        int h = nStrains*CELLHEIGHT+2*COORDHEIGHT;
        setPreferredSize(new Dimension(IMGWIDTH,h));
        img = new BufferedImage(IMGWIDTH, h, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = (Graphics2D)img.getGraphics();
        // Paint eader background.
        Area hdr = new Area();
        Rectangle2D rect = new Rectangle.Float(RADIUS, 0, IMGWIDTH-RADIUS, h-1);         // Body rectangle
        Ellipse2D circt = new Ellipse2D.Float(1, 0, 2*RADIUS, 2*RADIUS);            // Top circle
        Ellipse2D circb = new Ellipse2D.Float(1, h-2*RADIUS-1, 2*RADIUS, 2*RADIUS);   // bottom circle
        Rectangle2D rl = new Rectangle.Float(1, RADIUS, RADIUS, h-2*RADIUS);        // Rectangle left
        hdr.add(new Area(rect));
        hdr.add(new Area(circt));
        hdr.add(new Area(circb));
        hdr.add(new Area(rl));
        GradientPaint hdrgrad = new GradientPaint(0, 0, new Color(255,255,255), 0, img.getHeight(), new Color(235,235,235));
        g.setPaint(hdrgrad);
        g.fill(hdr);
        g.setColor(new Color(205,190,190));
        g.draw(hdr);
        // Paint lines and hints.
        g.setColor(new Color(0,0,0,128));
        g.setFont(new Font("Times New Roman", Font.PLAIN, 12));
        int y = COORDHEIGHT;
        g.drawLine(3, y, IMGWIDTH-6, y);
        g.drawString("Gene position", 10, y-4);
        y += nStrains*CELLHEIGHT+3;
        g.drawLine(3, y, IMGWIDTH-6, y);
        g.drawString("CDS position", 10, y+12);
        // Paint strain names.
        g.setColor(Color.BLUE);
        g.setFont(new Font("Times New Roman", Font.PLAIN, 14));
        for(int i=0;i<nStrains;i++)
        {
            g.drawString(String.format("%-10s", ge.getStrainEntry(i).getStrainName()),
                            OFFSET_LEFT, COORDHEIGHT+(i+1)*CELLHEIGHT);
        }        
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

    @Override
    public void paint(Graphics g)
    {
        g.drawImage(img, 0, 0, null);
    }
}
