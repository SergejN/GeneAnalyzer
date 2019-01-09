/*
    File:
        WaitDialog.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Displays animation while the host-application performs a time-consuming
        operation.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import javax.swing.JComponent;

public class WaitDialog extends JComponent implements IWaitDialog
{
    private class TicksRenderer implements Runnable
    {
        private boolean bActive = true;

        public void inactivate()
        {
            bActive = false;
        }

        public void run()
        {
            while(bActive)
            {
                repaint();
                try
                {
                    Thread.sleep(100);
                }
                catch(Exception e){}
            }
        }
    };

    private class Tick
    {
        private Area area;
        private Color clFrom;
        private Color clTo;
        private int iStep;

        private Color color;

        public Tick()
        {
            this.clFrom = COLOR_FROM;
            this.clTo = COLOR_TO;
            this.iStep = COLOR_STEP;
            this.color = COLOR_TO;
        }

        public void setArea(Area a)
        {
            area = a;
        }

        public void setColor(Color c)
        {
            this.color = c;
        }

        public void paintActive(Graphics2D g)
        {
            color = clFrom;
            g.setColor(color);
            g.fill(area);
        }

        public void paintInactive(Graphics2D g)
        {
            // First, update the color if necessary.
            if(color.getRed()<clTo.getRed())
            {
                color = new Color(color.getRed()+iStep,
                                  color.getGreen()+iStep,
                                  color.getBlue()+iStep);
            }
            g.setColor(color);
            g.fill(area);
        }
    };

    private final int TICKERS_COUNT     = 16;
    private final Color COLOR_FROM      = new Color(75,75,75);
    private final Color COLOR_TO        = new Color(215,215,215);
    private final int   COLOR_STEP      = 35;

    private RenderingHints hints        = null;
    private int iActiveTick             = 0;
    private int nInstances              = 0;

    private Tick[] tickers              = null;

    private TicksRenderer tr            = null;

    private boolean bRunning            = false;

    public WaitDialog()
    {
        addMouseListener( new MouseAdapter() {} );
        addMouseMotionListener( new MouseMotionAdapter() {} );
        hints = new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        hints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        hints.put(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
        tickers = new Tick[TICKERS_COUNT];
        for(int i=0;i<TICKERS_COUNT;i++)
            tickers[i] = new Tick();
        iActiveTick = TICKERS_COUNT-1;
    }

    public void show(TYPE type)
    {
        if(!bRunning)
        {
            setVisible(true);
            tr = new TicksRenderer();
            new Thread(tr).start();
            bRunning = true;
        }
        nInstances++;
    }

    public void setText(String strMainText, String strHintText)
    {
    }

    @Override
    public void paint(Graphics g)
    {
        Graphics2D g2d = (Graphics2D)g;
        g2d.addRenderingHints(hints);
        // Paint background.
        g2d.setColor(new Color(225,225,225,170));
        g2d.fillRect(0, 0, getWidth(), getHeight());
        Area[] areas = createTickers();
        for(int i=0;i<TICKERS_COUNT;i++)
        {
            tickers[i].setArea(areas[i]);
            if(i==iActiveTick)
                tickers[i].paintActive(g2d);
            else
                tickers[i].paintInactive(g2d);
        }
        iActiveTick--;
        if(iActiveTick<0)
            iActiveTick = TICKERS_COUNT-1;
    }

    public void close()
    {
        if(nInstances==0)
            return;
        nInstances--;
        if(nInstances==0)
        {
            tr.inactivate();
            setVisible(false);
            bRunning = false;
        }
    }

    private Area[] createTickers()
    {
        Area prim = createTickPrimitive();
        Area[] areas = new Area[TICKERS_COUNT];
        // Image center.
        Point2D.Float center = new Point2D.Float(getWidth()*0.5f, getHeight()*0.5f);
        double angle = 2.0*Math.PI/((double)TICKERS_COUNT);
        double th = getHeight()*0.1;        // Tick height.
        for(int i=0;i<TICKERS_COUNT;i++)
        {
            Area tick = (Area)prim.clone();
            AffineTransform at_center = AffineTransform.getTranslateInstance(center.getX(), center.getY());
            AffineTransform at_border = AffineTransform.getTranslateInstance(th, -10.0);
            AffineTransform at_circle = AffineTransform.getRotateInstance(-i*angle, center.getX(), center.getY());
            AffineTransform at_wheel = new AffineTransform();
            at_wheel.concatenate(at_center);
            at_wheel.concatenate(at_border);
            tick.transform(at_wheel);
            tick.transform(at_circle);
            areas[i] = tick;
        }
        return areas;
    }

    /**
     *  Creates a shape representing one tick.
     *  The tick is a rounded rectangle with the length  of 10% of the parent
     *  height and the width of 20% of the length.
     *
     *  @return
     */
    private Area createTickPrimitive()
    {
        Area tick = new Area();
        int l = (int)(getHeight()*0.1f);                                // Total tick length.
        int h = (int)(l*0.2f);                                          // Tick height.
        Rectangle2D rect = new Rectangle.Float(h*0.5f, 0.0f, l-h, h);   // Rectangular body
        tick.add(new Area(rect));
        Ellipse2D circl = new Ellipse2D.Float(0.0f, 0.0f, h, h);        // Left half-circle
        tick.add(new Area(circl));
        Ellipse2D circr = new Ellipse2D.Float(l-h, 0.0f, h, h);         // Right half-circle
        tick.add(new Area(circr));
        return tick;
    }
}
