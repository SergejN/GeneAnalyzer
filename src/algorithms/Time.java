/*
    File:
        Time.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Represents time with days, hours, minutes, seconds and milliseconds.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;


public class Time
{
    private int nDays       = 0;
    private int nHours      = 0;
    private int nMinutes    = 0;
    private int nSeconds    = 0;
    private int nMillis     = 0;

    /**
     *  Make the constructor private to avoid object creation by user.
     */
    private Time()
    {}

    public int getDays()
    {
        return nDays;
    }

    public int getHours()
    {
        return nHours;
    }

    public int getMinutes()
    {
        return nMinutes;
    }

    public int getSeconds()
    {
        return nSeconds;
    }

    public int getMillis()
    {
        return nMillis;
    }

    @Override
    public String toString()
    {
        StringBuffer sb = new StringBuffer();
        if(nDays>0)
            sb.append((nDays==1) ? "1 day, " : String.format("%d days, ", nDays));
        if(nHours>0)
            sb.append((nHours==1) ? "1 hour, " : String.format("%d hours, ", nHours));
        if(nMinutes>0)
            sb.append((nMinutes==1) ? "1 minute, " : String.format("%d minutes, ", nMinutes));
        if(nSeconds>0)
            sb.append((nSeconds==1) ? "1 second, " : String.format("%d seconds, ", nSeconds));
        if(nMillis>0)
            sb.append((nMillis==1) ? "1 millisecond, " : String.format("%d milliseconds, ", nMillis));
        if(sb.toString().endsWith(", "))
            sb.delete(sb.length()-2, sb.length());
        return sb.toString();
    }

    public static Time constructFromMillis(long millis)
    {
        Time time = new Time();
        // Days.
        time.nDays = (int)(millis/(1000*3600*24));
        millis -= time.nDays*(1000*3600*24);
        // Hours.
        time.nHours = (int)(millis/(1000*3600));
        millis -= time.nHours*(1000*3600);
        // Minutes.
        time.nMinutes = (int)(millis/(1000*60));
        millis -= time.nMinutes*(1000*60);
        // Seconds.
        time.nSeconds = (int)(millis/1000);
        millis -= time.nSeconds*1000;
        // Milliseconds.
        time.nMillis = (int)millis;
        return time;
    }
}
